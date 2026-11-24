"""コマンドラインインターフェース。

使い方の例::

    python -m func_sqrt "sinh(x)" -n 11
    python -m func_sqrt --derivs 0,1,0,1,0,1 --format latex
    python -m func_sqrt "exp(x)-1" -n 8 --order 3
    python -m func_sqrt "2*x + x**2" -n 6 --iterate 3/2
    python -m func_sqrt            # 対話モード
"""

from __future__ import annotations

import argparse
import sys
from typing import Callable, Sequence

import sympy as sp

from . import __version__
from .core import FunctionalRootError, fractional_iterate, functional_root, verify
from .desmos import DesmosError, build_state, create_link
from .formatting import format_json, format_series, format_table
from .series import Series, from_derivatives, from_expr, parse_number

DEFAULT_DEGREE = 10

FORMATS = ("text", "series", "latex", "table", "json")


def _parse_list(text: str) -> list[sp.Expr]:
    items = [s for s in text.replace(";", ",").split(",")]
    if any(not s.strip() for s in items):
        raise ValueError(f"空の要素があります: {text!r}")
    return [parse_number(s) for s in items]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="func-sqrt",
        description=(
            "f(0)=0, f'(0)≠0 を満たす関数 f について、関数的平方根 g∘g = f "
            "(および関数的 m 乗根・分数回反復) のテイラー展開を厳密計算します。"
        ),
        epilog=(
            "例: func-sqrt \"sinh(x)\" -n 11 / "
            "func-sqrt --derivs 0,1,0,1,0,1 / "
            "func-sqrt \"exp(x)-1\" -n 8 --order 3"
        ),
    )
    src = p.add_argument_group("入力 (いずれか 1 つ。何も指定しなければ対話モード)")
    src.add_argument(
        "expr",
        nargs="?",
        help='f(x) の式 (sympy 記法)。例: "sinh(x)", "exp(x)-1", "2*x + x**2"',
    )
    src.add_argument(
        "-d", "--derivs",
        metavar="LIST",
        help="f^(0)(0), f^(1)(0), ..., f^(N)(0) をカンマ区切りで (例: 0,1,0,1)",
    )
    src.add_argument(
        "-c", "--coeffs",
        metavar="LIST",
        help="テイラー係数 a_0, a_1, ..., a_N をカンマ区切りで (例: 0,1,0,1/6)",
    )
    src.add_argument("-i", "--interactive", action="store_true", help="対話モードで入力する")

    calc = p.add_argument_group("計算")
    calc.add_argument(
        "-n", "--degree",
        type=int,
        default=None,
        help=f"展開する最大次数 N (式で入力したとき。既定: {DEFAULT_DEGREE})",
    )
    mode = calc.add_mutually_exclusive_group()
    mode.add_argument(
        "-m", "--order",
        type=int,
        default=2,
        help="関数的 m 乗根 g^[m] = f を求める (既定: 2 = 関数的平方根)",
    )
    mode.add_argument(
        "-t", "--iterate",
        metavar="P/Q",
        help="分数回反復 f^[P/Q] を求める (例: 1/2, 3/2, -1)",
    )
    calc.add_argument(
        "--multiplier",
        metavar="VALUE",
        help="g'(0) に使う値 (既定: f'(0) の主値の m 乗根)。例: -1",
    )
    calc.add_argument("--var", default="x", help="変数名 (既定: x)")

    out = p.add_argument_group("出力")
    out.add_argument(
        "-f", "--format",
        choices=FORMATS,
        default="text",
        help="出力形式 (既定: text = 表 + 近似式)",
    )
    out.add_argument(
        "--digits",
        type=int,
        metavar="K",
        help="表に K 桁の数値近似を併記する",
    )
    out.add_argument(
        "--desmos",
        action="store_true",
        help="f と g を描いた Desmos グラフの共有リンクを作る (インターネット接続が必要)",
    )
    out.add_argument(
        "--check",
        action="store_true",
        help="結果を m 回合成して f に一致するか検証する",
    )
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
    return p


# ---------------------------------------------------------------------------
# 入力
# ---------------------------------------------------------------------------

def _ask(prompt: str, convert: Callable, input_fn: Callable[[str], str]):
    while True:
        raw = input_fn(prompt)
        try:
            return convert(raw)
        except (ValueError, TypeError, sp.SympifyError) as exc:
            print(f"    入力形式が正しくありません ({exc})。もう一度入力してください。")


def _interactive_input(var: str, input_fn: Callable[[str], str] = input) -> Series:
    print("=== 関数的平方根 計算ツール (厳密計算) ===")
    print("g(g(x)) = f(x) となる g(x) を、x=0 でのテイラー展開として求めます。")
    print("(前提: f(0)=0, f'(0)≠0)\n")

    text = input_fn(f"f({var}) の式を入力 (例: sinh({var}))。空のまま Enter で微分係数を入力: ")
    text = text.strip()

    def to_degree(s: str) -> int:
        n = int(s)
        if n < 1:
            raise ValueError("N は 1 以上にしてください")
        return n

    if text:
        N = _ask("最大次数 N: ", to_degree, input_fn)
        return from_expr(text, N, var)

    N = _ask("最大次数 N (0階〜N階の微分係数を使います): ", to_degree, input_fn)
    print(f"\nf^(0)(0) から f^({N})(0) までを順に入力してください。")
    print("整数・小数・分数 (例: -1/2)・根号 (例: sqrt(2)) が使えます。\n")
    derivs = [_ask(f"  f^({n})(0) = ", parse_number, input_fn) for n in range(N + 1)]
    return from_derivatives(derivs)


def read_input(args: argparse.Namespace, input_fn: Callable[[str], str] = input) -> Series:
    """入力を読んでテイラー係数を返す。式で入力されたときは args.f_expr にその式を入れる。"""
    args.f_expr = None
    given = [x for x in (args.expr, args.derivs, args.coeffs) if x is not None]
    if len(given) > 1 or (given and args.interactive):
        raise ValueError("入力は 式 / --derivs / --coeffs / --interactive のうち 1 つだけ指定してください")
    if args.derivs is not None:
        return from_derivatives(_parse_list(args.derivs))
    if args.coeffs is not None:
        return _parse_list(args.coeffs)
    if args.expr is not None:
        N = DEFAULT_DEGREE if args.degree is None else args.degree
        if N < 1:
            raise ValueError("--degree は 1 以上にしてください")
        args.f_expr = sp.sympify(args.expr, locals={args.var: sp.Symbol(args.var)})
        return from_expr(args.f_expr, N, args.var)
    return _interactive_input(args.var, input_fn)


# ---------------------------------------------------------------------------
# 本体
# ---------------------------------------------------------------------------

def _describe(args: argparse.Namespace) -> str:
    v = args.var
    if args.iterate is not None:
        return f"f^[{args.iterate}]({v})"
    if args.order == 2:
        return f"g({v})  (g∘g = f)"
    return f"g({v})  (g^[{args.order}] = f)"


def run(args: argparse.Namespace, input_fn: Callable[[str], str] = input) -> int:
    a = read_input(args, input_fn)
    multiplier = None if args.multiplier is None else parse_number(args.multiplier)

    if args.iterate is not None:
        t = parse_number(args.iterate)
        if not t.is_Rational:
            raise ValueError(f"--iterate には有理数を指定してください: {args.iterate}")
        b = fractional_iterate(a, t, multiplier=multiplier)
        meta = {"iterate": str(t)}
    else:
        b = functional_root(a, args.order, multiplier=multiplier)
        meta = {"order": args.order}

    fmt = args.format
    if fmt == "json":
        meta["input_coefficients"] = [str(c) for c in a]
        print(format_json(b, **meta))
    elif fmt == "series":
        print(format_series(b, args.var))
    elif fmt == "latex":
        print(format_series(b, args.var, style="latex"))
    elif fmt == "table":
        print(format_table(b, digits=args.digits))
    else:
        print(f"\nf({args.var}) ~ {format_series(a, args.var)}")
        print(f"\n--- 結果: {_describe(args)} ---")
        print(format_table(b, digits=args.digits))
        print(f"\ng({args.var}) ~ {format_series(b, args.var)}")

    if args.desmos:
        try:
            state = build_state(
                a, b, f_expr=getattr(args, "f_expr", None),
                order=None if args.iterate is not None else args.order,
                iterate=args.iterate,
            )
            print(f"\nDesmos: {create_link(state)}")
        except DesmosError as exc:
            print(f"\nDesmos: リンクを作れませんでした ({exc})", file=sys.stderr)

    if args.check:
        if args.iterate is not None:
            t = sp.Rational(parse_number(args.iterate))
            ok = t.p == 1 and all(c == 0 for c in verify(a, b, int(t.q)))
            if t.p != 1:
                print("\n検証: --iterate では分子が 1 (= 関数的 q 乗根) のときのみ検証できます", file=sys.stderr)
                return 0
        else:
            ok = all(c == 0 for c in verify(a, b, args.order))
        print("\n検証: " + ("OK (合成すると f に一致)" if ok else "NG (f に一致しません)"))
        return 0 if ok else 1
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (FunctionalRootError, ValueError, ZeroDivisionError, sp.SympifyError) as exc:
        print(f"エラー: {exc}", file=sys.stderr)
        return 2
    except (KeyboardInterrupt, EOFError):
        print("\n中断しました", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
