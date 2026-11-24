"""切り捨て形式的冪級数 (truncated formal power series) の基本演算。

級数は「係数のリスト」で表す::

    a = [a_0, a_1, ..., a_N]   <->   f(x) = a_0 + a_1 x + ... + a_N x^N + O(x^(N+1))

係数はすべて sympy の厳密な数 (Rational や根号を含む式) として扱うので、
浮動小数点誤差は入らない。
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

import sympy as sp

Series = List[sp.Expr]

__all__ = [
    "Series",
    "clean",
    "to_series",
    "truncate",
    "identity",
    "mul",
    "compose",
    "iterate",
    "inverse",
    "from_derivatives",
    "to_derivatives",
    "from_expr",
    "parse_number",
]


# ---------------------------------------------------------------------------
# 数値の整形
# ---------------------------------------------------------------------------

def clean(value: sp.Expr) -> sp.Expr:
    """係数を正規形に近づける。

    有理数ならそのまま返す (高速)。根号などを含む場合は展開・有理化して
    式が膨らむのを防ぐ。
    """
    value = sp.sympify(value)
    if value.is_Rational:
        return value
    return sp.radsimp(sp.expand(value))


def parse_number(text: str) -> sp.Expr:
    """``'1'``, ``'-1/2'``, ``'0.5'``, ``'sqrt(2)'`` などを厳密な sympy の数に変換する。"""
    text = text.strip()
    if not text:
        raise ValueError("空の値は数として解釈できません")
    try:
        value = sp.nsimplify(text, rational=True)
    except (sp.SympifyError, TypeError, SyntaxError) as exc:
        raise ValueError(f"数として解釈できません: {text!r}") from exc
    if value.free_symbols:
        raise ValueError(f"変数を含む値は指定できません: {text!r}")
    return value


# ---------------------------------------------------------------------------
# 生成・変換
# ---------------------------------------------------------------------------

def to_series(coeffs: Iterable) -> Series:
    """任意の数値の並びを sympy の係数リストに変換する。"""
    return [sp.sympify(c) for c in coeffs]


def truncate(a: Sequence, N: int) -> Series:
    """次数 N までに切り詰める (足りなければ 0 で埋める)。"""
    out = [sp.sympify(c) for c in a[: N + 1]]
    out.extend([sp.Integer(0)] * (N + 1 - len(out)))
    return out


def identity(N: int) -> Series:
    """恒等写像 x の級数。"""
    out = [sp.Integer(0)] * (N + 1)
    if N >= 1:
        out[1] = sp.Integer(1)
    return out


def from_derivatives(derivs: Sequence) -> Series:
    """f^(0)(0), ..., f^(N)(0) からテイラー係数 a_n = f^(n)(0) / n! を作る。"""
    return [clean(sp.sympify(d) / sp.factorial(n)) for n, d in enumerate(derivs)]


def to_derivatives(coeffs: Sequence) -> Series:
    """テイラー係数から n 階微分係数 f^(n)(0) = n! a_n を作る。"""
    return [clean(sp.sympify(c) * sp.factorial(n)) for n, c in enumerate(coeffs)]


def from_expr(expr, N: int, var: sp.Symbol | str = "x") -> Series:
    """sympy の式 (または文字列) を x=0 のまわりで N 次までテイラー展開する。

    >>> from_expr("sinh(x)", 5)
    [0, 1, 0, 1/6, 0, 1/120]
    """
    x = sp.Symbol(var) if isinstance(var, str) else var
    if isinstance(expr, str):
        expr = sp.sympify(expr, locals={x.name: x})
    expr = sp.sympify(expr)
    extra = expr.free_symbols - {x}
    if extra:
        names = ", ".join(sorted(s.name for s in extra))
        raise ValueError(f"変数 {x} 以外の記号が含まれています: {names}")
    poly = sp.series(expr, x, 0, N + 1).removeO()
    poly = sp.expand(poly)
    # x^x = 1 + x log x + ... や sqrt(x) のように x=0 で解析的でない関数は
    # 展開に log x や分数べきが残るので、テイラー展開できないとしてはじく
    if poly.has(sp.Order) or not poly.is_polynomial(x):
        raise ValueError(
            f"{str(expr).replace('**', '^')} は x=0 でテイラー展開できません "
            "(log x や x の分数乗などが現れます)"
        )
    return [clean(poly.coeff(x, n)) for n in range(N + 1)]


# ---------------------------------------------------------------------------
# 演算
# ---------------------------------------------------------------------------

def mul(a: Sequence, b: Sequence, N: int) -> Series:
    """積 a(x) * b(x) を次数 N で切り捨てる。"""
    out = [sp.Integer(0)] * (N + 1)
    for i, ai in enumerate(a[: N + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: N + 1 - i]):
            if bj == 0:
                continue
            out[i + j] += ai * bj
    return [clean(c) for c in out]


def compose(outer: Sequence, inner: Sequence, N: int | None = None) -> Series:
    """合成 outer(inner(x)) を次数 N で切り捨てる (Horner 法)。

    inner(0) = 0 であることが必要 (そうでないと形式的に合成できない)。
    """
    if N is None:
        N = min(len(outer), len(inner)) - 1
    if inner and sp.sympify(inner[0]) != 0:
        raise ValueError("内側の級数の定数項が 0 ではないため合成できません")
    outer = truncate(outer, N)
    inner = truncate(inner, N)
    result = [sp.Integer(0)] * (N + 1)
    for k in range(N, -1, -1):
        result = mul(result, inner, N)
        result[0] = clean(result[0] + outer[k])
    return result


def iterate(a: Sequence, k: int, N: int | None = None) -> Series:
    """k 回反復 f^[k] = f∘f∘...∘f を求める。

    k = 0 なら恒等写像、k < 0 なら逆関数の |k| 回反復。
    """
    if N is None:
        N = len(a) - 1
    if k < 0:
        return iterate(inverse(a, N), -k, N)
    result = identity(N)
    for _ in range(k):
        result = compose(a, result, N)
    return result


def inverse(a: Sequence, N: int | None = None) -> Series:
    """合成に関する逆関数 (級数の反転) を求める。f(0)=0, f'(0)≠0 が必要。"""
    if N is None:
        N = len(a) - 1
    a = truncate(a, N)
    if a[0] != 0:
        raise ValueError("f(0) ≠ 0 なので逆関数の級数は作れません")
    if N >= 1 and sp.simplify(a[1]) == 0:
        raise ValueError("f'(0) = 0 なので逆関数の級数は作れません")
    c = [sp.Integer(0)] * (N + 1)
    if N == 0:
        return c
    c[1] = clean(1 / a[1])
    # [x^n] f(c(x)) は c_n について a_1 * c_n + (c_1..c_{n-1} のみに依存する項)
    for n in range(2, N + 1):
        h = compose(a, c[: n + 1], n)
        c[n] = clean(-h[n] / a[1])
    return c
