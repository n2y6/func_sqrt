"""Web API 用のリクエスト処理 (フレームワーク非依存)。

:func:`handle` にクエリパラメータの辞書を渡すと ``(HTTP ステータス, JSON 用 dict)``
を返す。CGI などの薄いラッパーから呼び出して使う (``deploy/cgi/`` 参照)。

外部からの入力を扱うため、式は許可した記号・関数名だけからなる場合に限り
パースする (sympy の ``sympify`` は内部で ``eval`` を使うので直接渡さない)。
"""

from __future__ import annotations

import re
from typing import Mapping

import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication,
    parse_expr,
    standard_transformations,
)

from . import __version__
from .core import FunctionalRootError, functional_root, principal_root
from .desmos import DesmosError, build_state, create_link
from .formatting import format_series
from .series import from_derivatives, from_expr, iterate, to_derivatives

__all__ = ["ApiError", "LIMITS", "safe_parse", "handle", "usage"]

LIMITS = {
    "max_degree": 20,
    "max_order": 6,
    "max_iterate_part": 6,
    "max_input_length": 300,
    "max_list_items": 21,
}

X = sp.Symbol("x")

_FUNCTIONS = {
    name: getattr(sp, name)
    for name in (
        "sin cos tan sinh cosh tanh asin acos atan asinh acosh atanh "
        "exp log sqrt"
    ).split()
}
_CONSTANTS = {"pi": sp.pi, "E": sp.E, "e": sp.E, "I": sp.I}
_ALIASES = {
    "ln": sp.log,
    "arcsin": sp.asin, "arccos": sp.acos, "arctan": sp.atan,
    "arcsinh": sp.asinh, "arccosh": sp.acosh, "arctanh": sp.atanh,
}
_NORMALIZE = str.maketrans({"−": "-", "–": "-", "×": "*", "÷": "/", "・": "*"})

_ALLOWED_CHARS = re.compile(r"[0-9A-Za-z+\-*/^().\s]*")
_IDENT = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")


class ApiError(ValueError):
    """クライアント側の入力に問題がある (HTTP 400)。"""


def safe_parse(text: str, *, allow_x: bool = True) -> sp.Expr:
    """許可された関数・定数・変数 x だけからなる式を sympy の式に変換する。"""
    if text is None:
        raise ApiError("値がありません")
    text = text.translate(_NORMALIZE).strip()
    if not text:
        raise ApiError("空の式です")
    if len(text) > LIMITS["max_input_length"]:
        raise ApiError(f"式が長すぎます (最大 {LIMITS['max_input_length']} 文字)")
    if not _ALLOWED_CHARS.fullmatch(text):
        raise ApiError("使えない文字が含まれています")
    allowed = dict(_FUNCTIONS)
    allowed.update(_ALIASES)
    allowed.update(_CONSTANTS)
    if allow_x:
        allowed["x"] = X
    for name in _IDENT.findall(text):
        if name not in allowed:
            raise ApiError(f"使えない名前です: {name}")
    global_dict = {
        "__builtins__": {},
        "Integer": sp.Integer,
        "Float": sp.Float,
        "Rational": sp.Rational,
        "Symbol": sp.Symbol,
    }
    global_dict.update(allowed)
    try:
        expr = parse_expr(
            text.replace("^", "**"),
            local_dict=allowed,
            global_dict=global_dict,
            # 2x, 3sin(x) のような掛け算の省略も受け付ける
            transformations=standard_transformations + (implicit_multiplication,),
        )
    except Exception as exc:  # 構文エラーなど
        raise ApiError(f"式を解釈できません: {text}") from exc
    if not isinstance(expr, sp.Basic):
        raise ApiError(f"式を解釈できません: {text}")
    if expr.free_symbols - {X}:
        raise ApiError("x 以外の変数は使えません")
    return expr


def _number(text: str) -> sp.Expr:
    value = safe_parse(text, allow_x=False)
    if isinstance(value, sp.Float):
        value = sp.nsimplify(value, rational=True)
    return value


def _number_list(text: str) -> list[sp.Expr]:
    items = text.split(",")
    if len(items) > LIMITS["max_list_items"]:
        raise ApiError(f"値が多すぎます (最大 {LIMITS['max_list_items']} 個)")
    return [_number(s) for s in items]


def _int_param(params: Mapping[str, str], name: str, default: int, lo: int, hi: int) -> int:
    raw = params.get(name)
    if raw is None or raw == "":
        return default
    try:
        value = int(raw)
    except ValueError as exc:
        raise ApiError(f"{name} は整数で指定してください") from exc
    if not lo <= value <= hi:
        raise ApiError(f"{name} は {lo} 以上 {hi} 以下で指定してください")
    return value


def _branch_multiplier(params: Mapping[str, str], a, q: int, multiplier):
    """branch=-1 なら主値の m 乗根の符号を反転したものを g'(0) に使う。"""
    if multiplier is not None or params.get("branch", "") not in ("-1", "neg"):
        return multiplier
    if len(a) < 2 or a[1] == 0:
        return None
    if q % 2 == 1:
        raise ApiError(f"{q} 乗根は実数では 1 つだけなので、g'(0) の負の枝は選べません")
    return -principal_root(a[1], q)


def _floats(coeffs) -> list | None:
    try:
        return [float(sp.N(c)) for c in coeffs]
    except TypeError:  # 虚数を含む
        return None


SAMPLE_RANGE = (-5.0, 5.0, 801)


def _samples(expr) -> dict | None:
    """グラフ用に f(x) を等間隔の点で数値計算する (定義されない点は null)。"""
    import math

    lo, hi, count = SAMPLE_RANGE
    try:
        fn = sp.lambdify(X, expr, modules="math")
    except Exception:
        return None
    dx = (hi - lo) / (count - 1)
    ys = []
    for i in range(count):
        try:
            y = float(fn(lo + i * dx))
            ys.append(round(y, 8) if math.isfinite(y) and abs(y) < 1e9 else None)
        except (ValueError, ZeroDivisionError, OverflowError, TypeError):
            ys.append(None)
    return {"x0": lo, "dx": dx, "y": ys}


def usage() -> dict:
    """API の説明 (/func_sqrt にパラメータなしでアクセスしたときに返す)。"""
    return {
        "name": "func_sqrt API",
        "version": __version__,
        "description": "g(g(x)) = f(x) となる g (関数的平方根) などのテイラー展開を厳密計算します。",
        "endpoint": "/func_sqrt",
        "parameters": {
            "f": "f(x) の式。例: sinh(x), exp(x)-1, 2*x+x^2",
            "derivs": "f^(0)(0),...,f^(N)(0) をカンマ区切りで (f の代わり)",
            "coeffs": "テイラー係数 a_0,...,a_N をカンマ区切りで (f の代わり)",
            "n": f"展開次数 N (f を使うとき。既定 10、最大 {LIMITS['max_degree']})",
            "order": f"関数的 m 乗根の m (既定 2、最大 {LIMITS['max_order']})",
            "iterate": "分数回反復 p/q (order の代わり。例: 1/2, 3/2, -1)",
            "multiplier": "g'(0) に使う値 (省略時は主値)",
            "branch": "-1 にすると g'(0) に主値の符号を反転したものを使う (multiplier を省略したとき)",
            "samples": "1 にするとグラフ用に f(x) の値 (x = -5〜5) を f_samples に入れて返す",
            "desmos": "1 にすると f と g を描いた Desmos グラフのリンクを desmos に入れて返す",
        },
        "functions": sorted(_FUNCTIONS) + sorted(_ALIASES) + sorted(_CONSTANTS),
        "example": "/func_sqrt?f=sinh(x)&n=11",
        "source": "https://github.com/n2y6/func_sqrt",
    }


def handle(params: Mapping[str, str], *, desmos_link=create_link) -> tuple[int, dict]:
    """クエリパラメータを受け取り、(HTTP ステータス, レスポンス dict) を返す。

    desmos_link は Desmos のリンクを作る関数 (テスト時に差し替えられる)。
    """
    try:
        return 200, _compute(params, desmos_link)
    except (ApiError, FunctionalRootError, ValueError, ZeroDivisionError) as exc:
        return 400, {"ok": False, "error": str(exc)}


def _compute(params: Mapping[str, str], desmos_link=create_link) -> dict:
    sources = [k for k in ("f", "derivs", "coeffs") if params.get(k)]
    if len(sources) != 1:
        raise ApiError("f, derivs, coeffs のうちどれか 1 つを指定してください")
    source = sources[0]
    max_n = LIMITS["max_degree"]

    f_expr = None
    if source == "f":
        n = _int_param(params, "n", 10, 1, max_n)
        f_expr = safe_parse(params["f"])
        a = from_expr(f_expr, n, X)
    else:
        values = _number_list(params[source])
        if len(values) < 2:
            raise ApiError("値は 2 個以上必要です")
        a = from_derivatives(values) if source == "derivs" else values

    multiplier = _number(params["multiplier"]) if params.get("multiplier") else None

    result: dict = {"ok": True, "input": {source: params[source]}}
    if params.get("iterate"):
        if params.get("order"):
            raise ApiError("order と iterate は同時に指定できません")
        t = _number(params["iterate"])
        if not t.is_Rational:
            raise ApiError("iterate は有理数 p/q で指定してください")
        t = sp.Rational(t)
        lim = LIMITS["max_iterate_part"]
        if abs(t.p) > lim or t.q > lim:
            raise ApiError(f"iterate の分子・分母は絶対値 {lim} 以下にしてください")
        q = int(t.q)
        multiplier = _branch_multiplier(params, a, q, multiplier)
        root = a if q == 1 else functional_root(a, q, multiplier=multiplier)
        b = iterate(root, int(t.p), len(a) - 1)
        result["iterate"] = str(t)
    else:
        q = _int_param(params, "order", 2, 1, LIMITS["max_order"])
        multiplier = _branch_multiplier(params, a, q, multiplier)
        root = b = functional_root(a, q, multiplier=multiplier)
        result["order"] = q

    result.update(
        {
            "f_series": format_series(a),
            "input_coefficients": [str(c) for c in a],
            "input_derivatives": [str(d) for d in to_derivatives(a)],
            "series": format_series(b),
            "latex": format_series(b, style="latex"),
            "coefficients": [str(c) for c in b],
            "derivatives": [str(d) for d in to_derivatives(b)],
            # グラフ・表示用の近似値 (虚数を含むときは null)
            "numeric": {
                "f": _floats(a),
                "result": _floats(b),
                # 検算用: q 乗根 (これを q 回合成すると f に戻る)
                "root": _floats(root),
                "root_order": q,
            },
        }
    )
    if f_expr is not None and params.get("samples", "").lower() in ("1", "true", "yes", "on"):
        result["f_samples"] = _samples(f_expr)

    if params.get("desmos", "").lower() in ("1", "true", "yes", "on"):
        try:
            state = build_state(
                a, b, f_expr=f_expr,
                order=result.get("order"), iterate=result.get("iterate"),
            )
            result["desmos"] = desmos_link(state)
        except DesmosError as exc:
            # 計算自体は成功しているので、リンク作成の失敗はエラー扱いにしない
            result["desmos"] = None
            result["desmos_error"] = str(exc)
    return result
