"""計算結果を人間・他プログラム向けの文字列にする。"""

from __future__ import annotations

import json
from typing import Sequence

import sympy as sp

from .series import to_derivatives

__all__ = ["format_series", "format_table", "format_json"]


def _split_sign(c: sp.Expr) -> tuple[str, sp.Expr]:
    c = sp.sympify(c)
    if c.could_extract_minus_sign():
        return "-", -c
    return "+", c


def _plain_term(c: sp.Expr, n: int, var: str) -> str:
    power = "" if n == 0 else (var if n == 1 else f"{var}^{n}")
    if not power:
        return str(c)
    if c == 1:
        return power
    if c.is_Integer or c.is_Symbol or isinstance(c, (sp.Pow, sp.Function)) and not c.is_Rational:
        return f"{c}{power}"
    return f"({c}){power}"


def _latex_term(c: sp.Expr, n: int, var: str) -> str:
    power = "" if n == 0 else (var if n == 1 else f"{var}^{{{n}}}")
    if not power:
        return sp.latex(c)
    if c == 1:
        return power
    body = sp.latex(c)
    if isinstance(c, sp.Add):
        body = rf"\left({body}\right)"
    return f"{body} {power}"


def format_series(
    coeffs: Sequence,
    var: str = "x",
    *,
    style: str = "plain",
    big_o: bool = True,
) -> str:
    """係数リストを ``x + (1/12)x^3 - (1/160)x^5 + O(x^6)`` のような式にする。

    style には ``"plain"`` (既定) か ``"latex"`` を指定する。
    """
    if style not in ("plain", "latex"):
        raise ValueError(f"未対応の style です: {style!r}")
    term = _plain_term if style == "plain" else _latex_term
    parts: list[str] = []
    for n, c in enumerate(coeffs):
        c = sp.sympify(c)
        if c == 0:
            continue
        sign, mag = _split_sign(c)
        text = term(mag, n, var)
        if not parts:
            parts.append(text if sign == "+" else f"-{text}")
        else:
            parts.append(f"{sign} {text}")
    if big_o:
        N = len(coeffs)
        order = f"{var}^{N}" if style == "plain" else f"{var}^{{{N}}}"
        parts.append(f"+ O({order})" if parts else f"O({order})")
    return " ".join(parts) if parts else "0"


def format_table(coeffs: Sequence, *, digits: int | None = None, var: str = "g") -> str:
    """n, b_n (テイラー係数), g^(n)(0) (微分係数) の表を作る。

    digits を指定すると数値近似の列を追加する。
    """
    derivs = to_derivatives(coeffs)
    rows = [("n", "b_n", f"{var}^(n)(0)")]
    if digits is not None:
        rows[0] += ("b_n (approx)",)
    for n in range(1, len(coeffs)):
        row = (str(n), str(coeffs[n]), str(derivs[n]))
        if digits is not None:
            row += (str(sp.N(coeffs[n], digits)),)
        rows.append(row)
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = []
    for k, row in enumerate(rows):
        lines.append("  ".join(cell.rjust(w) for cell, w in zip(row, widths)).rstrip())
        if k == 0:
            lines.append("  ".join("-" * w for w in widths))
    return "\n".join(lines)


def format_json(coeffs: Sequence, **meta) -> str:
    """係数と微分係数を JSON 文字列にする (値は sympy の文字列表記)。"""
    data = dict(meta)
    data["coefficients"] = [str(c) for c in coeffs]
    data["derivatives"] = [str(d) for d in to_derivatives(coeffs)]
    data["series"] = format_series(coeffs)
    data["latex"] = format_series(coeffs, style="latex")
    return json.dumps(data, ensure_ascii=False, indent=2)
