"""計算結果を Desmos のグラフ (共有リンク) にする。

f(x) と g(x)、検算用の g(g(x)) (m 乗根なら g を m 回合成したもの) とその誤差 |g(g(x)) - f(x)| を描いた
Desmos グラフを保存し、``https://www.desmos.com/calculator/<id>`` を返す。

注意: 保存には Desmos のサイト内部で使われている非公開 API
(``/api/v1/calculator-shared/save``) を使っている。公式に提供されたものでは
ないので、仕様変更で動かなくなる可能性がある。
"""

from __future__ import annotations

import base64
import json
import secrets
import string
import struct
import urllib.error
import urllib.request
import zlib
from typing import Sequence

import sympy as sp

from .formatting import format_series

__all__ = ["DesmosError", "series_latex", "expr_latex", "build_state", "create_link"]

SAVE_URL = "https://www.desmos.com/api/v1/calculator-shared/save"
CALC_URL = "https://www.desmos.com/calculator/"

COLORS = {"f": "#c74440", "g": "#2d70b3", "check": "#388c46", "error": "#6042a6"}


class DesmosError(RuntimeError):
    """Desmos への保存に失敗した。"""


# ---------------------------------------------------------------------------
# LaTeX (Desmos 向け)
# ---------------------------------------------------------------------------

_DESMOS_REPLACEMENTS = {
    r"\operatorname{asinh}": r"\operatorname{arcsinh}",
    r"\operatorname{acosh}": r"\operatorname{arccosh}",
    r"\operatorname{atanh}": r"\operatorname{arctanh}",
    r"\operatorname{arsinh}": r"\operatorname{arcsinh}",
    r"\operatorname{arcosh}": r"\operatorname{arccosh}",
    r"\operatorname{artanh}": r"\operatorname{arctanh}",
}


def _desmos_fix(latex: str) -> str:
    # sympy は \sin{\left(x \right)} のように書くが、Desmos は外側の {} を嫌う
    latex = latex.replace(r"{\left(", r"\left(").replace(r"\right)}", r"\right)")
    for old, new in _DESMOS_REPLACEMENTS.items():
        latex = latex.replace(old, new)
    return latex


def _check_real(coeffs: Sequence) -> None:
    for c in coeffs:
        if sp.sympify(c).has(sp.I):
            raise DesmosError("係数に虚数が含まれるため Desmos では描けません")


def series_latex(coeffs: Sequence) -> str:
    """テイラー係数を Desmos で使える多項式の LaTeX にする (O(x^n) は付けない)。"""
    _check_real(coeffs)
    return _desmos_fix(format_series(coeffs, style="latex", big_o=False))


def expr_latex(expr) -> str:
    """sympy の式を Desmos で使える LaTeX にする (log は自然対数 \\ln)。"""
    expr = sp.sympify(expr)
    if expr.has(sp.I):
        raise DesmosError("式に虚数が含まれるため Desmos では描けません")
    return _desmos_fix(sp.latex(expr, ln_notation=True, inv_trig_style="full"))


# ---------------------------------------------------------------------------
# グラフ状態
# ---------------------------------------------------------------------------

def build_state(
    a: Sequence,
    b: Sequence,
    *,
    f_expr=None,
    order: int | None = 2,
    iterate: str | None = None,
) -> dict:
    """Desmos の state (JSON 化する前の dict) を作る。

    a: f のテイラー係数 / b: 結果 g のテイラー係数
    f_expr: f の式 (あれば f(x) はそのまま描き、なければテイラー多項式で描く)
    order: 関数的 m 乗根の m (g を m 回合成した検算用の曲線を描く)
    iterate: 分数回反復のとき 'p/q' (検算用の曲線は描かない)
    """
    f_latex = expr_latex(f_expr) if f_expr is not None else series_latex(a)
    g_name = "h" if iterate is not None else "g"
    items = []

    items.append(
        {"type": "expression", "id": "f", "color": COLORS["f"], "latex": f"f(x)={f_latex}"}
    )
    items.append(
        {
            "type": "expression",
            "id": "g",
            "color": COLORS["g"],
            "latex": f"{g_name}(x)={series_latex(b)}",
        }
    )
    if iterate is None and order and order >= 2:
        nested = "x"
        for _ in range(order):
            nested = f"g({nested})"
        items.append(
            {
                "type": "expression",
                "id": "check",
                "color": COLORS["check"],
                "lineStyle": "DASHED",
                "latex": nested,
            }
        )
        # 検算の誤差 |g^[m](x) - f(x)| (0 に近いほど正確)
        items.append(
            {
                "type": "expression",
                "id": "error",
                "color": COLORS["error"],
                "latex": rf"\left|{nested}-f(x)\right|",
            }
        )

    return {
        "version": 11,
        "graph": {"viewport": {"xmin": -2.5, "ymin": -2.5, "xmax": 2.5, "ymax": 2.5}},
        "expressions": {"list": items},
    }


# ---------------------------------------------------------------------------
# 保存
# ---------------------------------------------------------------------------

def _blank_png(width: int = 1, height: int = 1) -> bytes:
    """白一色の PNG (保存 API がサムネイル画像を必須とするため)。"""
    raw = b"".join(b"\x00" + b"\xff" * (width * 3) for _ in range(height))

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    header = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


def _new_hash() -> str:
    alphabet = string.ascii_lowercase + string.digits
    return "".join(secrets.choice(alphabet) for _ in range(10))


def create_link(state: dict, *, timeout: float = 15.0) -> str:
    """state を Desmos に保存し、共有リンクの URL を返す。"""
    graph_hash = _new_hash()
    thumb = "data:image/png;base64," + base64.b64encode(_blank_png()).decode()
    body = {
        "thumbnailData": thumb,
        "hash": graph_hash,
        "type": "quick-link",
        "state": json.dumps(state, ensure_ascii=False),
        "product": "graphing",
        "lang": "ja",
    }
    req = urllib.request.Request(
        SAVE_URL,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "func_sqrt",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as res:
            data = json.loads(res.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")[:200]
        raise DesmosError(f"Desmos への保存に失敗しました (HTTP {exc.code}): {detail}") from exc
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        raise DesmosError(f"Desmos への保存に失敗しました: {exc}") from exc
    return CALC_URL + data.get("hash", graph_hash)
