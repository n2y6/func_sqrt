"""関数的平方根・関数的 m 乗根・分数回反復の計算。

f(0) = 0 (原点が不動点) かつ f'(0) ≠ 0 の関数 f について、

* 関数的平方根      g∘g = f
* 関数的 m 乗根     g^[m] = g∘g∘...∘g (m 回) = f
* 分数回反復        f^[p/q] = (f の関数的 q 乗根)^[p]

の x=0 におけるテイラー係数を、係数比較によって低次から順に厳密に求める。

アルゴリズム
------------
g(x) = b_1 x + b_2 x^2 + ... とおく。b_1^m = a_1 となるように b_1 を選び、
n = 2, 3, ... の順に b_n を決める。g^[m] の x^n の係数は b_n について一次式で、

    [x^n] g^[m] = λ_n · b_n + (b_1, ..., b_{n-1} だけで決まる項)
    λ_n = Σ_{j=0}^{m-1} b_1^j · b_1^{n(m-1-j)}

となる (m=2 なら λ_n = b_1 + b_1^n)。そこで b_n = 0 として g^[m] を計算し、
その x^n の係数 h_n を使って b_n = (a_n - h_n) / λ_n と求める。
λ_n = 0 になる場合は解が一意に決まらない (または存在しない) のでエラーにする。
"""

from __future__ import annotations

from typing import Sequence

import sympy as sp

from .series import Series, clean, compose, iterate, to_series, truncate

__all__ = [
    "FunctionalRootError",
    "principal_root",
    "functional_root",
    "functional_sqrt",
    "fractional_iterate",
    "verify",
]


class FunctionalRootError(ValueError):
    """関数的累乗根が (この方法で) 一意に定まらないときに送出される。"""


def principal_root(value, m: int) -> sp.Expr:
    """value の m 乗根のうち「自然な」ものを返す。

    * value > 0 なら正の実数根
    * value < 0 かつ m が奇数なら負の実数根
    * それ以外は sympy の主値 (複素数になることがある)
    """
    value = sp.sympify(value)
    if m == 1:
        return value
    if value.is_extended_negative and m % 2 == 1:
        return clean(-sp.root(-value, m))
    return clean(sp.root(value, m))


def _linear_factor(b1: sp.Expr, n: int, m: int) -> sp.Expr:
    """λ_n = Σ_{j=0}^{m-1} b1^j · b1^{n(m-1-j)}  (b_n の係数)。"""
    return clean(sum(b1**j * b1 ** (n * (m - 1 - j)) for j in range(m)))


def _is_zero(value: sp.Expr) -> bool:
    value = sp.sympify(value)
    if value.is_Rational:
        return value == 0
    return sp.simplify(value) == 0


def functional_root(a: Sequence, m: int = 2, *, multiplier=None) -> Series:
    """f の関数的 m 乗根 g (g^[m] = f) のテイラー係数を返す。

    Parameters
    ----------
    a:
        f のテイラー係数 ``[a_0, a_1, ..., a_N]`` (a_0 = 0, a_1 ≠ 0)。
    m:
        何乗根か (2 なら関数的平方根)。
    multiplier:
        g'(0) = b_1 として使う値。省略時は :func:`principal_root` で a_1 の
        m 乗根を選ぶ。``multiplier**m == a_1`` を満たす必要がある。
        例えば f(x) = 4x + ... に対し ``multiplier=-2`` とすると
        g'(0) = -2 の枝の解が得られる。

    Returns
    -------
    g のテイラー係数 ``[0, b_1, ..., b_N]``。
    """
    if not isinstance(m, int) or m < 1:
        raise ValueError(f"m は 1 以上の整数で指定してください (m={m!r})")
    a = to_series(a)
    N = len(a) - 1
    if N < 0:
        raise ValueError("係数が 1 つもありません")
    if not _is_zero(a[0]):
        raise FunctionalRootError(
            "f(0) が 0 ではありません。x=0 を不動点とする関数のみ対応しています。"
        )
    b: Series = [sp.Integer(0)] * (N + 1)
    if N == 0:
        return b
    if _is_zero(a[1]):
        raise FunctionalRootError("f'(0) = 0 の場合、この方法では適切な展開ができません。")

    if multiplier is None:
        b1 = principal_root(a[1], m)
    else:
        b1 = clean(sp.sympify(multiplier))
        if not _is_zero(b1**m - a[1]):
            raise ValueError(
                f"multiplier={b1} は条件 multiplier^{m} = f'(0) = {a[1]} を満たしません"
            )
    b[1] = b1

    if all(sp.sympify(c).is_Rational for c in a):
        # 有理数係数なら代数体 Q(b_1) の要素のまま計算する (速く、式も膨らまない)
        algebraic = _root_in_number_field(a, m, b1)
        if algebraic is not None:
            return algebraic

    for n in range(2, N + 1):
        lam = _linear_factor(b1, n, m)
        if _is_zero(lam):
            raise FunctionalRootError(
                f"n={n} で b_{n} の係数 λ_{n} が 0 になり、解が一意に決まりません "
                f"(g'(0) = {b1})。"
            )
        # b_n = 0 のまま g^[m] を n 次まで計算し、x^n の係数を見る
        h = iterate(b[: n + 1], m, n)
        b[n] = clean((a[n] - h[n]) / lam)
    return b


def _root_in_number_field(a: Series, m: int, b1: sp.Expr) -> Series | None:
    """有理数係数の f について、Q(b_1) (b_1 が有理数なら Q) の中で高速に厳密計算する。

    sympy の汎用の式 (Expr) ではなく、代数体の要素 (QQ<b_1>) のまま演算するので、
    2^(1/3) のような無理数が出ても式が膨らまない。b_1 が代数的数でない場合は None。
    """
    try:
        K = sp.QQ if b1.is_Rational else sp.QQ.algebraic_field(b1)
    except (NotImplementedError, ValueError, sp.polys.polyerrors.PolynomialError):
        return None
    N = len(a) - 1
    A = [K.convert(sp.Rational(c)) for c in a]
    zero, one = K.zero, K.one
    t = K.from_sympy(b1)

    def mul(P, Q, L):
        R = [zero] * L
        for i, p in enumerate(P[:L]):
            if not p:
                continue
            for j, q in enumerate(Q[: L - i]):
                if q:
                    R[i + j] += p * q
        return R

    def compose(outer, inner, L):
        R = [zero] * L
        for k in range(min(len(outer), L) - 1, -1, -1):
            R = mul(R, inner, L)
            R[0] += outer[k]
        return R

    b = [zero] * (N + 1)
    b[1] = t
    for n in range(2, N + 1):
        lam = sum((t**j * t ** (n * (m - 1 - j)) for j in range(m)), zero)
        if not lam:
            raise FunctionalRootError(
                f"n={n} で b_{n} の係数 λ_{n} が 0 になり、解が一意に決まりません "
                f"(g'(0) = {b1})。"
            )
        # b_n = 0 のまま g^[m] を n 次まで計算し、x^n の係数を見る
        h = [zero, one] + [zero] * (n - 1)
        for _ in range(m):
            h = compose(b[: n + 1], h, n + 1)
        b[n] = (A[n] - h[n]) / lam
    return [sp.expand(K.to_sympy(c)) for c in b]


def functional_sqrt(a: Sequence, *, multiplier=None) -> Series:
    """関数的平方根 g∘g = f のテイラー係数 (``functional_root(a, 2)`` と同じ)。"""
    return functional_root(a, 2, multiplier=multiplier)


def fractional_iterate(a: Sequence, t, *, multiplier=None) -> Series:
    """分数回反復 f^[t] (t = p/q) のテイラー係数を返す。

    f の関数的 q 乗根 g を求め、それを p 回合成する (p < 0 なら逆関数を使う)。
    t = 1/2 なら関数的平方根、t = -1 なら逆関数、t = 3/2 なら f∘(√f) になる。

    ``multiplier`` は q 乗根 g の g'(0) を指定する (:func:`functional_root` 参照)。
    """
    t = sp.Rational(sp.nsimplify(t, rational=True))
    p, q = int(t.p), int(t.q)
    a = to_series(a)
    N = len(a) - 1
    root = a if q == 1 else functional_root(a, q, multiplier=multiplier)
    return iterate(root, p, N)


def verify(a: Sequence, b: Sequence, m: int = 2) -> Series:
    """g^[m] - f の各係数 (誤差) を返す。すべて 0 なら g は正しい解。"""
    a = to_series(a)
    N = len(a) - 1
    composed = iterate(truncate(b, N), m, N)
    return [clean(sp.simplify(c - ai)) for c, ai in zip(composed, a)]
