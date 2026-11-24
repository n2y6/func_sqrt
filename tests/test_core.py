import pytest
import sympy as sp
from sympy import Rational as R

from func_sqrt import (
    FunctionalRootError,
    compose,
    fractional_iterate,
    from_derivatives,
    from_expr,
    functional_root,
    functional_sqrt,
    inverse,
    iterate,
    verify,
)

SINH_SQRT = [0, 1, 0, R(1, 12), 0, R(-1, 160), 0, R(53, 40320), 0,
             R(-23, 71680), 0, R(92713, 1277337600)]


def test_sinh_matches_readme_example():
    assert functional_sqrt(from_expr("sinh(x)", 11)) == SINH_SQRT


def test_derivative_input_equals_expression_input():
    derivs = [0, 1] * 6
    assert from_derivatives(derivs) == from_expr("sinh(x)", 11)


def test_exp_minus_one_known_half_iterate():
    # 半反復 of e^x - 1 : x + x^2/4 + x^3/48 + 0 x^4 + x^5/3840 - 7x^6/92160
    b = functional_sqrt(from_expr("exp(x) - 1", 6))
    assert b == [0, 1, R(1, 4), R(1, 48), 0, R(1, 3840), R(-7, 92160)]


@pytest.mark.parametrize("expr", ["sinh(x)", "sin(x)", "exp(x)-1", "x + x**2", "4*x + x**3", "log(1+x)"])
@pytest.mark.parametrize("m", [2, 3])
def test_roundtrip(expr, m):
    a = from_expr(expr, 7)
    b = functional_root(a, m)
    assert all(c == 0 for c in verify(a, b, m))


def test_irrational_multiplier_roundtrip():
    a = from_expr("2*x + x**2", 5)
    b = functional_sqrt(a)
    assert b[1] == sp.sqrt(2)
    assert all(c == 0 for c in verify(a, b))


def test_negative_real_cube_root():
    a = from_expr("-8*x + x**2", 4)
    b = functional_root(a, 3)
    assert b[1] == -2
    assert all(c == 0 for c in verify(a, b, 3))


def test_custom_multiplier_branch():
    a = from_expr("4*x + x**2", 5)
    b = functional_sqrt(a, multiplier=-2)
    assert b[1] == -2
    assert all(c == 0 for c in verify(a, b))


def test_order_one_is_identity_operation():
    a = from_expr("sin(x)", 6)
    assert functional_root(a, 1) == a


def test_errors():
    with pytest.raises(FunctionalRootError):
        functional_sqrt(from_expr("cos(x)", 4))
    with pytest.raises(FunctionalRootError):
        functional_sqrt(from_expr("x**2", 4))
    with pytest.raises(ValueError):
        functional_sqrt(from_expr("4*x", 4), multiplier=3)
    with pytest.raises(FunctionalRootError, match="λ_2"):
        functional_sqrt(from_expr("x + x**2", 4), multiplier=-1)
    with pytest.raises(ValueError):
        functional_root([0, 1], 0)


def test_inverse_and_iterate():
    a = from_expr("x + x**2", 6)
    inv = inverse(a)
    assert inv == [0, 1, -1, 2, -5, 14, -42]  # Catalan 数 (符号付き)
    assert compose(a, inv) == [0, 1, 0, 0, 0, 0, 0]
    assert iterate(a, -1) == inv
    assert iterate(a, 0) == [0, 1, 0, 0, 0, 0, 0]
    assert iterate(a, 2) == compose(a, a)


def test_fractional_iterate():
    a = from_expr("x + x**2", 6)
    half = functional_sqrt(a)
    assert fractional_iterate(a, R(1, 2)) == half
    assert fractional_iterate(a, "1/2") == half
    assert fractional_iterate(a, 1) == a
    assert fractional_iterate(a, -1) == inverse(a)
    assert fractional_iterate(a, R(3, 2)) == compose(a, half)
    third = fractional_iterate(a, R(1, 3))
    assert fractional_iterate(a, R(2, 3)) == compose(third, third)
