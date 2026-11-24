import pytest
import sympy as sp

from func_sqrt import from_expr, functional_root, functional_sqrt
from func_sqrt.desmos import DesmosError, _blank_png, build_state, expr_latex, series_latex
from func_sqrt.web import handle


def test_latex_for_desmos():
    x = sp.Symbol("x")
    assert expr_latex(sp.sinh(x)) == r"\sinh\left(x \right)"
    assert expr_latex(sp.log(1 + x)) == r"\ln\left(x + 1 \right)"
    assert "arcsinh" in expr_latex(sp.asinh(x))
    assert series_latex([0, 1, 0, sp.Rational(1, 12)]) == r"x + \frac{1}{12} x^{3}"
    with pytest.raises(DesmosError):
        series_latex([0, sp.I])


def test_build_state_square_root():
    a = from_expr("sinh(x)", 5)
    state = build_state(a, functional_sqrt(a), f_expr=sp.sinh(sp.Symbol("x")))
    latex = [e.get("latex") for e in state["expressions"]["list"] if e["type"] == "expression"]
    assert latex[0] == r"f(x)=\sinh\left(x \right)"
    assert latex[1].startswith(r"g(x)=x + \frac{1}{12} x^{3}")
    assert latex[2] == "g(g(x))"


def test_build_state_cube_root_and_iterate():
    a = from_expr("exp(x)-1", 4)
    state = build_state(a, functional_root(a, 3), order=3)
    assert any(e.get("latex") == "g(g(g(x)))" for e in state["expressions"]["list"])
    state = build_state(a, a, order=None, iterate="3/2")
    latex = [e.get("latex", "") for e in state["expressions"]["list"]]
    assert any(s.startswith("h(x)=") for s in latex)
    assert not any(s.startswith("g(") for s in latex)


def test_blank_png_signature():
    assert _blank_png().startswith(b"\x89PNG\r\n\x1a\n")


def test_web_desmos_param():
    seen = {}

    def fake_link(state):
        seen["state"] = state
        return "https://www.desmos.com/calculator/test"

    status, body = handle({"f": "sinh(x)", "n": "5", "desmos": "1"}, desmos_link=fake_link)
    assert status == 200 and body["desmos"] == "https://www.desmos.com/calculator/test"
    assert seen["state"]["expressions"]["list"][1]["latex"].startswith("f(x)=")

    status, body = handle({"f": "sinh(x)", "n": "5"}, desmos_link=fake_link)
    assert "desmos" not in body


def test_web_desmos_failure_is_not_fatal():
    def broken(state):
        raise DesmosError("down")

    status, body = handle({"f": "sinh(x)", "n": "5", "desmos": "1"}, desmos_link=broken)
    assert status == 200 and body["ok"] and body["desmos"] is None and body["desmos_error"] == "down"
