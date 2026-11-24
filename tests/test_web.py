import pytest
import sympy as sp

from func_sqrt.web import ApiError, handle, safe_parse


def test_expression():
    status, body = handle({"f": "sinh(x)", "n": "7"})
    assert status == 200 and body["ok"]
    assert body["series"] == "x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 + O(x^8)"
    assert body["coefficients"][3] == "1/12"
    assert body["derivatives"][3] == "1/2"


def test_caret_power_and_constants():
    assert safe_parse("2*x + x^2") == 2 * sp.Symbol("x") + sp.Symbol("x") ** 2
    assert safe_parse("sqrt(2)*sin(pi*x)").has(sp.pi)


def test_derivs_order_iterate():
    status, body = handle({"derivs": "0,1,0,1,0,1"})
    assert status == 200 and body["series"].startswith("x + (1/12)x^3")
    status, body = handle({"f": "x + x^2", "n": "4", "iterate": "-1"})
    assert body["series"] == "x - x^2 + 2x^3 - 5x^4 + O(x^5)"
    status, body = handle({"f": "exp(x)-1", "n": "4", "order": "3"})
    assert status == 200 and body["order"] == 3


@pytest.mark.parametrize(
    "evil",
    [
        "__import__('os').system('id')",
        "x.__class__",
        "Symbol('y')",
        "open('/etc/passwd')",
        "y + x",
        "lambda: 1",
        "x; 1",
        "eval(1)",
    ],
)
def test_rejects_unsafe_input(evil):
    with pytest.raises(ApiError):
        safe_parse(evil)
    status, body = handle({"f": evil})
    assert status == 400 and not body["ok"]


@pytest.mark.parametrize(
    "params",
    [
        {},
        {"f": "sin(x)", "derivs": "0,1"},
        {"f": "sin(x)", "n": "999"},
        {"f": "sin(x)", "order": "0"},
        {"f": "cos(x)"},
        {"f": "sin(x)", "iterate": "1/100"},
        {"f": "sin(x)", "order": "2", "iterate": "1/2"},
        {"derivs": "0,x"},
        {"f": "x" * 400},
    ],
)
def test_bad_requests(params):
    status, body = handle(params)
    assert status == 400 and body["ok"] is False and body["error"]
