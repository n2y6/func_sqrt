import json

import pytest
import sympy as sp
from sympy import Rational as R

from func_sqrt import format_series, parse_number
from func_sqrt.cli import build_parser, main, run

README_SERIES = (
    "x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 - (23/71680)x^9"
    " + (92713/1277337600)x^11 + O(x^12)"
)


def test_format_series():
    assert format_series([0, 1, 0, R(1, 12), -1, 3]) == "x + (1/12)x^3 - x^4 + 3x^5 + O(x^6)"
    assert format_series([0, -2], big_o=False) == "-2x"
    assert format_series([0, sp.sqrt(2)], big_o=False) == "sqrt(2)x"
    assert format_series([0, 0], big_o=False) == "0"
    assert format_series([0, 1, R(-1, 2)], style="latex") == r"x - \frac{1}{2} x^{2} + O(x^{3})"


def test_parse_number():
    assert parse_number(" -1/2 ") == R(-1, 2)
    assert parse_number("0.25") == R(1, 4)
    assert parse_number("sqrt(2)") == sp.sqrt(2)
    for bad in ["", "abc", "1/"]:
        with pytest.raises(ValueError):
            parse_number(bad)


def test_cli_expression(capsys):
    assert main(["sinh(x)", "-n", "11", "-f", "series"]) == 0
    assert capsys.readouterr().out.strip() == README_SERIES


def test_cli_derivs_and_check(capsys):
    assert main(["--derivs", ",".join(["0", "1"] * 6), "--check"]) == 0
    out = capsys.readouterr().out
    assert README_SERIES in out
    assert "検証: OK" in out


def test_cli_json(capsys):
    assert main(["--coeffs", "0,1,0,1/6", "-f", "json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["coefficients"] == ["0", "1", "0", "1/12"]
    assert data["derivatives"] == ["0", "1", "0", "1/2"]
    assert data["order"] == 2


def test_cli_iterate(capsys):
    assert main(["x + x**2", "-n", "4", "-t", "-1", "-f", "series"]) == 0
    assert capsys.readouterr().out.strip() == "x - x^2 + 2x^3 - 5x^4 + O(x^5)"


def test_cli_errors(capsys):
    assert main(["cos(x)", "-n", "3"]) == 2
    assert "f(0) が 0 ではありません" in capsys.readouterr().err
    assert main(["sinh(x)", "--derivs", "0,1"]) == 2
    assert main(["--derivs", "0,1,,2"]) == 2


def test_interactive_derivatives(capsys):
    answers = iter(["", "x", "3", "0", "1", "0", "1/2"])  # 'x' は不正な N → 再入力
    args = build_parser().parse_args([])
    assert run(args, input_fn=lambda _: next(answers)) == 0
    assert "x + (1/24)x^3 + O(x^4)" in capsys.readouterr().out


def test_interactive_expression(capsys):
    answers = iter(["sinh(x)", "5"])
    args = build_parser().parse_args(["-i", "-f", "series"])
    assert run(args, input_fn=lambda _: next(answers)) == 0
    out = capsys.readouterr().out.strip()
    assert out.splitlines()[-1] == "x + (1/12)x^3 - (1/160)x^5 + O(x^6)"

