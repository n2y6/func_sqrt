"""func_sqrt — 関数的平方根 (g∘g = f) のテイラー展開を厳密計算するライブラリ。

基本的な使い方::

    >>> from func_sqrt import from_expr, functional_sqrt, format_series
    >>> a = from_expr("sinh(x)", 7)
    >>> format_series(functional_sqrt(a))
    'x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 + O(x^8)'
"""

__version__ = "0.2.0"

from .core import (
    FunctionalRootError,
    fractional_iterate,
    functional_root,
    functional_sqrt,
    principal_root,
    verify,
)
from .formatting import format_json, format_series, format_table
from .series import (
    compose,
    from_derivatives,
    from_expr,
    inverse,
    iterate,
    parse_number,
    to_derivatives,
)

__all__ = [
    "__version__",
    "FunctionalRootError",
    "functional_sqrt",
    "functional_root",
    "fractional_iterate",
    "principal_root",
    "verify",
    "compose",
    "iterate",
    "inverse",
    "from_expr",
    "from_derivatives",
    "to_derivatives",
    "parse_number",
    "format_series",
    "format_table",
    "format_json",
]
