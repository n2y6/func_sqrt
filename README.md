# func_sqrt

関数的平方根 `g(g(x)) = f(x)` のテイラー展開を、分数と根号のまま厳密に求める

[naisu.moe/func_sqrt](https://naisu.moe/func_sqrt/)

```text
$ python -m func_sqrt "sinh(x)" -n 11 -f series
x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 - (23/71680)x^9 + (92713/1277337600)x^11 + O(x^12)
```

| | 式 | オプション |
| --- | --- | --- |
| 関数的平方根 | `g∘g = f` | |
| 関数的 m 乗根 | `g^[m] = f` | `--order m` |
| 分数回反復 | `f^[p/q]` | `--iterate p/q` |

`f(0) = 0` かつ `f'(0) ≠ 0` の関数に対応

## インストール

Python 3.9+ と [sympy](https://www.sympy.org/)

```bash
pip install -e .
```

## コマンドライン

```bash
python -m func_sqrt "sinh(x)" -n 11            # 式で
python -m func_sqrt --derivs 0,1,0,1,0,1        # 微分係数 f^(n)(0) で
python -m func_sqrt --coeffs 0,1,0,1/6          # テイラー係数 a_n で
python -m func_sqrt                             # 対話モード
```

```text
f(x) ~ x + (1/6)x^3 + (1/120)x^5 + (1/5040)x^7 + (1/362880)x^9 + (1/39916800)x^11 + O(x^12)

--- 結果: g(x)  (g∘g = f) ---
 n               b_n  g^(n)(0)
--  ----------------  --------
 1                 1         1
 2                 0         0
 3              1/12       1/2
 ...
11  92713/1277337600  92713/32

g(x) ~ x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 - (23/71680)x^9 + (92713/1277337600)x^11 + O(x^12)
```

```bash
python -m func_sqrt "exp(x)-1" -n 6 --order 3        # 3 乗根
python -m func_sqrt "x + x**2" -n 6 --iterate 3/2    # f^[3/2]
python -m func_sqrt "x + x**2" -n 6 --iterate -1     # 逆関数
python -m func_sqrt "4*x + x**2" -n 5 --multiplier -2  # g'(0) = -2 の枝
```

| オプション | |
| --- | --- |
| `-f series` / `latex` / `table` / `json` | 出力形式 |
| `--digits K` | 数値近似を併記 |
| `--check` | 合成して f に戻るか検算 |
| `--desmos` | Desmos のグラフを作る |

式が `-` で始まるときは `--` の後に置く（`python -m func_sqrt -- "-8*x + x**2"`）

## Python

```python
from func_sqrt import from_expr, functional_sqrt, functional_root, fractional_iterate, format_series

a = from_expr("sinh(x)", 11)
print(format_series(functional_sqrt(a)))   # x + (1/12)x^3 - ...

functional_root(a, 3)
fractional_iterate(a, "3/2")
```

級数の演算（`compose`, `iterate`, `inverse` など）は `func_sqrt.series` に

## Web API

```bash
curl "https://api.naisu.moe/func_sqrt?f=sinh(x)&n=11"
```

```json
{
  "ok": true,
  "series": "x + (1/12)x^3 - (1/160)x^5 + ...",
  "latex": "x + \\frac{1}{12} x^{3} - ...",
  "coefficients": ["0", "1", "0", "1/12", ...],
  "derivatives": ["0", "1", "0", "1/2", ...]
}
```

| パラメータ | |
| --- | --- |
| `f` | 式（`derivs` / `coeffs` でも可） |
| `n` | 次数（最大 20） |
| `order` | m 乗根の m（最大 6） |
| `iterate` | `p/q` |
| `multiplier` | g'(0) |
| `desmos` | `1` で Desmos のリンクを付ける |

`https://api.naisu.moe/` で使い方を返す

Desmos のリンク作成は非公式の API を使っているので、そのうち動かなくなるかもしれない

## 仕組み

`g(x) = b_1 x + b_2 x^2 + …` とおいて、`g^[m]` と `f` の係数を低次から比べる

`b_1` は `f'(0)` の m 乗根で、`n ≥ 2` では `x^n` の係数が `b_n` の一次式

```
λ_n b_n + (b_1, …, b_{n-1} で決まる項),   λ_n = Σ_{j=0}^{m-1} b_1^j b_1^{n(m-1-j)}
```

になるので、`b_n` が一つに決まる

`λ_n = 0`（`f'(0) = 1` で `g'(0) = -1` など）のときは決まらない

## 開発

```bash
pip install -e ".[dev]"
pytest
```

## ライセンス

MIT
