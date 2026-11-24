# func_sqrt — 関数的平方根の厳密計算

`g(g(x)) = f(x)` を満たす関数 **g（f の関数的平方根）** のテイラー展開を、
`x = 0` のまわりで **厳密に**（有理数・根号のまま）求めるツールです。

さらに次のような一般化にも対応しています。

| 求めるもの | 条件 | オプション |
| --- | --- | --- |
| 関数的平方根 | `g∘g = f` | （既定） |
| 関数的 m 乗根 | `g∘g∘…∘g`（m 回）`= f` | `--order m` |
| 分数回反復 | `f^[p/q]`（例: `f^[3/2]`, `f^[-1]` = 逆関数） | `--iterate p/q` |

```text
$ python -m func_sqrt "sinh(x)" -n 11 -f series
x + (1/12)x^3 - (1/160)x^5 + (53/40320)x^7 - (23/71680)x^9 + (92713/1277337600)x^11 + O(x^12)
```

## 前提条件

- `f(0) = 0`（`x = 0` が不動点）
- `f'(0) ≠ 0`

`f'(0) = 0` の場合や、途中で解が一意に決まらない場合はエラーになります（[仕組み](#仕組み) 参照）。

## インストール

Python 3.9 以上と [sympy](https://www.sympy.org/) が必要です。

```bash
pip install sympy
```

リポジトリ直下でそのまま `python -m func_sqrt` として使えます。
`func-sqrt` コマンドとしてインストールする場合:

```bash
pip install -e .
```

## 使い方（コマンドライン）

### 1. 式で入力する

```bash
python -m func_sqrt "sinh(x)" -n 11
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

`b_n` はテイラー係数、`g^(n)(0) = n!·b_n` は n 階微分係数です。
式は sympy の記法で書きます（`exp(x)-1`, `log(1+x)`, `2*x + x**2` など）。
式が `-` で始まる場合は `--` の後に置いてください: `python -m func_sqrt -n 6 -m 3 -- "-8*x + x**2"`

### 2. 微分係数で入力する

`f(0), f'(0), f''(0), …, f^(N)(0)` をカンマ区切りで与えます（N は個数から決まります）。

```bash
python -m func_sqrt --derivs 0,1,0,1,0,1,0,1,0,1,0,1
```

テイラー係数 `a_n = f^(n)(0)/n!` で与える場合は `--coeffs` を使います。
値には整数・小数・分数（`-1/2`）・根号（`sqrt(2)`）が使えます。

### 3. 対話モード

引数なしで起動すると、旧 `sqrt.py` と同じように 1 つずつ入力できます。

```bash
python -m func_sqrt
```

### 一般化

```bash
# 関数的 3 乗根: g(g(g(x))) = e^x - 1
python -m func_sqrt "exp(x)-1" -n 6 --order 3

# 3/2 回反復: f^[3/2] = f∘√f
python -m func_sqrt "x + x**2" -n 6 --iterate 3/2

# 逆関数
python -m func_sqrt "x + x**2" -n 6 --iterate -1

# g'(0) の枝を指定する (4x + … に対して g'(0) = -2 の解)
python -m func_sqrt "4*x + x**2" -n 5 --multiplier -2
```

### 出力オプション

| オプション | 内容 |
| --- | --- |
| `-f text` | 表 + 近似式（既定） |
| `-f series` | 近似式のみ |
| `-f latex` | LaTeX 形式の近似式 |
| `-f table` | 係数表のみ |
| `-f json` | 係数・微分係数・式を JSON で |
| `--digits K` | 表に K 桁の数値近似を併記 |
| `--check` | 結果を合成して f に戻るか検証 |
| `--desmos` | f・g・g(g(x)) を描いた Desmos グラフの共有リンクを作る |

すべてのオプションは `python -m func_sqrt --help` で確認できます。

## 使い方（Python ライブラリ）

```python
from func_sqrt import (
    from_expr, from_derivatives, functional_sqrt, functional_root,
    fractional_iterate, verify, format_series,
)

a = from_expr("sinh(x)", 11)          # f のテイラー係数 [a_0, ..., a_11]
b = functional_sqrt(a)                 # g のテイラー係数
print(format_series(b))                # x + (1/12)x^3 - ...
print(format_series(b, style="latex")) # LaTeX
assert all(c == 0 for c in verify(a, b))

functional_root(a, 3)                  # 関数的 3 乗根
fractional_iterate(a, "3/2")           # f^[3/2]
from_derivatives([0, 1, 0, 1])         # 微分係数から係数へ
```

級数そのものの演算（`compose`, `iterate`, `inverse`, `to_derivatives` など）は
`func_sqrt.series` にあります。

## Web API

`https://api.naisu.moe/func_sqrt` で同じ計算を JSON として返します（CORS 許可済み）。

```bash
curl "https://api.naisu.moe/func_sqrt?f=sinh(x)&n=11"
```

```json
{
  "ok": true,
  "input": { "f": "sinh(x)" },
  "order": 2,
  "series": "x + (1/12)x^3 - (1/160)x^5 + ...",
  "latex": "x + \\frac{1}{12} x^{3} - ...",
  "coefficients": ["0", "1", "0", "1/12", ...],
  "derivatives": ["0", "1", "0", "1/2", ...]
}
```

| パラメータ | 内容 |
| --- | --- |
| `f` | f(x) の式（`sin cos tan sinh cosh tanh asin acos atan asinh acosh atanh exp log sqrt pi E I` と `x` が使えます。`^` も累乗） |
| `derivs` / `coeffs` | `f` の代わりに微分係数・テイラー係数をカンマ区切りで |
| `n` | 展開次数（既定 10、最大 20） |
| `order` | 関数的 m 乗根の m（既定 2、最大 6） |
| `iterate` | 分数回反復 `p/q`（`order` の代わり） |
| `multiplier` | g'(0) に使う値 |
| `desmos` | `1` にすると Desmos グラフのリンクを `desmos` に入れて返す |

パラメータなしで `https://api.naisu.moe/` を開くと使い方を JSON で返します。
エラー時は HTTP 400 と `{"ok": false, "error": "..."}`、計算が 10 秒を超えると 504 を返します。
外部入力は許可した関数名・記号だけを受け付けるので、任意コードは実行されません。

### Desmos グラフ

`--desmos`（CLI）や `desmos=1`（API）を付けると、f(x)・g(x)・検算用の g(g(x)) を描いた
Desmos グラフを保存して共有リンクを返します。

```bash
curl "https://api.naisu.moe/func_sqrt?f=sinh(x)&n=9&desmos=1"
# → "desmos": "https://www.desmos.com/calculator/xxxxxxxxxx"
```

保存には Desmos のサイト内部で使われている非公開の API を使っています。
公式に提供されたものではないため、仕様変更で動かなくなる可能性があります
（その場合も計算結果は返り、`desmos_error` に理由が入ります）。

## 仕組み

`g(x) = b_1 x + b_2 x^2 + …` とおき、`g^[m]`（m 回合成）の係数を `f` と比較して
低次から順に決めます。

1. `b_1^m = f'(0)` から `b_1` を決める（既定は主値。`--multiplier` で変更可）。
2. `n = 2, 3, …` について、`g^[m]` の `x^n` の係数は `b_n` の一次式
   `λ_n·b_n + (b_1 … b_{n-1} のみで決まる項)` になる。ここで

   ```
   λ_n = Σ_{j=0}^{m-1} b_1^j · b_1^{n(m-1-j)}     (m = 2 なら λ_n = b_1 + b_1^n)
   ```

3. `b_n = 0` として `g^[m]` を計算し、その `x^n` の係数 `h_n` から
   `b_n = (a_n - h_n) / λ_n` と求める。

`λ_n = 0` となる場合（例: `f'(0) = 1` で `g'(0) = -1` を選んだとき）は解が一意に
定まらないためエラーにします。計算はすべて sympy の厳密演算なので、
本来 0 になる係数が浮動小数点誤差で汚れることはありません。

## 開発

```bash
pip install -e ".[dev]"
pytest
```

```
func_sqrt/
├── series.py      # 形式的冪級数の演算 (積・合成・反復・逆関数・テイラー展開)
├── core.py        # 関数的平方根 / m 乗根 / 分数回反復
├── formatting.py  # 式・表・JSON への整形
├── cli.py         # コマンドラインインターフェース
├── desmos.py      # Desmos グラフの共有リンク作成
└── web.py         # Web API のリクエスト処理
deploy/cgi/        # api.naisu.moe 用 CGI と .htaccess
tests/             # pytest
```

## ライセンス

MIT
