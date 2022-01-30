---
title: "不完全Cholesky分解前処理付き共役勾配法"
emoji: "🌊"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["数学", "cpp"]
published: true
---

# 線型連立方程式を解く

線形な連立方程式は，行列 $A \in \mathbb{R}^{N\times N}$ とベクトル $\boldsymbol{x}, \boldsymbol{b} \in \mathbb{R}^N$ を用いて，

$$
A\boldsymbol{x} = \boldsymbol{b}
$$

と表されます．

これを数値的に解くことは，次元 $N$ が大きくなると計算量の大きな処理になります．

例えば，ガウスの消去法では， $\mathcal{O}(N^3)$ 程度の計算量になります．

線形連立方程式 $A\boldsymbol{x} = \boldsymbol{b}$ を数値的に解く作業は，行列 $A$ の具体的な形に依存して，その難易度も変わります．

例えば，極端な例を挙げると， $A$ が単位行列 $I$ の場合は，何もしなくても良いです．

## 前処理

したがって， 連立方程式 $A\boldsymbol{x} = \boldsymbol{b}$ を解く前に， この方程式に左から，行列 $P_L$ をかけ，また， $\boldsymbol{x} = P_R \boldsymbol{y}$ として，

$$
P_LAP_R \boldsymbol{y} = P_L\boldsymbol{b}
$$

としておいて， 行列 $P_L A P_R$ を数値的に解きやすい行列にしておくと，計算効率が劇的によくなることがあります．

このような処理のことを，前処理と呼びます．

一方で，前処理も計算を必要とする処理である以上，前処理の計算量は小さいものでないと意味がありません．

そこで，前処理の目指すところは，ざっくり言うと以下の2点と言うことになります．

- 前処理後の線形連立方程式 $P_LAP_R \boldsymbol{y} = P_L\boldsymbol{b}$ は解きやすい
- 前処理の計算量は小さい


# 前処理付き共役勾配法

この記事では，前処理一般について議論するのではなく，前処理を施した共役勾配法について述べます．

すなわち，

1. 前処理を施した連立線型方程式， $P_L A P_R \boldsymbol{y} = P_L\boldsymbol{b}$ を解くアルゴリズムの導出．
2. 前処理行列 $P_L, P_R$ の具体形の紹介.

について述べます．


## 解くべき方程式

今回は共役勾配法に対する前処理を考えます．

そのため $A$ は正定値対称行列として，前処理後の形も，正定値対称行列であるものを探します．すなわち，

$$
\begin{align*}
S^\mathrm{T} A S \boldsymbol{y} &= S^\mathrm{T} \boldsymbol{b} \\
x &= S\boldsymbol{y}
\end{align*}
$$

となる，行列 $S$ を求め，共役勾配法を実行します．

この記事では，まず，$S$ が与えられたとして，共役勾配法のアルゴリズムを示し，
その後に $S$ として良いものを紹介する，と言う順番で述べていきます．

## 共役勾配法のアルゴリズム

ここで，共役勾配法のアルゴリズムを書いておきます．

これは，$A\boldsymbol{x}=\boldsymbol{b}$ を満たす解 $\boldsymbol{x} = \boldsymbol{x}^*$ を求めるアルゴリズムです．

Input: $\boldsymbol{x}_0$ and $\epsilon$

Set $\boldsymbol{p}_0 = \boldsymbol{r}_0 = \boldsymbol{b} - A \boldsymbol{x}_0$

while $|| \boldsymbol{r}_j || > \epsilon$ do;

1. $\alpha_j \coloneqq || \boldsymbol{r}_j ||_2^2/\langle A\boldsymbol{p}_j, \boldsymbol{p}_j\rangle_2$
2. $\boldsymbol{x}_{j+1} \coloneqq \boldsymbol{x}_j + \alpha_j \boldsymbol{p}_j$
3. $\boldsymbol{r}_{j+1} \coloneqq \boldsymbol{r}_j - \alpha_j A\boldsymbol{p}_j$
4. $\beta_{j+1} \coloneqq ||\boldsymbol{r}_{j+1}||_2^2/||\boldsymbol{r}_j||_2^2$
5. $\boldsymbol{p}_{j+1}  \coloneqq \boldsymbol{r}_{j+1} + \beta_{j+1} \boldsymbol{p}_j$
6. $j \coloneqq j+1$

done

## 前処理つき共役勾配法のアルゴリズムの導出

まずは，前処理行列 $S$ が与えられているとして，前処理付き共役勾配法のアルゴリズムを導出します．

前処理後の方程式，

$$
S^\mathrm{T} A S \boldsymbol{y} = S^\mathrm{T} \boldsymbol{b}
$$ 

に対して，上記の共役勾配法を適用すると，

Input: $y_0$ and $\epsilon$

Set $\boldsymbol{p}_0 = \boldsymbol{r}_0 = S^\mathrm{T} (\boldsymbol{b} - AS \boldsymbol{y}_0)$

while $|| \boldsymbol{r}_j || > \epsilon$ do;

1. $\alpha_j \coloneqq || \boldsymbol{r}_j ||_2^2/\langle S^\mathrm{T}AS \boldsymbol{p}_j, \boldsymbol{p}_j\rangle_2$
2. $\boldsymbol{y}_{j+1} \coloneqq \boldsymbol{y}_j + \alpha_j \boldsymbol{p}_j$
3. $\boldsymbol{r}_{j+1} \coloneqq \boldsymbol{r}_j - \alpha_j S^\mathrm{T}AS \boldsymbol{p}_j$
4. $\beta_{j+1} \coloneqq ||\boldsymbol{r}_{j+1}||_2^2/||\boldsymbol{r}_j||_2^2$
5. $\boldsymbol{p}_{j+1}  \coloneqq \boldsymbol{r}_{j+1} + \beta_{j+1} \boldsymbol{p}_j$
6. $j \coloneqq j+1$

done

となる．この反復法によって得られた $\boldsymbol{y}$ に対し，

$$
\boldsymbol{x} = S\boldsymbol{y}
$$ 

とすることで，求めたい解 $\boldsymbol{x}$ が求まる．

これをもう少し変形して，直接，解 $\boldsymbol{x}$ を求めるアルゴリズムに書き換えておこう．


## 記号の導入と式変形

まずは，以下のように記号を定義する．

$$
\begin{align*}
\boldsymbol{p}^\prime_j &\coloneqq S\boldsymbol{p}_j \\
\boldsymbol{z}_j &\coloneqq S\boldsymbol{r}_j \\
\boldsymbol{r}^\prime_j &\coloneqq (S^\mathrm{T})^{-1} \boldsymbol{r}_j \\
P &\coloneqq SS^\mathrm{T}
\end{align*}
$$

今，$\boldsymbol{r}^\prime_j = (S^\mathrm{T})^{-1} \boldsymbol{r}_j = (S^\mathrm{T})^{-1} S^{-1} \boldsymbol{z}_j = P^\mathrm{T} \boldsymbol{z}_j$ である．
これらを使って，

$$
\begin{align*}
\langle \boldsymbol{r}_j, \boldsymbol{r}_j \rangle_2 &= \langle S^\mathrm{T} \boldsymbol{r}^\prime_j, S^{-1}\boldsymbol{z}_j \rangle_2 = \langle \boldsymbol{r}^\prime_j , SS^{-1}\boldsymbol{z}_j \rangle_2 = \langle \boldsymbol{r}_j^\prime, \boldsymbol{z}_j \rangle_2\\
\langle S^\mathrm{T}AS \boldsymbol{p}_j, \boldsymbol{p}_j \rangle_2 &= \langle S^\mathrm{T}A \boldsymbol{p}^\prime_j, S^{-1}\boldsymbol{p}_j^\prime\rangle_2 = \langle A \boldsymbol{p}^\prime_j, \boldsymbol{p}^\prime_j \rangle_2
\end{align*}
$$

となることがわかる．したがって，

$$
\begin{align*}
\alpha_j &= \frac{ \langle \boldsymbol{r}_j^\prime, \boldsymbol{z}_j \rangle_2}{\langle A\boldsymbol{p}_j^\prime, \boldsymbol{p}_j^\prime \rangle_2}\\
\beta_j &= \frac{\langle \boldsymbol{r}_{j+1}^\prime, \boldsymbol{z}_{j+1} \rangle_2}{\langle \boldsymbol{r}_j^\prime, \boldsymbol{z}_j \rangle_2}
\end{align*}
$$

となる

## 前処理付き共役勾配法のアルゴリズム

したがって，前処理付き共役勾配法のアルゴリズムは以下のようになる．
（注意：上までは $\prime$ がついていたものの $\prime$ を省略している．）

Input: $\boldsymbol{x}_0$ and $\epsilon$

Set $\boldsymbol{r}_0 = \boldsymbol{b} - A \boldsymbol{x}_0$ and $\boldsymbol{p}_0 = \boldsymbol{z}_0 = P\boldsymbol{r}_0$

while $|| \boldsymbol{r}_j || > \epsilon$ do;

1. $\alpha_j \coloneqq \langle \boldsymbol{r}_j, \boldsymbol{z}_j \rangle_2/\langle A\boldsymbol{p}_j, \boldsymbol{p}_j \rangle_2$
2. $\boldsymbol{x}_{j+1} \coloneqq \boldsymbol{x}_j + \alpha_j \boldsymbol{p}_j$
3. $\boldsymbol{r}_{j+1} \coloneqq \boldsymbol{r}_j - \alpha_j A\boldsymbol{p}_j$
4. $\beta_j \coloneqq \langle \boldsymbol{r}_{j+1}, \boldsymbol{z}_{j+1} \rangle_2/\langle \boldsymbol{r}_j, \boldsymbol{z}_j \rangle_2$
5. $\boldsymbol{z}_{j+1} \coloneqq P \boldsymbol{r}_{j+1}$
6. $\boldsymbol{p}_{j+1}  \coloneqq \boldsymbol{z}_{j+1} + \beta_{j+1} \boldsymbol{p}_j$
7. $j \coloneqq j+1$

done

# 前処理行列を決める

ここまでは，前処理行列 $S$ あるいは，$P = S^\mathrm{T} S$ が与えれていると言う前提で，前処理付き共役勾配法をアルゴリズムを示しました．

ここでは，前処理行列として，どのようなものを選べば良いかを述べ，具体例を挙げます．

最初にも述べましたが，前処理と，線形連立方程式を解くことは表裏一体で，
 $P_L = A^{-1}$ という前処理行列を選べば，それは，線形連立方程式 $A\boldsymbol{x} = \boldsymbol{b}$ を解くことそのものです．

今回は，共役勾配法に話を絞りましょう．

## 良い前処理とは

共役勾配法では， 行列 $A$ の相異なる固有値の個数 $d$ が最大の反復回数となります．

つまり，$d$ が小さいほど，連立方程式の解は求めやすいということなります．

よって，前処理を施して，相異なる固有値の個数を小さくできれば，前処理を施して良かったと言えます．

一方で，前処理には時間をさきたくありません．

行列 $A$ に前処理行列 $P_L$ や $P_R$ をかけるのが前処理だったので，行列 $A$ が密行列であると，処理時間の少ない前処理を施すことはなかなか大変です．

今回は， $A$ が疎行列と仮定することにします．

## 不完全Cholesky分解

前処理の例として， 不完全Cholesky分解を紹介します．


### Cholesky分解

その前にCholesky分解について簡単に説明します．

正定値対称行列 $A$ は下三角行列 $L$ と対角行列 $D$ を用いて，以下のように書ける．

$$
A = LDL^\mathrm{T}
$$

ここで，下三角行列 $L = (l_{ij})_{ij}$ は，

$$
l_{ij} = 
\begin{cases}
0 & \mathrm{for} \quad i < j \\
1 & \mathrm{for } \quad i = j \\
l_{ij} & \mathrm{for} \quad i > j
\end{cases}
$$

であるものとすると，このような分解は一意である．これをCholesky分解と呼ぶ．

これは，$S = LD^{1/2}$ とすれば，

$$
A = SS^\mathrm{T}
$$

となって，$S$ は下三角行列，$S^\mathrm{T}$ は上三角行列であるから，LU分解の対称行列の場合の特殊な表記と見ることができる．

$A = SS^\mathrm{T}$ に対して，$A = (a_{ij})_{ij}$, $S = (s_{ij})_{ij}$ とすると，$i \leq j$ に対して，

$$
a_{ij} = \sum_{k=0}^{j} s_{ik} s_{jk} = \sum_{k=0}^{j-1} s_{ik}s_{jk} + s_{ij}s_{jj}
$$

となるが，これは，

$$
\begin{align*}
s_{ii} &= \sqrt{a_{ii} - \sum_{k=1}^{i-1} s_{ik} s_{ik} }\\
s_{ij} &= \frac{1}{s_{jj}} \left[a_{ij} - \sum_{k=1}^{j-1} s_{ik}s_{jk} \right]
\end{align*}
$$

となって，これによって，$S$ の各成分 $s_{ij}$ を逐次的に求めることができる．

いま，$S$ と $L$, $D$ の関係は，$L = (l_{ij})_{ij}$, $D = \mathrm{diag}(d_1, \dots, d_n)$ とすると，

$$
s_{ij} = 
\begin{cases}
0 & \mathrm{for} \quad i < j \\
d_i & \mathrm{for}\quad i = j \\
l_{ij} & \mathrm{for}\quad i > j
\end{cases}
$$

であるから，結局，

$$
\begin{align*}
l_{ij} &= \frac{1}{d_j} \left[a_{ij} - \sum_{k=1}^{j-1} l_{ik}l_{jk} \right]\\
d_i  &= \sqrt{a_{ii} - \sum_{k=1}^{i-1} l_{ik} l_{ik}}
\end{align*}
$$

となる．

### Cholesky分解のアルゴリズム

Cholesky分解のアルゴリズムは以下のようになる．

Input: $A \in \mathbb{R}^{n\times n}$ symmetric and positive definite

1. $l_{11} = a_{11}$
2. for $i =2$ to $n$ do;
    1. for j = 1 to i-1 do;
        1. $l_{ij} \coloneqq a_{ij}$
        2. for $k = 1$ to $j-1$ do;
             - $l_{ij} \coloneqq l_{ij} - l_{ik}l_{jk}$
        3. $l_{ij} \coloneqq l_{ij}/d_j$ 
    2. $d_i \coloneqq a_{ii}$
    3. for $k=1$ to $i-1$ do;
       - $d_i \coloneqq d_i - l_{ik} l_{ik}$
    4. $d_i \coloneqq \sqrt{d_i}$

### 不完全Cholesky分解（0埋め）

行列 $A$ の不完全コレスキー分解とは，コレスキー分解を行う時に，

1. $\mathrm{M}$ に属する添字の場合には，コレスキー分解の計算を行わず，その要素を$0$とする．
2. $\mathrm{M}$ に属さない添字の場合には，通常通りコレスキー分解の計算を行い，要素の値を決定する．

という手順で，下三角行列 $\tilde{L}$ と対角行列 $\tilde{D}$ を決定する方法のことである．


つまり，$A$ の成分が$0$ のときは，$L$ の成分も $0$ にするということである．

$A$ が疎行列のとき，つまり，成分のほどんどが $0$ のときは，ほとんどの処理がスキップされるため，計算量は少ない．

しかし，これによって得られる分解，

$$
\tilde{L} \tilde{D} \tilde{L}^\mathrm{T}
$$

は，$LDL^\mathrm{T}$ ではない．つまり， $A$ ではない．が，とてもおしくはある．

# 不完全Cholesky分解前処理つき共役勾配法

不完全コレスキー分解によって得られる行列 $\tilde{L}$，$\tilde{D}$ をつかって前処理行列を構成するのが，不完全コレスキー分解前処理であり，これを用いた前処理付き共役勾配法が，不完全コレスキー分解前処理付き共役勾配法である．

すなわち，

$$
S = ((\tilde{L} \tilde{D}^{1/2})^\mathrm{T})^{-1}
$$

とする．

ただし， $D^{1/2}$ は 対角成分をそれぞれ $1/2$ 乗したもの， $D^{1/2} \equiv \mathrm{diag} ( \sqrt{d_1}, \dots, \sqrt{d_N})$ である．

また， このような $S$ を選ぶことは，

$$
P =  ((\tilde{L} \tilde{D}^{1/2})^\mathrm{T})^{-1}  (\tilde{L} \tilde{D}^{1/2}))^{-1} = (\tilde{L}^\mathrm{T})^{-1} \tilde{D}^{-1} \tilde{L}^{-1} = (\tilde{L} \tilde{D} \tilde{L}^\mathrm{T})^{-1}
$$

と選ぶことに対応する．

## これはよい毎処理か？

これが，確かに良い前処理であると納得するために，
一旦前処理として，不完全でない，完全な，Cholesky分解としてみよう．

つまり，

$$
A = LDL^\mathrm{T}
$$

として，

$$
S = ((LD^{1/2}))^{^1}
$$

を前処理行列に選んでみよう．

すると，

$$
A = LD^{1/2} D^{1/2}L^\mathrm{T} = (LD^{1/2})(LD^{1/2})^\mathrm{T}
$$

であるから，

$$
S^\mathrm{T} A S = (LD^{1/2})^{-1} A ((LD^{1/2})^\mathrm{T})^{-1} = I
$$

となる． $I$ は単位行列である．

もし，Cholesky分解が速やかに実行できて，$L$ と $D$ を得ることができるとすると，
前処理行列 $S$ として，$S = ((LD^{1/2})^\mathrm{T})^{-1}$ を選べば，
連立方程式 $A\boldsymbol{x} = \boldsymbol{b}$ の解は，

$$
\boldsymbol{x} = SS^Tb = (LDL^\mathrm{T})^{-1}b
$$

と求めることができる．

一方で，Cholesky分解の計算量は小さくない．

これを，不完全Cholesky分解に変えてみよう．
すなわち，

$$
S \equiv ((\tilde{L}\tilde{D}^{1/2})^\mathrm{T})^{-1}
$$

とする．すると， $S^\mathrm{T} A S$ は単位行列 $I$ にはならない．
しかし，それに近いもの，もう少し厳密に言うと， $S^\mathrm{T} A S$ の相異なる固有値の個数は，もとの行列 $A$ のそれ良いも小さくなっていると期待できる．


したがって，不完全Cholesky分解によって，前処理を施すと，よりはやく解を得ることができる．

# アルゴリズムを実装

## 不完全コレスキー分解前処理付き共役勾配法のアルゴリズム

不完全コレスキー分解前処理付き共役勾配法のアルゴリズムをまとめておく．

Input: $\boldsymbol{x}_0$ and $\epsilon$

Set $\boldsymbol{r}_0 = \boldsymbol{b} - A \boldsymbol{x}_0$ and $\boldsymbol{p}_0 = \boldsymbol{z}_0 = (\tilde{L} \tilde{D} \tilde{L}^\mathrm{T})^{-1} \boldsymbol{r}_0$

while $|| \boldsymbol{r}_j || > \epsilon$ do;

1. $\alpha_j \coloneqq \langle \boldsymbol{r}_j, \boldsymbol{z}_j \rangle_2/\langle A\boldsymbol{p}_j, \boldsymbol{p}_j \rangle_2$
2. $\boldsymbol{x}_{j+1} \coloneqq \boldsymbol{x}_j + \alpha_j \boldsymbol{p}_j$
3. $\boldsymbol{r}_{j+1} \coloneqq \boldsymbol{r}_j - \alpha_j A \boldsymbol{p}_j$
4. $\boldsymbol{z}_{j+1} \coloneqq (\tilde{L} \tilde{D} \tilde{L}^\mathrm{T})^{-1} \boldsymbol{r}_{j+1}$
5. $\beta_{j+1} \coloneqq \langle \boldsymbol{r}_{j+1}, \boldsymbol{z}_{j+1} \rangle_2/\langle \boldsymbol{r}_j, \boldsymbol{z}_j \rangle_2$
6. $\boldsymbol{p}_{j+1}  \coloneqq \boldsymbol{z}_{j+1} + \beta_{j+1} \boldsymbol{p}_j$
7. $j \coloneqq j+1$

done

ここで，$\boldsymbol{z}_j = (\tilde{L} \tilde{D} \tilde{L}^\mathrm{T})^{-1} \boldsymbol{r}_j \,\,$  ($j = 0, 1, \dots$) の計算は，$\tilde{L}\tilde{D}^{1/2}(\tilde{L}\tilde{D}^{1/2})^\mathrm{T}) \boldsymbol{z}_j = \boldsymbol{r}_j$ としておいて，

$$
\begin{align*}
\tilde{L}\tilde{D}^{1/2} \boldsymbol{w}_j = \boldsymbol{r}_j\\
(\tilde{L}\tilde{D}^{1/2})^\mathrm{T}\boldsymbol{z}_j = \boldsymbol{w}_j
\end{align*}
$$

と2段階で解く．第1式は前進代入，第2式は後退代入で解くことができる．

## C++による実装

以下の実装では，独自に実装した，行列演算のテンプレートクラス，`matrix<T>` とベクトルの演算を定義した，`vector` コンテナの拡張を用いているが，詳細は，[github](https://github.com/ksttr/linear-algebra/tree/main/solver/src)を参照してください．

### 前進代入・後退代入

以下は，前進代入と後退代入の実装である．

```cpp:sub.h
template <typename T>
vector<T> forward_substitution(const matrix<T> &A, const vector<T> &b)
{
    vector<T> x(b);

    for (size_t row = 0; row < A.row; row++)
    {
        for (size_t i = 0; i < row; i++)
        {
            x[row] -= A[row][i] * x[i];
        }
        x[row] *= 1 / A[row][row];
    }
    return x;
}

template <typename T>
vector<T> back_substitution(const matrix<T> &A, const vector<T> &b)
{
    vector<T> x(b);
    size_t N = b.size();

    x[N - 1] = b[N - 1] / A[N - 1][N - 1];
    for (int row = N - 2; row >= 0; row--)
    {
        for (size_t i = row + 1; i < N; i++)
        {
            x[row] -= A[row][i] * x[i];
        }
        x[row] *= 1 / A[row][row];
    }
    return x;
}
```

### 不完全Cholesky分解

以下は，不完全Cholesky分解の実装である．

```cpp:ic.h
template <typename T>
matrix<T> incomplete_cholesky_factorization(const matrix<T> &A, double precision)
{
  matrix<T> L(A.row, A.column, 0);

  L[0][0] = A[0][0];
  for (size_t i = 1; i < L.row; i++)
  {
      for (size_t j = 0; j < i; j++)
      {
          if (fabs(A[i][j]) < precision)
              continue;

          L[i][j] = A[i][j];
          for (size_t k = 0; k < j; k++)
          {
              L[i][j] -= L[i][k] * L[j][k];
          }
          L[i][j] *= 1 / L[j][j];
      }

      L[i][i] = A[i][i];
      for (size_t k = 0; k < i; k++)
      {
          L[i][i] -= L[i][k] * L[i][k];
      }
      L[i][i] = sqrt(L[i][i]);
  }

  return L;
}
```

### 不完全Cholesky分解前処理つき共役勾配法(ICCG)

以下は，不完全Cholesky分解前処理つき共役勾配法の実装である．

```cpp:iccg.h
template <typename T>
vector<T> incomplete_cholesky_factorization_conjugate_gradient(const matrix<T>& A, const vector<T>& b, const vector<T> &initial_guess, const T &precision)
{
  auto x = initial_guess;

  matrix<T> L = incomplete_cholesky_factorization(A);
  matrix<T> L_T = L.transpose();

  vector<T> r(N);
  auto k = A * x;
  r = b - (A * x);

  vector<T> z(N);
  auto w = forward_substitution(L, r);
  z = back_substitution(L_T, w);

  vector<T> p(z);

  for (size_t i = 0; i < N; i++)
  {
      auto t = A * p;
      auto dot_r_z = r * z;
      auto alpha = dot_r_z / (t * p);

      x += alpha * p;
      r -= alpha * t;

      if (r * r < precision)
          break;

      w = forward_substitution(L, r);
      z = back_substitution(L_T, w);

      auto beta = r * z / dot_r_z;
      p = z + beta * p;
  }

  return x;
}
```


# 参考文献

1. Wendland, H. (2017). Numerical Linear Algebra: An Introduction (Cambridge Texts in Applied Mathematics). Cambridge: Cambridge University Press. doi:10.1017/9781316544938
2. 森 正武, 共立数学講座12 数値解析 (第2版), 共立出版 (2002)