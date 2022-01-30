---
title: "共役勾配法"
emoji: "⛳"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["数学", "cpp"]
published: true
---

# 共役勾配法とは

線形連立方程式，

$$
Ax = b
$$

を数値的に解く方法．

ここで，行列 $A$ は，正定値対称行列でないといけない．

共役勾配法は，上記の連立方程式を，

$$
f(\boldsymbol{x}) = \boldsymbol{x}^\mathrm{T} A \boldsymbol{x} - \boldsymbol{x}^\mathrm{T} \boldsymbol{b}
$$

の極小値として求める方法です．


## 対称行列

行列 $A$ が対称とは，

$$
A^\mathrm{T} = A
$$

となることである．

## 正定値行列

対称行列 $A$ が正定値とは，任意のベクトル $\boldsymbol{x}$ に対して，

$$
\boldsymbol{x}^\mathrm{T} A \boldsymbol{x} > 0
$$

となることである．

- 正定値対称行列の固有値はすべて正である．
- 正定値対称行列は正則行列で，$A^{-1}$ も正定値対称行列である．


# 共役勾配法の詳細

ここで，共役勾配法の詳細を少し書いておきます．

共役勾配法とは，

$$
\min_{\boldsymbol{x} \in \mathbb{R}^N} f(x), \quad \mathrm{where} \quad f(\boldsymbol{x}) = \boldsymbol{x}^{\mathrm{T}} A \boldsymbol{x} - \boldsymbol{x}^{\mathrm{T}} \boldsymbol{b}
$$

を解く方法であって，求まった解$x = x^*$は，線形連立方程式

$$
A\boldsymbol{x}^* = \boldsymbol{b}
$$

の解でもある．

具体的には，反復法を用いる．

つまり，

$$
f(\boldsymbol{x}_{j+1}) \leq f(\boldsymbol{x}_j), \quad j = 1,2,\dots
$$

を満たす点列 $\boldsymbol{x}_1, \boldsymbol{x}_2, \dots$ を構成する．

そして，$f(\boldsymbol{x}_{M+1}) = f(\boldsymbol{x}_M)$ となる$M$があったとき，$x_M$ が求めたかった 解 $x^*$ となる．

実際には，数値計算で厳密な解を求められることはないので，前もって求めたい精度 $\epsilon$ を決めておいて，

$$
||\boldsymbol{b} - A\boldsymbol{x}^*||_2 < \epsilon
$$

を満たす $\boldsymbol{x}^*$ を求める．

## 最小化問題と線形連立方程式

$f(\boldsymbol{x}) = \boldsymbol{x}^\mathrm{T} A \boldsymbol{x} - \boldsymbol{x}^\mathrm{T} \boldsymbol{b}$ の極小値が線形連立方程式 $A\boldsymbol{x} = \boldsymbol{b}$ の解になっていることを確かめよう．

$f(x)$ を微分してみると，

$$
\begin{align*}
\nabla f(\boldsymbol{x}) &= A\boldsymbol{x} - \boldsymbol{b} \\
H f(\boldsymbol{x}) &= A
\end{align*}
$$

である．
ただし， $\boldsymbol{x} = (x_1, \dots, x_{N})^\mathrm{T}$ に対して， $\nabla f = (\partial f/ \partial x_1, \dots, \partial f/ \partial x_N)^{\mathrm{T}}$， $Hf = (\partial^2 f / \partial x_i \partial x_j)_{ij}$ のことである．つまり勾配とヘッセ行列．

これを使って，$f(\boldsymbol{x})$ を $\boldsymbol{y}$ の周りでテイラー展開してみると，

$$
\begin{align*}
f(\boldsymbol{x}) &= f(\boldsymbol{y}) + (\boldsymbol{x}-\boldsymbol{y})^{\mathrm{T}}\nabla f(\boldsymbol{y}) + \frac{1}{2} (\boldsymbol{x}-\boldsymbol{y})^{\mathrm{T}} Hf(\boldsymbol{y}) (\boldsymbol{x}-\boldsymbol{y}) \\
       &= f(\boldsymbol{y}) +  (\boldsymbol{x}-\boldsymbol{y})^{\mathrm{T}}(A\boldsymbol{x}-\boldsymbol{b}) + \frac{1}{2} (\boldsymbol{x}-\boldsymbol{y})^{\mathrm{T}} A (\boldsymbol{x}-\boldsymbol{y}) \\
      & \geq f(\boldsymbol{y}) +  (\boldsymbol{x}-\boldsymbol{y})^{\mathrm{T}}(A\boldsymbol{x}-\boldsymbol{b})
\end{align*}
$$

となって，$\boldsymbol{y}$ が $A\boldsymbol{x}=\boldsymbol{b}$ の解 $\boldsymbol{x}^*$ の時に $f$ が最小となることがわかる．
ここで，最後の不等式は， $A$ が正定値だからである．

つまり，最小化問題

$$
\min_{\boldsymbol{x} \in \mathbb{R}^N} f(\boldsymbol{x})
$$

の解 $\boldsymbol{x}^*$ と，線形連立方程式

$$
A\boldsymbol{x}=\boldsymbol{b}
$$

の解 $\boldsymbol{x}^*$ は同じである．


## 反復法の具体的な構成

以上をまとめると，

連立方程式 $A\boldsymbol{x} = b$ の解 $\boldsymbol{x} = \boldsymbol{x}^*$ を求めるには，

$f(\boldsymbol{x}) = \boldsymbol{x}^\mathrm{T} A \boldsymbol{x} - \boldsymbol{x}^\mathrm{T} \boldsymbol{b}$ に対して，

$$
f(\boldsymbol{x}_{j+1}) \leq f(\boldsymbol{x}_{j}), \quad j = 1,2, \dots
$$

となる点列 $\{\boldsymbol{x}_j\}$ を構成すればよい．

ここからは，この点列 $\{\boldsymbol{x}_j\}$ を具体的に構成していきましょう．


$$
\boldsymbol{x}_{j+1} = \boldsymbol{x}_j + \alpha_j \boldsymbol{p}_j
$$

と表しておいて，$\alpha_j$ と $\boldsymbol{p}_j$ を決めます．

### $\alpha_j$ を求める

さて，$f(\boldsymbol{x})$の形を思いだすと，$f(\boldsymbol{x})$は，下に凸であるから，

$$
\phi(\alpha) \equiv f(\boldsymbol{x}_j  + \alpha \boldsymbol{p}_j)
$$

と表しておくと，

$$
\phi^\prime (\alpha) = 0
$$

を満たす $\alpha$ を $\alpha_j$とすれば，$\boldsymbol{x}_j$ と $\boldsymbol{p}_j$ が与えられているときに，もっとも小さい $f(x_{j+1})$ を得ることができる．

したがって，

$$
0 = \phi^\prime(\alpha) = \mathrm{p}_j^{\mathrm{T}} \nabla f(\boldsymbol{x}_j + \alpha \boldsymbol{p}_j) =  \boldsymbol{p}_j^{\mathrm{T}} \left[A(\boldsymbol{x}_j + \alpha \boldsymbol{p}_j)- \boldsymbol{b}\right]
$$

すなわち，

$$
\alpha_j = \frac{\boldsymbol{p}_j^{\mathrm{T}} (\boldsymbol{b}-A\boldsymbol{x}_j)}{\boldsymbol{p}_j^{\mathrm{T}} A \boldsymbol{p}_j}
$$

となる．

ここで，いくつか記号を導入しておく．

$$
\boldsymbol{r}_j \equiv \boldsymbol{b} - A\boldsymbol{x}_j
$$

は，解 $\boldsymbol{x}^*$ との誤差を表すベクトルである．

また，任意の $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^N$ に対して， $A$ に関する内積を，

$$
\langle \boldsymbol{x}, \boldsymbol{y} \rangle_A \equiv \boldsymbol{y}^{\mathrm{T}} A \boldsymbol{x}
$$

特に，$A$ が単位行列 $I$ のときの内積を $\langle \boldsymbol{x}, \boldsymbol{y} \rangle_2$ と表し，

$$
\langle \boldsymbol{x}, \boldsymbol{y} \rangle_I \equiv \langle \boldsymbol{x}, \boldsymbol{y} \rangle_2 \equiv \boldsymbol{x}^\mathrm{T} \boldsymbol{y} = \sum_{i=1}^{N} x_i y_i.
$$

これは，いつもの2乗して足していく内積である．

また，

$$
\langle x,y \rangle_A = \langle A\boldsymbol{x}, \boldsymbol{y}\rangle_2 = \langle \boldsymbol{x}, A^{\mathrm{T}} \boldsymbol{y} \rangle_2
$$

が成り立つ．

これらの記号をつかうと， $\alpha_j$ は，

$$
\alpha_j = \frac{\langle \boldsymbol{r}_j, \boldsymbol{p}_j \rangle_2}{\langle \boldsymbol{p}_j, \boldsymbol{p}_j \rangle_A}
$$

と表せる．

### $\boldsymbol{p}_j$ を求める

次に，探索方向を表す $\boldsymbol{p}_j$ を求めよう．


#### $A$ 共役なベクトルの組 $\{\boldsymbol{p}_j\}$

これは，$A$ 共役であるベクトルの組 $\boldsymbol{p}_0, \dots, \boldsymbol{p}_{N-1}$ が選ばれる．

ここで，$\boldsymbol{p}_0, \dots, \boldsymbol{p}_{N-1}$ が $A$ 共役であるとは， $0 \leq i \not= j \leq N-1$ に対して，

$$
\langle \boldsymbol{p}_j, \boldsymbol{p}_k \rangle_A = 0
$$

であることを言う．

つまり， $A$ に関する内積による直交基底を選べ，と言うことである．

なぜ，$A$ 共役を選ぶかというと，上で選んだ $\alpha_j$に対して，解 $\boldsymbol{x}^*$ が，

$$
x^* = x_0 + \sum_{j=0}^{N-1} \alpha_j p_j
$$

と表すことができるからである．

つまり，解 $\boldsymbol{x}^*$ を基底 $\{\boldsymbol{p_j}\}$ の線型結合で表した時の係数がちょうど $\alpha_j$ になるのである．

よって，このように$\alpha_j, \boldsymbol{p}_j$ を選んでおくと，最悪でも，$N$ 回で解に辿り着くことができる．

#### $\boldsymbol{p}_j$ の具体的な構成方法

では， $\boldsymbol{p}_j$ は具体的にはどのように選べば良いでしょうか．

$A$ 共役に選べば，最大で $N$ 回の探索で十分でしたが，より少ない回数の反復で終わらせたいです．

また，$A$ 共役なベクトルの組 $p_0, \dots, p_{N-1}$ は，最初にすべて計算しておくのではなく，反復するごとに，その都度計算した方が，計算効率がよいです．

そこで，ここでは，$p_0, \dots, p_i$ から $p_{j+1}$ を構成し，また，
$x_{j}$が終了条件を満たさなかった時にのみ，次の探索方向$p_{j+1}$を決めることとします．

この場合，使用できるのは，

- これまでの探索方向，$\boldsymbol{p}_0,\dots, \boldsymbol{p}_{j-1}$，
- これまでの解の候補 $\boldsymbol{x}_0, \dots, \boldsymbol{x}_j$，
- 残差ベクトル $\boldsymbol{r}_0, \dots, \boldsymbol{r}_{j}$

です．

天下り的ですが，以下のように表してみましょう．

$$
\boldsymbol{p}_{j} = \boldsymbol{r}_j + \sum_{k=0}^{j-1} \beta_j^{(k)} \boldsymbol{p}_k
$$

$\boldsymbol{p}_0, \dots, \boldsymbol{p}_{N-1}$ は $A$ 共役であるので，$0 \leq l \leq j-1$ に対して，

$$
0 = \langle \boldsymbol{p}_l, \boldsymbol{p}_{j} \rangle_A = \langle \boldsymbol{p}_l, \boldsymbol{r}_j \rangle_A + \sum_{k=0}^{j-1} \beta_j^{(k)}  \langle \boldsymbol{p}_l, \boldsymbol{p}_k \rangle_A =  \langle \boldsymbol{p}_l, \boldsymbol{r}_j \rangle_A + \beta_j^{(l)} \langle \boldsymbol{p}_l, \boldsymbol{p}_l \rangle_A
$$

であるから，

$$
 \beta_j^{(l)} = -\frac{\langle \boldsymbol{p}_l, \boldsymbol{r}_j \rangle_A}{\langle \boldsymbol{p}_l, \boldsymbol{p}_l \rangle_A} \quad \mathrm{for} \quad l = 0,\dots, j-1
$$

このように，$\boldsymbol{p}_j$ を構成するメリットは，実は， $\beta_j^{(0)} = \cdots = \beta_j^{(j-2)} = 0$ となることである． つまり，$\beta_j = \beta_j^{(j-1)}$ とすると，

$$
p_j = r_j + \beta_j p_{j-1}
$$

となって，これはメモリー使用量的にもとても良い形である．

さらに，この形は，他の探索方向の候補よりも効率がよい．次はこれをより正確に述べよう．

### Krylov空間

まずは，より正確に述べるために，いくつか記号を定義する．

Krylov空間とは，ある行列 $A \in M_n(\mathbb{R})$ とベクトル $\boldsymbol{p} \in \mathbb{R}^N$ に対して， $i$ 個のベクトルの組 $\boldsymbol{p}, A\boldsymbol{p}, \dots, A^{i-1}\boldsymbol{p}$ が生成する部分空間

$$
\mathcal{K}_i(A,\boldsymbol{p}) \equiv \mathrm{span}\{\boldsymbol{p}, A\boldsymbol{p}, \dots, A^{i-1} \boldsymbol{p}\}
$$

である．

ただし， $n$ 個のベクトル $\boldsymbol{v}_1,\dots, \boldsymbol{v}_n$ が生成する空間 $\mathrm{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_n\}$ とは，$\boldsymbol{v}_1, \dots, \boldsymbol{v}_n$ の線型結合で表されるベクトル全体の集合のことである．

上記で構成した，$\boldsymbol{r}_0, \dots, \boldsymbol{r}_{N-1}$ と $\boldsymbol{p}_0, \dots, \boldsymbol{p}_{N-1}$ に関して，

$$
\mathcal{K}_i(A, \boldsymbol{r}_0) = \mathrm{span} \{\boldsymbol{r}_0, A\boldsymbol{r}_0, \dots, A^{i-1} \boldsymbol{r}_0\} = \mathrm{span}\{\boldsymbol{r}_0, \dots, \boldsymbol{r}_{i-1}\}= \mathrm{span}\{\boldsymbol{p}_0, \dots, \boldsymbol{p}_{i-1}\}
$$

が成り立ち， $\dim \mathcal{K}_i(A, \boldsymbol{r}_0) = i$ である．

また， $\boldsymbol{x}-\boldsymbol{y} \in \mathcal{K}_i(A, \boldsymbol{p})$ のことを， $\boldsymbol{x} \in \boldsymbol{y} + \mathcal{K}_i(A, \boldsymbol{p})$ と書く．

これらの記号を使うと， $j$ 番目の解候補 $\boldsymbol{x}_j$ は，$\boldsymbol{x}_0 + \mathcal{K}_j(A,\boldsymbol{r}_0)$ の元であると言える．

したがって，共役勾配法では，解候補 $\boldsymbol{x}_j$ を， $\boldsymbol{x}_0 + \mathcal{K}_j(A, \boldsymbol{r}_0)$ から探索していることになる．

そして，共役勾配法で得られる解候補 $\boldsymbol{x}_j$ は，$\boldsymbol{x}_0 + \mathcal{K}_j(A, \boldsymbol{r}_0)$ の他の元よりも，解 $\boldsymbol{x}^*$ に近い．すなわち，


$$
|| \boldsymbol{x}^* - \boldsymbol{x}_j||_A \leq || \boldsymbol{x}^* - \boldsymbol{x} ||_A, \quad \boldsymbol{x} \in \boldsymbol{x}_0 + \mathcal{K}_j(A, \boldsymbol{r}_0)
$$

がなりたつ．ただし， $|| \boldsymbol{a}||_A = \langle \boldsymbol{a}, \boldsymbol{a} \rangle_A$ である．

さらに，反復によって，解に毎回近づいていく．

$$
|| \boldsymbol{x}^* - \boldsymbol{x}_{j+1} ||_A \leq || \boldsymbol{x}^* - \boldsymbol{x}_j||_A
$$

### 最大の反復回数

反復の最大回数は， $A$ の相異なる固有値の個数 $d$ である．

これは，以下の不等式からわかる．

$\pi_i$ を最大次数が $i$ 次の多項式全体の集合とし，

$P(t) = \sum_{j=0}^{i-1} \gamma_j t^{j} \in \pi_{i-1}$ に対して，
$P(A) = \sum_{j=0}^{i-1} \gamma_j A^{j}$ と定める．

$A$ の固有値 $\lambda$ の固有ベクトルを $\boldsymbol{x}$ とすると，

$$
P(A)\boldsymbol{x} = P(\lambda) \boldsymbol{x}
$$

である．

このとき， $A$の固有値を，$\lambda_1, \dots, \lambda_N$， $A\boldsymbol{x} = \boldsymbol{b}$ の解を $\boldsymbol{x}^*$ ，共役勾配法の $j$ ステップ目の回候補を， $\boldsymbol{x}_j$ とすると，

$$
|| \boldsymbol{x}^* - \boldsymbol{x}_j ||_A \leq \min_{P \in \pi_{i}, P(0)=1} \max_{1\leq j \leq N} |P(\lambda_j)| ||\boldsymbol{x}^* - \boldsymbol{x}_0||_A
$$

が成り立つ．

この不等式で， $A$ の相異なる固有値の個数を $d$ に対して， $i \geq d$ とすると，

$$
\min_{P \in \pi_{i}, P(0)=1} \max_{1\leq j \leq N} |P(\lambda_j)|  = 0 \quad \mathrm{for} \quad i \geq d
$$

であるから ( $d$ 次多項式以上だと， $P(\lambda_j) = 0$ をすべて満たすようなものが存在する)，最大の反復回数は $d$ ということがわかる．


ここらへんの証明は，たとえば，参考文献1を参照してください．

## 若干の変形

最後に，数値計算しやすいように，若干の式変形をする．

ここで，数値計算しやすい形とは，

- 行列とベクトルの積を含む内積，$\langle \cdot \rangle_A$ よりも $\langle \cdot \rangle_2$ が現れるような形
- 内積よりもノルムで計算できるような形
 
という感じである．

今，

$$
\begin{align*}
\alpha_j &= \frac{\langle \boldsymbol{r}_j, \boldsymbol{p}_j \rangle_2}{\langle \boldsymbol{p}_j, \boldsymbol{p}_j \rangle_A} \\
\beta_{j+1} &= -\frac{\langle \boldsymbol{r}_{j+1}, p_j \rangle_A }{\langle \boldsymbol{p}_j, \boldsymbol{p}_j \rangle_A}
\end{align*}
$$

であって，これらを変形しておく．

$$
\langle \boldsymbol{r}_j, \boldsymbol{p}_j \rangle_2 = \langle \boldsymbol{r}_j, \boldsymbol{r}_j + \beta_j \boldsymbol{p}_{j-1} \rangle_2 = || \boldsymbol{r}_j ||_2^2 + \beta_j \langle \boldsymbol{r}_j, \boldsymbol{p}_{j-1} \rangle_2 = || \boldsymbol{r}_j ||_2^2
$$

最後の等号は$\alpha_j$を決めるときの等式($\phi^\prime(\alpha) = 0$)による．

また，

$$
A\boldsymbol{p}_j = \frac{1}{\alpha_j} (\boldsymbol{x}_{j+1} - \boldsymbol{x}_{j}) =  \frac{1}{\alpha_j}(\boldsymbol{r}_j - \boldsymbol{r}_{j+1})
$$

より，

$$
\begin{align*}
\alpha_j &= \frac{||\boldsymbol{r}_j||_2}{\langle A\boldsymbol{p}_j,\boldsymbol{p}_j \rangle_2} \\
\beta_{j+1} &= \frac{||\boldsymbol{r}_{j+1}||_2^2}{||\boldsymbol{r}_j||_2^2}
\end{align*}
$$

となる．

# アルゴリズムと実装

## 共役勾配法のアルゴリズム

ここまでをアルゴリズムにすると，以下のようになる．

Input: $x_0$ and $\epsilon$

Set $p_0 = r_0 = b - A x_0$

while $|| r_j || > \epsilon$ do;

1. $\alpha_j \coloneqq \frac{|| r_j ||_2^2}{\langle Ap_j, p_j\rangle_2}$
2. $x_{j+1} \coloneqq x_j + \alpha_j p_j$
3. $r_{j+1} \coloneqq r_j - \alpha_j Ap_j$
4. $\beta_{j+1} \coloneqq \frac{||r_{j+1}||_2^2}{||r_j||_2^2}$
5. $p_{j+1}  \coloneqq r_{j+1} + \beta_{j+1} p_j$
6. $j \coloneqq j+1$

done

## C++での実装

```cpp:cg.h
   template <typename T>
    vector<T> conjugate_gradient(const matrix<T>& A, const vector<T> &b, const vector<T> &initial_guess, const T &precision)
    {
        auto x = initial_guess;

        vector<T> p(N);
        auto t = A * x;
        p = b - t;

        vector<T> r(p);
        auto r_square = r * r;

        for (size_t i = 0; i < N; i++)
        {
            t = A * p;

            auto alpha = r_square / (t * p);

            x += alpha * p;
            r -= alpha * t;

            auto r_square_new = r * r;
            if (r_square_new < precision)
                break;

            auto beta = r_square_new / r_square;

            p = r + beta * p;

            r_square = r_square_new;
        }

        return x;
    }
```

ここで，テンプレートクラス `matrix` と，ベクトルの演算を定義した `vector` クラスの拡張を別途定義している．詳細は[github](https://github.com/ksttr/linear-algebra/tree/main/solver/src)に公開している．

# この後の展開

この後の展開としては，以下がある．

## 一般の行列に対する反復解法

行列 $A$ は対称正定値行列に限られていたので，これを一般の行列$A$にまで拡張させようというもの

### CGNR法

$A\boldsymbol{x} = \boldsymbol{b}$ の両辺に転置 $A^\mathrm{T}$ をかけて

$$
A^{\mathrm{T}} A \boldsymbol{x} = A^{\mathrm{T}}\boldsymbol{b}
$$

をとすれば，CG法が使えるぞっという方法．

$$
|| \boldsymbol{x}^* - \boldsymbol{x}_j||_{A^\mathrm{T}A} = || \boldsymbol{b} - A\boldsymbol{x}_j||_2
$$

であるから，結局，残差ベクトル $\boldsymbol{r}_j$ を最小化させていることとなる．
これより，Conjugate Gradient on th Normal equation to minimise the Residual と呼ばれる．

### GMRES and MINRES

上の発想から，はじめから

$$
\min\{||\boldsymbol{b}-A\boldsymbol{x}||_2 \, | \, \boldsymbol{x} \in \boldsymbol{x}_0 + \mathcal{K}_j(A,\boldsymbol{r}_0)\}
$$

を考えようというものが，GMRESとMINRES.

## 制約付きの最小化問題

制約がついている場合にも対応可能．

### リーマン多様体上の共役勾配法

制約を $\mathbb{R}^N$ に埋め込まれた多様体として表現しておいて，この多様体上で共役勾配法をすれば，無制約の最小化問題に帰着できるというもの．

などなど．

少しづづ記事にしていく予定ではある．

# 参考文献

1. Wendland, H. (2017). Numerical Linear Algebra: An Introduction (Cambridge Texts in Applied Mathematics). Cambridge: Cambridge University Press. doi:10.1017/9781316544938
2. 森 正武, 共立数学講座12 数値解析 (第2版), 共立出版 (2002)
