---
title: "直交多項式のゼロ点の数値計算（行列の固有値問題）"
emoji: "🌟"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["数学", "rust"]
published: false
---

# 行列の固有値問題として求める

ルジャンドル多項式やラゲール多項式などの直交多項式は三項漸化式をを持つ。この性質を利用して、直交多項式のゼロ点を行列の固有値問題に帰着させて求める。

直交多項式 $P_n$ は一般に以下のような三項漸化式を持つ。

$$
a_n P_{n+1}(x) + b_n P_n(x) + c_n P_{n-1}(x) = x P_n(x), \quad n = 0,1,\dots
$$

ここで、$a_n, b_n, c_n$ は $x$ に依らない定数であり、

$$
P_{-1}(x) = 0, \quad P_{0}(x) = 1
$$

である。これを行列の形に表すと、

$$
a_{N-1} P_{N}(x)
\begin{pmatrix}
    0 \\
    0 \\
    0 \\
    \vdots \\
    0 \\
    1
\end{pmatrix}
+
\begin{pmatrix}
b_{0} & a_{0} & 0 & 0 & \cdots & 0 & 0 & 0 \\
c_1 & b_1 & a_1 & 0 & \cdots & 0 & 0 & 0 \\
0 & c_2 & b_2 & a_2 & \cdots & 0 & 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \\
0 & 0 & 0 & 0 & \cdots & c_{N-2} & b_{N-2} & a_{N-2}\\
0 & 0 & 0 & 0 & \cdots & 0 & c_{N-1} & b_{N-1}
\end{pmatrix}
\begin{pmatrix}
    P_0(x) \\
    P_1(x) \\
    P_2(x) \\
    \vdots \\
    P_{N-2}(x) \\
    P_{N-1}(x)
\end{pmatrix}
 = x
\begin{pmatrix}
    P_0(x) \\
    P_1(x) \\
    P_2(x) \\
    \vdots \\
    P_{N-2}(x) \\
    P_{N-1}(x)
\end{pmatrix}
$$

となる。これを以降は、

$$
a_{N-1} P_N(x) \boldsymbol{e}_N + \boldsymbol{\mathrm{J}}_N \boldsymbol{P}_N(x) = x \boldsymbol{P}_N(x)
$$

と表記する。

$N$ 次の直交多項式 $P_{N}(x)$ のゼロ点を $x_0$ とすると、 $P_N(x_0) = 0$ であるから、

$$
\boldsymbol{\mathrm{J}}_N \boldsymbol{P}(x_0) = x_0 \boldsymbol{P}_N(x_0)
$$

となり、ゼロ点 $x_0$ は、行列 $\boldsymbol{\mathrm{J}}_N$ の固有値であり、ベクトル $\boldsymbol{P}_N(x_0)$ はその時の固有ベクトルであることがわかる。ここで、$\boldsymbol{P}_{N}(x_0)$ は、$N-1$ 次以下の直交多項式の $x=x_0$ での値を集めたものである。

したがって、直交多項式のゼロ点を求める問題は、行列の固有値問題に帰着された。

さらに、対称行列の固有値問題まで帰着させることもできる。

$$
a_n P_{n+1}(x) + b_n P_n(x) + c_n P_{n-1}(x) = x P_n(x), \quad n = 0,1,\dots
$$

両辺に、適当な定数 $\lambda_n$ をかけると

$$
a_n \lambda_n P_{n+1}(x) + b_n \lambda_n P_n(x) + c_n \lambda_n P_{n-1}(x) = x \lambda_n P_n(x), \quad n = 0,1,\dots
$$

$\tilde{P}_{n}(x) = \lambda_n P_n(x)$ とすると、

$$
a_n \frac{\lambda_n}{\lambda_{n+1}} \tilde{P}_{n+1}(x) + b_n \tilde{P}_n(x) + c_n \frac{\lambda_n}{\lambda_{n-1}} \tilde{P}_{n-1}(x) = x \tilde{P}_n(x), \quad n = 0,1,\dots
$$

となる。対称行列となるには、

$$
a_n \frac{\lambda_n}{\lambda_{n+1}} = c_{n+1} \frac{\lambda_{n+1}}{\lambda_{n}}
$$

となればよい。すなわち、

$$
\frac{\lambda_{n+1}}{\lambda_n} = \sqrt{\frac{a_n}{c_{n+1}}}
$$

よって、

$$
\begin{pmatrix}
b_{0} & \sqrt{c_1 a_0} & 0 & 0 & \cdots & 0 & 0 & 0 \\
\sqrt{c_1 a_0} & b_1 & \sqrt{c_2 a_1} & 0 & \cdots & 0 & 0 & 0 \\
0 & \sqrt{c_2 a_1} & b_2 & \sqrt{c_3 a_2} & \cdots & 0 & 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \\
0 & 0 & 0 & 0 & \cdots & \sqrt{c_{N-2} a_{N-3}} & b_{N-2} & \sqrt{c_{N-1} a_{N-2}}\\
0 & 0 & 0 & 0 & \cdots & 0 & \sqrt{c_{N-1} a_{N-2}} & b_{N-1}
\end{pmatrix}
\begin{pmatrix}
    \tilde{P}_0(x_0) \\
    \tilde{P}_1(x_0) \\
    \tilde{P}_2(x_0) \\
    \vdots \\
    \tilde{P}_{N-2}(x_0) \\
    \tilde{P}_{N-1}(x_0)
\end{pmatrix}
 = x_0
\begin{pmatrix}
    \tilde{P}_0(x_0) \\
    \tilde{P}_1(x_0) \\
    \tilde{P}_2(x_0) \\
    \vdots \\
    \tilde{P}_{N-2}(x_0) \\
    \tilde{P}_{N-1}(x_0)
\end{pmatrix}
$$

となる。

# Rustによる実装

Rustで実装する。固有値問題を解く部分は、線形代数ライブラリ[nalgebra](https://nalgebra.org/)を使用する。

## 例1: ラゲール多項式

(一般化された)ラゲール多項式 $L_n^{(\alpha)}(x)$ は、以下の漸化式によって定義する。

$$
-(n+1) L_{n+1}(x) + (2n+ \alpha + 1) L_{n}(x) - (n + \alpha) L_{n-1}(x) = x L_{n}(x)
$$


# 参考文献

1. A. Gil, J. Segura, and, N.M. Temme, Numerical Methods for Special Functions. Society for Industrial and Applied Mathematics, 2007.

2. R. クーラント(著), D. ヒルベルト(著), 藤田宏(訳), 高見頴郎(訳), 石村直之(訳), 数理物理学の方法 上, 丸善出版 (2013)
