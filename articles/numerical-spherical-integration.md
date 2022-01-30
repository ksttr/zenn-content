---
title: "球面上積分の数値計算法"
emoji: "🗂"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["c","fft","nfft", "数学"]
published: false
---

# 球面上積分

球面上の積分

$$
\int_{\mathbb{S}^2} f(\boldsymbol{x}) ~d\boldsymbol{x}
$$

を数値的に実行する方法についての記事です．

$M$個の代表点，

$$
\{ \boldsymbol{x}_{0}, \dots, \boldsymbol{x}_{M-1}\},
$$

が与えられると，数値積分は，

$$
\int_{\mathbb{S}^2} f(\boldsymbol{x}) ~d\boldsymbol{x} \simeq \sum_{i=0}^{M-1} w_i f(\boldsymbol{x}_i)
$$

となる重み$w_i$を探すという問題になります．

ここで，$w_i \,\,(i = 0,\dots,M-1)$ は，以下の性質を満たします．

1. $w_i > 0$ for $\forall i$
2. $\sum_{i=0}^{M-1} w_i = 4\pi$

1.については，ベクトル $\boldsymbol{w} \equiv (w_i) > 0$ と表すことにします．

この記事では，

1. 代表点を自由に決められる場合
2. 代表点があらかじめ与えられている場合

の２つの場合それぞれについて説明します．

また，ここからは，簡単のため，$\mathbb{S}^2$は単位球面であるとします．
つまり，球面上の点$\boldsymbol{x}$は，

$$
\boldsymbol{x} = (\cos\theta\sin\varphi, \sin\theta\sin\varphi, \cos\varphi)
$$

と $(\theta, \varphi)$でかける場合を考えます．一般の半径$r$の球面上積分に対しては，$\boldsymbol{x} \rightarrow \boldsymbol{x}/r \equiv \boldsymbol{z}$と変数変換すると，同じ問題に帰着します．


# 代表点を自由に決められる場合

## Gauss-Legendre求積法

## 例

```math
f(x) = \frac{1}{1+x^2}
```

# 代表点があらかじめ与えられている場合

例えば，実験や観測をする場合，機材の設置場所の関係によって，得られるデータの位置が決まっている場合には，こちらの方法を使う必要があります．

この場合，代表点$\{\boldsymbol{x}_0, \dots, \boldsymbol{x}_{M-1}\}$は，

1. 等間隔でないかもしれない．
2. 一様に分布しているわけではないかもしれない．

という，なかなか厳しい場合も含まれます．

## 定式化

具体的な数値計算（プログラミング）に入る前に，数値計算の理論的な側面をざっくりと説明します．
この場合でも，問題は，

$$
\int_{\mathbb{S}^2} f(\boldsymbol{x}) ~d\boldsymbol{x} \simeq \sum_{i=0}^{M-1} w_i f(\boldsymbol{x}_i)
$$

を満たす（最も良い近似となるような）$w_i$を求めよ，ということに変わりありません．

しがしながら，代表点 $\boldsymbol{x}_i, \,\, (i=0,\dots,M-1)$ はすでに与えれらていて，こちらで，適当なものを選ぶことができません．

したがって，代表点の取り方は決まっているものとして，重み $w_i$ を決めるという問題になります．

### 球面調和関数による展開

この場合，上記のように勝手に代表点を決めることができないので難しいです．
そこで，関数$f(x)$を球面調和関数

$$
Y_k^n(\boldsymbol{x}) = Y_k^n(\theta, \varphi) \equiv \sqrt{\frac{2n+1}{4\pi}} P^n_{|k|}(\cos\theta) e^{ik\varphi}
$$

でできる限り展開します．ここで，$P^n_{|k|}$はLegendre関数です．
つまり，

$$
f(\boldsymbol{x}) \simeq f_N(\boldsymbol{x}) = \sum_{k=0}^{N} \sum_{n=-k}^{n=k} a^n_k Y_k^n(\boldsymbol{x})
$$

と近似します．$N$が大きければ大きいほどよい近似になります．ここで，展開係数$a_k^n$は，

$$
a_k^n = \int_{\mathbb{S}^2} f(\boldsymbol{x})\overline{Y_k^n(\boldsymbol{x})} ~d\boldsymbol{x}
$$

です．

### 重みの計算

これらから，重み$w_i$は，

$$
\boldsymbol{Y}^* \boldsymbol{w} = \sqrt{4\pi} \boldsymbol{e}_0
$$

を満たすことがわかります．ここで，

$$
\begin{align*}
\boldsymbol{Y} &\equiv (Y_k^n(\boldsymbol{x}_i)_{i=0,\dots, M-1; n=0,\dots N, |k| \leq n} \in \mathbb{C}^{M\times(N+1)^2} \\
\boldsymbol{e}_0 &\equiv (1,0,\dots, 0)^T \in \mathbb{R}^{(N+1)^2} \\
\boldsymbol{w} &\equiv (w_i)_{i=0,\dots,M-1} \in \mathbb{R}^M
\end{align*}
$$

で，$\boldsymbol{Y}^*$は，$\boldsymbol{Y}$の随伴行列です．

したがって，上記の連立方程式を満たす$\boldsymbol{w}$を求めれば良いです．連立方程式の解がない場合もありうるので，実際には，

$$
\min_{\boldsymbol{w} \geq 0} \left| \boldsymbol{Y}^* \boldsymbol{w} - \sqrt{4\pi} \boldsymbol{e}_0 \right|
$$

を求めることになります．

### （数値計算しやすい）同値な問題

この最小化問題は以下と同値です．

$$
\min_{\boldsymbol{w} \geq 0} \left| \boldsymbol{A} \boldsymbol{w} - \sqrt{4\pi}\boldsymbol{e}_0 \right| 
$$

ここで，

$$
\begin{align*}
\boldsymbol{A} &\equiv \begin{pmatrix} \boldsymbol{A}_1^T \\ \boldsymbol{A}_2^T \end{pmatrix} \in \mathbb{R}^{(N+1)^2\times M} \\
\boldsymbol{A}_1 &\equiv \mathrm{Re} (Y_k^n(\boldsymbol{x}_i))_{i=0,\dots, M-1, n=0,\dots,N, 0\leq k\leq n} \in \mathbb{R}^{M \times \frac{(N+1)(N+2)}{2}} \\
\boldsymbol{A}_1 &\equiv \mathrm{Im} (Y_k^n(\boldsymbol{x}_i))_{i=0,\dots, M-1, n=1,\dots,N, -n\leq k\leq -1} \in \mathbb{R}^{M \times \frac{N(N+1)}{2}}
\end{align*}
$$

以下の数値計算では，この最小化問題を解いていくことになります．

## 数値計算方法

上記の最小化問題を計算するときに，以下の2つの困難があります．

1. 巨大な行列とベクトルの積，$\boldsymbol{A} \boldsymbol{w}$を計算しなければならない．
2. 最小化問題を解かなければならない．

1.に関しては，離散フーリエ変換（DFT）の問題と同じですが，通常のDFTでは，代表点$\{\boldsymbol{x}_0, \dots, \boldsymbol{x}_{M-1}\}$は等間隔に並んでいますが，今回の問題では，等間隔とは限りません．
このような，不等間隔のフーリエ変換を実行するのはコストの大きな問題でしたが，最近では，$\mathcal{O}(N^2\log^2 N  + M)$で解けるアルゴリズムがあって，[C言語ライブラル](https://www-user.tu-chemnitz.de/~potts/nfft/index.php)も作られています．

2.に関しては，様々な方法が知られていますが，今回は共役勾配法の一種であるCGNR法を用います．ここで問題になるのは，重み$\boldsymbol{w}$は正定値でなければならないということです．普通にCGNR法によって最小化問題を解いても，正定値性は担保してくれません．

そこで，探索範囲を $\boldsymbol{w} > 0$ に限ったCGNR法を行なって，正定値性を担保します．

### NFFTによる計算

まず，NFFTによって， $\boldsymbol{Aw}$ を計算します．今回は，[NFFTライブラリ](https://www-user.tu-chemnitz.de/~potts/nfft/index.php)を使います．NFFTライブラリの特に，球面上フーリエ変換の方法については，[別記事](https://zenn.dev/ksttr/articles/nfft-tutorial-nfsft)を参照して下さい．



### CGNR法による最小化問題の解法

## 例

# 参考文献

1. Manuel Gr$\ddot{\mathrm{a}}$f, Stefan Kunis, and, Daniel Potts. On the computation of nonegative quadrature weights on the sphere. Appl. Comput. Harm. Anal. 27(1), 2009, doi:10.1016/j.acha.2008.12.003
2. Kerstin Hesse, Ian H. Sloan, and Robert S. Womersley. Numerical Integration on the Sphere. In W. Freeden et al. (eds.), Handbook of Geomathematics, Springer: Berlin 2015, doi:10.1007/978-3-642-54551-1_4