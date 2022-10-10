---
title: "直交多項式のゼロ点の数値計算（行列の固有値問題）"
emoji: "🌟"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["数学", "rust", "nalgebra"]
published: true
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

となる。対称行列になるには、

$$
a_n \frac{\lambda_n}{\lambda_{n+1}} = c_{n+1} \frac{\lambda_{n+1}}{\lambda_{n}}
$$

を満たせば良い。すなわち、

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

となる。したがって，

$$
\boldsymbol{\mathrm{S}}_N \equiv
\begin{pmatrix}
b_{0} & \sqrt{c_1 a_0} & 0 & 0 & \cdots & 0 & 0 & 0 \\
\sqrt{c_1 a_0} & b_1 & \sqrt{c_2 a_1} & 0 & \cdots & 0 & 0 & 0 \\
0 & \sqrt{c_2 a_1} & b_2 & \sqrt{c_3 a_2} & \cdots & 0 & 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \\
0 & 0 & 0 & 0 & \cdots & \sqrt{c_{N-2} a_{N-3}} & b_{N-2} & \sqrt{c_{N-1} a_{N-2}}\\
0 & 0 & 0 & 0 & \cdots & 0 & \sqrt{c_{N-1} a_{N-2}} & b_{N-1}
\end{pmatrix}
$$

の固有値が $N$ 次の直交多項式のゼロ点となる．

# Rust による実装

Rust で実装する。固有値問題を解く部分は、線形代数ライブラリ[nalgebra](https://nalgebra.org/)を使用する。

まずは，任意の直交多項式のゼロ点を計算する trait `Zeros` を実装する．

```rust
use nalgebra::{DMatrix, SymmetricEigen};

trait Zeros {
    // 直交多項式の自由度 N を求める
    fn degree(&self) -> usize;

    // 三項漸化式の係数 a_n, b_n, c_n を求める
    fn coefficients(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>);

    // 対称三重対角行列 S_N を求める
    // 返り値の型は，nalgebra::DMatrix<f64>
    fn symmetric_jacobi(&self) -> DMatrix<f64> {
        let n = self.degree();
        let (a, b, c) = self.coefficients();
        let mut symmetric_coeff = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            symmetric_coeff.push(libm::sqrt(c[i + 1] * a[i]));
        }
        DMatrix::<f64>::from_fn(n, n, |r, c| {
            if r == c {
                b[r]
                // ex. (2,1), (3,2), etc..
            } else if r == c + 1 {
                symmetric_coeff[c]
                // ex. (1,2), (2,3), etc..
            } else if r == c - 1 {
                symmetric_coeff[r]
            } else {
                0.
            }
        })
    }

    // ゼロ点を求める
    // 返り値はVec<f64>
    fn zeros(&self) -> Vec<f64> {
        let sym_j = self.symmetric_jacobi();
        let eigen = SymmetricEigen::new(sym_j);
        eigen.eigenvalues.data.into()
    }
}
```

## 例: ラゲール多項式

次に，具体的な直交多項式に対するゼロ点を計算する構造体を実装する．
今回は例として，[（一般化された）ラゲール多項式](https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials)のゼロ点を計算する．

(一般化された)ラゲール多項式 $L_n^{(\alpha)}(x)$ は、以下の三項間漸化式を持つ．

$$
-(n+1) L_{n+1}(x) + (2n+ \alpha + 1) L_{n}(x) - (n + \alpha) L_{n-1}(x) = x L_{n}(x)
$$

今回は，`Laguerre` 構造体を定義し，この `Laguerre` 構造体に対して上記で実装した `Zeros` trait を実装する．

```rust
// Laguerre 構造体
// n は自由度
// alpha は一般化されたラゲール多項式の係数 alpha > -1
pub struct Laguerre {
    n: usize,
    alpha: f64,
}

impl Laguerre {
    // Laguerre 構造体を生成する関連関数
    pub fn new(n: usize, alpha: f64) -> Self {
        if alpha <= -1. {
            panic!("can't define Laguerre polynomials for alpha <= -1")
        }
        Self { n, alpha }
    }
}

// Laguerre 構造体に対するZeros traitの実装
impl Zeros for Laguerre {
    // 自由度 N を取得する
    fn degree(&self) -> usize {
        self.n
    }
    fn coefficients(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut a = Vec::with_capacity(self.n);
        let mut b = Vec::with_capacity(self.n);
        let mut c = Vec::with_capacity(self.n);

        for i in 0..self.n {
            // a_i = -(i + 1)
            a.push(-(i as f64 + 1.));
            // b_i = 2i + alpha + 1
            b.push(2. * i as f64 + self.alpha + 1.);
            // c_i = -(i + alpha)
            c.push(-(i as f64 + self.alpha));
        }

        (a, b, c)
    }
}
```

これで，ラゲール多項式のゼロ点を求めることができる．

```rust
#[cfg(test)]
mod tests {
    use super::{Laguerre, Zeros};

    #[test]
    fn test_laguerre() {
        let l = Laguerre::new(2, 0.);
        let zeros = l.zeros();
        let zeros_rounding: Vec<f64> = zeros
            .iter()
            .map(|z| (z * 10000.).floor() / 10000.)
            .collect();
        let want = vec![3.4142, 0.5857];

        println!("{:?}", zeros);
        assert_eq!(zeros_rounding, want);
    }

```

```shell
running 1 test
[3.414213562373095, 0.5857864376269049]
test special::tests::test_laguerre ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 9 filtered out; finished in 0.00s
```

# 参考文献

1. A. Gil, J. Segura, and, N.M. Temme, Numerical Methods for Special Functions. Society for Industrial and Applied Mathematics, 2007.
