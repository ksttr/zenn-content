---
title: "NFFT逆変換"
emoji: "🕌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["c","cpp","fft", "nfft"]
published: false
---

# NFFTとは，

NFFTとは，離散フーリエ変換，

$$
f_j = f(\boldsymbol{x}_j) = \sum_{\boldsymbol{k}} \hat{f}_{\boldsymbol{k}} e^{-2\pi i \boldsymbol{k}\cdot\boldsymbol{x}_j}
$$

に対して，

- フーリエ級数 $\hat{f}_{\boldsymbol{k}}$ が$2N+1$個の波数$\boldsymbol{k}$で与えられている
- ノード$\bm{x}_j \, (j=0,\dots, M-1)$，が，**不等間隔**に与えれている

ときに，元の関数$f_j$を求めるアルゴリズムのことです．これを行列ベクトルで表すと，


$$
\boldsymbol{f} = \boldsymbol{A} \boldsymbol{\hat{f}},
$$

$$
\boldsymbol{f} = (f_j)_{j=0,\dots, M-1}, \quad \boldsymbol{\hat{f}} = (\hat{f}_{\boldsymbol{k}})_{\boldsymbol{k}}, \quad \boldsymbol{A} = \left(e^{-2\pi i \boldsymbol{k}\cdot\boldsymbol{x}_j } \right)_{j=0,\dots, M-1; \boldsymbol{k}}, 
$$

となります．

NFFTとそのC言語ライブラリである[NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/)の簡単な説明・使い方は以下の記事を参照してください，

- [NFFTライブラリの使い方](https://zenn.dev/ksttr/articles/nfft-tutorial)

# NFFT逆変換

今回の記事では，NFFTの逆変換（INFFT）の説明をします．

NFFTの逆変換（INFFT）とは，ノードとその上の値$(\boldsymbol{x}_j, y_j), \, (j=0,\dots,M-1)$が与えれた時に，フーリエ級数$\boldsymbol{\hat{f}}$を求めることです．つまり，以下の線形連立方程式の解$\boldsymbol{\hat{f}}$を求めることになります．

$$
\boldsymbol{A\hat{f}} = \boldsymbol{y}
$$

等間隔にノード$\boldsymbol{x}_j$が与れれている通常のFFTとは異なり，フーリエ逆変換が厳密に存在しているとは限りません．**フーリエ行列$\boldsymbol{A}$の随伴行列$\boldsymbol{A}^*$を左から掛けかけると求まるわけではありません．**

ノードの個数がフーリエ級数の個数個数より多い場合と少ない場合で，それぞれ計算法が異なります．

## ノードの個数が過剰な場合

ノードの個数が過剰な場合，フーリエ級数$\hat{f}$が元の関数$f$を表現できるほど十分ないということになります．

この場合は，最も良いものを最小二乗法で求めます．

実際には，

$$
\min_{\boldsymbol{\hat{f}}} \left(\sum_{j=0}^{M-1} w_j |y_j - f(\boldsymbol{x}_j)|^2\right)^{1/2}, \quad f(\boldsymbol{x}_j) = (\boldsymbol{A\hat{f}})_{j}
$$

を解いているようですが，この重み$w_j$がどのように決められているのかは，[ドキュメント](https://www-user.tu-chemnitz.de/~potts/nfft/guide/paper-nfft33.pdf)にも書いてありませんでした．判明したら追記するかもしれません．


## ノードの個数が不足な場合

ノードの個数が不足している場合は，フーリエ級数$\hat{f}$が余分に余っている状態です．つまり，適当な$\hat{f}$が複数個存在する場合があるということです．

この場合は，サンプル$y$の個数を内挿によって増やしてやる必要があります．
これはすなわち，複数個ある解$\hat{f}$の中から，最も良さそうなものを選ぶということです．

最も良いものを選ぶ基準は，

$$
\min_{\boldsymbol{\hat{f}}} \left(\sum_{\boldsymbol{k}}\frac{|\hat{f}_{\boldsymbol{k}|}}{\hat{w}_{\boldsymbol{k}}}
\right)^{1/2}, \quad\mathrm{subject \,\, to} \quad \boldsymbol{A\hat{f}} = \boldsymbol{y}
$$

と書いてありますが，重み$\hat{w}_{\boldsymbol{k}}$の詳細は書いてありませんでした．詳細がわかったら追記するかもしれません．

# サンプルコード

今回は，コサイン，

$$
\cos(x) = \cos(2\pi x), \quad x \in \left[-\frac{1}{2}, \frac{1}{2}\right]
$$

のフーリエ級数，

$$
\begin{align*}
\hat{f}(k) = 
\begin{cases}
1/2 & \mathrm{for} \quad k = \pm 1 \\
0 & \mathrm{otherwise}
\end{cases}
\end{align*}
$$

を求めるコードを示します．

**(注意) ドキュメントが3.3でとまっていて，かつ，逆変換のコードは3.4以降でかなり加筆・修正されているようです．ここで示すものも3.5以外では使えなくなっている可能性が高いです．**


```c
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <nfft3.h>

int main(void)
{
    nfft_plan plan;
    solver_plan_complex inverse_plan;
    int N = 4;
    int M = 16;
    int iter = 4;

    nfft_init_1d(&plan, N, M);

    solver_init_complex(&inverse_plan, (nfft_mv_plan_complex *)(&plan));

    nfft_vrand_shifted_unit_double(plan.x, plan.M_total);

    for (int i = 0; i < plan.M_total; i++)
    {
        inverse_plan.y[i] = cos(2 * M_PI * plan.x[i]);
    }

    // initial guess for f_hat
    for (int k = 0; k < plan.N_total; k++)
    {
        inverse_plan.f_hat_iter[k] = 0.0;
    }

    solver_before_loop_complex(&inverse_plan);

    for (int l = 0; l < iter; l++)
    {
        printf("----- %d iteration -----\n", l + 1);

        solver_loop_one_step_complex(&inverse_plan);

        for (int k = 0; k < plan.N_total; k++)
        {
            printf("%e %e\n", creal(inverse_plan.f_hat_iter[k]), cimag(inverse_plan.f_hat_iter[k]));
        }

        printf("\n Residual r=%e\n-----------------------\n", inverse_plan.dot_r_iter);
    }

    solver_finalize_complex(&inverse_plan);
    nfft_finalize(&plan);

    return 0;
}
```


この結果は以下のようになります，

```shell
----- 1 iteration -----
2.967704e-02 1.008778e-01
4.712111e-01 4.425935e-02
4.601135e-03 0.000000e+00
4.712111e-01 -4.425935e-02

 Residual r=2.235322e-01
-----------------------
----- 2 iteration -----
5.654262e-03 8.838023e-03
4.935071e-01 1.039638e-02
-2.747496e-02 -5.693536e-03
5.010603e-01 -1.039638e-02

 Residual r=1.301747e-02
-----------------------
----- 3 iteration -----
3.159839e-04 -1.878513e-03
4.988973e-01 2.774002e-03
1.477890e-04 -4.532437e-04
5.009242e-01 -2.774002e-03

 Residual r=3.316033e-04
-----------------------
----- 4 iteration -----
-1.651681e-17 1.323644e-18
5.000000e-01 -8.873240e-18
-2.119126e-17 -5.392640e-19
5.000000e-01 7.110964e-18

 Residual r=2.474691e-32
-----------------------
```


コサインだと4回の反復で十分な精度のフーリエ係数が求まります．

以下はコードの詳細な説明です．


## 初期化

```c
    nfft_plan plan;
    solver_plan_complex inverse_plan;
    int N = 4;
    int M = 16;
    int iter = 4;

    nfft_init_1d(&plan, N, M);

    solver_init_complex(&inverse_plan, (nfft_mv_plan_complex *)(&plan));
```

- 内部で，フーリエ変換と擬逆変換(上記説明の逆変換解法)を繰り返すため，
逆変換用のプラン`solver_plan_complex`に加え，`nfft_plan`も必要です．
- 逆変換用プラン`inverse_plan`の初期化には，`solver_init_complex`を使います．引数として，逆変換用プラン`solver_plan`とNFFTプラン`nfft_plan`を渡しますが，型としては`nfft_mv_plan_complex`に直してから渡す必要があるようです．


## 値の代入

```c
    nfft_vrand_shifted_unit_double(plan.x, plan.M_total);

    for (int i = 0; i < plan.M_total; i++)
    {
        inverse_plan.y[i] = cos(2 * M_PI * plan.x[i]);
    }

    // initial guess for f_hat
    for (int k = 0; k < plan.N_total; k++)
    {
        inverse_plan.f_hat_iter[k] = 0.0;
    }
```

- ここでは，必要な値を代入します．
- `plan.x`は，ノード$x$の値を入れます．今回はランダムに入れています．
- `inverse_plan.y`はノード$x$に対応する関数$f$の値です．今回はコサインの値を入れています．
- `inverse_plan.f_hat_iter`はこれから求める$\hat{f}$の値の初期値を入れます．反復解法で求めるので，このイニシャルゲスが良ければ，良い値が見つかることになります．


## 反復前処理

```c
    solver_before_loop_complex(&inverse_plan);
```

- `solver_before_loop_complex`は反復前の処理を実行しています（多分）．詳細はよくわかりません．

## 反復解法


```c
    for (int l = 0; l < iter; l++)
    {
        solver_loop_one_step_complex(&inverse_plan);
        ...
    }
```

- `solver_loop_one_step_complex`は，実際に反復解法によって，逆変換を求める関数です．１ループごとにこの関数を呼ぶ必要があります．


## 終了処理

```c
    solver_finalize_complex(&inverse_plan);
    nfft_finalize(&plan);
```

- 2つのプランを終了します．


以上です．