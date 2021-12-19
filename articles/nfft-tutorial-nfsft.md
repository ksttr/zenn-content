---
title: "不等間隔ノードに対する高速球面上フーリエ変換（NFSFT）"
emoji: "🕌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["c","cpp","fft","nfft"]
published: true
---


# 球面上フーリエ変換

球面上フーリエ変換とは，球面調和関数による展開のことです:

$$
f(\theta, \varphi) = \sum_{k,n} \hat{f}^n_k Y^n_k(\theta,\varphi)
$$

ここで，$f$は単位球面上で定義された関数で，$Y^n_k$は球面調和関数です．
球面調和関数の詳細は，[Wikipedia](https://ja.wikipedia.org/wiki/%E7%90%83%E9%9D%A2%E8%AA%BF%E5%92%8C%E9%96%A2%E6%95%B0)等を参照してください．

今回は，この展開の数値計算を実行する[C言語ライブラリNFFT](https://www-user.tu-chemnitz.de/~potts/nfft/)について説明します．

NFFTでは，任意に与えられた，$M$ 個のノード，$\boldsymbol{x}_0, \dots, \boldsymbol{x}_{M-1}$，に対して，$N$ 次までの球面調和関数展開を計算します．つまり，

$$
f_j = f(\theta_j, \varphi_j) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}^n_k Y^n_k(\theta_j,\varphi_j), \quad (j = 0,\dots,M-1)
$$

を計算します．ベクトルと行列で表すと，


$$
\boldsymbol{f} = \boldsymbol{Y \hat{f}}
$$

$$
\begin{align*}
\boldsymbol{f} &= (f_j)_0^{M-1} \in \mathbb{C}^{M}, \\
\boldsymbol{Y} &= (Y_n^k(\theta_j, \varphi_j))_{j=0,\dots,M-1; k=0,\dots,N, -k \leq n \leq k} \in \mathbb{C}^{M \times(N+1)^2}\\
\boldsymbol{\hat{f}} &= (\hat{f}_n^k)_{k=0,\dots,N, -k \leq n \leq k} \in \mathbb{C}^{(N+1)^2}
\end{align*}
$$

となって，結局，巨大な疎行列 $\boldsymbol{Y}$ と巨大なベクトル $\boldsymbol{\hat{f}}$ の積を計算する問題に帰着します．これは，通常の離散フーリエ変換と全く同じ問題ということになります．

今回は，アルゴリズムの詳細には立ち入らず，[C言語ライブラリNFFT](https://www-user.tu-chemnitz.de/~potts/nfft/)による計算方法だけ説明します．


# 計算量

計算量は，

$$
O(N^2\log^2N + M)
$$

です．

# サンプルコード

実際のサンプルを示しながら説明します．

今回は，

$$
\hat{f}_k^n =
\begin{cases}
\mp \sqrt{\frac{2\pi}{3}} & \mathrm{for} \quad (k,n) = (1,\pm 1) \\
0 & \mathrm{otherwise}
\end{cases}
$$

に対して，ランダムにノード $\boldsymbol{x_0}, \dots, \boldsymbol{x}_{M-1}$ を与えて計算していきます．計算結果は，

$$
f_j = \sin\theta_j \cos\varphi_j \, \, (= x_j)
$$

となります．

以下が，サンプルコードです．

```c
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <nfft3.h>

int main(void)
{

    const int N = 2;
    const int M = 4;
    nfsft_plan plan;

    nfsft_precompute(N, 1000.0, 0U, 0U);

    nfsft_init(&plan, N, M);

    for (int i = 0; i < plan.M_total; i++)
    {
        plan.x[2 * i] = nfft_drand48() - 0.5;
        plan.x[2 * i + 1] = nfft_drand48() * 0.5;
    }

    nfsft_precompute_x(&plan);

    for (int k = 0; k < plan.N; k++)
    {
        for (int n = -k; n <= k; n++)
        {
            if (k == 1 && (n == 1 || n == -1))
            {
                plan.f_hat[NFSFT_INDEX(k, n, &plan)] = - sqrt(2*M_PI/3)*k;
            }
            else
            {

                plan.f_hat[NFSFT_INDEX(k, n, &plan)] = 0;
            }
        }
    }

    nfsft_trafo(&plan);

    nfsft_finalize(&plan);
    nfsft_forget();

    return 0;
}
```

## 宣言

```c
    const int N = 2;
    const int M = 4;
    nfsft_plan plan;
```

- NFSFTに必要な変数を宣言します．
- `N` はフーリエ級数の個数
- `M` はノード点数
- `nfsft_plan`はNFSFTに必要な変数をメンバに持つ構造体

## 予備計算1

```c
    nfsft_precompute(N, 1000.0, 0U, 0U);
```

- NFSFTを計算するために，事前にグローバル変数を計算する関数です．
- `N`は以下で計算する変換における，自由度の最大値です．
- `1000`は，内部で実行される多項式変換(FTP)の上限です．
- 2つの`0U`はNFSFTとFTPのフラッグです．詳細は分かりません．

よくわからない場合は，上記の設定にしておくと良いようです．


## 初期化

```c
    nfsft_init(&plan, N, M);
```

- 初期化です．
- 内部で必要な変数のメモリーを動的確保しています．


## ノードの設定

```c
    for (int i = 0; i < plan.M_total; i++)
    {
        plan.x[2 * i] = nfft_drand48() - 0.5;
        plan.x[2 * i + 1] = nfft_drand48() * 0.5;
    }
```

- ノードを$\boldsymbol{x}_j$与えます．
- $\boldsymbol{x}_j = (\theta_j, \varphi_j)$ を直接与えるのではなく，$(\tilde{\theta_j},\tilde{\varphi}_j) \in ([0,1/2] \times [-1/2,1/2)$ を与えます．つまり，

$$
\tilde{\varphi} = 
\begin{cases} 
\varphi/2\pi, & \mathrm{if} \quad 0 \leq \varphi < \pi \\
\varphi/2\pi -1, & \mathrm{if} \quad \pi \leq \varphi < 2\pi \\
\end{cases}
, \quad\quad

\tilde{\theta} = \theta/2\pi
$$
を与えます．


## 予備計算2

```c
    nfsft_precompute_x(&plan);
```

- 不等間隔フーリエ変換に必要な予備計算を行う関数です．
- これに関する詳細は，[NFFTの記事](https://zenn.dev/ksttr/articles/nfft-tutorial)をご覧ください．

## フーリエ波数

```c
    for (int k = 0; k < plan.N; k++)
    {
        for (int n = -k; n <= k; n++)
        {
            if (k == 1 && (n == 1 || n == -1))
            {
                plan.f_hat[NFSFT_INDEX(k, n, &plan)] = - sqrt(2*M_PI/3)*k;
            }
            else
            {

                plan.f_hat[NFSFT_INDEX(k, n, &plan)] = 0;
            }
        }
    }
```

- フーリエ波数を与えます．
- 内部では，配列の順番が異なるようなので，補助的なマクロ`NFSFT_INDEX`使いましょう．

## 球面フーリエ変換


```c
    nfsft_trafo(&plan);
```

- ここまで来れば変換は簡単です．


## 終了処理

```c
    nfsft_finalize(&plan);
    nfsft_forget();
```

- 終了処理には，２段階あります．
- 1つ目の`nfsft_finalize`で`plan`で確保したメモリーを破棄します．
- 2つ目の`nfsft_forget`で`nfsft_precompute`で確保したメモリーを破棄します．

以上です．