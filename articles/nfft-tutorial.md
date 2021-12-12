---
title: "NFFTライブラリの使い方"
emoji: "🕌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["c","cpp", "fft", "nfft"]
published: false
---


# はじめに

この記事では，[C言語ライブラリNFFT](https://www-user.tu-chemnitz.de/~potts/nfft/index.php)の使用方法について説明します．

# NFFTとは

まずは，とても簡単にNFFTについて説明します． 

離散フーリエ変換

$$
f_j = f(\boldsymbol{x}_j) = \sum_{\boldsymbol{k}} \hat{f}_{\boldsymbol{k}} e^{-2\pi i \boldsymbol{k}\cdot\boldsymbol{x}_j}
$$

のノード$\boldsymbol{x}_j$が不等間隔である場合を，不等間隔離散フーリエ変換（NDFT）と言います．

フーリエ係数$\hat{f}_{\boldsymbol{k}}$が与えれた時に，元の関数$f_j$を高速に求めるアルゴリズムを不等間隔高速フーリエ変換（NFFT）と言います．

行列で書くと，

$$
\boldsymbol{f} = \boldsymbol{A} \boldsymbol{\hat{f}},
$$

$$
\boldsymbol{f} = (f_j)_{j=0,\dots, M-1}, \quad \boldsymbol{\hat{f}} = (f_{\boldsymbol{k}})_{\boldsymbol{k}}, \quad \boldsymbol{A} = \left(e^{-2\pi i \boldsymbol{k}\cdot\boldsymbol{x}_j } \right)_{j=0,\dots, M-1; \boldsymbol{k}}, 
$$

となり，NFFTとは，巨大な疎行列$\boldsymbol{A}$とベクトル$\boldsymbol{\hat{f}}$の積を効率よく求めるアルゴリズム，ということになります．


通常のFFTは，ノード$\boldsymbol{x}_0, \dots, \boldsymbol{x}_{M-1}$が等間隔に並んでいる場合にのみ適用できますが，今回のNFFTでは，不等間隔に並んでいる場合でも実行できるアルゴリズムです．

# 計算量

1次元NFFTの計算量は，

$$
\mathcal{O}\left(N\log N + \log \left(1/\varepsilon\right) M \right)
$$

です．　ここで，　$N$はフーリエ波数$\boldsymbol{k}$の個数，$M$はノード$\boldsymbol{x}_j$の個数，$\varepsilon$は精度を表します．

$d$次元の場合は，これの$d$乗です．

# NFFTライブラリ

NFFTのC言語ライブラリの使い方について説明します．

## インストール
[ここ](https://www-user.tu-chemnitz.de/~potts/nfft/installation.php)の手順通りにインストールしましょう．

**(注意) 事前にFFTWをインストールしておく必要があります．**

1. ダウンロードし，解凍後，ディレクトリに入る．

```shell 
wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.5.2.tar.gz
tar -zxf nfft-3.5.2.tar.gz
cd nfft-3.5.2
```

2. コンパイルする．　

```shell
./configure --enable-all --enable-openmp && make
```

3. インストール

```shell
sudo make install
```

## 使い方

簡単な例を示します．
今回は， コサイン

$$
f(x) = \cos(2\pi x),  \quad x \in \left[-\frac{\pi}{2}, \frac{\pi}{2}\right]
$$

を，そのフーリエ級数

$$
\begin{align}
\hat{f}(k) = 
\begin{cases}
 1/2 & \mathrm{for} \quad k = \pm 1 \\
 0 & \mathrm{otherwise}
\end{cases}
\end{align}
$$

から求める，コードを示します．

ノード$x$はランダムに用意します．

```c
#include <complex.h>
#include <nfft3.h>

int main(void)
{
    nfft_plan my_plan;
    int N = 4;　// 波数の個数
    int M = 32;　// ノードxの個数

    // planの初期化
    nfft_init_1d(&my_plan, N, M);

    // ランダムにxの値を決めてくれる関数
    nfft_vrand_shifted_unit_double(my_plan.x, my_plan.M_total);

　　　　　　　　// 窓関数の計算
    if (my_plan.flags & PRE_ONE_PSI)
        nfft_precompute_one_psi(&my_plan);
　　　　　　　　
　　　　　　　　// フーリエ級数
    // 今回はcos(2*pi*x)のフーリエ変換を行う
    for (int i = 0; i < my_plan.N_total; i++)
    {
        my_plan.f_hat[i] = 0;
    }
    my_plan.f_hat[my_plan.N_total / 2 - 1] = 0.5;　 // 波数の入れ方に注意
    my_plan.f_hat[my_plan.N_total / 2 + 1] = 0.5; // 波数の入れ方に注意

　　　　　　　　// NFFTの実行
    nfft_trafo(&my_plan);

　　　　　　　　// planの終了
    nfft_finalize(&my_plan);

    return 0;
}
```

1つずつ説明します．

### 宣言

まずは，必要な変数を宣言します．

```c
    nfft_plan my_plan;
    int N = 4;　// 波数の個数
    int M = 32;　// ノードxの個数
```

- `nfft_plan`は，FFTWのときと同じようなplanの宣言です． このplan変数に，NFFTの実行に必要なものをすべて詰め込んでいきます．
- `N`は波数の個数です．　`N=4`なので，$k=-2, -1, 0, 1$となります．
- `M`はノード$x$の個数です．　ノード$x$は以下で設定します．

### 初期化

プランを初期化します．

```c
    // planの初期化
    nfft_init_1d(&my_plan, N, M);
```

プランの初期化をしています．

構造体nfft_plan型のmy_planのメンバには，


| my_plan.N_total | ノード$x$の個数       |
|     ----        | ---- |
| my_plan.M_total | フーリエ波数$k$の個数  |
| my_plan.x       | ノード$x$を表すポインタ |
| my_plan.f_hat   | フーリエ級数$\hat{f}$を表すポインタ |
| my_plan.f       | 求めたい関数$f$を表すポインタ|


などがいます．


## ノードを決める

ノード$x$を決めます．
これは，$[-1/2, 1/2]$の範囲であればなんでも良いです（**等間隔である必要もないし，小さい方から順番に格納する必要もない！！**）．

```c
    // ランダムにxの値を決めてくれる関数
    nfft_vrand_shifted_unit_double(my_plan.x, my_plan.M_total);
```

今回は，NFFTライブラリが提供する，ランダムに格納してくれる関数`nfft_vrand_shifted_unit_double`を使います．

### 窓関数の計算

今回理論の説明で端折ってしまいましたが，NFFTはフーリエ変換を窓関数の和で近似します．
そのため，窓関数を決め，あらかじめ計算するというステップが入ります．

```c
　　　　　　　　// 窓関数の計算
    if (my_plan.flags & PRE_ONE_PSI)
        nfft_precompute_one_psi(&my_plan);
```

今回は，`nfft_precompute_one_psi`という関数を使って，計算します．
窓関数の詳細な説明や設定は，記事を分けようと思います．

### フーリエ級数の計算

フーリエ級数の計算をします．

```c
　　　　　　　　// フーリエ級数
    // 今回はcos(2*pi*x)のフーリエ変換を行う
    for (int i = 0; i < my_plan.N_total; i++)
    {
        my_plan.f_hat[i] = 0;
    }
    my_plan.f_hat[my_plan.N_total / 2 - 1] = 0.5;　 // 波数の入れ方に注意
    my_plan.f_hat[my_plan.N_total / 2 + 1] = 0.5; // 波数の入れ方に注意
```

今回はコサインなので簡単です．

### NFFTの実行

NFFTを実行しましょう．

```c
　　　　　　　　// NFFTの実行
    nfft_trafo(&my_plan);
```

- `nfft_trafo`でNFFTを実行します．
- `my_plan.f`に実行結果が格納されます．

また，符号が逆の変換，

$$
\hat{f}_k = \sum_{j=0}^{M-1} f_j e^{+2\pi i kx_j}
$$

を求めるときは， `nfft_adjoint`を使います．

### 終了処理

最後にプランを破棄して終了です．

```c
　　　　　　　　// planの終了
    nfft_finalize(&my_plan);
```



# おわりに

今回は，NFFTを実行するC言語ライブラリNFFTの説明をしました．
簡単な使用例にとどめましたが，これ以外にも，

- 逆フーリエ変換
- 2次元以上の高次元フーリエ変換
- 高速サイン・コサイン変換（NFCT,NFST)
- 球面上の変換（球面調和関数展開）NSFST

など，様々な機能があります．これらについても今後記事にしていこうと思います．

