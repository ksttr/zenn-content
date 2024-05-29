---
title: "[Fortran] 疎行列をCOOフォーマットからCSRフォーマットに変換する"
emoji: "💾"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["Fortran", "mkl"]
published: false
---

Intel MKL の inspector-excutor Sparse Blas は疎行列の作成・解析・行列演算などを行うことができるライブラリである．

Inspector-excutor Sparse Blas には任意の疎行列フォーマットを CSR へ変換する `mkl_sparse_convert_csr` が用意されているが，
これを使って COO から CSR へ変換した疎行列を OpenMP で並列化しようとしたところうまくできなかった．
一方で，直接作成した CSR は期待通り並列化された．
そこで，COO を CSR へ変換するサブルーチンを自作することにした．

# 疎行列フォーマット

はじめに簡単に疎行列のフォーマットについて説明しよう．
[Intel の解説ページ](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2024-1/sparse-matrix-storage-formats.html#MKL_APPA_SMSF_3)を全面的に参考にしている．

## COO

[COO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2024-1/sparse-blas-coordinate-matrix-storage-format.html) は Coordinate Matrix Strage の略で，
0 でない要素の行番号(`rows`)・列番号(`columns`)及びその値(`values`)の 3 つの変数で疎行列を表現する方法である．

例えば，

$$
C =
\begin{pmatrix}
0 & 2 & -5 & 0 & 0 \\
1 & 0 & 0  & 4 & 0 \\
0 & 0 & 0  & 0 & 2 \\
0 & -3 & 8  & 0 & 0 \\
5 & 7 & 0  & 3 & 0 \\
\end{pmatrix}
$$

という行列を COO で表現すると，

```python
values  = (2, -5, 1, 4, 2, -3, 8, 5, 7, 3)
rows    = (1,  1, 2, 2, 3,  4, 4, 5, 5, 5)
columns = (2,  3, 1, 4, 5,  2, 3, 1, 2, 4)
```

となる．(今回は Fortran で実装するのでインデックスはすべて**1 始まり**である．)

COO は作成が簡単であるのがメリットであるので，
例えば筆者のようなズボラな人間はとりあえず COO で作成してその後で CSR へ変換しようというふうに使う．

しかし，Intel MKL ライブラリでは COO は OpenMP による並列化をサポートしていない (そもそも並列化に不向きなフォーマットである)ので，並列化が必要なときは次の CSR を用いる必要がある．

## CSR

CSR は Compressed Sparse Row の略で，COO の行の表現 (`rows`) を更に圧縮する．
COO では列番号の配列 `columns` と行番号の配列 `rows` の両方を指定して非ゼロ要素の場所を特定したが，
列番号の配列があれば，この配列の要素のどこで行が切り替わるかさえわかれば非ゼロ要素の場所は特定できる．

COO の例と同じ行列

$$
C =
\begin{pmatrix}
0 & 2 & -5 & 0 & 0 \\
1 & 0 & 0  & 4 & 0 \\
0 & 0 & 0  & 0 & 2 \\
0 & -3 & 8  & 0 & 0 \\
5 & 7 & 0  & 3 & 0 \\
\end{pmatrix}
$$

で具体的に考えよう．
この場合，行が切り替わるときの列番号 `columns`のインデックスは `3, 5, 6, 8` である (しつこいようだがインデックスは 1 から始まっている)．
CSR では，これに 1 行目の最初の要素のインデックス`1`を加えた配列を使って，行列 $C$ を

```python
values   = (2, -5, 1, 4, 2, -3, 8, 5, 7, 3)
columns  = (2,  3, 1, 4, 5,  2, 3, 1, 2, 4)
pointerB = (1,  3, 5, 6, 8)
```

と表現する．
`pointerB`は各行の先頭へのポインタとみなせる．
したがって，この表現では行へのアクセスが高速である．
例えば，2 行目にアクセスするには`pointerB[2]`番目から`pointerB[3]`までの `columns`のスライスを取得すれば良い．

Intel MKL では `pointerB`に加えて行末尾の次の要素へのポインタ `pointerE`も用意する．
すなわち，`pointer[j]-1`は`j`行目の末尾へのポインタでなる．

```python
pointerE = (3, 5, 6, 8, 11)
```

# 変換方法

上記の CSR の作成方法がそのまま COO を CSR へ変換する方法になっているので，
ここまでくれば実装は簡単である．

1. COO の行番号配列 `rows` の値が変わるインデックスを記録し配列 `p` に格納する
2. 配列 `p` の先頭に `1`，末尾に `rows`の要素数に 1 加えた値を追加する．
3. `pointerB = p[1:size(p)-1]`
4. `pointerE = p[2:size(p)]`

# 実装

Fortran での実装例を示す．
コード内の変数名などは [Intel のリファレンス](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2024-1/inspector-executor-sparse-blas-routines.html)と同じものを使っている．

```fortran
subroutine compress_coo_rows(rows, row_indx, rows_start, rows_end)
    integer, intent(in) :: rows ! number of rows of matrix
    integer, dimension(:), intent(in) :: row_indx
    integer, dimension(:), intent(out) :: rows_start, rows_end

    integer :: p(rows + 1)

    integer :: rid
    integer :: pid
    integer :: k, i

    rid = 1
    pid = 1

    p(pid) = 1
    pid = pid + 1

    do k = 1, size(row_indx)
        if (rid < row_indx(k)) then
            p(pid) = k
            rid = row_indx(k)
            pid = pid + 1
        end if
    end do

    do i = pid, rows + 1
        p(i) = size(row_indx) + 1
    end do

    rows_start = p(1:size(p) - 1)
    rows_end = p(2:size(p))

end subroutine compress_coo_rows
```
