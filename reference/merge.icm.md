# Merge Data across Stochastic Individual Contact Model Simulations

Merges epidemiological data from two independent simulations of
stochastic individual contact models from
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
# S3 method for class 'icm'
merge(x, y, ...)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md).

- y:

  Another `EpiModel` object of class
  [`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md), with
  the identical model parameterization as `x`.

- ...:

  Additional merge arguments (not used).

## Value

An `EpiModel` object of class
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md)
containing the data from both `x` and `y`.

## Details

This merge function combines the results of two independent simulations
of [`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md) class
models, simulated under separate function calls. The model
parameterization between the two calls must be exactly the same, except
for the number of simulations in each call. This allows for manual
parallelization of model simulations.

This merge function does not work the same as the default merge, which
allows for a combined object where the structure differs between the
input elements. Instead, the function checks that objects are identical
in model parameterization in every respect (except number of
simulations) and binds the results.

## Examples

``` r
param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
init <- init.icm(s.num = 1000, i.num = 100)
control <- control.icm(
  type = "SI", nsteps = 10,
  nsims = 3, verbose = FALSE
)
x <- icm(param, init, control)

control <- control.icm(
  type = "SI", nsteps = 10,
  nsims = 1, verbose = FALSE
)
y <- icm(param, init, control)

z <- merge(x, y)

# Examine separate and merged data
as.data.frame(x)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   986   114 1100      14
#> 3    1    3   968   132 1100      18
#> 4    1    4   956   144 1100      12
#> 5    1    5   937   163 1100      19
#> 6    1    6   912   188 1100      25
#> 7    1    7   887   213 1100      25
#> 8    1    8   859   241 1100      28
#> 9    1    9   829   271 1100      30
#> 10   1   10   788   312 1100      41
#> 11   2    1  1000   100 1100       0
#> 12   2    2   989   111 1100      11
#> 13   2    3   978   122 1100      11
#> 14   2    4   959   141 1100      19
#> 15   2    5   931   169 1100      28
#> 16   2    6   905   195 1100      26
#> 17   2    7   876   224 1100      29
#> 18   2    8   848   252 1100      28
#> 19   2    9   828   272 1100      20
#> 20   2   10   798   302 1100      30
#> 21   3    1  1000   100 1100       0
#> 22   3    2   992   108 1100       8
#> 23   3    3   980   120 1100      12
#> 24   3    4   960   140 1100      20
#> 25   3    5   944   156 1100      16
#> 26   3    6   927   173 1100      17
#> 27   3    7   903   197 1100      24
#> 28   3    8   876   224 1100      27
#> 29   3    9   851   249 1100      25
#> 30   3   10   823   277 1100      28
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   986   114 1100      14
#> 3    1    3   964   136 1100      22
#> 4    1    4   940   160 1100      24
#> 5    1    5   917   183 1100      23
#> 6    1    6   890   210 1100      27
#> 7    1    7   868   232 1100      22
#> 8    1    8   840   260 1100      28
#> 9    1    9   803   297 1100      37
#> 10   1   10   765   335 1100      38
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   986   114 1100      14
#> 3    1    3   968   132 1100      18
#> 4    1    4   956   144 1100      12
#> 5    1    5   937   163 1100      19
#> 6    1    6   912   188 1100      25
#> 7    1    7   887   213 1100      25
#> 8    1    8   859   241 1100      28
#> 9    1    9   829   271 1100      30
#> 10   1   10   788   312 1100      41
#> 11   2    1  1000   100 1100       0
#> 12   2    2   989   111 1100      11
#> 13   2    3   978   122 1100      11
#> 14   2    4   959   141 1100      19
#> 15   2    5   931   169 1100      28
#> 16   2    6   905   195 1100      26
#> 17   2    7   876   224 1100      29
#> 18   2    8   848   252 1100      28
#> 19   2    9   828   272 1100      20
#> 20   2   10   798   302 1100      30
#> 21   3    1  1000   100 1100       0
#> 22   3    2   992   108 1100       8
#> 23   3    3   980   120 1100      12
#> 24   3    4   960   140 1100      20
#> 25   3    5   944   156 1100      16
#> 26   3    6   927   173 1100      17
#> 27   3    7   903   197 1100      24
#> 28   3    8   876   224 1100      27
#> 29   3    9   851   249 1100      25
#> 30   3   10   823   277 1100      28
#> 31   4    1  1000   100 1100       0
#> 32   4    2   986   114 1100      14
#> 33   4    3   964   136 1100      22
#> 34   4    4   940   160 1100      24
#> 35   4    5   917   183 1100      23
#> 36   4    6   890   210 1100      27
#> 37   4    7   868   232 1100      22
#> 38   4    8   840   260 1100      28
#> 39   4    9   803   297 1100      37
#> 40   4   10   765   335 1100      38
```
