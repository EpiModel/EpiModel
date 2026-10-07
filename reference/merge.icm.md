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
control <- control.icm(type = "SI", nsteps = 10,
                       nsims = 3, verbose = FALSE)
x <- icm(param, init, control)

control <- control.icm(type = "SI", nsteps = 10,
                       nsims = 1, verbose = FALSE)
y <- icm(param, init, control)

z <- merge(x, y)

# Examine separate and merged data
as.data.frame(x)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   985   115 1100      15
#> 3    1    3   969   131 1100      16
#> 4    1    4   944   156 1100      25
#> 5    1    5   921   179 1100      23
#> 6    1    6   897   203 1100      24
#> 7    1    7   874   226 1100      23
#> 8    1    8   853   247 1100      21
#> 9    1    9   822   278 1100      31
#> 10   1   10   788   312 1100      34
#> 11   2    1  1000   100 1100       0
#> 12   2    2   986   114 1100      14
#> 13   2    3   970   130 1100      16
#> 14   2    4   952   148 1100      18
#> 15   2    5   933   167 1100      19
#> 16   2    6   909   191 1100      24
#> 17   2    7   876   224 1100      33
#> 18   2    8   848   252 1100      28
#> 19   2    9   819   281 1100      29
#> 20   2   10   785   315 1100      34
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   971   129 1100      13
#> 24   3    4   960   140 1100      11
#> 25   3    5   941   159 1100      19
#> 26   3    6   918   182 1100      23
#> 27   3    7   895   205 1100      23
#> 28   3    8   867   233 1100      28
#> 29   3    9   844   256 1100      23
#> 30   3   10   817   283 1100      27
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   986   114 1100      14
#> 3    1    3   967   133 1100      19
#> 4    1    4   947   153 1100      20
#> 5    1    5   930   170 1100      17
#> 6    1    6   903   197 1100      27
#> 7    1    7   879   221 1100      24
#> 8    1    8   846   254 1100      33
#> 9    1    9   817   283 1100      29
#> 10   1   10   789   311 1100      28
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   985   115 1100      15
#> 3    1    3   969   131 1100      16
#> 4    1    4   944   156 1100      25
#> 5    1    5   921   179 1100      23
#> 6    1    6   897   203 1100      24
#> 7    1    7   874   226 1100      23
#> 8    1    8   853   247 1100      21
#> 9    1    9   822   278 1100      31
#> 10   1   10   788   312 1100      34
#> 11   2    1  1000   100 1100       0
#> 12   2    2   986   114 1100      14
#> 13   2    3   970   130 1100      16
#> 14   2    4   952   148 1100      18
#> 15   2    5   933   167 1100      19
#> 16   2    6   909   191 1100      24
#> 17   2    7   876   224 1100      33
#> 18   2    8   848   252 1100      28
#> 19   2    9   819   281 1100      29
#> 20   2   10   785   315 1100      34
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   971   129 1100      13
#> 24   3    4   960   140 1100      11
#> 25   3    5   941   159 1100      19
#> 26   3    6   918   182 1100      23
#> 27   3    7   895   205 1100      23
#> 28   3    8   867   233 1100      28
#> 29   3    9   844   256 1100      23
#> 30   3   10   817   283 1100      27
#> 31   4    1  1000   100 1100       0
#> 32   4    2   986   114 1100      14
#> 33   4    3   967   133 1100      19
#> 34   4    4   947   153 1100      20
#> 35   4    5   930   170 1100      17
#> 36   4    6   903   197 1100      27
#> 37   4    7   879   221 1100      24
#> 38   4    8   846   254 1100      33
#> 39   4    9   817   283 1100      29
#> 40   4   10   789   311 1100      28
```
