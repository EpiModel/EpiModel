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
#> 2    1    2   992   108 1100       8
#> 3    1    3   976   124 1100      16
#> 4    1    4   963   137 1100      13
#> 5    1    5   946   154 1100      17
#> 6    1    6   919   181 1100      27
#> 7    1    7   894   206 1100      25
#> 8    1    8   864   236 1100      30
#> 9    1    9   834   266 1100      30
#> 10   1   10   804   296 1100      30
#> 11   2    1  1000   100 1100       0
#> 12   2    2   983   117 1100      17
#> 13   2    3   966   134 1100      17
#> 14   2    4   947   153 1100      19
#> 15   2    5   930   170 1100      17
#> 16   2    6   908   192 1100      22
#> 17   2    7   881   219 1100      27
#> 18   2    8   852   248 1100      29
#> 19   2    9   819   281 1100      33
#> 20   2   10   780   320 1100      39
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   967   133 1100      17
#> 24   3    4   951   149 1100      16
#> 25   3    5   933   167 1100      18
#> 26   3    6   917   183 1100      16
#> 27   3    7   889   211 1100      28
#> 28   3    8   861   239 1100      28
#> 29   3    9   832   268 1100      29
#> 30   3   10   805   295 1100      27
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   985   115 1100      15
#> 3    1    3   960   140 1100      25
#> 4    1    4   941   159 1100      19
#> 5    1    5   921   179 1100      20
#> 6    1    6   895   205 1100      26
#> 7    1    7   868   232 1100      27
#> 8    1    8   839   261 1100      29
#> 9    1    9   814   286 1100      25
#> 10   1   10   788   312 1100      26
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   992   108 1100       8
#> 3    1    3   976   124 1100      16
#> 4    1    4   963   137 1100      13
#> 5    1    5   946   154 1100      17
#> 6    1    6   919   181 1100      27
#> 7    1    7   894   206 1100      25
#> 8    1    8   864   236 1100      30
#> 9    1    9   834   266 1100      30
#> 10   1   10   804   296 1100      30
#> 11   2    1  1000   100 1100       0
#> 12   2    2   983   117 1100      17
#> 13   2    3   966   134 1100      17
#> 14   2    4   947   153 1100      19
#> 15   2    5   930   170 1100      17
#> 16   2    6   908   192 1100      22
#> 17   2    7   881   219 1100      27
#> 18   2    8   852   248 1100      29
#> 19   2    9   819   281 1100      33
#> 20   2   10   780   320 1100      39
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   967   133 1100      17
#> 24   3    4   951   149 1100      16
#> 25   3    5   933   167 1100      18
#> 26   3    6   917   183 1100      16
#> 27   3    7   889   211 1100      28
#> 28   3    8   861   239 1100      28
#> 29   3    9   832   268 1100      29
#> 30   3   10   805   295 1100      27
#> 31   4    1  1000   100 1100       0
#> 32   4    2   985   115 1100      15
#> 33   4    3   960   140 1100      25
#> 34   4    4   941   159 1100      19
#> 35   4    5   921   179 1100      20
#> 36   4    6   895   205 1100      26
#> 37   4    7   868   232 1100      27
#> 38   4    8   839   261 1100      29
#> 39   4    9   814   286 1100      25
#> 40   4   10   788   312 1100      26
```
