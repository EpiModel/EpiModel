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
#> 2    1    2   982   118 1100      18
#> 3    1    3   961   139 1100      21
#> 4    1    4   939   161 1100      22
#> 5    1    5   916   184 1100      23
#> 6    1    6   899   201 1100      17
#> 7    1    7   866   234 1100      33
#> 8    1    8   832   268 1100      34
#> 9    1    9   803   297 1100      29
#> 10   1   10   771   329 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   986   114 1100      14
#> 13   2    3   969   131 1100      17
#> 14   2    4   951   149 1100      18
#> 15   2    5   934   166 1100      17
#> 16   2    6   915   185 1100      19
#> 17   2    7   898   202 1100      17
#> 18   2    8   873   227 1100      25
#> 19   2    9   849   251 1100      24
#> 20   2   10   816   284 1100      33
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   967   133 1100      17
#> 24   3    4   951   149 1100      16
#> 25   3    5   928   172 1100      23
#> 26   3    6   910   190 1100      18
#> 27   3    7   884   216 1100      26
#> 28   3    8   856   244 1100      28
#> 29   3    9   828   272 1100      28
#> 30   3   10   797   303 1100      31
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   982   118 1100      18
#> 3    1    3   966   134 1100      16
#> 4    1    4   950   150 1100      16
#> 5    1    5   927   173 1100      23
#> 6    1    6   893   207 1100      34
#> 7    1    7   859   241 1100      34
#> 8    1    8   827   273 1100      32
#> 9    1    9   795   305 1100      32
#> 10   1   10   765   335 1100      30
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   982   118 1100      18
#> 3    1    3   961   139 1100      21
#> 4    1    4   939   161 1100      22
#> 5    1    5   916   184 1100      23
#> 6    1    6   899   201 1100      17
#> 7    1    7   866   234 1100      33
#> 8    1    8   832   268 1100      34
#> 9    1    9   803   297 1100      29
#> 10   1   10   771   329 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   986   114 1100      14
#> 13   2    3   969   131 1100      17
#> 14   2    4   951   149 1100      18
#> 15   2    5   934   166 1100      17
#> 16   2    6   915   185 1100      19
#> 17   2    7   898   202 1100      17
#> 18   2    8   873   227 1100      25
#> 19   2    9   849   251 1100      24
#> 20   2   10   816   284 1100      33
#> 21   3    1  1000   100 1100       0
#> 22   3    2   984   116 1100      16
#> 23   3    3   967   133 1100      17
#> 24   3    4   951   149 1100      16
#> 25   3    5   928   172 1100      23
#> 26   3    6   910   190 1100      18
#> 27   3    7   884   216 1100      26
#> 28   3    8   856   244 1100      28
#> 29   3    9   828   272 1100      28
#> 30   3   10   797   303 1100      31
#> 31   4    1  1000   100 1100       0
#> 32   4    2   982   118 1100      18
#> 33   4    3   966   134 1100      16
#> 34   4    4   950   150 1100      16
#> 35   4    5   927   173 1100      23
#> 36   4    6   893   207 1100      34
#> 37   4    7   859   241 1100      34
#> 38   4    8   827   273 1100      32
#> 39   4    9   795   305 1100      32
#> 40   4   10   765   335 1100      30
```
