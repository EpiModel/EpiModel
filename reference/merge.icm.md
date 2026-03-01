# Merge Data across Stochastic Individual Contact Model Simulations

Merges epidemiological data from two independent simulations of
stochastic individual contact models from
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

## Usage

``` r
# S3 method for class 'icm'
merge(x, y, ...)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md).

- y:

  Another `EpiModel` object of class
  [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md), with
  the identical model parameterization as `x`.

- ...:

  Additional merge arguments (not used).

## Value

An `EpiModel` object of class
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md)
containing the data from both `x` and `y`.

## Details

This merge function combines the results of two independent simulations
of [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md) class
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
#> 3    1    3   978   122 1100      14
#> 4    1    4   963   137 1100      15
#> 5    1    5   947   153 1100      16
#> 6    1    6   934   166 1100      13
#> 7    1    7   913   187 1100      21
#> 8    1    8   887   213 1100      26
#> 9    1    9   853   247 1100      34
#> 10   1   10   821   279 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   983   117 1100      17
#> 13   2    3   963   137 1100      20
#> 14   2    4   948   152 1100      15
#> 15   2    5   931   169 1100      17
#> 16   2    6   906   194 1100      25
#> 17   2    7   880   220 1100      26
#> 18   2    8   851   249 1100      29
#> 19   2    9   818   282 1100      33
#> 20   2   10   781   319 1100      37
#> 21   3    1  1000   100 1100       0
#> 22   3    2   983   117 1100      17
#> 23   3    3   965   135 1100      18
#> 24   3    4   946   154 1100      19
#> 25   3    5   931   169 1100      15
#> 26   3    6   908   192 1100      23
#> 27   3    7   884   216 1100      24
#> 28   3    8   844   256 1100      40
#> 29   3    9   812   288 1100      32
#> 30   3   10   785   315 1100      27
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   991   109 1100       9
#> 3    1    3   980   120 1100      11
#> 4    1    4   969   131 1100      11
#> 5    1    5   946   154 1100      23
#> 6    1    6   922   178 1100      24
#> 7    1    7   892   208 1100      30
#> 8    1    8   873   227 1100      19
#> 9    1    9   845   255 1100      28
#> 10   1   10   815   285 1100      30
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   992   108 1100       8
#> 3    1    3   978   122 1100      14
#> 4    1    4   963   137 1100      15
#> 5    1    5   947   153 1100      16
#> 6    1    6   934   166 1100      13
#> 7    1    7   913   187 1100      21
#> 8    1    8   887   213 1100      26
#> 9    1    9   853   247 1100      34
#> 10   1   10   821   279 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   983   117 1100      17
#> 13   2    3   963   137 1100      20
#> 14   2    4   948   152 1100      15
#> 15   2    5   931   169 1100      17
#> 16   2    6   906   194 1100      25
#> 17   2    7   880   220 1100      26
#> 18   2    8   851   249 1100      29
#> 19   2    9   818   282 1100      33
#> 20   2   10   781   319 1100      37
#> 21   3    1  1000   100 1100       0
#> 22   3    2   983   117 1100      17
#> 23   3    3   965   135 1100      18
#> 24   3    4   946   154 1100      19
#> 25   3    5   931   169 1100      15
#> 26   3    6   908   192 1100      23
#> 27   3    7   884   216 1100      24
#> 28   3    8   844   256 1100      40
#> 29   3    9   812   288 1100      32
#> 30   3   10   785   315 1100      27
#> 31   4    1  1000   100 1100       0
#> 32   4    2   991   109 1100       9
#> 33   4    3   980   120 1100      11
#> 34   4    4   969   131 1100      11
#> 35   4    5   946   154 1100      23
#> 36   4    6   922   178 1100      24
#> 37   4    7   892   208 1100      30
#> 38   4    8   873   227 1100      19
#> 39   4    9   845   255 1100      28
#> 40   4   10   815   285 1100      30
```
