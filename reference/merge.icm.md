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
#> 2    1    2   993   107 1100       7
#> 3    1    3   979   121 1100      14
#> 4    1    4   964   136 1100      15
#> 5    1    5   948   152 1100      16
#> 6    1    6   936   164 1100      12
#> 7    1    7   914   186 1100      22
#> 8    1    8   890   210 1100      24
#> 9    1    9   856   244 1100      34
#> 10   1   10   824   276 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   982   118 1100      18
#> 13   2    3   962   138 1100      20
#> 14   2    4   947   153 1100      15
#> 15   2    5   928   172 1100      19
#> 16   2    6   902   198 1100      26
#> 17   2    7   876   224 1100      26
#> 18   2    8   840   260 1100      36
#> 19   2    9   804   296 1100      36
#> 20   2   10   766   334 1100      38
#> 21   3    1  1000   100 1100       0
#> 22   3    2   985   115 1100      15
#> 23   3    3   960   140 1100      25
#> 24   3    4   943   157 1100      17
#> 25   3    5   919   181 1100      24
#> 26   3    6   893   207 1100      26
#> 27   3    7   857   243 1100      36
#> 28   3    8   830   270 1100      27
#> 29   3    9   808   292 1100      22
#> 30   3   10   773   327 1100      35
as.data.frame(y)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   987   113 1100      13
#> 3    1    3   969   131 1100      18
#> 4    1    4   944   156 1100      25
#> 5    1    5   925   175 1100      19
#> 6    1    6   901   199 1100      24
#> 7    1    7   875   225 1100      26
#> 8    1    8   848   252 1100      27
#> 9    1    9   818   282 1100      30
#> 10   1   10   784   316 1100      34
as.data.frame(z)
#>    sim time s.num i.num  num si.flow
#> 1    1    1  1000   100 1100       0
#> 2    1    2   993   107 1100       7
#> 3    1    3   979   121 1100      14
#> 4    1    4   964   136 1100      15
#> 5    1    5   948   152 1100      16
#> 6    1    6   936   164 1100      12
#> 7    1    7   914   186 1100      22
#> 8    1    8   890   210 1100      24
#> 9    1    9   856   244 1100      34
#> 10   1   10   824   276 1100      32
#> 11   2    1  1000   100 1100       0
#> 12   2    2   982   118 1100      18
#> 13   2    3   962   138 1100      20
#> 14   2    4   947   153 1100      15
#> 15   2    5   928   172 1100      19
#> 16   2    6   902   198 1100      26
#> 17   2    7   876   224 1100      26
#> 18   2    8   840   260 1100      36
#> 19   2    9   804   296 1100      36
#> 20   2   10   766   334 1100      38
#> 21   3    1  1000   100 1100       0
#> 22   3    2   985   115 1100      15
#> 23   3    3   960   140 1100      25
#> 24   3    4   943   157 1100      17
#> 25   3    5   919   181 1100      24
#> 26   3    6   893   207 1100      26
#> 27   3    7   857   243 1100      36
#> 28   3    8   830   270 1100      27
#> 29   3    9   808   292 1100      22
#> 30   3   10   773   327 1100      35
#> 31   4    1  1000   100 1100       0
#> 32   4    2   987   113 1100      13
#> 33   4    3   969   131 1100      18
#> 34   4    4   944   156 1100      25
#> 35   4    5   925   175 1100      19
#> 36   4    6   901   199 1100      24
#> 37   4    7   875   225 1100      26
#> 38   4    8   848   252 1100      27
#> 39   4    9   818   282 1100      30
#> 40   4   10   784   316 1100      34
```
