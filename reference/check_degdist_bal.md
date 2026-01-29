# Check Degree Distribution for Balance in Target Statistics

Checks for consistency in the implied network statistics of a two-group
network in which the group size and group-specific degree distributions
are specified.

## Usage

``` r
check_degdist_bal(num.g1, num.g2, deg.dist.g1, deg.dist.g2)
```

## Arguments

- num.g1:

  Number of nodes in group 1.

- num.g2:

  Number of nodes in group 2.

- deg.dist.g1:

  Vector with fractional degree distribution for group 1.

- deg.dist.g2:

  Vector with fractional degree distribution for group 2.

## Details

This function outputs the number of nodes of degree 0 to g, where g is
the length of a fractional degree distribution vector, given that vector
and the size of the group. This utility is used to check for balance in
implied degree given that fractional distribution within two-group
network simulations, in which the degree-constrained counts must be
equal across groups.

## Examples

``` r
# An unbalanced distribution
check_degdist_bal(num.g1 = 500, num.g2 = 500,
                  deg.dist.g2 = c(0.40, 0.55, 0.03, 0.02),
                  deg.dist.g1 = c(0.48, 0.41, 0.08, 0.03))
#> Degree Distribution Check
#> =============================================
#>         g1.dist   g1.cnt   g2.dist   g2.cnt
#> Deg0       0.48      240      0.40      200
#> Deg1       0.41      205      0.55      275
#> Deg2       0.08       40      0.03       15
#> Deg3       0.03       15      0.02       10
#> Edges      1.00      330      1.00      335
#> =============================================
#> ** Group 1 Edges < Group 2 Edges: -0.015 Rel Diff 

# A balanced distribution
check_degdist_bal(num.g1 = 500, num.g2 = 500,
                  deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
                  deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03))
#> Degree Distribution Check
#> =============================================
#>         g1.dist   g1.cnt   g2.dist   g2.cnt
#> Deg0       0.40      200      0.48      240
#> Deg1       0.55      275      0.41      205
#> Deg2       0.04       20      0.08       40
#> Deg3       0.01        5      0.03       15
#> Edges      1.00      330      1.00      330
#> =============================================
#> ** Edges balanced ** 
```
