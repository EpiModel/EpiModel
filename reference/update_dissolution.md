# Adjust Dissolution Component of Network Model Fit

Adjusts the dissolution component of a dynamic ERGM fit using the
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)
function with the edges dissolution approximation method.

## Usage

``` r
update_dissolution(old.netest, new.coef.diss, nested.edapprox = TRUE)
```

## Arguments

- old.netest:

  An object of class `netest`, from the
  [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)
  function.

- new.coef.diss:

  An object of class `disscoef`, from the
  [`dissolution_coefs()`](http://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
  function.

- nested.edapprox:

  Logical. If `edapprox = TRUE` the dissolution model is an initial
  segment of the formation model (see details in
  [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)).

## Value

An updated network model object of class `netest`.

## Details

Fitting an ERGM is a computationally intensive process when the model
includes dyad dependent terms. With the edges dissolution approximation
method of Carnegie et al, the coefficients for a temporal ERGM are
approximated by fitting a static ERGM and adjusting the formation
coefficients to account for edge dissolution. This function provides a
very efficient method to adjust the coefficients of that model when one
wants to use a different dissolution model; a typical use case may be to
fit several different models with different average edge durations as
targets. The example below exhibits that case.

## Examples

``` r
if (FALSE) { # \dontrun{
nw <- network_initialize(n = 1000)

# Two dissolutions: an average duration of 300 versus 200
diss.300 <- dissolution_coefs(~offset(edges), 300, 0.001)
diss.200 <- dissolution_coefs(~offset(edges), 200, 0.001)

# Fit the two reference models
est300 <- netest(nw = nw,
                formation = ~edges,
                target.stats = c(500),
                coef.diss = diss.300)

est200 <- netest(nw = nw,
                formation = ~edges,
                target.stats = c(500),
                coef.diss = diss.200)

# Alternatively, update the 300 model with the 200 coefficients
est200.compare <- update_dissolution(est300, diss.200)

identical(est200$coef.form, est200.compare$coef.form)
} # }
```
