# Summary for Network Model Fit

Prints the summary model fit statistics for an ERGM or STERGM fit.

## Usage

``` r
# S3 method for class 'netest'
summary(object, ...)
```

## Arguments

- object:

  An `EpiModel` object of class `netest`.

- ...:

  Additional summary function arguments.

## Details

This function is simply a wrapper function for `summary.ergm`.
Additionally, if the edges dissolution approximation was used to fit the
temporal ERGM, then the dissolution coefficient information will be
printed.

If the `fit` object is attached to the `netest` object, then
`summary.netest` will call `summary` on `fit` using the `...` passed to
`summary.netest`. Otherwise, `summary.netest` will print the stored
summary of the fit generated in the original `netest` call, using the
`...` passed to `netest`.
