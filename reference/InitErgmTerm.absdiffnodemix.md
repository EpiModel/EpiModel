# Definition for absdiffnodemix ERGM Term

This function defines and initializes the absdiffnodemix ERGM term that
allows for targeting homophily based on a non-binary attribute (e.g.,
age) by combinations of a binary attribute (e.g., race).

## Usage

``` r
InitErgmTerm.absdiffnodemix(nw, arglist, ...)
```

## Arguments

- nw:

  An object of class `network`.

- arglist:

  A list of arguments as specified in the `ergm.userterms` package
  framework.

- ...:

  Additional data passed into the function as specified in the
  `ergm.userterms` package framework.

## Details

This ERGM user term was written to allow for age-based homophily in
partnership formation that is heterogeneous by race. The `absdiff`
component targets the distribution of age mixing on that continuous
variable, and the `nodemix` component differentiates this for
black-black, black-white, and white-white couples.
