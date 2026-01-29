# Definition for absdiffby ERGM Term

This function defines and initializes the absdiffby ERGM term that
allows for representing homophily with respect to a non-binary attribute
(e.g., age) differentially by a binary attribute (e.g., sex).

## Usage

``` r
InitErgmTerm.absdiffby(nw, arglist, ...)
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
partnership formation that is asymmetric by sex. The `absdiff` component
targets age-based homophily while the `by` component allows that to be
structured by a binary attribute such as "male", in order to enforce an
offset in the average difference. This allows, for example, a average
age difference in partnerships, but with males (on average) older than
females.
