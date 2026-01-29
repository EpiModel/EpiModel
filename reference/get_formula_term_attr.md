# Output ERGM Formula Attributes into a Character Vector

Given a formation formula for a network model, outputs a character
vector of vertex attributes to be used in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
simulations.

## Usage

``` r
get_formula_term_attr(form, nw)
```

## Arguments

- form:

  An ERGM model formula.

- nw:

  A network object.

## Value

A character vector of vertex attributes.
