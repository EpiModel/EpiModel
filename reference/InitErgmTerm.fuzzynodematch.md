# Definition for fuzzynodematch ERGM Term

This function defines and initializes the fuzzynodematch ERGM term that
allows for generalized homophily.

## Usage

``` r
InitErgmTerm.fuzzynodematch(nw, arglist, ...)
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

This ERGM user term was written to allow for generalized homophily.The
`attr` term argument should specify a character vertex attribute
encoding the "venues" associated to each node. The `split` argument
should specify a string that separates different "venues" in the
attribute value for each node, as handled by `strsplit` with
`fixed = TRUE`. For example, if `split` is `"|"` (the default), and the
attribute value for a given node is `"a12|b476"`, then the associated
venues for this node are `"a12"` and `"b476"`. The empty string `""` is
interpreted as "no venues".

If the `binary` term argument is `FALSE` (the default), the change
statistic for an on-toggle is the number of unique venues associated to
both nodes (informally speaking, this could be described as the number
of venues on which the two nodes "match"); if `binary` is `TRUE`, the
change statistic for an on-toggle is `1` if any venue is associated to
both nodes, and `0` otherwise.
