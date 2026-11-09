# Shift Vertex IDs Down Past Deleted IDs

Subtracts from each vertex ID in `orig` the number of deleted IDs in
`del` that are less than or equal to it, so that the IDs left after a
deletion are consecutive again. Both vectors must be sorted in
increasing order.

## Usage

``` r
shiftVec(orig, del)
```

## Arguments

- orig:

  Sorted vector of vertex IDs.

- del:

  Sorted vector of deleted vertex IDs.

## Value

An integer vector the length of `orig`.

## Details

For example, `orig = 1:10` and `del = c(4, 6, 8, 9)` return
`c(1, 2, 3, 3, 4, 4, 5, 5, 5, 6)`, so the IDs kept (1, 2, 3, 5, 7, 10)
become 1 to 6.
