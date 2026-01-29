# Grow a Vector to a Given Size, Padding it With Empty Elements

Grow a vector to a given size, padding it with `NULL` if `orig` is a
`list` and with `NA` otherwise

## Usage

``` r
padded_vector(orig, size)
```

## Arguments

- orig:

  A vector to grow.

- size:

  The final size of the vector.

## Value

A vector of size `size` padded with `NULL`s or `NA`s at the end.
