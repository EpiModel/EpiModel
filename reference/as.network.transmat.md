# Convert transmat Infection Tree into a network Object

Converts a transmission matrix from the `get_transmat` function into a
[`network::network`](https://rdrr.io/pkg/network/man/network.html) class
object.

## Usage

``` r
# S3 method for class 'transmat'
as.network(x, ...)
```

## Arguments

- x:

  An object of class `transmat` to be converted into a `network` class
  object.

- ...:

  Unused.

## Value

A [`network::network`](https://rdrr.io/pkg/network/man/network.html)
object.

## Details

When converting from a `transmat` to a `network` object, this functions
copies the edge attributes within the transmission matrix (`'at'`,
`'infDur'`, `'transProb'`, `'actRate'`, and `'finalProb'`) into edge
attributes on the network.
