# Convert an object to a `cumulative_edgelist`

Convert an object to a `cumulative_edgelist`

## Usage

``` r
as_cumulative_edgelist(x)
```

## Arguments

- x:

  An object to be converted to a cumulative edgelist

## Value

A `cumulative_edgelist` object, a `data.frame` with at least the
following columns: `head`, `tail`, `start`, `stop`.

## Details

The edges are active from time `start` to time `stop` included. If stop
is `NA`, the edge was not disolved in the simulation that generated the
list.
