# Returns all the node connected directly or indirectly to a set of nodes

Returns all the node connected directly or indirectly to a set of nodes

## Usage

``` r
get_connected_nodes(adj_list, nodes)
```

## Arguments

- adj_list:

  The network represented as an adjacency list

- nodes:

  A set of nodes

## Value

A vector of nodes indexes that are connected together with the ones
provided in the `nodes` argument. The `nodes` themselves are not listed
in this output
