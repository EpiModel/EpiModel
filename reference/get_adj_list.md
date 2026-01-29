# Returns an adjacency list from an edge list

Returns an adjacency list from an edge list

## Usage

``` r
get_adj_list(el, n_nodes)
```

## Arguments

- el:

  An edge list as a data.frame with columns `head` and `tail`

- n_nodes:

  The size number of node in the network

## Value

An adjacency list for the network

## Details

The adjacency list is a `list` of length `n_nodes`. The entry for each
node is a integer vector containing the index of all the nodes connected
to it. This layout makes it directly subsetable in O(1) at the expanse
of memory usage. To get all connections to the nodes 10 and 15 :
`unlist(adj_list[c(10, 15)]`
