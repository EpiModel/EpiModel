# Return an adjacency list of subnets

Return an adjacency list of subnets

## Usage

``` r
get_subnet_adj_list(adj_list)
```

## Arguments

- adj_list:

  A normal adjacency list

## Value

An adjacency list where only the first node of a subnet contains the
subnet and all other contain only the first node

## Details

A graph with 4 components: 1, 2, 3, 4, and 5 and 6, 7, 8 would yield a
list like so: 1: 2, 3, 4 2: 1 3: 1 4: 1 5: numeric(0) 6: 7, 8 7: 6, 8: 6

This format speeds up the construction of reachable sets on dense
networks
