# Remove Edges That Include Specified Vertices

Given a current two-column matrix of edges and a vector of vertex IDs,
this function removes any rows of the edgelist in which the IDs are
present.

## Usage

``` r
delete_edges(el, vid)
```

## Arguments

- el:

  A two-column matrix of current edges (edgelist).

- vid:

  A vector of vertex IDs whose edges are to be deleted from the
  edgelist.

## Value

Returns an updated edgelist object, with any edges including the
specified vertices removed.
