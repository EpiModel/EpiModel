# Return the Historical Contacts (Partners) of a Set of Index Nodes

From a full cumulative edgelist that contains the history of contacts
(both persistent and one-time), this function returns a data frame
containing details of the index (head) and partner (tail) nodes, along
with start and stop time steps for the partnership and the network
location.

## Usage

``` r
get_partners(
  dat,
  index_posit_ids,
  networks = NULL,
  truncate = Inf,
  only.active.nodes = FALSE
)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- index_posit_ids:

  The positional IDs of the indexes of interest.

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multi-layer networks.) If `NULL`, extract
  from all networks.

- truncate:

  After how many time steps a partnership that is no longer active
  should be removed from the output.

- only.active.nodes:

  If `TRUE`, then inactive (e.g., deceased) partners will be removed
  from the output.

## Value

A `data.frame` with 5 columns:

- `index`: the unique IDs of the indexes.

- `partner`: the unique IDs of the partners/contacts.

- `start`: the time step at which the edge started.

- `stop`: the time step in which the edge stopped; if ongoing, then `NA`
  is returned.

- `network`: the numerical index for the network on which the
  partnership/contact is located.

## Details

Note that `get_partners` takes as input the positional IDs of the
indexes of interest but returns the unique IDs. That is by design,
because while `get_partners` would be expected to be called for active
nodes, some partners (contacts) of nodes may be inactive in the network
history. Therefore, both index and partner IDs are returned as unique
IDs for consistency. To convert between a positional to a unique ID, you
may use
[`get_posit_ids`](http://epimodel.github.io/EpiModel/reference/unique_id-tools.md);
to convert between a unique ID to a positional ID, you may use
[`get_unique_ids`](http://epimodel.github.io/EpiModel/reference/unique_id-tools.md).
