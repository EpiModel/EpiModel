# Get Discordant Edgelist Based on Specified Status Variable

This function returns a `data.frame` with a discordant edgelist, defined
as the set of edges for which the status attribute of interest is
discordant between the two partners.

## Usage

``` r
get_discordant_edgelist(
  dat,
  status.attr,
  head.status,
  tail.status,
  networks = NULL
)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- status.attr:

  The name of the status attribute of interest.

- head.status:

  The value(s) of `status.attr` for which to look for the head of the
  edge. Can be a single value or a vector.

- tail.status:

  The value(s) of `status.attr` for which to look for the tail of the
  edge. Can be a single value or a vector.

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multiple overlapping networks.) If
  `NULL`, extract from all networks.

## Value

A `data.frame` with the following columns:

- `head`: Positional ID of the head node.

- `tail`: Positional ID of the tail node.

- `head_status`: Status of the head node.

- `tail_status`: Status of the tail node.

- `network`: The numerical index of the network on which the partnership
  is located.

## Details

This is a generalized version of the `discord_edgelist` function. It
creates an edgelist of current partnerships in which the status
attribute of interest (as specified by the parameter `status.attr`) of
one partner matches the value (or one of the values) of the
`head.status` parameter while the corresponding status attribute of the
other partner matches the value (or one of the values) of the
`tail.status` parameter.

## See also

[`discord_edgelist()`](http://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
