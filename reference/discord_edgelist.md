# Discordant Edgelist

This function returns a `data.frame` with a discordant edgelist, defined
as the set of edges in which the status of the two partners is one
susceptible and one infected.

## Usage

``` r
discord_edgelist(dat, at, network = 1, infstat = "i", include.network = FALSE)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

- network:

  In case of models with multiple networks, the network to pull the
  current edgelist from. Default of `network = 1`.

- infstat:

  Character vector of disease status values that are considered
  infectious, defining the SI pairs.

- include.network:

  Should the `network` value be included as the final column of the
  discordant edgelist?

## Value

This function returns a `data.frame` with the following columns:

- **time:** time step queried.

- **sus:** ID number for the susceptible partner.

- **inf:** ID number for the infectious partner.

The output from this function is added to the transmission `data.frame`
object that is requested as output in `netsim` simulations with the
`save.trans=TRUE` argument.

## Details

This internal function works within the parent
[`infection.net()`](http://epimodel.github.io/EpiModel/reference/infection.net.md)
function to pull the current edgelist from the dynamic network object,
look up the disease status of the head and tails on the edge, and subset
the list to those edges with one susceptible and one infected node.

EpiModel v2.0.3 extended the function by allowing flexibility in the
definition what disease status counts as infectious, with the `infstat`
parameter. For extension models with multiple infectious states, this
can be a vector of length greater than 1: `infstat = c("i", "a")`.

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md),
[`infection.net()`](http://epimodel.github.io/EpiModel/reference/infection.net.md)
