# Primary Infection Module for netsim

This function simulates the main infection process given the current
state of the partnerships and disease in the system.

## Usage

``` r
infection.net(dat, at)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- at:

  Current time step.

## Value

The updated `netsim_dat` main list object.

## Details

The main steps in this infection module are as follows:

1.  Get IDs for current infected and susceptible nodes given the current
    disease status.

2.  Call
    [`discord_edgelist()`](http://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
    to get the current discordant edgelist given step 1.

3.  Determine the transmission rates (e.g., as a function of group).

4.  Pull the number of acts per partnership in a time step from the
    `act.rate` parameter.

5.  Calculate the final transmission probabilities given the
    transmission rates and act rates.

6.  Randomly transmit on the discordant edgelist.

7.  Conduct bookkeeping for new infections to update status on the nodes
    and calculate disease incidence.

## See also

[`discord_edgelist()`](http://epimodel.github.io/EpiModel/reference/discord_edgelist.md)
is used within `infection.net` to obtain a discordant edgelist.
