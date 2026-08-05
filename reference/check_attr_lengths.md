# Check that All Attributes in the Main `netsim_dat` Object are of Equal Length

Every nodal attribute must hold exactly one value per node, the number
of nodes being the length of the `active` attribute. This is checked at
the end of each time step of
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).
An attribute of the wrong length is usually created by an
[`append_attr()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
call that does not cover every new node, or by direct assignment into
`dat$run$attr` in place of the accessors.

## Usage

``` r
check_attr_lengths(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

invisible(TRUE) if everything is correct; an error if not.
