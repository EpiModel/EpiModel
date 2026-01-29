# Create a TEA Variable for Infection Status for `ndtv` Animations

Creates a new color-named temporally-extended attribute (TEA) variable
in a `networkDynamic` object containing a disease status TEA in numeric
format.

## Usage

``` r
color_tea(
  nd,
  old.var = "testatus",
  old.sus = "s",
  old.inf = "i",
  old.rec = "r",
  new.var = "ndtvcol",
  new.sus = NULL,
  new.inf = NULL,
  new.rec = NULL,
  verbose = TRUE
)
```

## Arguments

- nd:

  An object of class `networkDynamic`.

- old.var:

  Old TEA variable name.

- old.sus:

  Status value for susceptible in old TEA variable.

- old.inf:

  Status value for infected in old TEA variable.

- old.rec:

  Status value for recovered in old TEA variable.

- new.var:

  New TEA variable name to be stored in `networkDynamic` object.

- new.sus:

  Status value for susceptible in new TEA variable.

- new.inf:

  Status value for infected in new TEA variable.

- new.rec:

  Status value for recovered in new TEA variable.

- verbose:

  If `TRUE`, print progress to console.

## Value

The updated object of class `networkDynamic`.

## Details

The `ndtv` package (<https://cran.r-project.org/package=ndtv>) produces
animated visuals for dynamic networks with evolving edge structures and
nodal attributes. Nodal attribute dynamics in `ndtv` movies require a
temporally extended attribute (TEA) containing a standard R color for
each node at each time step. By default, the `EpiModel` package uses
TEAs to store disease status history in network model simulations run in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).
But that status TEA is in numeric format (0, 1, 2). The `color_tea`
function transforms those numeric values of that disease status TEA into
a TEA with color values in order to visualize status changes in `ndtv`.

The convention in
[`plot.netsim()`](http://epimodel.github.io/EpiModel/reference/plot.netsim.md)
is to color the susceptible nodes as blue, infected nodes as red, and
recovered nodes as green. Alternate colors may be specified using the
`new.sus`, `new.inf`, and `new.rec` parameters, respectively.

Using the `color_tea` function with a `netsim` object requires that TEAs
for disease status be used and that the `networkDynamic` object be saved
in the output: `tergmListe` must be set to `FALSE` in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) and
the `ndtv` package documentation.
