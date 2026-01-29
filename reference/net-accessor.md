# Functions to Access and Edit the Main netsim_dat Object in Network Models

These `get_`, `set_`, `append_`, and `add_` functions allow a safe and
efficient way to retrieve and mutate the main `netsim_dat` class object
of network models (typical variable name `dat`).

## Usage

``` r
get_attr_list(dat, item = NULL)

get_attr(dat, item, posit_ids = NULL, override.null.error = FALSE)

add_attr(dat, item)

set_attr(dat, item, value, posit_ids = NULL, override.length.check = FALSE)

append_attr(dat, item, value, n.new)

remove_node_attr(dat, posit_ids)

get_epi_list(dat, item = NULL)

get_epi(dat, item, at = NULL, override.null.error = FALSE)

add_epi(dat, item)

set_epi(dat, item, at, value)

get_param_list(dat, item = NULL)

get_param(dat, item, override.null.error = FALSE)

add_param(dat, item)

set_param(dat, item, value)

get_control_list(dat, item = NULL)

get_control(dat, item, override.null.error = FALSE)

get_network_control(dat, network, item, override.null.error = FALSE)

add_control(dat, item)

set_control(dat, item, value)

get_init_list(dat, item = NULL)

get_init(dat, item, override.null.error = FALSE)

add_init(dat, item)

set_init(dat, item, value)

append_core_attr(dat, at, n.new)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  A character vector containing the name of the element to access (for
  `get_` functions), create (for `add_` functions), or edit (for `set_`
  and `append_` functions). Can be of length \> 1 for `get_*_list`
  functions.

- posit_ids:

  For `set_attr` and `get_attr`, a numeric vector of posit_ids to subset
  the desired `item`.

- override.null.error:

  If TRUE, `get_` will return NULL if the `item` does not exist instead
  of throwing an error. (default = FALSE).

- value:

  New value to be attributed in the `set_` and `append_` functions.

- override.length.check:

  If TRUE, `set_attr` allows the modification of the `item` size.
  (default = FALSE).

- n.new:

  For `append_core_attr`, the number of new nodes to initiate with core
  attributes; for `append_attr`, the number of new elements to append at
  the end of `item`.

- at:

  For `get_epi`, the timestep at which to access the specified `item`;
  for `set_epi`, the timestep at which to add the new value for the epi
  output `item`; for `append_core_attr`, the current time step.

- network:

  index of network for which to get control

## Value

A vector or a list of vectors for `get_` functions; the main list object
for `set_`, `append_`, and `add_` functions.

## Core Attribute

The `append_core_attr` function initializes the attributes necessary for
EpiModel to work (the four core attributes are: "active", "unique_id",
"entrTime", and "exitTime"). These attributes are used in the
initialization phase of the simulation, to create the nodes (see
[`initialize.net()`](http://epimodel.github.io/EpiModel/reference/initialize.net.md));
and also used when adding nodes during the simulation (see
[`arrivals.net()`](http://epimodel.github.io/EpiModel/reference/arrivals.net.md)).

## Mutability

The `set_`, `append_`, and `add_` functions DO NOT modify the
`netsim_dat` object in place. The result must be assigned back to `dat`
in order to be registered: `dat <- set_*(dat, item, value)`.

## `set_` and `append_` vs `add_`

The `set_` and `append_` functions edit a pre-existing element or create
a new one if it does not exist already by calling the `add_` functions
internally.

## Examples

``` r
dat <- create_dat_object(control = list(nsteps = 150))
dat <- append_core_attr(dat, 1, 100)

dat <- add_attr(dat, "age")
dat <- set_attr(dat, "age", runif(100))
dat <- set_attr(dat, "status", rbinom(100, 1, 0.9))
dat <- append_attr(dat, "status", 1, 10)
dat <- append_attr(dat, "age", NA, 10)
get_attr_list(dat)
#> $active
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $entrTime
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $exitTime
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> $unique_id
#>   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#>  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
#>  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
#>  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
#>  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
#>  [91]  91  92  93  94  95  96  97  98  99 100
#> 
#> $age
#>   [1] 0.472155346 0.827517861 0.997468345 0.171409699 0.553224359 0.235689790
#>   [7] 0.641282232 0.592536345 0.950813536 0.769726436 0.632953239 0.691117570
#>  [13] 0.719409832 0.856726790 0.448163730 0.252906574 0.027796780 0.490980056
#>  [19] 0.943569213 0.785268735 0.637415222 0.677191104 0.540372382 0.178387259
#>  [25] 0.798645704 0.800875106 0.626758761 0.193208164 0.888248340 0.287432256
#>  [31] 0.862127952 0.914751872 0.234589809 0.779132877 0.113940740 0.561120051
#>  [37] 0.403765588 0.998114243 0.525125574 0.077033859 0.248949234 0.252873791
#>  [43] 0.062847466 0.312056552 0.441575141 0.543781811 0.502625705 0.634648598
#>  [49] 0.191642583 0.599501718 0.629196057 0.014272752 0.255305569 0.571998335
#>  [55] 0.050432306 0.393724435 0.316647325 0.861556549 0.749570003 0.629912602
#>  [61] 0.515326896 0.967560790 0.694990798 0.966807820 0.022100390 0.148553824
#>  [67] 0.243629732 0.261067565 0.775994096 0.083854648 0.599293890 0.080073340
#>  [73] 0.135773981 0.578292555 0.003942028 0.657847548 0.372409806 0.042153687
#>  [79] 0.885667579 0.447681102 0.931397314 0.763073731 0.704600006 0.448510038
#>  [85] 0.633974422 0.765366134 0.032040004 0.826327592 0.854571128 0.078676174
#>  [91] 0.062692564 0.816444727 0.439895075 0.281943856 0.265352782 0.131272081
#>  [97] 0.297060508 0.568750000 0.540987363 0.946943816          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $status
#>   [1] 1 1 1 0 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr_list(dat, c("age", "active"))
#> $age
#>   [1] 0.472155346 0.827517861 0.997468345 0.171409699 0.553224359 0.235689790
#>   [7] 0.641282232 0.592536345 0.950813536 0.769726436 0.632953239 0.691117570
#>  [13] 0.719409832 0.856726790 0.448163730 0.252906574 0.027796780 0.490980056
#>  [19] 0.943569213 0.785268735 0.637415222 0.677191104 0.540372382 0.178387259
#>  [25] 0.798645704 0.800875106 0.626758761 0.193208164 0.888248340 0.287432256
#>  [31] 0.862127952 0.914751872 0.234589809 0.779132877 0.113940740 0.561120051
#>  [37] 0.403765588 0.998114243 0.525125574 0.077033859 0.248949234 0.252873791
#>  [43] 0.062847466 0.312056552 0.441575141 0.543781811 0.502625705 0.634648598
#>  [49] 0.191642583 0.599501718 0.629196057 0.014272752 0.255305569 0.571998335
#>  [55] 0.050432306 0.393724435 0.316647325 0.861556549 0.749570003 0.629912602
#>  [61] 0.515326896 0.967560790 0.694990798 0.966807820 0.022100390 0.148553824
#>  [67] 0.243629732 0.261067565 0.775994096 0.083854648 0.599293890 0.080073340
#>  [73] 0.135773981 0.578292555 0.003942028 0.657847548 0.372409806 0.042153687
#>  [79] 0.885667579 0.447681102 0.931397314 0.763073731 0.704600006 0.448510038
#>  [85] 0.633974422 0.765366134 0.032040004 0.826327592 0.854571128 0.078676174
#>  [91] 0.062692564 0.816444727 0.439895075 0.281943856 0.265352782 0.131272081
#>  [97] 0.297060508 0.568750000 0.540987363 0.946943816          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $active
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr(dat, "status")
#>   [1] 1 1 1 0 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1
get_attr(dat, "status", c(1, 4))
#> [1] 1 0

dat <- add_epi(dat, "i.num")
dat <- set_epi(dat, "i.num", 150, 10)
dat <- set_epi(dat, "s.num", 150, 90)
get_epi_list(dat)
#> $i.num
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 10
#> 
#> $s.num
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 90
#> 
get_epi_list(dat, c("i.num", "s.num"))
#> $i.num
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 10
#> 
#> $s.num
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 90
#> 
get_epi(dat, "i.num")
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 10
get_epi(dat, "i.num", c(1, 4))
#> [1] NA NA

dat <- add_param(dat, "x")
dat <- set_param(dat, "x", 0.4)
dat <- set_param(dat, "y", 0.8)
get_param_list(dat)
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_param_list(dat, c("x", "y"))
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_param(dat, "x")
#> [1] 0.4

dat <- add_init(dat, "x")
dat <- set_init(dat, "x", 0.4)
dat <- set_init(dat, "y", 0.8)
get_init_list(dat)
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_init_list(dat, c("x", "y"))
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_init(dat, "x")
#> [1] 0.4

dat <- add_control(dat, "x")
dat <- set_control(dat, "x", 0.4)
dat <- set_control(dat, "y", 0.8)
get_control_list(dat)
#> $nsteps
#> [1] 150
#> 
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_control_list(dat, c("x", "y"))
#> $x
#> [1] 0.4
#> 
#> $y
#> [1] 0.8
#> 
get_control(dat, "x")
#> [1] 0.4
```
