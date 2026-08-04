# Functions to Access and Edit the Main netsim_dat Object in Network Models

These `get_`, `set_`, `append_`, and `add_` functions allow a safe and
efficient way to retrieve and mutate the main `netsim_dat` class object
of network models (typical variable name `dat`). They are intended for
use *inside module functions* that run during a
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
simulation, not for editing the user-facing `param.net`, `init.net`, or
`control.net` inputs before a run or the output object returned by
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).
See the **Intended Usage** section below for details and alternatives.

This function returns an exhaustive named list of the attributes managed
by EpiModel itself. It can be used to check the validity of an
attributes list and of its types.

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

get_core_attributes()

append_core_attr(dat, at, n.new)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- item:

  A character vector containing the name of the element to access (for
  `get_` functions), create (for `add_` functions), or edit (for `set_`
  and `append_` functions). Can be of length greater than 1 for
  `get_*_list` functions.

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

## Intended Usage

These accessors operate on the live `netsim_dat` object that EpiModel
constructs internally when a simulation starts and passes through each
module at every time step. The expected calling context is **inside a
user-supplied module function** (e.g., a custom `infection.FUN`,
`recovery.FUN`, or arrivals/departures module) registered with
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).

They are *not* a general-purpose tool for editing the user-facing input
or output objects of
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md):

- Calling them on a
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
  [`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md),
  or
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  object before a run will appear to succeed (those objects are also
  list-like) but produces an object that is not a valid simulation
  input.

- Calling them on the object returned by
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  modifies only the saved input slot, not anything that would change the
  result of a re-run. Note also that
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  strips the `*.FUN` entries from its returned `control` slot, so
  feeding that slot back into a new simulation will fail.

For modifications outside a running simulation, use the following
instead:

- **Editing parameters before a run:**
  [`update_params()`](https://epimodel.github.io/EpiModel/reference/update_params.md)
  for a `param.net` object, or direct list assignment (e.g.,
  `p$inf.prob <- 0.5`).

- **Editing init or control before a run:** direct list assignment
  (e.g., `ctrl$nsteps <- 1000`), or rebuild with a fresh call to the
  [`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md)
  /
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  constructor.

- **Scheduled mid-simulation changes:** the `.param.updater.list`
  argument to
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md)
  and the `.control.updater.list` argument to
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).
  See the "model-parameters" vignette for the underlying updater module
  and the scenario API built on top of it.

## Core Attribute

The `append_core_attr` function initializes the attributes necessary for
EpiModel to work (the four core attributes are: "active", "unique_id",
"entrTime", and "exitTime"). These attributes are used in the
initialization phase of the simulation, to create the nodes (see
[`initialize.net()`](https://epimodel.github.io/EpiModel/reference/initialize.net.md));
and also used when adding nodes during the simulation (see
[`arrivals.net()`](https://epimodel.github.io/EpiModel/reference/arrivals.net.md)).

## Mutability

The `set_`, `append_`, and `add_` functions DO NOT modify the
`netsim_dat` object in place. The result must be assigned back to `dat`
in order to be registered: `dat <- set_*(dat, item, value)`.

## `set_` and `append_` vs `add_`

The `set_` and `append_` functions edit a pre-existing element or create
a new one if it does not exist already by calling the `add_` functions
internally.

## See also

[`update_params()`](https://epimodel.github.io/EpiModel/reference/update_params.md)
for editing a `param.net` object outside a simulation;
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
[`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md),
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
for input constructors;
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
for running a simulation.

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
#>   [1] 0.981842521 0.480152659 0.497900751 0.969944620 0.730625570 0.871249147
#>   [7] 0.991377368 0.664240878 0.713891583 0.143768221 0.896165183 0.622561361
#>  [13] 0.699182787 0.242552443 0.909468295 0.601935858 0.954815769 0.189696660
#>  [19] 0.806798475 0.614939201 0.337632246 0.695872342 0.945397054 0.588769985
#>  [25] 0.202820597 0.361969674 0.604657195 0.927379446 0.720764208 0.851808469
#>  [31] 0.995168352 0.352305684 0.441294787 0.484173094 0.743150527 0.175656549
#>  [37] 0.581930148 0.474467288 0.115034923 0.774499406 0.453538384 0.556696330
#>  [43] 0.926556017 0.443646785 0.008501619 0.062822039 0.744874349 0.065412350
#>  [49] 0.630897832 0.961829548 0.996465983 0.668148058 0.835421227 0.343280025
#>  [55] 0.725047796 0.494380954 0.245206871 0.914296501 0.903821214 0.523425133
#>  [61] 0.652749121 0.663911420 0.893451170 0.595002735 0.608075619 0.724174157
#>  [67] 0.113207193 0.139429949 0.112331380 0.462294281 0.042387796 0.247349081
#>  [73] 0.857424026 0.850155975 0.206830324 0.029709226 0.340550933 0.150765312
#>  [79] 0.881567638 0.001950755 0.231866038 0.250428017 0.326047915 0.943153746
#>  [85] 0.257344903 0.219001535 0.117412576 0.835769577 0.117100854 0.987048781
#>  [91] 0.950680589 0.167390233 0.288777680 0.303878748 0.937346140 0.060707961
#>  [97] 0.553356577 0.681577198 0.135135154 0.018200780          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $status
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1
#>  [75] 1 0 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr_list(dat, c("age", "active"))
#> $age
#>   [1] 0.981842521 0.480152659 0.497900751 0.969944620 0.730625570 0.871249147
#>   [7] 0.991377368 0.664240878 0.713891583 0.143768221 0.896165183 0.622561361
#>  [13] 0.699182787 0.242552443 0.909468295 0.601935858 0.954815769 0.189696660
#>  [19] 0.806798475 0.614939201 0.337632246 0.695872342 0.945397054 0.588769985
#>  [25] 0.202820597 0.361969674 0.604657195 0.927379446 0.720764208 0.851808469
#>  [31] 0.995168352 0.352305684 0.441294787 0.484173094 0.743150527 0.175656549
#>  [37] 0.581930148 0.474467288 0.115034923 0.774499406 0.453538384 0.556696330
#>  [43] 0.926556017 0.443646785 0.008501619 0.062822039 0.744874349 0.065412350
#>  [49] 0.630897832 0.961829548 0.996465983 0.668148058 0.835421227 0.343280025
#>  [55] 0.725047796 0.494380954 0.245206871 0.914296501 0.903821214 0.523425133
#>  [61] 0.652749121 0.663911420 0.893451170 0.595002735 0.608075619 0.724174157
#>  [67] 0.113207193 0.139429949 0.112331380 0.462294281 0.042387796 0.247349081
#>  [73] 0.857424026 0.850155975 0.206830324 0.029709226 0.340550933 0.150765312
#>  [79] 0.881567638 0.001950755 0.231866038 0.250428017 0.326047915 0.943153746
#>  [85] 0.257344903 0.219001535 0.117412576 0.835769577 0.117100854 0.987048781
#>  [91] 0.950680589 0.167390233 0.288777680 0.303878748 0.937346140 0.060707961
#>  [97] 0.553356577 0.681577198 0.135135154 0.018200780          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $active
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr(dat, "status")
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1
#>  [75] 1 0 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
get_attr(dat, "status", c(1, 4))
#> [1] 1 1

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
