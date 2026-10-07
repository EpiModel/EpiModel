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

## One Value Per Node

Every nodal attribute holds exactly one value per node, and
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
checks this at the end of each time step. When `set_attr` and
`append_attr` create an attribute, the nodes they do not write to are
filled with `NA` to preserve that: `set_attr` with `posit_ids` fills the
unselected nodes, and `append_attr` fills the nodes that predate the
ones being appended. This means `append_core_attr` must be called before
`append_attr` when adding new nodes, as it is what makes those nodes
known to the rest of the object.

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
#>   [1] 0.794723696 0.079009692 0.994147693 0.333935289 0.092798396 0.118813677
#>   [7] 0.814425972 0.902433855 0.626641747 0.848945338 0.528402529 0.148639226
#>  [13] 0.426021027 0.833002899 0.402158567 0.747636655 0.783656693 0.644317521
#>  [19] 0.232141421 0.517684197 0.838909269 0.593803046 0.246336028 0.240527395
#>  [25] 0.235886996 0.654152293 0.838277736 0.251922707 0.665510750 0.151081642
#>  [31] 0.225265846 0.094681246 0.823617215 0.175607729 0.230003665 0.304878073
#>  [37] 0.582349576 0.917549005 0.907157657 0.857565414 0.414115999 0.371066698
#>  [43] 0.675518590 0.005220825 0.456325933 0.436814286 0.422702053 0.257331764
#>  [49] 0.403097231 0.047443524 0.173852789 0.786053162 0.537042114 0.414580266
#>  [55] 0.511543171 0.520385041 0.677273601 0.752910193 0.750365018 0.246252101
#>  [61] 0.478001623 0.994308035 0.322863226 0.556724218 0.210162415 0.154073077
#>  [67] 0.090549808 0.525575908 0.240835476 0.341524981 0.605908218 0.317186866
#>  [73] 0.864370837 0.297693344 0.533690562 0.911750742 0.158485555 0.984418668
#>  [79] 0.213457823 0.095546803 0.432033776 0.102071220 0.936309566 0.286762381
#>  [85] 0.913731117 0.385140918 0.024500575 0.136184629 0.487673842 0.115797226
#>  [91] 0.866911953 0.692078316 0.792277329 0.621852847 0.837726101 0.634371946
#>  [97] 0.884076100 0.407352951 0.404075008 0.750559513          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $status
#>   [1] 1 1 0 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr_list(dat, c("age", "active"))
#> $age
#>   [1] 0.794723696 0.079009692 0.994147693 0.333935289 0.092798396 0.118813677
#>   [7] 0.814425972 0.902433855 0.626641747 0.848945338 0.528402529 0.148639226
#>  [13] 0.426021027 0.833002899 0.402158567 0.747636655 0.783656693 0.644317521
#>  [19] 0.232141421 0.517684197 0.838909269 0.593803046 0.246336028 0.240527395
#>  [25] 0.235886996 0.654152293 0.838277736 0.251922707 0.665510750 0.151081642
#>  [31] 0.225265846 0.094681246 0.823617215 0.175607729 0.230003665 0.304878073
#>  [37] 0.582349576 0.917549005 0.907157657 0.857565414 0.414115999 0.371066698
#>  [43] 0.675518590 0.005220825 0.456325933 0.436814286 0.422702053 0.257331764
#>  [49] 0.403097231 0.047443524 0.173852789 0.786053162 0.537042114 0.414580266
#>  [55] 0.511543171 0.520385041 0.677273601 0.752910193 0.750365018 0.246252101
#>  [61] 0.478001623 0.994308035 0.322863226 0.556724218 0.210162415 0.154073077
#>  [67] 0.090549808 0.525575908 0.240835476 0.341524981 0.605908218 0.317186866
#>  [73] 0.864370837 0.297693344 0.533690562 0.911750742 0.158485555 0.984418668
#>  [79] 0.213457823 0.095546803 0.432033776 0.102071220 0.936309566 0.286762381
#>  [85] 0.913731117 0.385140918 0.024500575 0.136184629 0.487673842 0.115797226
#>  [91] 0.866911953 0.692078316 0.792277329 0.621852847 0.837726101 0.634371946
#>  [97] 0.884076100 0.407352951 0.404075008 0.750559513          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $active
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr(dat, "status")
#>   [1] 1 1 0 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
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
