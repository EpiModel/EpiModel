# Functions to Access and Edit the Main netsim_dat Object in Network Models

These `get_`, `set_`, `append_`, and `add_` functions allow a safe and
efficient way to retrieve and mutate the main `netsim_dat` class object
of network models (typical variable name `dat`).

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
#>   [1] 0.704090369 0.006688182 0.817856184 0.214666902 0.491504228 0.559755405
#>   [7] 0.805175302 0.686336085 0.222169009 0.778594022 0.420193425 0.871153415
#>  [13] 0.608496110 0.078330556 0.444365631 0.162669908 0.760845326 0.722566311
#>  [19] 0.701210755 0.488646645 0.111995193 0.206219819 0.396012436 0.226479796
#>  [25] 0.550735588 0.823862592 0.369507225 0.907415280 0.988624807 0.848326449
#>  [31] 0.397229077 0.171014216 0.414827363 0.392055724 0.956480091 0.643862929
#>  [37] 0.750201265 0.761700642 0.478871147 0.203750858 0.141421030 0.118388054
#>  [43] 0.288137309 0.853066888 0.749456638 0.582563226 0.457053448 0.137151998
#>  [49] 0.486554086 0.866291023 0.783467444 0.811919667 0.681531546 0.016453522
#>  [55] 0.084276739 0.512968185 0.437201015 0.169162008 0.980899066 0.063738914
#>  [61] 0.923840689 0.331321878 0.619526952 0.438955051 0.139812627 0.455290131
#>  [67] 0.531221183 0.506365841 0.249873911 0.488281931 0.752215375 0.546390898
#>  [73] 0.882220431 0.986256457 0.503518561 0.725824599 0.386454617 0.014247679
#>  [79] 0.175181189 0.236391941 0.956912601 0.886869618 0.230511166 0.249074377
#>  [85] 0.720701747 0.905780326 0.654541427 0.644006039 0.936669014 0.558416004
#>  [91] 0.155095161 0.337097100 0.187390231 0.951736721 0.645011278 0.624031018
#>  [97] 0.756518194 0.466310136 0.173742835 0.473355886          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $status
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 0
#>  [75] 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1
#> 
get_attr_list(dat, c("age", "active"))
#> $age
#>   [1] 0.704090369 0.006688182 0.817856184 0.214666902 0.491504228 0.559755405
#>   [7] 0.805175302 0.686336085 0.222169009 0.778594022 0.420193425 0.871153415
#>  [13] 0.608496110 0.078330556 0.444365631 0.162669908 0.760845326 0.722566311
#>  [19] 0.701210755 0.488646645 0.111995193 0.206219819 0.396012436 0.226479796
#>  [25] 0.550735588 0.823862592 0.369507225 0.907415280 0.988624807 0.848326449
#>  [31] 0.397229077 0.171014216 0.414827363 0.392055724 0.956480091 0.643862929
#>  [37] 0.750201265 0.761700642 0.478871147 0.203750858 0.141421030 0.118388054
#>  [43] 0.288137309 0.853066888 0.749456638 0.582563226 0.457053448 0.137151998
#>  [49] 0.486554086 0.866291023 0.783467444 0.811919667 0.681531546 0.016453522
#>  [55] 0.084276739 0.512968185 0.437201015 0.169162008 0.980899066 0.063738914
#>  [61] 0.923840689 0.331321878 0.619526952 0.438955051 0.139812627 0.455290131
#>  [67] 0.531221183 0.506365841 0.249873911 0.488281931 0.752215375 0.546390898
#>  [73] 0.882220431 0.986256457 0.503518561 0.725824599 0.386454617 0.014247679
#>  [79] 0.175181189 0.236391941 0.956912601 0.886869618 0.230511166 0.249074377
#>  [85] 0.720701747 0.905780326 0.654541427 0.644006039 0.936669014 0.558416004
#>  [91] 0.155095161 0.337097100 0.187390231 0.951736721 0.645011278 0.624031018
#>  [97] 0.756518194 0.466310136 0.173742835 0.473355886          NA          NA
#> [103]          NA          NA          NA          NA          NA          NA
#> [109]          NA          NA
#> 
#> $active
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
get_attr(dat, "status")
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 0
#>  [75] 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1
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
