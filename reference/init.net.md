# Initial Conditions for Stochastic Network Models

Sets the initial conditions for stochastic network models simulated with
`netsim`.

## Usage

``` r
init.net(i.num, r.num, i.num.g2, r.num.g2, status.vector, infTime.vector, ...)
```

## Arguments

- i.num:

  Number of initial infected persons. For two-group models, this is the
  number of initial group 1 infected persons.

- r.num:

  Number of initial recovered persons. For two-group models, this is the
  number of initial group 1 recovered persons. This parameter is only
  used for the `SIR` model type.

- i.num.g2:

  Number of initial infected persons in group 2. This parameter is only
  used for two-group models.

- r.num.g2:

  Number of initial recovered persons in group 2. This parameter is only
  used for two-group `SIR` models.

- status.vector:

  A vector of length equal to the size of the input network, containing
  the status of each node. Setting status here overrides any inputs
  passed in the `.num` arguments.

- infTime.vector:

  A vector of length equal to the size of the input network, containing
  the (historical) time of infection for each of those nodes with a
  current status of `"i"`. Can only be used if `status.vector` is used,
  and must contain `NA` values for any nodes whose status is not `"i"`.

- ...:

  Additional initial conditions passed to model.

## Value

An `EpiModel` object of class `init.net`.

## Details

The initial conditions for a model solved with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
should be input into the `init.net` function. This function handles
initial conditions for both base models and new modules. For an overview
of specifying initial conditions across a variety of base network
models, consult the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) tutorials.

## See also

Use
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
to specify model parameters and
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
to specify the control settings. Run the parameterized model with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Examples

``` r
# Example of using status.vector and infTime.vector together
n <- 100
status <- sample(c("s", "i"), size = n, replace = TRUE, prob = c(0.8, 0.2))
infTime <- rep(NA, n)
infTime[which(status == "i")] <- -rgeom(sum(status == "i"), prob = 0.01) + 2

init.net(status.vector = status, infTime.vector = infTime)
#> Network Model Initial Conditions
#> =================================
#> status.vector = s s s s s s s i i s s s s s s i s s s s s i s s s s s s s s i s 
#> s s i s s s i i s s i s s s i s s s s s s s i s s s s s s s s s s i i s s s s s 
#> s s s s s s s s s s s s s s s s s s s s s s s s s s s s
#> infTime.vector = NA NA NA NA NA ...
```
