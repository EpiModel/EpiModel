# Utility Function for Printing netdx Object

Prints basic information and statistics from a `netdx` object.

## Usage

``` r
# S3 method for class 'netdx'
print(x, digits = 3, ...)
```

## Arguments

- x:

  an object of class `netdx`

- digits:

  number of digits to print in statistics tables

- ...:

  additional arguments (currently ignored)

## Details

Given a `netdx` object, `print.netdx` prints the diagnostic method
(static/dynamic), number of simulations, and (if dynamic) the number of
time steps per simulation used in generating the `netdx` object, as well
as printing the formation statistics table and (if present) the duration
and dissolution statistics tables. The statistics tables are interpreted
as follows.

Each row has the name of a particular network statistic. In the
formation table, these correspond to actual network statistics in the
obvious way. In the duration and dissolution tables, these correspond to
dissolution model dyad types: in a homogeneous dissolution model, all
dyads are of the `edges` type; in a heterogeneous dissolution model, a
dyad with a nonzero `nodematch` or `nodemix` change statistic in the
dissolution model has type equal to that statistic, and has type equal
to `edges` otherwise. The statistics of interest for the duration and
dissolution tables are, respectively, the mean age of extant edges and
the edge dissolution rate, broken down by dissolution model dyad type.
(The current convention is to treat the mean age and dissolution rate
for a particular dissolution dyad type as 0 on time steps with no edges
of that type; this behavior may be changed in the future.)

The columns are named `Target`, `Sim Mean`, `Pct Diff`, `Sim SE`,
`Z Score`, `SD(Sim Means)`, and `SD(Statistic)`. The `Sim Mean` column
refers to the mean statistic value, across all time steps in all
simulations in the dynamic case, and across all sampled networks in all
simulations in the static case. The `Sim SE` column refers to the
standard error in the mean, estimated using
[`coda::effectiveSize`](https://rdrr.io/pkg/coda/man/effectiveSize.html).
The `Target` column indicates the target value (if present) for the
statistic, and the `Pct Diff` column gives `(Sim Mean - Target)/Target`
when `Target` is present. The `Z Score` column gives
`(Sim Mean - Target)/(Sim SE)`. The `SD(Sim Means)` column gives the
empirical standard deviation across simulations of the mean statistic
value within simulation, and `SD(Statistic)` gives the empirical
standard deviation of the statistic value across all the simulated data.
