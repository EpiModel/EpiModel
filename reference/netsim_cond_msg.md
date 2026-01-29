# Message to Find in Which Module a `condition` Occurred

This function returns a formatted string describing when, where, and why
an error, message, or warning occurred.

## Usage

``` r
netsim_cond_msg(cond, module, at, msg)
```

## Arguments

- cond:

  The type of `condition` handled (message, warning, error).

- module:

  The name of the module where the `condition` occurred.

- at:

  The time step the `condition` occurred.

- msg:

  The `condition`'s message.

## Value

A formatted string describing where and when the `condition` occurred as
well as the `condition`'s message.
