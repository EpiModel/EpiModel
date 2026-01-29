# Update List `x` Using the Elements of List `new.x`.

Update List `x` Using the Elements of List `new.x`.

## Usage

``` r
update_list(x, new.x)
```

## Arguments

- x:

  A list.

- new.x:

  A list.

## Value

The full `x` list with the modifications added by `new.x`.

## Details

This function updates list `x` by name. If `x` and `new.x` elements are
not named, the function will not work properly. If a function is
provided to replace an element that was originally not a function, this
function will be applied to the original value.
