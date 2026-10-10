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

This function updates list `x` by name. An element of `new.x` that is
itself a named list is merged into the matching element of `x` the same
way, by name. An element that is a list without names, such as a
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
object, has no names to merge by, so it replaces the matching element of
`x` whole. If a function is provided to replace an element that was
originally not a function, this function will be applied to the original
value.
