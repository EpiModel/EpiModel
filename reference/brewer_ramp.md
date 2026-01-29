# RColorBrewer Color Ramp for EpiModel Plots

Returns a vector of colors consistent with a high-brightness set of
colors from an `RColorBrewer` palette.

## Usage

``` r
brewer_ramp(n, plt, delete.lights = TRUE)
```

## Arguments

- n:

  Number of colors to return.

- plt:

  `RColorBrewer` palette from
  [`RColorBrewer::brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).

- delete.lights:

  If TRUE, delete the lightest colors from the color palette; this helps
  with plotting in many high-contrast palettes.

## Value

A vector of length equal to `n` with a range of color values consistent
with an RColorBrewer color palette.

## Details

`RColorBrewer` provides easy access to helpful color palettes, but the
built-in palettes are limited to the set of colors in the existing
palette. This function expands the palette size to any number of colors
by filling in the gaps. Also, colors within the "div" and "seq" set of
palettes whose colors are very light (close to white) are deleted by
default for better visualization of plots.

## See also

[RColorBrewer::RColorBrewer](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)

## Examples

``` r
# Shows a 100-color ramp for 4 RColorBrewer palettes
par(mfrow = c(2, 2), mar=c(1, 1, 2, 1))
pals <- c("Spectral", "Greys", "Blues", "Set1")
for (i in seq_along(pals)) {
 plot(1:100, 1:100, type = "n", axes = FALSE, main = pals[i])
 abline(v = 1:100, lwd = 6, col = brewer_ramp(100, pals[i]))
}

```
