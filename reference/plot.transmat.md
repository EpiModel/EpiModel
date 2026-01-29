# Plot transmat Infection Tree in Three Styles

Plots the transmission matrix tree from from `get_transmat` in one of
three styles: a phylogram, a directed network, or a transmission
timeline.

## Usage

``` r
# S3 method for class 'transmat'
plot(x, style = c("phylo", "network", "transmissionTimeline"), ...)
```

## Arguments

- x:

  A
  [`transmat()`](http://epimodel.github.io/EpiModel/reference/get_transmat.md)
  object to be plotted.

- style:

  Character name of plot style. One of `"phylo"`, `"network"`, or
  `"transmissionTimeline"`.

- ...:

  Additional plot arguments to be passed to lower-level plot functions
  (`plot.network`, `plot.phylo`, or `transmissionTimeline`).

## Details

The `phylo` plot requires the `ape` package. The `transmissionTimeline`
plot requires that the `ndtv` package.

## See also

[`network::plot.network`](https://rdrr.io/pkg/network/man/plot.network.html),
[`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
[`ndtv::transmissionTimeline()`](https://rdrr.io/pkg/ndtv/man/transmissionTimeline.html).
