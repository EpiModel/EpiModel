# Assign Group Ids to an Existing Population

Assigns nodes to groups of a target size distribution, subject to a
composition rule that keeps dependent members (such as children) in
groups that contain an anchor member (such as an adult). This is the
reverse of
[`sample_groups()`](https://epimodel.github.io/EpiModel/reference/sample_groups.md):
the nodal attributes already exist, from a census age distribution for
example, and the group ids are built to match them.

## Usage

``` r
assign_groups(
  size.dist,
  role = NULL,
  anchor = NULL,
  dependent = NULL,
  n = NULL
)
```

## Arguments

- size.dist:

  Distribution of group sizes: a numeric vector of weights named by
  group size, such as
  `c("1" = 0.28, "2" = 0.34, "3" = 0.15, "4" = 0.14, "5" = 0.09)`. An
  unnamed vector is taken over sizes `1, 2, ...`. Weights are normalized
  to sum to one.

- role:

  Optional vector of node roles, one per node, such as an age group. Its
  length sets the population size.

- anchor:

  Values of `role` that can anchor a group. Every group that contains a
  `dependent` member also contains an anchor, as long as there are
  enough anchors.

- dependent:

  Values of `role` that must share a group with an anchor.

- n:

  Population size when `role` is not given.

## Value

An integer vector of group ids, one per node, with ids running from `1`
to the number of groups.

## Details

The algorithm has four steps:

1.  Group sizes are drawn from `size.dist` until they sum to at least
    `n`; the last group is trimmed so that the sizes sum to exactly `n`.

2.  One anchor is placed in each group, in random order, until either
    the groups or the anchors run out. When there are fewer anchors than
    groups, the groups without an anchor cannot receive dependents.

3.  Dependents are placed in the open slots of anchored groups, each
    open slot being equally likely, so that a group receives dependents
    in proportion to its remaining size. When the dependents outnumber
    those slots, the surplus is added to random anchored groups, which
    then exceed their drawn size; this is the only case in which the
    realized size distribution departs from `size.dist`, and a message
    reports it.

4.  All remaining nodes (further anchors, and roles that are neither)
    fill the remaining open slots in the same way.

Without `role`, or without `anchor`, nodes are assigned to the drawn
group sizes at random and no composition rule applies.

The rule guarantees cross-generational contact for every dependent. It
does not otherwise control group composition: the number of dependents
per group, the pairing of anchors, and whether a third role (older
adults, for example) lives alone or with others all follow from the size
distribution and the population shares of the roles. When the
composition itself is the quantity to control, build the population from
a table of group types with
[`sample_groups()`](https://epimodel.github.io/EpiModel/reference/sample_groups.md)
instead.

## See also

[`sample_groups()`](https://epimodel.github.io/EpiModel/reference/sample_groups.md)
and
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md).

## Examples

``` r
# Ages from a census-like distribution, then households around them
set.seed(1)
n <- 1000
age <- sample(c("child", "adult", "elderly"), n, replace = TRUE,
              prob = c(0.22, 0.60, 0.18))
size.dist <- c("1" = 0.28, "2" = 0.34, "3" = 0.15, "4" = 0.14, "5" = 0.09)
hh_id <- assign_groups(size.dist, role = age,
                       anchor = c("adult", "elderly"), dependent = "child")

table(tabulate(hh_id))                          # realized size distribution
#> 
#>   1   2   3   4   5 
#> 116 156  60  53  36 
has_anchor <- tapply(age %in% c("adult", "elderly"), hh_id, any)
all(has_anchor[as.character(hh_id[age == "child"])])  # every child has one
#> [1] TRUE

nw <- network_initialize(n)
nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "hh_id", hh_id)
est_hh <- netclique(nw, group.attr = "hh_id")
print(est_hh, by = "age")
#> EpiModel Clique Network Layer
#> =======================
#> Model class: netclique
#> Grouping attribute: hh_id
#> 
#> Layer Summary
#> -----------------------
#> Nodes: 1000
#> Groups: 421
#> Edges: 1014
#> Mean degree: 2.028
#> Isolates: 116
#> 
#> Group Size Distribution
#> -----------------------
#>                        
#> size     1   2  3  4  5
#> groups 116 156 60 53 36
#> 
#> Mean Degree by `age`
#> -----------------------
#>   adult   child elderly 
#>   1.931   2.452   1.864 
#> 
#> Arrivals: isolate
```
