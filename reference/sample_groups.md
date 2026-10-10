# Sample a Population of Groups from a Table of Group Types

Builds a population one group at a time from a table of group types,
such as household compositions, returning each node's group id and the
attributes its group type implies. The result supplies both the grouping
attribute for
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
and the nodal attributes for the dynamic layers of the same model.

## Usage

``` r
sample_groups(n, types, prob = NULL, attr.name = "member", sep = " ")
```

## Arguments

- n:

  Number of nodes in the population.

- types:

  The group types. Either a character vector in which each element lists
  the members of one type separated by `sep`, such as
  `c("adult", "adult adult", "adult adult child")`; a named numeric
  vector whose names are those templates and whose values are the type
  weights, in which case `prob` is taken from the values; or a list of
  `data.frame`s, each with one row per member and any number of
  attribute columns, for types that carry more than one attribute per
  member.

- prob:

  Sampling weights for the types, one per type. They are normalized to
  sum to one. The default draws every type with equal probability.

- attr.name:

  Name of the attribute column in the result when `types` is a character
  vector.

- sep:

  Separator between member labels within a character template.

## Value

A `data.frame` with one row per node and `n` rows, in group order:
`group` (integer group id, `1` to the number of groups), `type` (index
of the type the group was drawn from), and the attribute columns (one
column named `attr.name` for character templates, or the columns of the
`data.frame`s in `types`).

## Details

Types are drawn with replacement until the population reaches `n`. The
final group is drawn from the types whose size fits the remaining slots,
so that the population has exactly `n` nodes without cutting a group
short. When no type fits (no type has size one, and the remainder is
smaller than every type), the last group is truncated and a message says
so.

The realized distribution of group types, and so of group sizes and of
the attribute values, follows the `prob` weights up to sampling
variation; the expected mean group size is the weighted mean of the type
sizes. The expected share of nodes with a given attribute value is the
weighted share of that value among the members of all types. Tune the
type table to hit both.

## See also

[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
to turn the `group` column into a clique layer, and
[`assign_groups()`](https://epimodel.github.io/EpiModel/reference/assign_groups.md)
for the reverse problem of assigning group ids to a population whose
attributes already exist.

## Examples

``` r
# Household types by age group, weights chosen so that every child lives
# with at least one adult and about a quarter of older adults live alone
hh_types <- c("adult"                    = 0.15,
              "adult adult"              = 0.20,
              "elderly"                  = 0.12,
              "elderly elderly"          = 0.10,
              "adult elderly"            = 0.03,
              "adult child"              = 0.05,
              "adult adult child"        = 0.15,
              "adult adult child child"  = 0.15,
              "adult adult child elderly" = 0.05)
pop <- sample_groups(1000, hh_types, attr.name = "age")
head(pop, 10)
#>    group type   age
#> 1      1    1 adult
#> 2      2    8 adult
#> 3      2    8 adult
#> 4      2    8 child
#> 5      2    8 child
#> 6      3    8 adult
#> 7      3    8 adult
#> 8      3    8 child
#> 9      3    8 child
#> 10     4    8 adult
table(tabulate(pop$group))                 # household size distribution
#> 
#>   1   2   3   4 
#> 119 160  71  87 
prop.table(table(pop$age))                 # person-level age distribution
#> 
#>   adult   child elderly 
#>   0.582   0.240   0.178 

# Every child shares a household with an adult
has_adult <- tapply(pop$age == "adult", pop$group, any)
all(has_adult[as.character(pop$group[pop$age == "child"])])
#> [1] TRUE

# Types with more than one attribute per member
types <- list(
  data.frame(age = "adult", sex = "F"),
  data.frame(age = c("adult", "adult"), sex = c("F", "M")),
  data.frame(age = c("adult", "adult", "child"), sex = c("F", "M", "F"))
)
pop2 <- sample_groups(100, types, prob = c(0.3, 0.4, 0.3))
head(pop2)
#>   group type   age sex
#> 1     1    2 adult   F
#> 2     1    2 adult   M
#> 3     2    1 adult   F
#> 4     3    2 adult   F
#> 5     3    2 adult   M
#> 6     4    1 adult   F
```
