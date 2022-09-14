## EpiModel 2.3.1

### NEW FEATURES

- `get_transmat` adds a `deduplicate` argument to randomly select one transmitting act in the case that multiple potential transmissions occur within a time step to the newly infected person.

### BUG FIXES

### OTHER

- The `tracker.net` module was removed as an optional extension module, and this functionality is now  default for all network models (built-in and extension models). See the *Working with Attributes and Summary Statistics Vignette* for more details.
- Speed up the transmission matrix (`transmat`) storage by using `padded_vector` instead of repeated `rbind` calls. The result of `get_transmat(netsim_object)` is unchanged.
- Use `unflatten_params` in `generate_random_params` instead of bespoke code.

## EpiModel 2.3.0

### NEW FEATURES

- Network models can now use predefined scenarios with the `use_scenario` function. See the **Working with Model Parameters** vignette for details.
- The `trim_netest` function has been added to reduce `netest` object sizes in common use cases.
- The `networkLite` class now supports some additional functionalities.
- Duration = 1 time step models in `netsim` with `tergmLite = FALSE` now use `tergm` simulation to produce a `networkDynamic` class object, similar to the duration > 1 time step case.
- `netdx` now supports heterogeneous dissolution model diagnostics. See the example in `help("netdx")`.
- `plot.netsim` now supports `type = "duration","dissolution"` for homogeneous ("edges-only") dissolution models.
- Network model parameters can now be supplied to `param.net` via a `data.frame`. See the **Working with Model Parameters** vignette for details.
- A new term `fuzzynodematch` was added, to support modeling of more generalized notions of homophily.
- Network models can now be checkpointed to automate saving, recovering, and restarting simulations after interruption, as might occur with large-scale simulation jobs on high-performance computing (HPC) systems. See `help("control.net)`.

### BUG FIXES

- EpiModel's custom terms `absdiffby` and `absdiffnodemix` are now usable again (due to enabling symbol search).
- The correct formula is now used in `sim_nets_t1` when initializing the network with the non-nested edges dissolution approximation approach in `netest`.
- References to `network` class internals have been removed from `net.utils.R` to make the code function properly with `networkLite`s.
- `update_dissolution` now correctly handles duration 1 models.
- `get_edgelists` and `get_cumulative_edgelists` now throw informative errors when the `network` refers to non-existing networks. `get_partners` and `get_cumulative_edgelists_df` will throw errors as well as they call those functions internally.
- fix double memory allocation in `netsim` with `ncores > 1`

### OTHER

- The `updater.net` module was removed as an optional module but its functionalities are now the default behavior for all network models (built-in and custom).
- The parameters `param.updater.list` and `control.updater.list` in `param.net` and `control.net` respectively were renamed `.param.updater.list` and `.control.updater.list`. The leading dot indicates that these are built-in EpiModel elements not to be confused with the user-defined ones.
- Dissolution models using the `nodefactor` term have been deprecated.
- `stergm` controls have been deprecated in favor of `tergm` controls.
- Formula-style simulations are used consistently for both `ergm` and `tergm` simulation, requiring control arguments of class `control.simulate.formula` and `control.simulate.formula.tergm`.
- `ndtv` is added to "Suggests" (again) and `plot.transmat` now accepts the `transmissionTimeline` in the `style` argument (again).
- Systematic review and improvement of documentation across the package.
- `netest` and `netdx` now have default control arguments following R's default argument mechanism.
- Simulation names are now attached to `netsim` object fields produced via the `save.other` mechanism.
- Imputed durational corrections for onset-censored edges are now sampled from a geometric distribution rather than simply using the mean of that distribution.


## EpiModel 2.2.1

### NEW FEATURES

- Improved optional module `updater.net` allowing it to update the model controls as well as the parameters. See the vignette, "Working with model parameters."
- General updates to the names and content of the included vignettes.

### BUG FIXES

- Fix dissolution model statistics calculations for `netsim` in the case with a model with an "end horizon" (when the network is not resimulated at the end of the time series).
- Fix duplicate printing issues across `print.netsim`, `print.param.net`, and `print.control.net`.
- Fix use of `all.equal` in unit tests as requested by CRAN.
- Change defaults of newly introduced cumulative edgelist functionality to not store it (improves speed).
- Allow `set_attr` for `posit_ids` to return unchanged `dat` object.

### OTHER

- We have changed the names of arguments from the function `get_partners` newly introduced in EpiModel v2.2.0: `max.age` is renamed to `truncate` for consistency with the other cumulative edgelist functions; `only.active` is renamed `only.active.nodes` to clarify that this argument subsets by *nodes* and not by *partnerships*.

## EpiModel 2.2.0

### NEW FEATURES

- Developed a general approach to tracking and querying historical and contacts, called a cumulative edgelist. This may be used, for example, to query the recent but non-current contacts of newly infected nodes. See the vignette, "Working with network objects".
- A `create_dat_object` helper function was added to standardize the creation of the core `dat` object within `initialize.net`.
- The current timestep within `netsim` simulations is now stored in the `dat` object and accessible with `get_current_timestep`. This eliminates the need to explicitly pass `at` as a function argument, although that is still allowed.
- Addition of the `get_param_set` function that extracts from a `netsim` object the set of parameters used by each simulation. See the help page: `help("get_param_set")`.
- Developed a mechanism to store nodal attribute history over the course of a `netsim` simulation. See the vignette, "Working with attributes and summary statistics."
- Developed an optional module to define prevalence statistics (also called "epi stats") as functions to be passed to the model as control settings before each `netsim` simulation. This allows users to avoid updating the `prevalence.net` module. See the vignette, "Working with attributes and summary statistics."
- Developed an optional module allowing the update of the model controls and parameters over timesteps within `netsim` simulations (i.e., time-varying parameters). See the vignette, "Working with model parameters."
- Improved the random parameterization programming interface to allow correlation between parameters in each simulation (e.g., the ability to pass in a multivariate parameter set for each simulation). See the vignette, "Working with model parameters."

### BUG FIXES

- When calling `plot` on a `netsim` object, the arguments in the ellipsis (`...`) are now correctly passed to the `base::plot` call.
- When trying to use the built-in `group` attribute, `netsim` will now output a more explicit error if the values used are not only `1` and `2`.
- Fixed the names of the target formation statistics in `netdx` when `edapprox == FALSE` that were causing the plotting functions to misbehave.
- Simplification of the `set_transmat` function removing the assumption that `dat$stats$transmat` was to exist only if `at != 2` (thanks to @thednainus).
- More consistent formation and dissolution statistics print between `netdx` and `netsim`.
- Removed duplication in the printing of the parameters when a parameter was defined both as fixed and as random.
- When using custom `netsim` modules with `type == NULL`, some built-in modules no longer stop because they required `type` to be a string.
- Fixed issue with `Error`, `Warning` or `Message` in `netsim` printing twice.
- Fixed problem with unique ID counter not saved by `saveout.net`, resulting in the unique ids to start a 1 again when restarting a model from a previous simulation.

### OTHER

- The new home for EpiModel on Github is: <https://github.com/EpiModel/EpiModel>. It was previously located on the `statnet` organization on Github.

## EpiModel 2.1.0

### NEW FEATURES

- Summary network statistics for `netsim` class objects (epidemic simulations) are now available when `tergmLite` is used. Previously these network statistics were only available when `tergmLite = FALSE` but updates to the `networkLite` class in `tergmLite` made this possible. These network stats are output with `save.nwstats` and `nwstats.formula` arguments in `control.net`.
- Developed a general storage and printing mechanism for the recently developed random parameterization interface. See the help page `help("param.net")`.
- Cleaned up the handling of the initial network simulation in `initialize.net`, so that the user-facing code in that function is more readable, and the more complex code is put in `sim_nets_t1`.

### BUG FIXES

- Addressed errors in handling cross-sectional network simulations embedded within a dynamic epidemic model (e.g., contacts that last the length of one time step, specified by `duration = 1` in `dissolution_coefs`).

## EpiModel 2.0.5

- Minor updates for CRAN submission.

## EpiModel 2.0.4

### NEW FEATURES

- Added new `update_params` function to add new parameters to an existing list of network parameters specified in `param.net`. This aids in workflows that distinguish fixed parameters versus varying parameters that may change across scenarios or simulations.
- Added new general interface for random parameters in network models that allows randomly drawing a parameter value from a specified statistical distribution, where the distribution may either be a sampling of discrete values or a factory function for any of R's random statistical distribution functions. See help file for `generate_random_params` for examples.
- Implemented a standardized approach with helper functions for setting core attributes (those nodal attributes which should be present in any workflows) in network models. This functionality is specified with `append_core_attr` function in the initialization and arrival modules in any extension models. This includes a standardized implementation of persistent, unique IDs as an attribute that remains constant for nodes even with open population models.

### BUG FIXES

- With the use of a standardized core attribute framework that now correctly handles unique IDs in all models, now the transmission matrix objects output from any network model work consistently and correctly for both closed population and open population models.
- For DCM models with `dt < 1`, fix `NA` output for any `.flow` variables.
- Reduce complexity of some unit tests that were stochastically generating errors due to ERGM MCMC estimation problems.
- Fix problem with temporally extended status variables in network models (i.e., tracking of disease status history across time steps) by simplifying the general approach that works across built-in and extension model types.

### OTHER

- Reimplemented the handling of relational age diagnostics in `netdx`, with updated numerical summaries in `print.netdx` and visuals in `plot.netdx`. Because relational ages are left-censored for any edges that existed at time zero, this led to a misleading diagnostic that ages were lower than the targeted durations. Imputation of a start time for those edges was added, with the option in `plot.netdx` to visualize with imputed start times (default = TRUE) or not.
- More consistent approach to trimming unneeded environmental data from ERGM objects implemented with `statnet.common::trim_env()`, used in `netest`.
- Reimplemented the `netest` "edges dissolution approximation" for efficient estimation of a temporal ERGM via a cross-sectional ERGM estimation with adjustment of formation model coefficients (see `netest` help page). This new approach further reduces bias in the approximation method, plus now allows for non-nested dissolution models (i.e., dissolution formula does not need to be a subset of the formation formula).

## EpiModel 2.0.3

### NEW FEATURES

- Implemented an error catching approach for `netsim` so that epidemic modules with errors or warnings are clearly identified in the console.
- Allow saving the transmission matrix with `save.transmat` in `control.net` independently of using tergmLite methods (previously use of tergmLite did not allow for saving these data).
- Added an `infstat` parameter to the internal `discord_edgelist` function used in the infection module, to allow for arbitrary specification of which disease statuses are considered infectious for the purpose of dyad discordance.
- Added ability to vary node size in `plot.netsim` with `type = 'network'` with `vertex.cex` parameter.

### BUG FIXES

- Fix issue for `plot.netdx` when plot legend set to `FALSE`.
- Updated `print.netsim` that does not error when displaying new epidemic modules for extension models.
- Use appropriate tergmLite resimulation methods for `netsim` for networks with duration of 1 (i.e., one-time contacts handled with cross-sectional ERGMs).

### OTHER

- Further minor edits/updates to EpiModel 2.0 migration documentation (posted on <https://www.epimodel.org/>).

## EpiModel 2.0.2

- Minor updates for CRAN submission.

## EpiModel 2.0.1

- Minor updates for CRAN submission.

## EpiModel 2.0.0

- This release introduces a major update to the EpiModel package infrastructure and application programming interface for both built-in models (primarily used for teaching purposes) and extension models (primarily used for research purposes). The major substantive changes are summarized in a EpiModel 1.x to EpiModel 2.0 migration guide on our primary website: <https://www.epimodel.org/>.

## EpiModel 1.8.0

### BUG FIXES

- Improve error handling for inputs to `dissolution_coefs`.
- Skip dissolution diagnostics in `netdx` if a static ERGM is passed.

### NOTE

- This will be the last version of EpiModel 1.0 before major revisions to the package infrastructure and API to be released in EpiModel 2.0.

## EpiModel 1.7.5

### NEW FEATURES

- Add foundation updates to support the `tergmLite` package (to be released).
- Print network statistic diagnostics stored in `netsim` object with `print.netsim(x, formation.stats = TRUE)`.

### BUG FIXES

- Fix issue `plot.netsim` default colors when number of variables exceeds 3.

### OTHER

- Add example for differential homophily in a TERGM dissolution model in `dissolution_coefs` documentation.

## EpiModel 1.7.3

### NEW FEATURES

- Add a `skip.dissolution` argument to `netdx` to skip dissolution diagnostics for computational efficiency.

### BUG FIXES

- Fixes output formatting of network stats saved during `netsim` runs.
- Correctly errors when running dynamic network diagnostics with `netdx` on a cross-sectional ergm.
- Remove old unused utility functions.
- Enforce depend on ergm >= 3.10 package.

## EpiModel 1.7.2

### BUG FIXES

- Fix logical evaluations in `get_sims` and unit tests in `plot.transmat` found during CRAN submission for v1.7.1.

## EpiModel 1.7.1

### NEW FEATURES

- Two helper functions, `get_args` and `apportion_lr`, ported over from EpiModelHIV.
- Two custom ERGM terms, `absdiffby` and `absdiffnodemix`, ported over from EpiModelHIV.

### BUG FIXES

- Fix linked functions in embedded Shiny apps broken in v1.7.0.
- Update handling of parameter and module name changes related to births/deaths to arrivals/departures renames in v1.7.0.
- Reduce complexity of verbose output so that it can generalize across EpiModel extension packages.
- Fix bug in `as.data.frame.dcm` when `dt` control setting \< 1.

## EpiModel 1.7.0

### NEW FEATURES

- Add `keep.tnetwork` parameter to `netdx` to allow for retaining the full `networkDynamic` object during dynamic network simulations. Relatedly, add support for `get_network` to extract those networks from `netdx` objects.
- Change the default handling of `as.data.frame` function for processing model output for all three model classes (DCM, ICM, and Network Models) to generate a stacked data frame of all simulations (instead of row means across simulations). **This is a breaking change that may require updating old code.**
- New `as.data.frame.netdx` function extracts the timed edgelists directly from a `netdx` object.
- The `get_nwstats` function now extracts data frames of network statistics from both `netdx` and `netsim` objects.

### BUG FIXES

- Improve functionality and error handling of `get_sims`.
- Fix problems with color handling of network statistics plots in `plot.netdx` and `plot.netsim`.
- Enforce maximum number `ncores` argument in `netdx` to prevent over parallelization of simulations.
- Removed the redundant storage of the timed edgelist data in `netdx` objects.
- Fix errors in calculation of population sizes in verbose module that prints simulation output to the console.

### OTHER

- Add warning for input parameter with a name `act.rate.m2` for network models in `param.net`, as this is an unused parameter for built-in models.
- Updated parameters and documentation throughout EpiModel for vital dynamics parameters and processes to reflect a more general method of demographic in-flows and out-flows from the population. Previous terms were births and deaths; new terms are arrivals and departures. The default parameter for births was previously `b.rate`; it is now `a.rate`. Inputs of a `b.rate` parameter yield a message and will automatically set `a.rate` to the `b.rate` value. **This is a breaking change that may require updating old code.**
- Adapted y-axis limit calculation for all stochastic plots to depend on dynamic range of data displayed instead of full data range.
- Changed the default plot type for static diagnostics in `netdx` (when `dynamic` is set to `FALSE`) to smoothed rolling averages instead of the full MCMC trace. The trace plots may be turned back on with `sim.lines = TRUE`.

<br>

## EpiModel 1.6.5

### NEW FEATURES

- `netdx` now includes a new argument, `sequential`, for static diagnostics that mirrors the same argument from `ergm::simulate.ergm` to simulate from MCMC chains based on previous draws versus new draws.

### BUG FIXES

- Fix `mutate_epi` output when new variable is a constant.

### OTHER

- Move `ggplot2` from depend to import.
- References added for publication of Journal of Statistical Software methods paper on EpiModel: Jenness SM, Goodreau SM, Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. DOI: 10.18637/jss.v084.i08.

<br>

## EpiModel 1.6.1

### OTHER

- Fixed minor issue with unit tests using `identical` function causing some to fail under alternative BLAS/LAPACK implementations.

<br>

## EpiModel 1.6.0

### NEW FEATURES

- `as.data.frame` methods for `netsim` and `icm` classes now allow creation of a single data frame with epidemic outcomes across multiple simulations, where previous only single individual simulations would be output. This is specified with the `sim = "all"` parameter when `out = "vals"`. See the help page for examples. This "tidy" data format allows for easier integration with external plotting and analysis approaches, including ggplot2.
- `geom_bands` is a new "geom" for use by `ggplot2` to facilitate plotting of simulation intervals given a specified lower and upper quantile set. Examples of plotting ICM simulations are provided, and the same principle applies for network models. As a result of this, `ggplot2` was added as a depend.
- `truncate_sims` is a new utility function that takes truncates the time series of a `netsim` or `icm` class object at a specified time step. This truncation will remove all epidemic output before that time step, and reset the control settings to start at that time step. This is useful in our modeling workflows when we need to remove a pre-intervention burnin period from the model simulations.
- `init.net` allows you to pass in a vector of backwards-looking infection times for those initally infected at t_1 through the `infTime.vector` parameter. Combined with the `status.vector` parameter, this provides users maximal control over who is infected and for how long as initial conditions.

### BUG FIXES

- Fixed bug in DCM Shiny app related to plotting prevalence vs count outcomes.
- Removed unneeded and unused input parameters from `discord_edgelist` function.
- Fixed issue where SIS/SIR models with vital dynamics, and a low mortality rate relative to the recovery rate (which is typical) would get very long initial infection times assigned at t_1.

### OTHER

- Changed the title (actually, it's a subtitle) in the DESCRIPTION to: "Mathematical Modeling of Infectious Disease Dynamics".
- Deprecated the `status.rand` argument for `init.net` and `init.icm` that allowed users to specify a random number of initially infected. Support for this got too complex for a little (or never) used argument, and users interested in randomly setting the initial number infected may control this more flexibly with the `status.vector` parameter.

<br>

## EpiModel 1.5.0

### NEW FEATURES

- Add `grid` argument to plot functions to overlay a grid on line plots.

### BUG FIXES

- Fix bug in `plot.netdx` examples in help file.

### OTHER

- Reset the `verbose` default for network models to `TRUE` (reverts change in v1.3.0 specifically for network models).
- Rename `leg` argument name (to add default legends to plots) to `legend`. Note this is backwards-incompatible because of fuzzy matching with other function arguments starting `leg`; prior model code must be updated.
- Change default transparency level to 0.5 (if unspecified).

<br>

## EpiModel 1.3.0

### NEW FEATURES

- In `control.dcm`, `nsteps` may now be a vector of time steps or, as before, an integer containing the number of time steps within a DCM simulation. For example, `control.dcm(..., nsteps = seq(1980, 2015, 1/12), ...)` for solve for monthly outputs from a range of dates from 1980 to 2015.
- `mutate_epi` for adding new variables to a epidemic simulation object now works for all three model classes.

### BUG FIXES

- Outputs from `param`, `init`, and `control` functions are now dual-classed as lists as well as their native classes.
- When passing a `new.mod` into `control.dcm`, printing the `control.dcm` object no longer yields a warning and instead prints the function name.

### OTHER

- Update handling of transparent colors within `transco` to use the base `adjustcolor` function.
- Derivatives tracking a "flow" or the size of a transition between compartments for DCM simulations (e.g., disease incidence) often output `NA` for the final value, creating issues with analyzing those data. Those `NA`s are replaced with the penultimate value of that vector.
- Simplify printing of `dcm`, `icm`, and `netsim` objects to list "Variables" together instead of dividing them into compartments, flows, and other.
- Change the `popfrac` default for plotting `dcm`, `icm`, and `netsim` objects to `FALSE`. This avoids any problems when prevalences are already stored within the model simulation.
- Change the `verbose` default for control functions to `FALSE`.

<br>

## EpiModel 1.2.8

### NEW FEATURES

- Print simulation number and prevalence value for static network plots in `print.netsim` when `sims` is `mean`, `min`, or \`max.

### BUG FIXES

- Add new line at end of `print.coefdiss` output.
- Tighten the default ylim ranges for `plot.netsim`

### OTHER

- Include error check for duration \< 1 in `dissolution_coefs`.
- Update documentation in a number of places.

<br>

## EpiModel 1.2.7

### NEW FEATURES

- Add new `mutate_epi` function inspired by the `dplyr` package, to add post-hoc summary statistic calculations to completed network simulations. See the function help file for examples.
- Added a speedy `get_degree` function that returns a vector of current network degree for each person in a network.

### BUG FIXES

- Updated internal plot functions that calculate prevalences.
- Disable verbose output if running network models in parallel.

### OTHER

- Allow network simulations of 1 time step (mainly used for debugging and testing).

<br>

## EpiModel 1.2.6

### NEW FEATURES

- Updates to `as.phylo.transmat` to fix issues with vertex exit times and to now accept multiple seed vertices if multiple seeds are detected, returning a list of phylo objects of class `multiPhylo` following the convention of `ape::read.tree`.

### BUG FIXES

- Corrected an error governing the birth rate of 2-group, open-population deterministic compartmental models (DCMs).

### OTHER

- Updated license to GPL-3.

<br>

## EpiModel 1.2.5

### NEW FEATURES

- Added multicore functionality to simulating stochastic network models with `netsim`. This only supports single-node frameworks currently, using the `doParallel` package. Run models in parallel by using the `ncores` parameter in `control.net`.
- Modifications to the `as.phylo.transmat` function to construct the phylo tree with all network vertices as phylo-tips and all transmissions as phylo nodes.

### OTHER

- General code cleanup and improvement of package tests to increase coverage about 90%.

<br>

## EpiModel 1.2.3

### NEW FEATURES

- The stochastic network model Shiny application now features adaptive concurrency levels with ERGMs including that network statistic.

### BUG FIXES

- `plot.netsim` now correctly functions for diagnostic plots (`type = "formation"`) when summary statistics contain variable names with numeric values as suffixes.
- Avoided duplicate reinitialization of persistent IDs for network models started with a full STERGM fit (`edapprox = FALSE` in the `netest` function).
- Fixed error for stochastic network model simulations in `netest` when models were fit with the full STERGM method.
- Automatically set `depend` parameter to `TRUE` in `control.net` when user passes in any new birth or death modules.

<br>

## EpiModel 1.2.2

### NEW FEATURES

- New translation and plotting functions for temporal transmission chains measured in stochastic network models. These include a dendogram using methods from the `ape` package and a transmission timeline from the `ndtv` package. See the help files for the `as.phylo.transmat` and `plot.transmat` functions.
- Added a Shiny application for stochastic network models. This may be accessed from within the package with `epiweb(class = "net")`. It is also hosted online at shinyapps.
- `get_sims` function is used to extract individual simulations from larger `netsim` objects. This function has been updated to include a `var` argument that allows for automatic calculation of which simulation is closest to the mean across all simulations for extraction.
- Added a quantile extraction method for `as.data.frame` method for `icm` and `netsim` classes. This will provide a data frame of output corresponding to defined quantiles across all simulations contained within a model object.

### BUG FIXES

- Supress warnings with the lowess smoother in `plot.netsim` in cases where there are `NA` values in the epidemiological output.
- Removed error check for `control.net` when `type` is missing, and automatically sets `type` to `"SI"`. This will impact extensions to EpiModel in the case when the default transmission module is replaced.
- Fixed bug in `netdx` on calculating summary statistics from models with multiple structural zeros for target statistics.

### OTHER

- Changed the default of `status.rand`, which controls whether the number initially infected in stochastic epidemic models, to `FALSE`. This will ensure that exactly the number specified in `init.icm` and `init.net` are matched in each simulation.
- Fully removed the `netsim_parallel` function from the package. See the EpiModelHPC extension package at <https://github.com/EpiModel/EpiModelHPC> for running network simulations in parallel.

<br>

## EpiModel 1.2.1

### BUG FIXES

- `check_bip_degdist` now uses more tolerant checks of equality when comparing bipartite mode statistics.
- Fixes a formatting issue with output for DCMs run with the `dcm` function.
- Fixes a variable name collision problem for epidemic plotting functions.
- Removes long burn-ins from network model estimation in `netest` to improve performance of fitting models.
- In stochastic network models, one may now remove built-in modules, such as `births.FUN`, from the dynamic workflow by setting the argument value for that module to `NULL` in the `control.net` inputs.

<br>

## EpiModel 1.2.0

### NEW FEATURES

- The `calc_eql` function now returns test statistics invisibly.
- Major overhaul of plotting functions for stochastic model plots. `plot.netsim` is now a separate method for epidemic plots (it was previously a function call to `plot.icm`), with function arguments and default settings more consistent across plotting functions. There may be minor backwards incompatibility for some epidemic plots. Network statistic plots in `plot.netdx` and `plot.netsim` now use the same methods and share the same defaults. The defaults for these plots will be to plot smoothed quantile bands (the IQR) and means of simulations without the individual simulation lines. Any individual elements may be toggled on or off as before.
- Modules are now listed in the output for `param.icm` and `param.net` classes.
- Removed `dissolution` argument to `netest`. This argument specified the right- hand sided dissolution formula for temporal ERGMs. It was removed because this formula was already specified in the `dissolution_coefs` function, the output of which is passed to `netest`, thereby removing the duplication.

### BUG FIXES

- `as.data.frame` methods for stochastic models remove `NA` from individual simulations when calculating row means.
- Fixed bug in network birth module for assigning infection status for incoming nodes.
- The `verbose` parameter in `netest` now correctly controls the model fitting output level in the underlying `ergm` and `stergm` functions.
- `merge.netsim` now correctly checks elements of two objects to be merged when the classes of those elements may be of length greater than 1.

### OTHER

- Major updated internal package function testing for more reliable performance.
- Added `...` argument to `epiweb` to pass additional arguments to `shiny::runApp`.
- Importing the `graphics`, `grDevices`, `stats`, and `utils` packages as required by CRAN.

<br>

## EpiModel 1.1.6

### NEW FEATURES

- Built-in parallelization of stochastic network model simulations directly within the package with the `netsim_parallel` function has been deprecated. This functionality has been replaced with model simulation functions within the `EpiModelHPC` extension package: <https://github.com/EpiModel/EpiModelHPC>
- Cosmetic and functional updates to built-in Shiny applications accessible within the package via `epiweb`.
- New function, `calc_eql`, calculates whether a model of any class in EpiModel has reached an equilibrium state over a defined time series. Equilibrium is defined as the absolute value of the difference of the maximum prevalence and minimum prevalence over a specified time series falling below a specified threshold. For stochastic models, these values are calcualted based on the mean of the individual time series simulations.
- `netest` now includes a new argument, `nonconv.error`, that will send the function to an error state if the ERGM did not coverge after the specified number of interations. The default is to allow for a nonconverged model fit to be returned. Requiring an error may be helpful when running a number of models in batch mode.

### BUG FIXES

- Within the built-in deterministic compartmental models solved with the `dcm` function, there was an error in the calculation of flows (e.g., disease incidence or number of deaths per unit time) when the models were integrated with methods other than the "Euler" solution. Flows are now calculated correctly for all numerical integration methods supported via the `deSolve` package.
- Minor bugs in the default deaths module for stochastic network models were corrected.

### OTHER

- `netest` will now check to ensure that the formation and dissolution models are in allignment (terms specified in the same order) and that dissolution model is of proper forms (see v1.1.4 notes).

<br>

## EpiModel 1.1.4

### NEW FEATURES

- A limited set of heterogeneous dissolution models now allowed for network models (#184): edges + nodematch, nodemix, or nodefactor formulas now supported. See help file for `dissolution_coefs` for examples.
- Network models now feature more consistent and flexible use of persistent IDs for networkDynamic objects (#199). This involved adding a new control setting, `use.pids` in `control.net`. See help("persistent.ids") in the `networkDynamic` package for more background.
- Interventions are added to all model classes (#20). For DCMs, ICMs, and network models, there are new parameters, inter.eff and inter.start, for the efficacy and starting time of the intervention. This generic intervention has the effect of reducing the probability of transmission given a contact between a susceptible and infected person by the efficacy parameter.

### BUG FIXES

- Fixed error in `births.net` module that set the default `entrTime` and `exitTime` attributes twice for bipartite models (#205).
- Plotting for all model classes now allow setting `xlab` and `ylab` (#206).
- `get_sims` extraction now outputs correct data when object contains single simulation.

### OTHER

- More robust testing of functions.
- Updated hyperlinks within doc files to new github-based website.

<br>

## EpiModel 1.1.3

### NEW FEATURES

- The `skip.check` argument for `control.net` is even more flexible, to allow for passing different class elements into `netsim` with original models.
- New `param.error` argument for `merge.netsim` that allows bypassing the stop error if the parameters and control settings from the two merged objects are not identical.
- `control.net` has new `module.order` argument to provide control of the order in which modules are evaluated within each time step. The default ordering is maintained as explained in the updated help file.

### BUG FIXES

- `netsim_parallel` now returns the correct object if used for single simulations or on single cores.
- `plot.icm` removes NA values from the data when calculating `ylim` and the quantile bands.

### OTHER

- `netest` now implements an improved "Edges Dissolution Approximation" via the `edapprox` argument.
- Several documentation updates.

<br>

## EpiModel 1.1.2

### NEW FEATURES

- Implement `control.dcm` option `dede`, which if true allows for delayed differential equations to be passed into a new model solved with `dcm`.
- New option for `netdx` to simulate static diagnostics from an ERGM, rather than the temporal diagnostics (still the default). This will help better diagnose poor dynamic model fit when using the edges dissolution approximation (#175).
- Plot option added for `netdx`, with the `method` parameter, to plot boxplots of the simulations against the target statistics. The default is still the line plots (#191).
- Additional summary elements may now be plotted with `netdx` objects, similar to epidemic data plots: mean lines and quantile bands. Additional arguments added to allow toggling of these along with individual simulation lines and target lines.
- Print method for `netdx` is updated, along with a new statistic for the percent deviation between the simulation means and target statistics (#192).
- Added other epidemiological outcomes saved in user-defined modules to print output with `print.netsim` (#183).
- New function `get_sims` will subset and extract entire simulations from `netsim` objects with multiple simulations. A vector of simulation numbers may be specified, or if set as "mean", the simulation with the infected prevalence closest to the means across all simulations will be chosen.

### BUG FIXES

- Object elements saved in stochastic network models with the `save.other` parameter in `control.net` may now be merged with `merge.netsim` (#185).
- Quantile band is displayed in `plot` for ICMs and network models when the `y` argument is specified (#188).

### OTHER

- Package `deSolve` moved from import to depend (#194).

<br>

## EpiModel 1.1.1

### NEW FEATURES

- Added dissolution diagnostics in `netdx`, for the proportion of edges that dissolve per time step, as another diagnostic for the dissolution model (#53).
- Network plots with `plot.netsim` now allow specifying `"mean"`, `"min"`, or `"max"` to plot the network at with the most average, maximum, and minumum disease prevalence at the specified time step (#73).
- Network models may now use time-varying recovery rates, similar to the previous time-varying infection probabilities and act rates. The documentation for the `param.net` function has been updated with details (#65).
- New control setting for DCMs, `param.sens`, that allows bypassing the default behavior of evaluating parameters with length greater than 1 as sensitivity analyses. This should be used for single-run models if passing in parameters with arbitrary form.

### BUG FIXES

- Print functions for initial condition processing functions now handle list and data frame structures (#135).
- Fix bug for new DCMs in which the initial condition names include standard integrated initial condition names (#160).
- Several bugs fixes related to network diagnostics for models with offset terms in the formation model. Also related formation diagnostics plots in `plot.netsim` fixed (#164).
- Initialization of infection time for stochastic SIS/SIR models with two groups or modes now fixed (#102).
- Edges population size correction module, `edges_correct`, now runs for any dependent network simulations, not just if built-in vital dynamics modules are called (#141).

### OTHER

- The new website for the EpiModel project is <https://www.epimodel.org/>
- Added a new example of a SEIR Ebola DCM in the "Solving New DCMs with EpiModel" tutorial.
- The shiny apps now use the single file method (#155).
- Exported and added documentation for the `verbose.icm` function (#71).
- Other elements saved in network simulations with the `save.other` control setting in `control.net` are now printed as output in `print.netsim` (#174).

<br>

## EpiModel 1.1

### NEW FEATURES

- Added three new extraction functions for network models (`get_network`, `get_transmat`, and `get_nwstats`) which extract the network objects, transmission matrices, and data frame of network statistics from a completed `netsim` simulation. These functions also support extraction of network model simulations with multiple networks (see API note).
- The plot function for `netsim` objects now has an argument, network, for plotting network statistics and static networks (`type = "formation"` and `"network"`, respectively) in simulations with multiple networks.
- For stochastic epidemic plots, added an option `mean.smooth`. If `TRUE`, this uses a lowess smoother on the outcome variables of interest. This is helpful in visualization of low-count outcomes like disease incidence.
- Automatic parallelization of network models is now possible with the `netsim_parallel` function. Note that this is experimental and has not been tested extensively across platforms, so bug reports are welcome. Two parallel methods are supported: `doParallel` for multiple cores on a single node, and `doMPI` for multiple cores across multiple nodes. The latter requires an MPI installation on a linux-based cluster.
- Network diagnostics in `netdx` also accepts a new `ncores` argument, which will run the diagnostic simulations and calculations on those simulations in parallel on a specified number of cores (single node only).
- Added an argument, `skip.check`, for the control settings in both ICM and network model classes, which overrides the default error checking of parameters, initial conditions, and control settings. This should only be used for original models with new modules that may unnecessarily trigger a check error.
- Added an argument, `save.other`, for the control settings in network models, which is a character vector of other elements from the main data list, `dat`, to save out in the simulation.
- Added an argument, `start`, for the control settings in network models, which is a starting time step to resume simulations. In this case, the `x` argument in `netsim` is a previously saved `netsim` object rather than a `netest` object. The `start` argument should be one integer higher than the `nsteps` in that earlier `netsim` object. The `nsteps` argument should now be the final steps for the simulation. Note that this requires specifying `save.other = "attr"` in the control settings, as well as saving the networks.
- Added progress bars for `netdx` diagnostic simulations for computationally intensive parts of the simulations.
- Network model estimation with `netest` now provides an output argument. When using the edges dissolution approximation (`edapprox = TRUE`), one may set output to `"sim"` to save a static simulation network instead of the `ergm` object as an element of the `netest` output. This is mainly for file size efficiency.

### USER INTERFACE CHANGES

- The internal representation of disease status as an individual-level attribute in the stochastic ICM and network models has been changed from number `(0, 1, 2)` to character `("s", "i", "r")`. This changes little when running the integrated models, and has greater implications for the API when editing modules. But one change for integrated models is that the status vector passed into the initial conditions functions must now be in this new format. This also impacts the expansion of EpiModel for original models.
- The `zeromarg` argument has been removed from `plot.netsim` for static network plots (`type = "network"`) to reduce potential issues with setting default margins on plots. Now they must explicitly be set with standard par options.

### API CHANGES

- For ICMs and network models, the internal main data object has been renamed from `all` to `dat` to prevent function name conflicts. Additionally, all summary output is now stored within `dat$epi`, whereas the previous location was `all$out`.
- The ordering of built-in modules within a time step for network simulations has been changed such that the network resimulation module is run before the infection module. There should be no substantive differences in model results, but this provides a more logical consistency between edges toggled on at a time step and the infections that may occur over those edges.
- In network models, two preset functions have been changed to replaceable modules: `edges_correct` and `verbose.net`. The former performs the adjustment to the edges coefficient for network models with population size changes, in order to preserve the mean degree; for mass action epidemic models, for example, one would not want this adjustment, so the module should be set to NULL in `control.net`. The latter performs the printing of simulation results to the console. Both functions are now listed in the modules help file accessed by: `help(modules.net)`.
- Evaluation of parameters, initial conditions, and control settings in the core parameterization functions is now more stable, and also more flexible. Defaults for the fixed arguments are now included in the documentation.
- Users may now bypass the built-in `param`, `init`, and `control` functions altogether for original ICM and network models, because the definition of new and replacement modules occurs within the control functions themselves. The existing control functions should be used as a template if one is considering replacing these parameterization functions.
- Users may also bypass any of the built-in modules in network models (see the list in `control.net`) by setting the argument for that module to `NULL`. This may be replaced in the future by a user-defined ordered vector of modules.
- The `x` argument in netsim may now be a list of `netest` objects. This would be used only if supplying new simulation modules that know how to process that data structure. The motivation for this is to allow original models with multiple networks simulated (e.g., a main partnership network and a casual partnership network).

### BUG FIXES

- Minor updates and bug fixes to the two built-in Shiny applications (accessed via the `epiweb` function). These apps now benefit from the more stable parameterization functions.
- Updated the `print` method for the `param.net` class to handle parameters that are lists or data frames.
- `merge.netsim` now ignores any differences in the environment of the `nwstats.formula` control, previously preventing proper merging of some network model simulations.

### OTHER

- Added new test cases for running new DCMs, ICMs, and network models, following the vignette examples (see <https://www.epimodel.org/>).

<br>

## EpiModel 1.0.2

### INTERFACE CHANGES

- The trans.rate and trans.rate.g2/m2 parameters have been renamed to inf.prob and inf.prob.g2/m2 to better characterize that they are probabilities, rather than rates, and towards infection of persons in that group/mode.
- Added new documentation for newly exported utility functions for network models, mostly used in the birth/entry modules. Now users may directly edit these modules and use the utility functions without explicitly adding them to the global environment.
- Added warning message for network models in which there is a vertex attribute for status added to the network but not referenced in the formation formula, in which case the initial conditions for status will still be derived from the input for init.net. This does not apply to "serosorting models", which reference status in the formation formula, and which require setting status as a vertex attribute on the network before calling netest.

### BUG FIXES

- Network models with passed network attributes in the formation formula in open populations now do not generate an error for persistent ID numbers in the latest versions of the tergm and networkDynamic package.
- Fix bug in printing simulation progress in network and ICM class models when verbose is not specified in the control settings.
- Running netdx diagnostics with offset terms in the formation no longer generates an error.
- Simulating an SIR serosorting network epidemic model (status in the formation formula) no longer stops due to missing r.num in init.net.
- Fix bug when calling netdx for one simulation only with a network model fit with the full STERGM method (i.e., using the edapprox = FALSE in netest).

### OTHER

- The dissolution_coefs function now stores and prints the death/exit rate.
- Removed explicit parameters for xlim, ylim, and main from plot.dcm and plot.icm functions, although they still may be passed through the ... argument.
- Introduction vignette has been updated to list new EpiModel website.

<br>

## EpiModel 1.0.1

### NEW FEATURES

- Added coef.form argument to netest for network model formation formulas with offset terms.
- Allow edge duration of 1 in netest when using the edges dissolution approximation (handles one-off partnerships in network models when using the approximation).
- Death modules for network models are now contained in one function, deaths.net, to facilitate replacement death modules from users. This is also now consistent with the death module for ICMs.
- Automated plotting of target statistic lines to plot.netsim formation plots, matching the methods of plot.netdx formation plots.

### BUG FIXES

- netdx now simulates from a different starting network at the begining of each dynamic simulation, eliminating correlation at time 1 across simulations.
- Can now pass status.vector into init.net for bipartite simulations.
- Several plotting and printing bugs fixed.
- Fixed bug in network models for open populations in which an attribute was passed to the network in the formation formula (e.g., serostatus mixing models).

### OTHER

- Added internal test structure for build checking.
- Added a help file document for building ICM modules at ?modules.icm.
- Expanded and clarified tutorial documentation, available at: <https://www.epimodel.org/>

<br>

## EpiModel 1.0

- Model parameterization for all model classes has been substantially revised to improve organization and ability for expansion. Whereas previous models required input of parameters directly into the main functions (now: dcm, icm, and netsim), now the parameters are input into three parameter-processing functions: param, init, and control. The param function sets the core epidemic parameters, the init function sets the initial conditions, and the control function specifies other model settings. These functions are class-specific, so each function has a .dcm, .icm, or .net suffix.

- Modeling functions have been renamed for clarity and consistency:

 - dcm is now used for deterministic compartmental models (replaces epiDCM)
 - icm is now used for stochastic individual contact models (replaces epiICM)
 - netest is now used for network model estimation (replaces epiNet.est)
 - netsim is now used for network model simulation (replaces epiNet.simTrans)

- Network models with independence between epidemic/demographic processes and network structures (independent models) were previously first simulated with epiNet.simNet, and then those pre-simulated networks were input to epiNet.simTrans. Now the network model simulation is all handled within the simulation function, netsim.

- Network model diagnostics have been moved from within the network estimation process (netest) to their own function: netdx. The parameter names for running, printing, and plotting the results of these diagnostics have been updated for consistency. See ?netdx and related functions.

- Internal model functions have been significantly revised to improve efficiency.

- The dcm function can handle model functions, parameter sets, and initial conditions of arbitrary complexity. See the HTML vignette on this topic at: <http://www.epimodel.org/>

- Moved the package vignettes external to the package to reduce package size and build time.

<br>

## EpiModel 0.95

### INITIAL RELEASE

- The EpiModel package provides functions for building, solving, and plotting mathematical models of infectious disease.

- See the main package help function ?EpiModel-package, and the EpiModel tutorials online.
