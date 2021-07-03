# Command line options

There are a few command line options which are very useful when building statistical models with `TRExFitter`.
Take a look at the [list of options in the README](https://gitlab.cern.ch/TRExStats/TRExFitter#command-line-options).
We will outline a couple of them here.

## Split by region for n-step parallelization

The `Regions` command line option can be used to build workspaces only containing specific regions, or produce histograms only for a region or a list of regions.
This is very useful for histogram production in the `n` step in `TRExFitter`.
If you enable `SplitHistoFiles`, the `.root` file containing the histogram inputs for the workspace are split by region.
You can take advantage of that and produce these histograms in parallel, for example with one region per job:
```bash
trex-fitter n example.config Regions="region_1"
```
and in a second job
```bash
trex-fitter n example.config Regions="region_2"
```

When you have many regions, this parallelization can significantly speed up your `n` step.

It is possible to go one step further and also split systematics into different jobs.
This approach is described [here in the README](https://gitlab.cern.ch/TRExStats/TRExFitter#input-file-merging-with-hupdate).

## Useful options for workspace building

The command line options can also be useful when constructing a workspace.
You can use the `Exclude` option for example to exclude a region:
```bash
trex-fitter wf example.config Exclude="region_1"
```

Similarly, you can use it to exclude samples or systematics.
This can be convenient when you have many different signal models defined in your config, and want to build workspaces with one signal each.
