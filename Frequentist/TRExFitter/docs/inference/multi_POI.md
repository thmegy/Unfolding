# Multiple POIs

`TRExFitter` supports multiple POIs (parameter of interests) in a single config.
This option is very intuitive to use, just provide a comma separated list of names for the in the `POI` option of the `Job` setting.
However, you need to keep in mind that many settings are also tied to the POIs, e.g. `POIAsimov`.
The general rule for these settings is that you need to provide as many parameters e.g. for `POIAsimov` as you provide for `POI` and their order also has to match!
Do not worry, the code will complain when it encounters some inconsistencies.

So what is the benefit of using mutiple POIs?
The reasons is twofold: sometimes, you don't want to treat some NormalisationFactors as nuisance parameters, and you can also evaluate e.g. ranking on two POIs simultaneously (it's really simultaneous in the code).

!!! warning "Limit and significance"
    In case of limit and significance calculation, there is a strict requirement of using only single POI, as it is needed to construct the test statistics. This is done by specifying `POI` in the limit and significance blocks. This does not need to be set when only one POI is set in the `Job` block.

!!! hint "POIs and unfolding"
    Note, that in the case of `Unfolding` fit, normalisations of all truth bins are automatically added as POIs.

## Hands-on example
To demonstrate the usefulness of the multiple POIs, let us have a look at one of the jobs that are run in the CI.

Open the config file in `test/configs/FitExampleMorphing.config` and check the setup.
This config file uses two NormFactors, one is for the top width (which uses `Morphing`, but that is not relevant here), the second one is the ttbar cross-section.
Only the `topWidth` parameter is set as the POI.
Now, change the POI to include also the cross-section in the `Job` block:

```bash
  POI: "topWidth","ttbarXsec"
```

Now, run the usual `w` and `f` steps:

```bash
trex-fitter wf test/configs/FitExampleMorphing.config
```

As you can see, there is nothing special about the output, since the multiple POIs do not affect the simple fit.
Now, try to run the full ranking:

```bash
trex-fitter r test/configs/FitExampleMorphing.config
```

You should be able to see _two_ ranking files as an output, `Ranking_topWidth` and `Ranking_ttbarXsec`, each one showing the impact on the respective parameter.
This can be extended to any number of POIs!

!!! warning "Impact of POI1 on POI2"
    By default, when two POIs are set, ranking will _not_ show the impact of POI1 on POI2 and vice versa. However, you can enable this by setting `UsePOISinRanking: TRUE` in the `Fit` block.
