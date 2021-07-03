# Fixing and injecting NPs

In many analyses, it is important to check the fit model by injecting or fixing some of the NPs in the likelihood and repeating the fit.
`TRExFitter` provides multiple convenient options to help do this.

## Fixing NPs

To fix NPs in the model, use `FixNPs` option, the syntax is
```yaml
FixNPs: NPname1:value1,NPname2:value2
```

where `NPname` is the name of the NP, and `value` is the target value for the NP in "sigmas" (also negative value can be set).
Once you fix a NP, it no longer appears in the likelihood as a parameter that can be fitted and thus will not show up in the results.
The parameter will be held constant at whatever value you set it to.
You do not need to reproduce the workspace when using this option, you can directly run `f` (assuming a workspace already exists).

## Shift NPs to build pseudo-dataset

By default, the Asimov dataset is built with all NP values set to 0.
However, you may need to test the impact of some of the NPs by shifting them to build an alternative pseudo-dataset.
This can be done by using `NPValues` option in the `Fit` block.
The syntax is:
```yaml
NPValues: NPname1:value1,NPname2:value2
```

The dataset will be built from the prediction of your model when setting those parameters to the values provided, and the remaining parameters to their nominal values.
Note that when this is set, the automatic data/MC fit (based on `Region`->`DataType` option) is disabled.

You can also inject the values from a `.txt` file produced by the `TRExFitter` `f` step.
You need to use the `NPValuesFromFitResults` `Fit` option for this and pass a path to the `.txt` file.

!!! warning "Nuisance parameter names"
    When injecting NP values, the name to use is **NOT** the name you set for systematics in your config.
    Instead, you have to use the name used internally in HistFactory.
    For systematics, that means adding an `alpha_` prefix: the parameter attached to `Systematic: JER` is called `alpha_JER`.
    Normalization factors do not receive this prefix, and gammas keep their existing `gamma_` prefix without receiving another one.

!!! info "Example and pulls when fitting pseudo-dataset"
    When using the `NPValues` option, and then fitting the resulting pseudo-dataset, your fit will generally not pull the nuisance parameters to the same value you injected.
    We can give this a try with the example setup in `config/ttH_tutorial.config`.
    Start by adding `FitBlind: TRUE` to the `Fit` block and run the following:
    ```bash
    trex-fitter nwf ttH_tutorial.config
    ```

    As expected, when fitting the Asimov dataset there are no pulls present.
    Now inject a pull of the JER systematic via `NPValues: alpha_JER:1.5` and repeat the fit.
    You will see a slight pull of the JER NP (alongside a pull of the JES NP 1), however that pull is much smaller than 1.5.
    This behavior is normal.
    When fitting a dataset, the parameter setting maximizing the likelihood is the best-fit point.
    That means finding the optimal trade-off between stronger pulls (decreasing the likelihood) and better modeling of data (increasing the likelihood).
    In this case, a stronger pull does not improve the likelihood further.

## Injection of Global Observables

If you set `InjectGlobalObservables` to `TRUE` in the Fit block (default is `FALSE`), the NP shifts from `NPValues` or `NPValuesFromFitResults` will also shift the global observables in the likelihood.
Global observables represent the observed quantities in auxiliary measurements.
We model our constraints as auxiliary measurements, suitably normalized such that the observed data from them is usually 0 (for the normal types of systematics).
See the [HistFactory documentation](https://cds.cern.ch/record/1456844) for more information about this concept.
For example, when you implement a JER NP, then in the likelihood there will be a Gaussian constraint term $G(0 | \theta, \sigma=1)$.
$0$ here is the global observable: the nominal value coming from the b-tagging group, implemented in a simplified way in this Gaussian form.
$\theta$ is the nuisance parameter in the fit.
By using `InjectGlobalObservables`, this global observable will change, and hence the position at which the Gaussian is evaluated will change.
That means that the new maximum likelihood (when just considering this Gaussian by itself) is now at the value of $\theta$ matching the global observable.

In other words, if the NP is set to a non-zero value, fitting such dataset will result in a pull of the given NP that will come at a likelihood penalty.
Setting `InjectGlobalObservables` will remove this penalty.

When using `InjectGlobalObservables`, you turn the pseudo-dataset created with `NPValues` again into an Asimov dataset: the best-fit values for all NPs will maximize their associated constraint terms.

!!! info "Example continued"
    Continuing the example from above, add `InjectGlobalObservables: TRUE` to your fit block and re-run the fit.
    You will now recover the pull of 1.5 in the JER NP.
