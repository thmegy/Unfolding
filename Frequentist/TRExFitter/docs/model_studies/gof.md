# Goodness-of-fit calculation

## Saturated model

Whenever a fit procedure is used, it is important to check the Goodness-of-fit (GoF) status.
GoF is a metric that quantifies how well the fit model describes the observed data.
If the GoF assigns a very small probability, the model should be checked.

For a long time a very ad hoc GoF test was implemented in `TRExFitter`, which compared the likelihood values for the fit to data and a fit to Asimov dataset.
This is obviously not a proper test.
The proper test is provided by [saturated model](http://www.physics.ucla.edu/~cousins/stats/cousins_saturated.pdf).
In this test, the likelihoods of two fits are compared.
One is obtained by fitting the data with the nominal model.
The other one is obtained by fitting the same data with the so-called _saturated model_, a model that has enough freedom that it will fit the data perfectly, without requiring nuisance parameter pulls that would result in likelihood penalties due to the associated constraint terms.
In other words, the saturated model is a modified model with sufficient expressiveness to exactly describe data.
The ratio of these two likelihoods follows the $\chi^2$ distribution asymptotically (Wilks theorem), and thus can be used in a standard GoF test.
To further compare it to the $\chi^2$ test, the saturated model represents $\chi^2 = 0$ (perfect agreement).
It can also be viewed as the constant term $c$ that is removed from the full likelihood
$$
 -2\ln L(\mu) = \chi^2(\mu) + c
$$

### TRExFitter implementation

The saturated model as GoF test has been implemented in `TRExFitter` since tag `TtHFitter-00-04-05`.
Since tag `TRExFitter-00-04-08` it is the default option when `GetGoodnessOfFit` is set to `TRUE`.
Technically, it is implemented by assigning _shape factors_ — "normalisation factors per bin" — to each bin to allow the model to fit the data perfectly without the need to pull any constrained nuisance parameter.
The minimisation procedure is then run to calculate the the absolute likelihood value that is then used in the GoF calculation.
There is also a setting `SaturatedModel` to build workspaces with those shape factors independently of whether `GetGoodnessOfFit` is used.

Let us try to use this option now.
We will use a configuration file used in our CI tests.
First produce the inputs

```bash
trex-fitter n test/configs/FitExampleNtuple.config
```

And now, run the fit (and also create the workspace first)

```bash
trex-fitter wf test/configs/FitExampleNtuple.config
```

We get the fitted values and everything looks fine.
Now, modify the `Fit` block of the config file and add the following option: `GetGoodnessOfFit: TRUE` .
Since we are using a very recent version of `TRExFitter`, this will by default use the saturated model as the GoF test.
Now, run the workspace creation and fit again

```bash
trex-fitter wf test/configs/FitExampleNtuple.config
```

You should now see that the fit is run _twice_.
First the standard fit is run, then the fit with the saturated model is run.
You should see that in the second step no pulls are present (only the shape factors will be used to fit data).
Check the lines that print the likelihood values, e.g.:

```bash
=== INFO::FittingTool::FitPDF:    -> Reduced Final value of the NLL = -998511.89450029970612376928
=== INFO::FittingTool::FitPDF:    -> Final value of the NLL = 1488.105500
=== INFO::FittingTool::FitPDF:    -> Final value of offset = 1446.032751
=== INFO::FittingTool::FitPDF:    -> Final NLL - offset = 42.072749
```

```bash
=== INFO::FittingTool::FitPDF:    -> Reduced Final value of the NLL = -998516.03562775580212473869
=== INFO::FittingTool::FitPDF:    -> Final value of the NLL = 1483.964372
=== INFO::FittingTool::FitPDF:    -> Final value of offset = 1446.032751
=== INFO::FittingTool::FitPDF:    -> Final NLL - offset = 37.931622
```

And finally, the GoF value is printed

```bash
=== INFO::TRExFit::PerformFit: ----------------------- GOODNESS OF FIT EVALUATION -----------------------
=== INFO::TRExFit::PerformFit:   NLL0        = 1483.964372
=== INFO::TRExFit::PerformFit:   NLL         = 1488.105500
=== INFO::TRExFit::PerformFit:   ndof        = 10
=== INFO::TRExFit::PerformFit:   dNLL        = 4.141127
=== INFO::TRExFit::PerformFit:   2dNLL/nof   = 0.828225
=== INFO::TRExFit::PerformFit:   probability = 0.601288
=== INFO::TRExFit::PerformFit: ----------------------- -------------------------- -----------------------
```

This example uses 10 degrees of freedom in the p-value calculation.
The number is automatically determined, and given by the number of bins minus the number of unconstrained parameters in the nominal fit.
The example uses 12 bins and 2 normalization factors, resulting in 12-2=10 degrees of freedom.

!!! warning "Using workspaces with saturated model parameters"
    The HistFactory workspaces produced by `TRExFitter` when using `SaturatedModel: TRUE` or `GetGoodnessOfFit: TRUE` contain the saturated model shape factors.
    Those shape factors are **NOT** set to constant in the workspace, since this is not supported by HistFactory.
    When fitting such workspaces with tools other than `TRExFitter`, please take care to setting those parameters to constant if intended.
    The parameters are called `gamma_saturated_model_sf_*`.
