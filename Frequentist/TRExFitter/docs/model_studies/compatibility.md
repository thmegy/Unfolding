# Evaluating compatibility

Analyses frequently measure a parameter of interest (POI) using data in different channels / phase space regions.
If the best-fit results for a POI disagree between channels / regions, it can be useful to evaluate the compatibility.
This can be done by using multiple POIs in a fit, and letting the POIs act on different channels or regions.
For example, a fit measuring the ttH signal strength (measured cross-section divided by the SM prediction) may use data in the l+jets and dilepton channels, and ultimately report the combined result.
To evaluate the compatibility between channels, one normalization factor can be assigned to ttH in the l+jets regions and another normalization factor to ttH in the dilepton regions.
A likelihood ratio test can then be used, comparing the likelihood in the setup with two normalization factors to the likelihood in the setup with a single normalization factor acting on both channels.
The default setup with a single normalization factor corresponds to our null hypothesis, which states that data in both channels can be described with a single normalization factor.
If data in both channels prefers very different normalization factors, we can reject the null hypothesis.

The likelihood ratio (with a factor $-2$ multiplied in) asymptotically follows the $\chi^2$ distribution (Wilk's theorem) with a number of degrees of freedom equal to the difference in degrees of freedom between the two likelihoods, and assuming the null hypothesis is true.
To obtain a p-value quantifying compatibility with the null hypothesis, we thus need to calculate the likelihood for the two hypotheses (single normalization factor, two normalization factors).
This is described with an example below.

This method is not exclusive to signal strength measurements, but can be used more generally with free-floating parameters.
Another example can be found in this [ATLAS Higgs mass measurement paper](https://arxiv.org/abs/1307.1427) in section 7.
Instead of using two normalization factors scaling a signal strength like above, the problem can also be set up by using a single normalization factor scaling the signal strength in one channel, and then scaling another channel by the sum of the first normalization factor and a second one.
The second normalization factor describes the offset between channels, and its compatibility with zero can again be tested.

The method also works for evaluating the compatibility of more than two measurements: instead of comparing the likelihood between models with one or two normalization factors, more generally the likelihood can be compared for models with any number of normalization factors.
The difference in the number of normalization factors is the number of degrees of freedom to be used in the $\chi^2$ test.

## Example: ttH

We will use the example config file provided in the repository under `config/ttH_tutorial.config` to demonstrate the compatibility calculation.
There is a single normalization factor `mu_ttH` scaling the ttH signal strength in all regions.
This example has three regions: two l+jets regions `ljets_HThad_5j3b` and `ljets_HThad_ge6jge4b`, and a dilepton signal region `dilep_HThad`.

### Single signal strength

To start, we create a workspace and run a maximum likelihood fit:
```bash
trex-fitter nwf config/ttH_tutorial.config
```

Make sure your `DebugLevel` is set to a value of 1 or larger, and the following output is visible in the fit results:
```txt
=== INFO::FittingTool::FitPDF: ***********************************************************
=== INFO::FittingTool::FitPDF:   Final value of the NLL = -7877.695420
=== INFO::FittingTool::FitPDF: ***********************************************************
```

This is the NLL for our null hypothesis.
The best-fit signal strength is
```txt
mu_ttH  2.87996 +1.54354 -1.19859
```

### Two signal strengths

We now modify our fit model to have two normalization factors, one for l+jets channels and one for dilepton channels.
This can be achieved with the `Regions` option for normalization factors:
```yaml
NormFactor: "mu_ttH_lj"
  Title: "#it{#mu}(#it{t#bar{t}H}) l+jets"
  Nominal: 1
  Min: -10
  Max: 20
  Samples: ttH
  Regions: ljets_HThad_5j3b,ljets_HThad_ge6jge4b

NormFactor: "mu_ttH_dil"
  Title: "#it{#mu}(#it{t#bar{t}H}) dilepton"
  Nominal: 1
  Min: -10
  Max: 20
  Samples: ttH
  Regions: dilep_HThad
```

Update your config to have these two normalization factors, and also update the `POI` in the `Job` block at the top of the config.
It does not matter which of the normalization factors you put there.
Then run the fit after building a new workspace (required to have the normalization factors available in the model):
```txt
=== INFO::FittingTool::FitPDF: ***********************************************************
=== INFO::FittingTool::FitPDF:   Final value of the NLL = -7877.878194
=== INFO::FittingTool::FitPDF: ***********************************************************
```

The best-fit results for the two signal strengths are
```txt
mu_ttH_lj   1.60839 +2.61885 -3.00015
mu_ttH_dil  3.70343 +2.23008 -1.82930
```

### Compatibility between the results

`TRExFitter` provides a small script to help convert the difference in NLLs into a p-value.
It can be found at `util/CalculateCompatibility.py`.
The script expects the NLL for the fit with a signal normalization factor first, then the NLL for the second fit, and finally the number of degrees of freedom for the $\chi^2$ test.
In this case, the difference in degrees of freedom between the two fits is 1, so we use that value.
```bash
python util/CalculateCompatibility.py -7877.695420 -7877.878194 1
```
results in
```txt
The compatibility is 54.54410472689661%.
```

With a p-value this large, we cannot claim incompatibility between the two signal strengths.
This is consistent with a more naive calculation by considering the best-fit values and associated uncertainties: the two best-fit results considerably overlap within their associated uncertainties.
