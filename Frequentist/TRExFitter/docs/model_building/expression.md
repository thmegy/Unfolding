# Using the `Expression` config option

## Introduction

It is technically not possible to add correlations of the NPs to the likelihood by hand.
Correlations of the NPs are estimated from the fit during the likelihood minimisation step.
However, it is possible to correlate the normalisation factors with other normalisation factors.
This allows to modify the standard fits to non-standard ones, e.g. fraction fitting where the sum of some normalisation has to add up to one.
These correlations can be set via the `Expression` option.

!!! hint "Expression and morphing"
    The underlying functionality, the `AddPreprocess()` RooFit function, is used for the `Expression` part as well as for the fits with multiple templates.


## Example: W helicity
As an example use case for the `Expression` option is the W helicity measurement.
This measurement fits the fractions of the individual helicity templates (pure left-handed, pure longitudinal and pure right-handed), however, the fractions are not independent.
The fractions have to satisfy

$$
F_{L} + F_{0} + F_{R} = 1
$$

for all allowed values of the fractions.
This leaves only two independent fractions.
The fit setup then needs three NormFactors, one for each fraction, but one of the NormFactors needs to be fixed from the other ones.

Have a look at the config file in `test/configs/FitExampleExpression.config` and focus on the line withe the Expression in the `NormFactor: "norm_left"` block

```bash
 Expression: (1.-norm_long-norm_right):norm_long[0.687,0,1],norm_right[0.002,0,1]
```

In this line, we tell the code to replace the `norm_left` normalisation parameter with `1 - norm_long - norm_right` and then set the intial and minimum/maximum values for the normalisation factors.

!!! hint "Note"
    The formula part of the expression can be anything that can be read by `TFormula`.

Now produce the histograms and run the fit

```bash
trex-fitter hwf test/configs/FitExampleExpression.config
```

You will see that the `norm_left` disappeared from the results completely, but that is expected since this parameter no longer exists in the likelihood.

!!! hint "Flexibility"
    The `Expression` functionality provides high level of flexibility for different kind of measurements. Do not be afraid to experiment with it.

## Example: charge asymmetry

Imagine we want to measure the top quark charge asymmetry $A_C = \frac{\sigma^+ - \sigma^-}{\sigma^+ + \sigma^-}$.
We might have two samples `ttbar_plus` and `ttbar_minus` that are scaled by two normalisation factors `sigma_plus` and `sigma_minus`:
```yaml
NormFactor: "sigma_plus"
  Title: "sigma +"
  Nominal: 1
  Min: 0
  Max: 10
  Samples: ttbar_plus

NormFactor: "sigma_minus"
  Title: "sigma -"
  Nominal: 1
  Min: 0
  Max: 10
  Samples: ttbar_minus
```

We could use those in a fit, get the best-fit result for both normalisation factors and calculate the asymmetry, plus an approximation of its uncertainty via error propagation.
What if we would like to calculate the asymmetry directly, and get its uncertainty without any assumptions due to error propagation?
The asymmetry then needs to be a parameter in the fit.
We can re-parameterize the setup to have the asymmetry be one of the two degrees of freedom, and pick either `sigma_plus` or `sigma_minus` as the other degree of freedom in the fit.
The example below choses `sigma_minus`.
Since we still need to know the value of `sigma_plus`, which scales the `ttbar_plus` sample, we determine it via an expression:
```yaml
NormFactor: "asymmetry"
  Title: "asymmetry"
  Nominal: 0
  Min: -10
  Max: 10
  Exclude: ttbar_minus,ttbar_plus

NormFactor: "sigma_minus"
  Title: "sigma -"
  Nominal: 1
  Min: 0
  Max: 10
  Samples: ttbar_minus

NormFactor: "sigma_plus_from_expr"
  Title: "sigma +"
  Expression: minus_NF*(asymmetry+1)/(1-asymmetry):minus_NF[1,0,10],asymmetry[0,-10,10]
  Nominal: 1
  Min: 0
  Max: 10
  Samples: ttbar_plus
```

Note that `asymmetry` does not scale anything sample in this setup, but it is (together with `sigma_minus`) one of the two degrees of freedom in the fit.
It only acts indirectly via its effect on `ttbar_plus` through the `sigma_plus_from_expr` normalisation factor.
The `sigma_plus_from_expr` normalisation factor on the other hand is not a free parameter in the fit, but can be determined at every step during minimisation from the two other normalisation factors.

If we do not specify the `Exclude` setting, the `asymmetry` normalisation factor will by default to apply to all samples.
This is not the behavior we intend here, so we explicitly exclude samples (or we could exclude all regions) so that it does not scale anything directly.

The best-fit result for the `asymmetry` parameter in this parameterisation should be equivalent to the manual calculation of the parameter from the best-fit results of `sigma_plus` and `sigma_minus` in the alternative parameterisation (within minimisation tolerance).

!!! tip "Tip"
    If you would like to quote uncertainties for `sigma_plus`, `sigma_minus`, and `asymmetry` without error propagation, you run into a problem: you can only get two out of the three from a single fit.
    One possibility is fitting with two parameterisations.
    If you go that route, you should ensure that the two minima are sufficiently similar: if the minimisation ends up in different local minima, the numbers you will obtain from the two fits are not directly comparable!

## Example: ratios

We conclude with a brief example of a common scenario for which no expressions are needed.
If you want to measure the ratio of two normalisation factors, a suitable parameterisation allows you to implement everything without `Expression` usage.
Here we want to measure the ratio between a signal strength for di-muon processes and di-electron processes.
We could set this up by having one normalisation factor each scaling di-muon and di-electron samples, but then would have to resort to error propagation to get the error on the ratio.
Instead, we use the ratio of di-electron to di-muon normalisation as our first degree of freedom, and the di-muon normalisation factor as the second one.
Then we scale all di-muon samples by the di-muon normalisation factor.
The di-electron samples are scaled by both the ratio and the di-muon normalisation factor, which causes the di-muon contribution to cancel out.
We are thus left scaling the di-electron samples by the di-electron normalisation as desired:
```yaml
NormFactor: "dimuon_NF"
  Title: "di-muon NF"
  Samples: mumu_sample_1, mumu_sample_2, ee_sample_1, ee_sample_2

NormFactor: "ratio"
  Title: "ratio di-electron / di-muon"
  Samples: ee_sample_1, ee_sample_2
```

!!! question "Using expressions"
    How would the above setup look like if we want the di-electron to di-muon ratio to be one degree of freedom in the fit, and the di-electron normalisation factor to be the second degree of freedom?
    In this case you need to make use of expressions again.
