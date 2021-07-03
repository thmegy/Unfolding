# Minimisation of the likelihood

A vital part of the "fit" is the minimisation of the (negative logarithm of) the likelihood.
This is internally done by calling [Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html).
The minimisation and estimation of parameters is done in multiple steps.
The first step uses MIGRAD algorithm to find the minimum of the likelihood function.
The second step includes calling HESSE method that provides more precise estimate of the uncertainties than MIGRAD, but reports only symmetric uncertainties, by construction.
The last (optional) step runs MINOS algorithm to get the most precise estimate of the uncertainties of the parameters by scanning the likelihood, and generally results in asymmetric uncertainties.

As can be seen, the crucial point is to find the correct minimum by MIGRAD, as other steps depend on the MIGRAD estimate.
Minuit allows different minimisation strategies within MIGRAD, denoted by integer numbers 0, 1 and 2.
The smaller the number, the faster and less precise the algorithm is.
By default, in TRExFitter, the first strategy that is tried is strategy 1 and up to three retries are tested if the previous strategy failed.
If after the retries Minuit still reports convergence problems, the fit is marked as failed and reported properly.

## TRExFitter options
Several options in the minimisation, and likelihood creation in general, can be set via a config file.

The first option is `BinnedLikelihoodOptimisation` in the `Fit` block, which turns-on internal optimisation in RooFit that leads to significant (about 10x) faster minimisation compared to the default.
The option is set to `FALSE` by default, due to numerical instabilities reported in this [ticket](https://root-forum.cern.ch/t/unexpected-behaviour-of-turning-on-off-binned-likelihood-fit-in-the-context-of-reordering-samples/39983?u=philtk74).

Users can also choose the starting minimisation strategy via `FitStrategy` in the `Fit` block.
Note that when you choose strategy 2, only up to two retries will be tested.

Last option is `UseHesseBeforeMigrad` as well in the `Fit` block, which tells the code to run HESSE after a _failed_ MIGRAD minimisation, before the next retry is tested.
This can improve the next iteration of the minimisation, but it requires some time to run.
By default, this option is turned on.
