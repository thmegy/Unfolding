# Unfolding using profile likelihood

## Basic idea

Standard profile likelihood fitting can be used for unfolding.
The key is to transform the unfolding problem into a standard problem of fitting normalisation of distributions, and then we can use the standard profile-likelihood machinery.
This transformation can be done very easily.
Consider a truth distribution (can be either on parton on particle level) with N bins and a corresponding "response matrix" of size N x M.
The response matrix here represents the combination of the selection efficiency, migration matrix and acceptance.
The distribution on the detector level with M bins is simply obtained by a matrix multiplication of the truth distribution with the response matrix.
Thus, for each bin of the truth distribution we can obtain one distribution with M bins on the detector level.
The last missing piece is to realise that the normalisation (signal strength) of these _folded_ distribution on the detector level are identical to the normalisation of the truth distribution.
Thus, measuring the normalisation of each folded distribution on the detector level (which is a "standard" fit), one directly measures the normalisation of the truth level - which is the goal of the unfolding.
For more details, check this [presentation](https://indico.cern.ch/event/890060/contributions/3754199/attachments/1991168/3320058/Unfolding_with_TRExFitter.pdf).

## Config file

All of the machinery that builds the folded distribution that are then fitted is incorporated inside TRExFitter, provided all the needed input histograms are available.
There is one slight conceptual difference compared to the regular (non-unfolding) fits, the signal samples need to be treated differently due to the folding procedure, while background samples are treated in the same way as in the regular fits.

The relevant truth distribution and response matrix need to be provided for the signal sample.
Note that for each systematic uncertainty, the response matrix has to be provided to generate (fold) the templates needed for systematic variations.
Alternatively, per-bin selection efficiency (probability to pass the reco-level selection for an event in a given bin of the truth-level distribution),
migration matrix and per-bin acceptance (probability to pass the truth-level selection for an event in a given bin of the reco-level distribution)
can be provided instead of the full response matrix.

!!! tip "Definition of response matrix"
    The response matrix is defined as selection efficiency times migration matrix divided by acceptance. All bin contents must be $\leq 1$.

Now, let us look at how this is done in `TRExFitter`.
We will use a simple example from the ttZ analysis.
Let us have a look at the [config file](https://gitlab.cern.ch/TRExStats/TRExFitter/-/blob/master/config/ttZ3l_example.config).

### Job block

Firstly, let us look at what is defined in the Job block of the config file.
The first thing that you should notice is that in the `Job` part, new path settings appeared: `MigrationPath`, `SelectionEffPath` and `AcceptancePath`.
As the names suggest, these just set the paths for the input files/histograms
The standard `TRExFitter` logic is used here that `Path`, `File`, `Name` are combined to get the relevant histograms.
The settings in `Job` can be overwritten by the same setting in `Region` which can in turn be overwritten by settings from  `UnfoldingSample` or `Unfolding Systematics`.

### Fit block

Another important difference to the standard config file is the `FitType: UNFOLDING` setting in the `Fit` block.
This tells the code that the main result of this configuration is unfolding and not a regular fit.

### Unfolding block
This is a completely new block used only in unfolding jobs.
The first setting `MatrixOrientation` simply tells the code what your definition of horizontal and vertical axes is (which axis represents the truth distribution).

`NumberOfTruthBins` just tells the code how many truth bins there are, this is mainly to allow some automatic cross-checks in the code, however, this needs to be set.

An important option is `NominalTruthSample`, which tells the code which of the defined `TruthSample` (see below) is the nominal one that is used for folding (others cane be used just for plotting purposes).

`UnfoldingResultMin` and `UnfoldingResultMax` set range of the normalisation parameters for the truth bins.

All the other options are just cosmetic options for the final plot, e.g. plotting the final distribution with yields divided by the bin width.

### TruthSample

These blocks define the truth distribution for the signal that is used throughout the procedure.
Multiple distributions can be defined, these will then be plotted in the final unfolded distribution, however only one truth distribution can be defined that will be used for the unfolding (even when multiple regions are used!).
Note that one sample needs to be defined as nominal using the `NominalTruthSample` setting in the `Unfolding` block (see above).
File/histogram paths are defined in this block as well as some of the cosmetic options.

### Region

The standard definition of regions is extended to define paths for the migration matrix, selection efficiency and acceptance.
Importantly, `NumberOfRecoBins` is used to define the number of reco bins (detector level), mostly for the sake of automatic tests (but needs to be defined!).
As you can see in the provided example, the config file uses two regions.
One signal region, where the unfolding is to be performed and one control region used to constrain background normalisation.

### UnfoldingSample

A new block is provided for the signal samples that will be used in the unfolding.
Their idea is similar to the regular Sample block.
The relevant paths for the acceptance, selection efficiency and migration matrix are defined here.
The information form the TruthSample is combined with the information from UnfoldingSample to build the folded truth distribution that are then used in the fit.
Note, multiple unfolding samples can be used at the same time, this allows to combine multiple detector level regions for unfolding.
In the current example, the UnfoldingSample is set only for the signal region, and in the control region a regular Sample is added for the ttZ.
This means that the ttZ contribution in the control region is treated as _background_, which is not correct, but is a good approximation if the signal contribution in the control region is small.

### Sample and NormFactor

Standard configuration is used for the _background_ samples.
NormFactors can also be attached to the _background_ samples.
The NormFactors used for unfolding are automatically added and you do not need to worry about them.

### UnfoldingSystematic

This new block is used for defining systematic uncertainties that act on the signal and need the folding step.
Similarly to standard systematics, `MigrationNameUp` can be defined (same for selection efficiency and acceptance) which will replace the nominal migration matrix when this uncertainty is being folded.
You can correlate the standard `Systematic` and `UnfoldingSystematic` using the same `NuisanceParameter` name as is done in the example config with `bTagSF_DL1r_Light_0`.
This will tell the code to use only one parameter to describe the systematic shift in the (folded) signal distributions and the background distributions.

## Preparing the fit

Now with the basic parts of the config file explained, try to run the (new) `u` step.

```bash
trex-fitter u ttZ3l_example.config
```

First, you will see that some plots are created in the `ttZ3l_PtLepNonZ` folder.
These plots show the migration matrix as well as the response matrix.
The folded distributions can be found in `ttZ3l_PtLepNonZ/UnfoldingHistograms/`.
You can inspect them.

## Fitting
Now that all the hard work is done to prepare the inputs needed for unfolding.
You can now run the usual steps:

```bash
trex-fitter h ttZ3l_example.config
```

Which will read the input histograms.
It will also read the folded distributions from the previous step without you needing to change anything!
Have a look at the `ttZ3l_PtLepNonZ/Systematics/` folder

Now we can make the workspace and run the fit

```bash
trex-fitter wf ttZ3l_example.config
```

Apart from the usual plots made by `TRExFitter`, you should also see a new plot called `UnfoldedData` that shows the unfolded distributions as well as all truth distributions defined in the config.

You can also produce the standard pre/post fit plots running:

```bash
trex-fitter dp ttZ3l_example.config
```

!!! tip "Running toys"
    It may be very useful to run pseudoexperiments to test the unfolding procedure. You can do this using the standard config option in `Fit` block: `FitToys: XXX` where `XXX` is the number of toys to be used. Note that for toys, the Asimov prediction will be used.

## Normalized differential cross-section

You can tell the code to produce the normalized differential cross-section plots setting `UnfoldNormXSec: TRUE` in the `Unfolding` block.
This is done technically in three steps.
Firstly, a new NormFactor for the total cross-section is added, $\mu_{tot}$.
Then, all the NormFactors $\mu_{i}$ are replaced with $\mu_{i} = \mu_{rel, i}\mu_{tot}$.
Finally, the code will replace one of the free parameters for the truth bin normalization with a formula that will take into account the total normalization that will be used for the normalized cross-section fitting. By default, the last bin of the truth distribution will be replaced, but you can change this with the `UnfoldNormXSecBinN` option.
For $N$ truth bins, where the normalisation of the bin $j$ is used for the replacement, the formula is

$$
\mu_{rel,j} = (1-\sum_{i, i\neq j}^{N}\mu_{rel,i}(N_i/N)) / (N_j/N),
$$

where $N_{i}$ is the initial number of events in bin $i$ and $N$ is the total number of events (sum of all $N_{i}$).

The fitted number of events (or cross-section) in the excluded bin is equal to the total minus the sum of the other bins, see below

$$
\mu_{rel,j} N_{j} = ( 1 - \sum \mu_{rel,i} N_{i}/N ) / ( 1/N ),\\
\mu_{rel,j} N_{j} = N - \sum \mu_{rel,i} N_{i},\\
\mu_{tot} \mu_{rel,j} N_{j} = \mu_{tot} N - \sum \mu_{tot} \mu_{rel,i} N_{i}.
$$

## Regularization

Running the fit out-of-the-box is equivalent to un-regularized unfolding.
However, `TRExFitter` provides a way how to apply regularization.
This is technically done by adding additional terms to the likelihood.
The regularisation (and its strength) is controlled by the `Tau` config option in the `Unfolding` block.
It can be set to a single (float) value that will be applied to all bins, or it can be applied per-bin with the following syntax: e.g. `Tau: 1:2,3:1.7` that will apply the constraint term to the first and third bin.
There are two types of regularization currently available, controlled by `RegularizationType` option in the `Unfolding` block.
It can either be set to 0 or 1.
When set to 0 (default) a bin-by-bin constraint, of the form $e^{-(1/2) \tau^2 (NF-\mu)^2}$, where $\mu$ is the Nominal value for the NF is used.
Thus, this biases the result towards the nominal prediction, and the bias is stronger the larger the tau value is.
Second regularization option, `RegularizationType: 1` uses discretized second derivative constraint.
Try to add these regularization option to the example config and rerun the `w` and `f` steps and see what happens with the results!

!!! tip "Uncertainty estimate with regularization"
    Generally speaking, when regularization is used, the standard way how to estimate the uncertainty of the parameters from the fit, either using the Hessian matrix or from scanning the likelihood is not the appropriate one. One must run the toys to estimate the uncertainty of the parameters. Futher studies are needed for regularized fits.
