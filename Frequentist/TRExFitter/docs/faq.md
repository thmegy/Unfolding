# FAQ


## General questions

???+ question "The `n` step takes a very long time to run, how do I speed this up?"
    Run multiple jobs, with each job only processing all histograms for one region.
    For a region called `RegionA`, this is achieved by running

    ```bash
    trex-fitter n your.config Regions="RegionA"
    ```

    Run one job per region, then proceed with `w` only when all jobs are finished.
    Multiple regions per job can be processed via `Regions="RegionA,RegionB"`.

    It is possible to split up the histogram creation into even more jobs, splitting up systematics or samples.
    See section "Input File Merging with hupdate" in the README for more details, in this case an extra step is required before proceding with `w`.

    The runtime of the `n` step also depends on the amount of events in the ntuple(s) that are being read for each region.
    Events in the ntuple(s) that do not pass the region selection requirements further increase the execution time.
    A pre-selection or splitting of the input events across several files (e.g. one file per region) can decrease the execution time.

???+ question "Why is the E in TRExFitter capitalized?"
    `TRExFitter` stands for "Top Related Experiment Fitter".


## Statistical model building

???+ question "Can I just use more bins to gain sensitivity?"
    While more bins generally increase sensitivity, it is important to keep two things in mind.
    Limit the size of the MC statistical uncertainties ("gammas") to at most 20% per bin.
    Larger values can bias the signal extraction, see [these slides](https://indico.cern.ch/event/615262/contributions/2484815/) for a study and the corresponding recommendation.
    Furthermore, one aspect that is not considered in the fit is the statistical uncertainty in the templates that describe the ±1σ variations of a nuisance parameter.
    This uncertainty increases when adding more bins, and can lead to fluctuations in the templates.
    The nuisance parameter then does not describe physical effects, but rather fluctuations.
    It is very important to closely study all templates found in the `Systematics/` folder and verify that the input distributions to the fit (solid lines) look reasonable.
    These plots are obtained by enabling the settings `SystControlPlots` and `SystErrorBars` (both enabled by default).

???+ question "Can we insert correlations for NPs?"
    No, this is not possible as the fit finds the correlations as described in [Minimisation procedure](advanced_topics/likelihood.md#minimisation-procedure).
    Also this is problematic due to the split of the effects.
    All we can do it apply 100% correlation (fully correlating two NPs) - this can be done by providing same `NuisanceParameter` string for NPs that you want to fully correlate.

???+ question "Can we insert correlations for NFs?"
    This is possible!
    There is a difference between correlating NPs and NFs as NFs do not have an associated constraint term.
    Have a look at `Expression` option in the settings.
    This can also be used to calculate some other properties and not only signal strengths (like ratio of signal strengths, fitting fractions directly, etc.)

???+ question "How should I fix empty bins?"
    `TRExFitter` fixes empty bins for you by adding some arbitrary (small) values, otherwise the fit would not work at all.
    *But this is not a good solution!* You should find a proper solution - merging small backgrounds, rebinning to get rid of zero/negative bin contents.

???+ question "What should I do with systematic uncertainties that have both up and down shift in the same direction and I want to symmetrise them?"
    `TRExFitter` will print warnings when this happens.
    Problem is that the symmetrisation procedure does not make sense when this happen.
    We provide two alternative symmetrisation options: `Symmetrisation: ABSMEAN` or `Symmetrisation: MAXIMUM`.
    However, there is no technical way how to fix this problem and these options are not enough.
    So you need to adjust the inputs, unfortunately.


## Fit results

???+ question "There is some problem with the fit, what can I do?"
    Have a look at the [fit problems tutorial](fit_problems/problems.md).

???+ question "I am getting warnings about underconstrained nuisance parameters, should I worry?"
    Have a look at the relevant [fit problems tutorial](fit_problems/problems.md).
    Typically an underconstraint indicates an issue, though there can be cases where it is expected.
    See [this notebook](https://cernbox.cern.ch/index.php/s/DiPdvBlRBQPfHEy) for such an example.

???+ question "Why do I get different results every time I run identical setup?"
    This should not happen.
    It can happen only when you use `SetRandomInitialNPval` as the fit will start from a different point each time and thus the result can vary slightly (rule of thumb - the difference should not be larger than 0.01 on the POI).


## Algorithms

???+ question "How do the smoothing algorithms work?"
    Some info can be found in these [slides](https://indico.cern.ch/event/691683/contributions/2873279/attachments/1593521/2522846/PruningSmoothing.pdf).
    More information about `TTBARRESONANCE` is found in these [slides](https://indico.cern.ch/event/669913/contributions/2769795/attachments/1549339/2433688/ttres-fullunblind-smooth2-summary2.pdf).
    Further info regarding `TTBARRESONANCE`: The Smoothing parameter in the Systematics area can be set to 40 to treat the systematic uncertainty as correlated with the nominal (e.g. when obtained via reweighting) or 400 to treat it as uncorrelated with the nominal (e.g. for two-point systematics, then the statistical uncertainties on nominal sample and systematics variation are added in quadrature and compared to the smoothing threshold).

    A description of `TCHANNEL` can be found [here](https://gitlab.cern.ch/atlas-phys/exot/CommonSystSmoothingTool/uploads/466f09c1bb0272f5782bcd39ad3d9f7e/smoothing.pdf). This algorithm does not depend on the statistical uncertainties of the templates, which is useful when the degree of correlation of the statistical uncertainties between the two samples is not known. Additionally the algorithm is designed not to underestimate slopes at the start and end of the histograms. Shape effects that only affect one bin will be removed and finer binning leads to weaker smoothing.

???+ question "How do the automatic binning algorithms work?"
    The algorithms for `TransfoD` and `TransfoF` can be found in [these slides](https://indico.cern.ch/event/455289/contributions/1953694), `TransfoJ` is found in [slides here](https://indico.cern.ch/event/472696/contributions/1992693/). See also [this thesis](https://cds.cern.ch/record/2296985/), section 5.3.1.
    In practice, the `TransfoD` algorithm was found to work well, and a popular setting for the two parameters is to have both equal to the same integer.
    The amount of bins in the distribution is then equal to the sum, i.e. `"AutoBin","TransfoD",4,4` will create a distribution with 8 bins.


## Docker image

???+ question "How do I access files on EOS via the TRExFitter image?"
    You need to install `xrootd`, which is not available by default within the image.
    See [this question in the RECAST forum](https://atlas-talk.web.cern.ch/t/how-do-i-use-xrootd-with-the-trexfitter-image/104).
    See also [this gitlab issue](https://gitlab.cern.ch/atlas-amglab/atlstats/-/issues/16).

???+ question "Where can I learn more about docker?"
    There are many resources available online.
    See for example [this tutorial](https://www.docker.com/101-tutorial), the [official documentation](https://docs.docker.com/get-started/), and [an introduction that also uses gitlab CI](https://matthewfeickert.github.io/intro-to-docker/).


## Crashes

???+ question "Opening a workspace results in a segmentation fault without additional information, what can I do?"
    If you have a particularly complex workspace, you might run into a stack size limitation. Try running this command:

    ```bash
    ulimit -S -s unlimited
    ```

    It may solve the segmentation fault.

???+ question "TRExFitter crashes, what can I try to solve this?"
    Try re-compiling `TRExFitter` from scratch: delete the `build/` folder and re-compile (see [Compiling the code](setup.md#compiling-the-code)).
