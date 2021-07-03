Wrapper to perform profile likelihood unfolding with TRExFitter.

[[_TOC_]]

# Setup

Get the code and set it up with the following commands:
```
git clone ssh://git@gitlab.cern.ch:7999/thmegy/pl_unfolding.git
cd pl_unfolding
source install.sh <architecture>
```
where the architecture is for instance `x86_64-centos7-gcc8-opt`

<br>

# Run the code


## Input structure

The input folder is expected to have the following structure:
```
<input-folder>/<campaign>/<channel>/<variable>.root
```

The `root` files are expected to contain:
- the signal truth-level distribution
- the migration matrix
- the reco-level distribution of every background
- any data or pseudo-data distribution you plan on using (e.g. nominal and reweighted asimov datasets)
- a directory for every systematic variation, containing:
  - the systematically-varied migration matrix (if signal affected by the systematic)
  - the signal truth distribution, if different from nominal (if signal affected by the systematic)
  - the varied reco-level distribution of every background affected by the systematic
- a directory for every reference sample use for systematic estimation (e.g. signal sample generated with atlasfast), with the same content as for systematic variations

----

## The config file

The settings for the various steps of the data treatment are defined in a `json` config file.  
The following settings should be defined:
- luminosity corresponding to the MC campaigns
- background processes
- systematic uncertainties
- reference samples to compare systematic variations to (if different from nominal)
- options to be included in the TRExFitter configuration file

The default file defined in the different python scripts is [config/config.json](config/config.json).

----

## Rebin the input

If your input histograms are produced with a high granularity, you can rebin them with:
```
python scripts/rebin.py -i <path-to-input-folder> -v <variable> -c <channel> --campaign <campaign>
```
The desired binning, as well as the histograms you wish to rebin shoud be indicated in a json config file, whose default is [config/rebin_config.json](config/rebin_config.json). The path of this file should be given in the main config file described in the previous section.  

The nominal and systematically varied response matrices used in the unfolding are derived by this script from the provided migration matrices and truth distribution(s).  
The computation of the response matrices also depends on the systematics treatment chosen with the `-s` argument. There are two alternative approaches:
- The default approach is to account for both shape and normalisation uncertainties (`-s full`);
- Alternatively, by arguing that you are not interested in measuring the normalisation of the unfolded distribution, you can account only for shape uncertainties (`-s shape`). In this case it is ensured that for a given truth distribution, the nominal and systematically varied reco distributions have the same normalisation. This is achieved by reweighting the systematically varied response matrix.

This script also extracts the truth values of the spin parameter corresponding to the considered observable, which are needed for the linearity test and the mean calibration.

**Caveat**: If your input already has the desired binning and you want to skip this step, you then need to add to your input the nominal and systematically varied response matrices, which are otherwise computed by the rebin script.

----

## Run TRExFitter

Once your input is rebinned, your can generate the TRExFitter config file and run the unfolding with:
```
python scripts/runTRExFitter.py -i <path-to-rebinned-input> -v <variable> -c <channel(s)> --campaign <campaign> -d <dataset>
```
Several channels can be given simultaneously with the `-c` argument (e.g. `-c ee em mm`). The reponse matrix and the unfolded distribution can be plotted by adding the `-p` argument.  

Separate measurements have to be performed in order to get the unfolded distribution and the spin parameter with a correct error estimation. This is set by the `-m` argument, which can take the `param` (default) and `distrib` values.  
The output of the fit, and input for the EFT interpretation, is saved to a `json` file called `<variable>_unfolded_distribution.json` or `<variable>_unfolded_parameter.json`, depending on the type of measurement you perform.

----

## Linearity test and mean calibration

The linearity of the unfolding response can be assessed by running the unfolding with various signals reweighted to inject a chosen spin parameter. This can be performed with:
```
python scripts/runTRExFitter.py -i <path-to-rebinned-input> -v <variable> -c <channel(s)> --campaign <campaign> -d <list-of-datasets>
```

Once the previous step is done, the linearity test itself can be performed with:
```
python scripts/linearity_test.py -r <list-of-datasets> -v <variable> -c <channel(s)> --campaign <campaign>
```

The mean-calibration curve is also draw by the above script (True parameter vs truth-level parameter computed from binned histogram).

----

## Systematics impact ranking

The impact of individual systematic uncertainties on the measured spin parameter(s) can be measured with:
```
python scripts/makeSystImpact.py -i <path-to-rebinned-input> -v <variable> -c <channel(s)> --campaign <campaign> -d <dataset>
```
The impacts are also ranked and a summary plot is drawn. If for some reason the ranking job has not run, the ranking can be performed locally with:
```
python scripts/makeSystImpact.py -i <path-to-output> -v <variable> -c <channel(s)> --campaign <campaign> -d <dataset> --rank-plot
```
where the path to the output of the previous step should be given with the `-i` argument.

----

## Uncertainties breakdown

The total uncertainty on the measured spin parameter(s) can be broken down into contributions from different groups of uncertainties with:
```
python scripts/makeUncertaintyBreakdown.py -i <path-to-rebinned-input> -v <variable> -c <channel(s)> --campaign <campaign> -d <dataset> -b <path-to-breakdown-config>
```
The groups of uncertainties should be defined in a `json` config file, which is given with the `-b` argument. The default is [config/breakdown.json](config/breakdown.json).
  
The breakdown can be derived either for the bins of the unfolded distribution or for the spin parameter, depending on what you measured with `runTRExFitter.py`. This choice is controlled with the `-m` argument, which can take the `param` (default) and `distrib` values. The breakdown is written to the EFT-interpretation input file created with `runTRExFitter.py`.  

A summary table of the breakdown is also produced (spin parameter measurement only). If for some reason the job producing it has not run, it can be made locally with:
```
python scripts/makeUncertaintyBreakdown.py -i <path-to-output> -v <variable> -c <channel(s)> --campaign <campaign> -d <dataset> --dobreakdown
```
where the path to the output of the previous step should be given with the `-i` argument.

----

## Significance

if you wish to derive significance of a measurement with respect to a given null hypothesis, you should add a `significance` block to you config file such as:

```
    significance:{
        'DeltaPhi' : -0.3333333
    }
```

In this block, you should provide the name of the unfolded variable as well as the value that the corresponding parameter should take in the null hypothesis (e.g. -1/3 for the quantum entanglement measurement).  

This allows to run a second unfolding where the measured parameter of interest is fixed to its null hypothesis value, and as a consequence to derive the significance (asymptotic approximation).


<br>


# Minimal working example

A non-rebinned unfolding input can downloaded from [https://cernbox.cern.ch/index.php/s/D4QNiyNz2S8X7yD](https://cernbox.cern.ch/index.php/s/D4QNiyNz2S8X7yD). If you keep the original folder name `TRExFitter_inputUnfolding`, you can then rebin the input and run the unfolding with:

```
python scripts/rebin.py -i TRExFitter_inputUnfolding/ -v CorrKK -c em --campaign mc16a --config config/example/config.json
python scripts/runTRExFitter.py -i unfoldingInput/base/ -v CorrKK -c em --campaign mc16a --config config/example/config.json
```


<br>


# Differential unfolding

Technically speaking, a differential measurement (e.g. as a function of the ttbar invariant mass) is a 2-dimensional unfolding. In such a case, the input is in principle composed of 2D histograms for the truth- and reco-level distributions and a 4D array for the response. In order to be used in this framework, the input however needs to be flattened to 1D distributions and a 2D matrix for the response:

![](plots_readme/distrib_flattening.png?raw=true)
![](plots_readme/response_flattening.png?raw=true)

With such an input, the only differences compared to an inclusive measurement lies in the [config/rebin_config.json](config/rebin_config.json) file. The following arguments should be given there:

```
    differential_measurement : true
    diff_variable : 'name-of-differential-variable'
    diff_variable_unit : 'unit-of-differential-variable'
    diff_bins : [ list-of-bin-edges-for-each-differential-bin ]
``` 

In addition, the binning you want to rebin your input to should have the right dimensionality, i.e. one list of bin edges per differential bin. It is expected that, before rebining, the spin observable has the same binning with `nbin_input` bins in each differential bin.

## Minimal working example

A non-rebinned unfolding input can downloaded from [https://cernbox.cern.ch/index.php/s/D4QNiyNz2S8X7yD](https://cernbox.cern.ch/index.php/s/D4QNiyNz2S8X7yD). If you keep the original folder name `TRExFitter_inputUnfolding`, you can then rebin the input and run the unfolding with:

```
python scripts/rebin.py -i TRExFitter_inputUnfolding/ -v CorrKK -c em --campaign mc16a --config config/example/config_diff.json
python scripts/runTRExFitter.py -i unfoldingInput/base/ -v CorrKK -c em --campaign mc16a --config config/example/config_diff.json
```


<br>


# Unfolding with control regions

When unfolding several reco-level channels, it is possible to use some of these as control regions rather than signal regions. While in a signal region the same observable is used as reco- and truth-level, the reco-level observable in a control region is different from the unfolded one. Such a procedure can allow to constrain systematic uncertainties by measuring the shape of the control observable (e.g. constrain ttbar modelling with transverse momentum of ttbar system).   

The control regions are defined in the [config/rebin_config.json](config/rebin_config.json) file. The following arguments should be given there:

```
    is_control : true
    control_channels : ['ee', 'mm']
    control_variable : 'PTttbar'
    control_variable_unit : 'GeV'
    control_variable_bounds :{
        'Mttbar' : [340, 1200],
        'PTttbar' : [0, 300]
    }
    nbin_control : 15 // rebin from nbin_input to nbin_control
``` 

## Minimal working example

A non-rebinned unfolding input can downloaded from [https://cernbox.cern.ch/index.php/s/PjOIjnczMe7q6yW](https://cernbox.cern.ch/index.php/s/PjOIjnczMe7q6yW). If you keep the original folder name `TRExFitter_inputUnfolding_Mttbar_control`, you can then rebin the input and run the unfolding with:

```
python scripts/rebin.py -i TRExFitter_inputUnfolding_Mttbar_control/ -o unfoldingInput/Mttbar_control/ -v CorrKK -c em ee mm --campaign mc16a --config config/example/config_control.json
python scripts/runTRExFitter.py -i unfoldingInput/Mttbar_control/ -o results/Mttbar_control/ -v CorrKK -c em ee mm --campaign mc16a --config config/example/config_control.json
```


<br>


# Quantum entanglement measurement

The quantum entanglement and the corresponding significance for different ttbar-invariant-mass thresholds can be derived in an automised manner with:
```
python scripts/measureQuantumEntanglement.py  -i <path-to-input-folder> -o <path-to-output> -or <path-to-rebinned-input> -c <channel(s)> --campaign <campaign> -l <ttbar-mass-thresholds>
```
where the list of tested mass thresholds is given with the `-l` argument. It is assumed that the original input has the following structure:
```
<input-folder>/<mass-threshold>/<campaign>/<channel>/<variable>.root
```

For instance, the command could look like:
```
python scripts/measureQuantumEntanglement.py  -i TRExFitter_inputUnfolding_PTttbar_control -o results/PTttbar_control -or unfoldingInput/PTttbar_control -c em ee mm --campaign all -l 395 400 405 450
```

<br>


# Evaluate uncertainty from limited MC statistics

Limited MC statistics can be taken into account in the fit via so-called gamma nuisance parameters. These are handled separately for background and signal. They are included for background by adding `UseMCstat : TRUE` in the 'bkgsample' block of the config file. The gammas are included by default for signal, but the can be removed by adding `Gammas : DISABLED` to the `unfoldingsample` block of the config.

As a cross-check, the uncertainty from limited MC statistics (for signal) can also be derived from an ensemble test, where each bin of the response matrix is varied according to its statistical uncertainty. This is achieved with the following command:
```
python scripts/evaluateMCstatUncertainty.py -i <path-to-rebinned-input> -v <variable> -c <channel(s)> --campaign <campaign> --ntoys <number-of-pseudo-experiments> -m <method-to-draw-toys>
```
The standard deviation of the distribution of the spin parameter measured for each pseudo-experiments is the estimate for the uncertainty. The type of p.d.f. used to draw toy response matrices can be chosen with the `-m` argument. It can take the values `Gaus` to use a gaussian p.d.f. Gaus(N_MC|sigma_MC), or 'MCstat' to use a Poisson p.d.f. N_MC * Poisson(N_eff) / N_eff, where N_eff = (N_MC/sigma_MC)^2.  
The histogram of the distribution should be produced. If for some reason the job producing it has not run, it can be made locally with:
```
python scripts/evaluateMCstatUncertainty.py -i <path-to-output> -v <variable> -c <channel(s)> --campaign <campaign> --ntoys <number-of-pseudo-experiments> --plot-ensemble
```


<br>


# Binning optimisation

The binning used to unfold a given observable can be optimised. For this, all possible binnings for a given number of bins can be generated with:
```
python scripts/generateBinnings.py -i <path-to-input-folder> -v <variable> -c <channel(s)> --campaign <campaign> -d <list-of-datasets> -n <list-of-bin-numbers>
```
This script will in addition apply a filter to keep only configuration for which all bins have a width which is not smaller than the resolution (-15% to release a bit the constraint (hard-coded)). A default `json` file containing the resolution of the spin observables is provided in `config/resolution.json`.  
It will then run the unfolding for each selected binning and perform the linearity test. The `-n` argument takes a list of bin numbers for which binnings should be generated (default = 6, 8, 10).

Once the previous step is done, the best binnings, i.e. binnings performing the best at the linearity test, can be picked with:
```
python scripts/pickBinnings.py -i <path-to-generateBinnings-output> -v <variable> -c <channel(s)> --campaign <campaign>
```
The number of picked binnings for each bin number can be set with the argument `-n` (default=2).  
In addition, an ensemble test can be performed in order to check the robustness of the selected best binnings against data statistical fluctuations. For this, pseudo-experiments can be run for the selected best binnings by adding the `-p --ntoys <number-of-pseudo-experiments>` arguments.  
The pseudo-experiments are run in parallel on the batch system. Once they are done running, the results of the ensemble test can be extracted with:
```
python scripts/pickBinnings.py -i <path-to-generateBinnings-output> -v <variable> -c <channel(s)> --campaign <campaign> --ntoys <number-of-pseudo-experiments> --plot-ensemble
```


<br>


# Measure correlations between spin observables

The correlations between the spin observables should be check in order to determine whether the observables can safely be unfolded independantly from one another.  
To achieve this, one should provide an input containing **bootstraps** maintaining correlations between observables. A bootstrap is pseudo-dataset, derived by varying the signal distribution and constructing the corresponding varied asimov dataset. It should be generated by assigning each signal MC event a weight drawn from a Poisson p.d.f. of mean 1. The same weights should be used in the distributions of all the spin observables, thus maintaining the correlations.  
For each bootstrap, the 15 spin observables are unfolded and the corresponding spin parameters are measured. A sufficient number of bootstraps allows then to extract the correlations between the spin parameters.

The name pattern and the number of bootstraps should be indicated in [config/rebin_config.json](config/rebin_config.json). After rebinning the input, the unfolding step can be achieved with:
```
python scripts/measureObservableCorrelations.py -i <path-to-rebinned-input> -v <list-of-variables> -c <channel(s)> --campaign <campaign>
```
where the default for the `-v` argument is the full list of the 15 spin observables.  
The correlation matrix between the spin parameters can then be drawn with:
```
python scripts/measureObservableCorrelations.py -i <path-to-output> -v <list-of-variables> -c <channel(s)> --campaign <campaign> --plot-correlation
```
where the path to the output of the previous step should be given with the `-i` argument.
