# Smooth Systematic histograms

## Setup ATLAS:  
    > $export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase  
    > $source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh    

## Build package
* build
    > $git clone https://:@gitlab.cern.ch:8443/atlas-phys/exot/CommonSystSmoothingTool.git    
    > $asetup AnalysisBase,21.2.20,here  
    > $mkdir build    
    > $cd build  
    > $ cmake ../   
    > $ make   
    > $ source x86_64-slc6-gcc62-opt/setup.sh    
    > $ cd ../

* test
    > $mkdir runSmooth/    
    > $ cd runSmooth/   
    > $ln -s ../CommonSystSmoothingTool/share/13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root  
    > $python ../CommonSystSmoothingTool/python/testSmooth.py  --methods smoothRebinParabolic smoothRebinMonotonic  



# Use as a library

[Integration with other tools](doc/use_as_library.md)




## Standalone usage  

Assume you have a file histograms.root which contains nominal and systematics histograms.  
You can look at `util/testInputFile.cxx` file 

### Example code
you use this packege in you `mycode.cxx` as bellow.

#### Create an instance of Smoothing tool  

```c++  
#include "SmoothSystematics/SmoothingTool.h"
.....  
string inFile = "histograms.root";
string outFile = "outputHist.root";
SmoothingTool ST; //create an instance
```  
- inFile contains nominal and systematics histograms.
- outFile will contain nominal, systematic and smoothed histograms. It is same as inFile but with Smoothed hists in it.

#### Set input and outout root files.

```c++
ST.addInputFileNom(inFile, "Nominal-TDir");  
ST.addInputFileSys(inFile, "Systematics-TDir");   
ST.addOutputFileSys(outFile); // Copy results to new file.  
```
- addInputFileNom(inputFileName, TDirectory) fuction set the input root file.  
- Here "Nominal-TDir" and "Systematics-TDir" are the names ofTDirectory object. If you don't have them leave these filds empty ST.addInputFileNom(inFile).   
- If Systematics hists are in different file it can be specified with ST.addInputFileSys(inFile, "Systematics-TDir").  
- ST.addOutputFileSys(outFile) will create a new file and copy nominal and systematic histograms to "outFile". Smoothed Histograms are also saved to new file.  

#### Skip list

```c++  
ST.skipHists({"data",   "ggZll", "SysEL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR"}); // Skip list. will not be smoothed.  

```
- If some of the systematics should not be smoothed, add their name (or part of the names) to skip list.

#### add list of hists for smoothing

```c++  
ST.generateList(10); //Test first 10 hists.
```
- generateList(10) will get first 10 hists for smoothing.  generateList() will get all hists except the hists in the skip list.

#### Alternative method to add list of hist for smoothng
```c++
ST.addHist("NominalHistName", {sys1_up, sys1_down, sys2_up, sys2_down,...})
```  
This is very useful when you want to smooth few specific histograms.  
`NOTE:` Don't use ST.addHist() togather with ST.generateList().

#### Which smoothing method to use?  
This tool provides several smoothing methods
- Predefined methods: 
    - From `WSM` and `RF`
        - [x] smoothRebinParabolic 
        - [x] smoothRebinMonotonic
        - [x] smoothRatioUniformKernel    
        - [x] smoothDeltaUniformKernel
        - [x] smoothRatioGaussKernel
        - [x] smoothDeltaGaussKernel
    - From `TRExFitter`
        - [x] smoothTtresDependent   
        - [x] smoothTtresIndependent
        - [x] smoothTRExDefault   


- Configure smoothing tool: Chose one of the smoothing methods  

```c++ 
ST.setSmoothingOption("smoothRebinParabolic");   
```

#### Configuration of smoothing
Different treatments for the error of the smoothed histogram can be configured with
the an enum `SmoothedHistErrors`. Three options are available to either set the error
to zero (`SmoothedHistErrors::NoErrors`, which is the default), to take the errors from the original histogram
(`SmoothedHistErrors::Original`), or to propagate the errors according to the
smoothing algorithm (`SmoothedHistErrors::Propagated`). The last of the three
is not supported for all smoothing algorithms yet. 
```c++ 
ST.setSmoothedHistErrors(SmoothedHistErrors::NoErrors);   
```

#### Plotting with this tool!  

```c++
ST.makePlot(1);  // enable plotting
```
Set plot option to true to make plots of nominal, original and smoothed systematics.    
Both roo and pdf version of plots are saved under plots/ directory. 

```c++
ST.setRebinFactorPlot(0); // auto rebin plots.  
```  
Rebin histograms befor plot: 0 means auto rebin plots, -1 means no rebinning and  
N is rebin hists: histogram.Rebin(N)

#### Run smoothing!
```c++
ST.runSmoothing();  
```
It will return true if everything is OK.  
This is not very useful for now. Implemented for future.

