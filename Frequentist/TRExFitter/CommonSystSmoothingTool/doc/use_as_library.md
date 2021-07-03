# Use this package as a library!  
> This package provide functions that can apply to smooth systematics histograms.
> Package includes methods from `WSM`, `RF` and TRExFitter. However user may not aware of
> all settings related to different methods. Therefore configured methods are available 
> from a common interface fuction just passing predefined names as an argument.

## Example of using  as a library  
Assume you have a nominal histogram `hnom` and a systematic histogram `hsys`.    
`SmoothSystematics/SmoothHist.h` provide a fuction `Smooth(TH1* hnom, TH1* hsys, string Option, SmoothedHistErrors smoothedErrors)`    
to smooth `hsys`.  

### example code  
```c++
//Creat Smoothing Alg
  SmoothHist smoothTool; 
  
  //Set some parameters. These parameters are already set by default.
  // call them only if you need different values.
  //smoothTool.setStatErrThreshold(0.05);
  //smoothTool.setTRExTolerance(0.08); 
  //smoothTool.setTRExNbins(1);
  //smoothTool.setSmoothedHistErrors(SmoothedHistErrors::NoErrors);
  
  TH1* hsmooth = smoothTool.Smooth(hnom, hsys, "smoothRebinParabolic");

```

The `Smooth()` fuction will modify `hsys` which means after this fuction is called  
the `hsys` is become a smoothed histogram. The `hsmooth` is also a smoothed histogram.  
In this case `hsys` and `hsmooth` poiting to same histogram.

### If you want to save the original `hsys` you can to this:
```c++
TH1* hsmooth = smoothTool.Smooth(hnom, (TH1*)hsys->Clone("hsmooth"), "smoothRebinParabolic");
```  
Here `hsmooth` is a smoothed histogram.

### Preconfigured methods

- Names of preconfigured methods: 
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

