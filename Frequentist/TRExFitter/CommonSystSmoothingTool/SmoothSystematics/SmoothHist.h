#ifndef __SMOOTH_HIST_H__
#define __SMOOTH_HIST_H__
/*
 *Add some description here
 *
 *
 */

#include "TH1.h"
//#include "SmoothSystematics/TRExTools.h"
#include <iostream>

enum class SmoothType { SmoothDelta, SmoothRatio };
enum class SmoothedHistErrors { NoErrors,  Original, Propagated };

namespace TREx{
  TH1* Smooth_Ttres(TH1 *hnom, TH1 *hsys, bool independentVar, SmoothedHistErrors smoothedHistErrors);  
  TH1* Smooth_maxVariations(TH1* hnom, TH1* hsys, int nbins, double tolerance, SmoothedHistErrors smoothedHistErrors);
}

class SmoothHist {
public:
  
  /************* This method is used in RF*********************//**
   * This method will merge bins based on stat>0.05 (default).\n
   * Better to clone sys hist. otherwise original hist might be modified.\n
   * hnom is pointer to nominal historgram.\n
   * hsys is a pointer to systematic histogram.\n
   * */
  TH1 *smoothHistogram(TH1 *hnom, TH1 *hsys, bool smooth=true);

  
  /************** This method is used in WSM ******************//**
   * Alternative method: use kernel fuction to smooth histograms.\n
   * Better to clone sys hist. otherwise original hist might be modified.\n
   */
  TH1 *smoothWithKernel(TH1 *const hnom, TH1 *const hsys);

   /** 
   * Set kernel fuction. kOption is s string.\n
   * There are two kernel fuction available: "box" and "normal".\n
   * "box"-> take average of several bins.\n
   * "normal"-> Gaussian fuction.\n
   * */
  void setKernelOption(std::string kOption) {
    kernel = kOption;
    if ((kernel != "box") && (kernel != "normal")) {
      std::cerr << "ERROR: available kernel options are {box, normal}." << std::endl;
      std::exit(-1);
    }
  }

  /** Set max number of maximum extrema.\n
   * Here it set it to 1 (default). This means total number of
   * extrema will be 3 (1 max extrema + 2 min extrema).
   * */
  void setNmax(int NumberOfMax = 1) { Nmax = NumberOfMax; }

  
  /** 
   * This is apply to smoothing with kernel.\n
   * There are two ways to apply smoothing:\n
   * 1) smoothe the ratio:  Sys/nominal.\n
   * 2) Smooth the delta: Sys-Nominal\n
   * */
  void setSmoothType(std::string sType) {
    if (sType == "delta")
      DeltaOrRatio = SmoothType::SmoothDelta;
    else if (sType == "ratio")
      DeltaOrRatio = SmoothType::SmoothRatio;
    else {
      std::cerr << "ERROR: available smoothType {delta, ratio}." << std::endl;
      std::exit(-1);
    }
  };

  /** Configure how errors are to be handled on the smoothed hist.\n
      There are three options:
      SmoothedHistErrors::NoErrors: errors are set to 0
      SmoothedHistErrors::Original: errors before smoothing are used
      SmoothedHistErrors::Propagated: errors are propagated
   * */
  void setSmoothedHistError(SmoothedHistErrors smoothedErrors = SmoothedHistErrors::NoErrors){
    smoothedHistErrors = smoothedErrors;
  }

  
  /**
   * This is apply to smoothing with kernel.\n
   * This is a list of radii for kernel fuction.\n
   * If you dont know what to set then dont call this fuction.\n
   * The default list will be used.*/
  void setSpans(std::vector<double> sListOfSpans) { spans = sListOfSpans; }

  /*
   * This is apply to smoothing with kernel.\n
   * In the default span list, maximum span is 2.\n
   * If you bin width is larger then 2, you should call this fuction.\n
   * */
  //void useRelativeSpans(bool Rel = true) { relative = Rel; }

  /**
   * Rebin histograms befor smoothing.\n
   * -1 means dont rebin.
   *  WARNING: if you set this and use it in a loop over systematic hists
   *  you may rebin nominal  histogram several times.\n
   *  So do not use this if you don know what are you doing.*/
  void setRebinFactor(int RebinF = -1) { rebinFactor = RebinF; }

  /** This applys to smoothing with kernel.\n
   * Set threshold for stat error.\n
   * if stat error of bins are larger then this threshold
   * merge bins between two extrema untill number of extrema reach
   * the value set by setNmax(N).\n
   * If total stat error of nominal hist larger then the threshold,
   * skip smoothing will appy.
   * */
  void setStatErrThreshold(float ErrFraction) { StatErr = ErrFraction; }

  /***
   * This is from TRExFitter\n
   * Merge bins with using new criteria.\n
   *
   * check if one bin has:
   * | (S(i)-N(i))/N(i) - (S(i-1)-N(i-1))/N(i-1)| <
   * sqrt[ dS(i)^2/N(i)^2 + S(i)^2 dN(i)^2/N(i)^4
   *       - 2 dS(i)/N(i) S(i) dN(i)/N(i)^2
   *     dS(i+1)^2/N(i+1)^2 + S(i+1)^2 dN(i+1)^2/N(i+1)^4
   *     - 2 dS(i+1)/N(i+1) S(i+1) dN(i+1)/N(i+1)^2 ]
   * */  
  TH1* Smooth_Ttres(TH1 *hnom, TH1 *hsys);

  /**
   *  Set independentVar for Ttres method.
   * */
  void setIndependentVar(bool inVar=false){independentVar=inVar;}

  void setTRExNbins(int nBins){TRExNbins=nBins;}
  void setTRExTolerance(double tole){TREx_tolerance=tole;}
  TH1* Smooth_maxVariations(TH1* hnom, TH1* hsys);

  std::vector<int> findPeaks(TH1 *h);
  int testBins(TH1 *hnom, TH1 *hsys);
  double kdist(TH1 *hnom, TH1 *hsys);
  TH1* tchannelSmooth(TH1 *hnom, TH1 *hsys);

  /**
   * This is a common interface to all methods.
   * */

  TH1* Smooth(TH1* hnom, TH1* hsys, std::string SmoothingMethod,
              SmoothedHistErrors smoothedHistErrors = SmoothedHistErrors::NoErrors);

private:

  /**
   * Find positon of local extrema in historgram
   * hnom is nominal historgram 
   * hsys is systematic histogram
   * */
  std::vector<int> getLocalExtremaBinning(TH1 *hnom, TH1 *hsys,
                                          unsigned int nmax = 1);
  
  TH1 *smoothRemovePoint(TH1 *const hnom, TH1 *const hsys, float span_up,
                         double &chi2_up, const int pointRemoved = -1);

  double getBestSpan(TH1 *const hnom, TH1 *const hsys);
  std::vector<double> spans = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
                               1.0, 1.2, 1.4, 1.6, 1.8, 2.0};


  unsigned int Nmax = 1;
  int rebinFactor = -1;
  bool relative = true;
  SmoothType DeltaOrRatio = SmoothType::SmoothRatio;
  SmoothedHistErrors smoothedHistErrors = SmoothedHistErrors::NoErrors;
  
  std::string kernel = "box";

  // stat uncertainty threshold;
  float StatErr = 0.05;
  int independentVar = false;
  double TREx_tolerance=0.08;
  int TRExNbins=1;
};
#endif
