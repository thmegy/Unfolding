#ifndef __SMOTTHINGTOOL_H__
#define __SMOTTHINGTOOL_H__

#include <iostream>
#include <map>
#include <string>

#include "SmoothSystematics/PlotHist.h"
#include "SmoothSystematics/SmoothHist.h"
#include "TDirectory.h"
#include "TFile.h"
//#include "TLegend.h"
#include "TObject.h"

//class PlotHist;

using namespace std;

class SmoothingTool {
public:
  SmoothingTool(){};
  ~SmoothingTool();

  /**
   * Add input root file for nominal histograms.\n
   * This fuction take a file name as an argument.\n
   * The second argument is name of a TDirectory. 
   * */
  void addInputFileNom(const string NomFileName, const string TdirName = "");
    
  /**
   * Add input root file for systematics histograms.
   * This fuction take a file name as an argument.\n
   * The second argument is name of a TDirectory. 
   * */
  void addInputFileSys(const string SysFileName, const string TdirName = "");

  /** Set out put file: This file is a copy of input file + smoothed hists.
   * Smoothed hists are at the same place as the original systematic hists.
   * If you want to write smoothed hists to original file, don't call this
   * fuction.*/
  void addOutputFileSys(const string SysFileName);

  /** Add name of nominal hist and its corresponding systematics
  *hsysNameList is a list of sytematics variations {sys_upName, sys_downName,
  * ...}
  * */
  void addHist(string hnomName, vector<string> hsysNameList);
  
  /*Smoothing option: this tool support four type of smoothing\n
   * 1) Merge bins then smooth: this is used in ResonanceFinder (RF)\n
   *    To use this option simpliy call this fuction with "merge".\n
   * 2) Kernel based smoothing: This is used WorlSpaceMaker (WSM).\n
   *    To use this method call this fuction with "kernel" \n
   * 3) Ttress method from TrexFitter \n
   *    To use this method call this fuction with "Ttres".\n
   * 4) TRExDefault method from TrexFitter \n
   *    o use this method call this fuction with "TRExDefault".\n
   */
  void setSmoothingOption(string sOption) { SmoothingOption = sOption; }

  /** Configure how errors are to be handled on the smoothed hist.\n
      There are three options:
      SmoothedHistErrors::NoErrors: errors are set to 0
      SmoothedHistErrors::Original: errors before smoothing are used
      SmoothedHistErrors::Propagated: errors are propagated
    * */
  void setSmoothedHistErrors(SmoothedHistErrors smoothedErrors = SmoothedHistErrors::NoErrors){
    smoothedHistErrors = smoothedErrors;
  }

  
  /** If some of the systematics should not be smoothed, add their name (or part of the names)
   * to skip list*/
  void skipHists(vector<string> sList) { SkipList = sList; }

  
  /** Smooth only first N histograms.\n
   * This is very useful when you want to do some tests with small numbe of hists.\n
   * NOTE: ST.generateList() will tall all histograms except the hists in the
   * skip list.
   */
  void generateList(int FirstN = 0);
  
  /** There is two type of kernel:\n
   *  1) "box": take average of of several bins.\n
   *  2) "normal": normal gaussian kernel\n
   *  */
  void setKernelOption(string kOption) { SmoothAlg.setKernelOption(kOption); }
  
  void setSmoothType(string sType){SmoothAlg.setSmoothType(sType);}
  void setSmoothedHistName(string sName) { SmoothedHistName = sName; }
  void setSpans(vector<double> sList){ SmoothAlg.setSpans(sList);}
  //void useRelativeSpans(bool rSpans = true){SmoothAlg.useRalativeSpans(rSpans);}
  void setRebinFactor(int RBF = -1){/*SmoothAlg.setRebinFactor(RBF)*/ rebinFactor = RBF;}
  void setStatErrThreshold(float rErrThreshold){SmoothAlg.setStatErrThreshold(rErrThreshold);}
  void setNmax(int nOfMax = 1){SmoothAlg.setNmax(nOfMax);}
  void setSmoothHist(SmoothHist &sAlg) { SmoothAlg = sAlg; }
  bool skip(string checkSkip);
  void setTRExTolerance(double tole=0.08){SmoothAlg.setTRExTolerance(tole);}
  void setTREx_nbins(int nBins=3){TREx_nbins=nBins;}
  bool runSmoothing();

  void makePlot(bool oPlot = true) { 
    MakePlot = oPlot; 
    PlotAlg = new PlotHist();
  }
  void normalizePlot(bool normPlot = true){PlotAlg->normalize(normPlot);}
  void setRebinFactorPlot(int rBinPlot){PlotAlg->setRebinFactor(rBinPlot);}
  void SetIndependentVar(bool InVar=false){independentVar = InVar;}
private:
  string InputFileNomName;
  string InputFileSysName;

  string TdirectorySys = "";
  string TdirectoryNom = "";
  string OutputFileSysName;

  int rebinFactor = -1;
  // List of histograms for smoothing.
  map<string, vector<string>> ListOfHistToSmooth;
  vector<string> SkipList = {""};

  string SmoothingOption = "merge"; //"merge", "kernel", "Ttres", "TRExDefault"
  SmoothedHistErrors smoothedHistErrors = SmoothedHistErrors::NoErrors; //"merge", "kernel", "Ttres", "TRExDefault"  
  // This will added to smoothed hist name;
  string SmoothedHistName = "Smooth";
  bool MakePlot = false;
  
  bool independentVar = false;
  int TREx_nbins=3;
  // SmoothHist Alg.
  SmoothHist SmoothAlg;
  // plot hist
  PlotHist* PlotAlg = nullptr;

  // Copy input file content to output file
  void copyHist(TDirectory *td);
};

#endif
