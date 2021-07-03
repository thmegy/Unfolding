#include <iostream>

//#include "TCanvas.h"
//#include "TColor.h"
//#include "TDirectory.h"
//#include "TF1.h"
//#include "TFile.h"
//#include "TH1F.h"
//#include "TLegend.h"
//#include "TRandom3.h"
#include "TStyle.h"
/*This is header file for Smoothing Algorithm*/
#include "SmoothSystematics/SmoothingTool.h"
int main() {
  gStyle->SetOptStat(0); //turn off Stat box in plots.

  /*This is input file: containig nominal and systematic histograms*/
  string inFile = "13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root";
  /*This is an output file: will contain nominal systimatic and smoothed systematic
   *histograms. Output directory structure is same as the input file.*/
  string outFile ="13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaSmoothed.root";

  // Test SmoothingTool
  cout << "Configuring Smoothing tool ..." << endl;
  SmoothingTool ST; //create instance 
  
  /*Configure tool*/
  ST.addInputFileNom(inFile); //Add input file

  /*Add input file for systematics histograms. Systematics hists are stored in
   * a TDirectory named "Systematics"*/
  ST.addInputFileSys(inFile, "Systematics");
  
  /*Set out put file: This file is a copy of input file + smoothed hists.
   * Smoothed hists are at the same place as the original systematic hists.
   * If you want to write smoothed hists to original file, don't call this
   * fuction.*/
  ST.addOutputFileSys(outFile);

  /*Smoothing option: this tool support two type of smoothing
   * 1) Merge bins then smooth: this is used in ResonanceFinder (RF)
   *    To use this option simpliy call this fuction with "merge".
   * 2) Kernel based smoothing: This is used WorlSpaceMaker (WSM).
   *    To use this method call this fuction with "kernel" 
   * 3) Ttres smoothing from TrexFitter.
   *    To use this method call this fuction with "Ttres"
   */
  //ST.setSmoothingOption("merge");
  //ST.setSmoothingOption("Ttres");
  ST.setSmoothingOption("smoothRebinMonotonic");

  /*How to configure kernel based smoothing?*/
  //ST.setSmoothingOption("kernel"); //use kernel based smoothing
  
  /*There is two type of kernel:
   *  1) "box": take average of of several bins.
   *  2) "normal": normal gaussian kernel*/
  ST.setKernelOption("box");  // set kernel type "box" or "normal"
  
  /*Rebin histograms befor smoothing.
   * smoothed hists are saved after rebinning. But nominal and
   * original Systematics hists are saved with original binning.*/
  ST.setRebinFactor(20);


  /*If some of the systematics should not be smoothed, add their name (or part of the names)
   * to skip list*/
  ST.skipHists({"data",   "ggZll",     "stop",     "Wcc",   "Wl",  "Wcl",
                "Zl",     "Zcc",       "Zcl",      "ttbar", "ttH", "Wbc",
                "Wb",     "Wbb",       "WlvH125",  "WZ",    "Zbc", "Zbl",
                "MadZee", "qqZllH125", "MadZmumu", "ZZ",
                "SysEL_EFF_Reco_TOTAyL_1NPCOR_PLUS_UNCOR", "SysEL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
                "SysEL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR", "SysEL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
                "SysMUON_EFF_SYS","SysMUON_EFF_STAT", "SysMUON_ISO_SYS","SysMUON_ISO_STAT","SysMUON_TTVA_SYS",
                "SysMUON_TTVA_STAT","SysFT_EFF_extrapolation_from_charm","SysFT_EFF_extrapolation","SysFT_EFF_Eigen_Light_0",
                "SysFT_EFF_Eigen_Light_1","SysFT_EFF_Eigen_Light_2","SysFT_EFF_Eigen_Light_3","SysFT_EFF_Eigen_Light_4",
                "SysFT_EFF_Eigen_Light_5","SysFT_EFF_Eigen_Light_6","SysFT_EFF_Eigen_Light_7","SysFT_EFF_Eigen_Light_8",
                "SysFT_EFF_Eigen_Light_9","SysFT_EFF_Eigen_Light_10","SysFT_EFF_Eigen_Light_11","SysFT_EFF_Eigen_Light_12",
                "SysFT_EFF_Eigen_Light_13","SysFT_EFF_Eigen_C_0","SysFT_EFF_Eigen_C_1","SysFT_EFF_Eigen_C_2","SysFT_EFF_Eigen_C_3",
                "SysFT_EFF_Eigen_B_0","SysFT_EFF_Eigen_B_1","SysFT_EFF_Eigen_B_2","SysFT_EFF_Eigen_B_3","SysFT_EFF_Eigen_B_4",
                "SysFT_EFF_Eigen_B_5","SysVVMbbME","SysVVPTVME","SysVHQCDscalePTV","SysVHQCDscaleMbb","SysVHPDFPTV",
                "SysVHQCDscalePTV_ggZH","SysVHQCDscaleMbb_ggZH","SysVHPDFPTV_ggZH","SysVHUEPSPTV",
                "SysVHUEPSMbb","SysVHNLOEWK","SysStoptPTV","SysStoptMBB","SysStopWtPTV","SysStopWtMBB"
              });
  
  /*Smooth only first N histograms.
   * This is very useful when you want to do some tests with small numbe of hists.
   * NOTE: ST.generateList() will tall all histograms except the hists in the
   * skip list.
   */
  ST.generateList(10); 
  
  /*If you want to smoothed few specific histograms, use 
   * ST.addHist("NominalHist", {sys1_up, sys1_down, sys2_up, sys2_down,...})
   * In this case only this sys1/2 will be smoothed.
   * NOTE: Don't use ST.addHist() togather with ST.generateList().
   *        Choose one of them.
   */
  // ST.addHist("Zbb", {"Zbb_SysJET_21NP_JET_EffectiveNP_1__1down"});
  // "Zbb_SysMUON_SAGITTA_RESBIAS__1up"});
  // "Zbb_SysMET_SoftTrk_ResoPara__1up",
  // "Zbb_SysJET_21NP_JET_EffectiveNP_3__1down"});
  
  cout << "running smoothing tool ..." << endl;
  
  /*This tool can plot nominal, sys and smoothed sys.
   * plots are saved onder plots/ directory. Both root and pdf vestions of 
   * plots are saved.
   * */
  ST.makePlot(1); //Turn on plotting.
  
  //Normalized plot.
  ST.normalizePlot(false);
  
  /*Rebin histograms to get nice plots: 
   * 0 is auto rebin, plot in 10 bins.
   * -1 is plot original bins.
   *  N is rebin hists: histogram.Rebin(N)
   */
  ST.setRebinFactorPlot(-1);

  /*run smoothing. It will return ture if every thing is OK.
   * This is note very useful for now. Impelemented is for future.
   */
  bool result = ST.runSmoothing();

  if (result) {
    cout << "Smoothing succeed." << endl;
  } else {
    cout << "Smoothing failed." << endl;
  }

  cout << "exiting ..." << endl;
}

