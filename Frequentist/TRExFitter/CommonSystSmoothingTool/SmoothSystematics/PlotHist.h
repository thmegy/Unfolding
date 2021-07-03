#ifndef __PLOTHIST_H_
#define __PLOTHIST_H_

#include <iostream>
//#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
//#include "TPad.h"
//#include "TLegend.h"

//class TLegend;
using namespace std;

class PlotHist {
public:
  PlotHist(std::string pFile = "");
  ~PlotHist();
  void plot(TH1 *hnom, TH1 *hsys, TH1 *hsysSmooth);
  void plotWithRatio(TH1 *hnom, TH1 *hsys, TH1 *hsysSmooth);
  void normalize(bool norm=true) { Normalize = norm; }
  void setRebinFactor(int rBin) { RebinFactor = rBin; }
  TH1* getPull(TH1* hNOM, TH1* hSYS, bool fillErr=false);
private:
  TFile *OutputPlotFile;

  bool Normalize=false;
  int RebinFactor = -1;
};

#endif
