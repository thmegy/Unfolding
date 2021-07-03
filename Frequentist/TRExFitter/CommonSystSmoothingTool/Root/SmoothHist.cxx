#include <SmoothSystematics/SmoothHist.h>
#include "SmoothSystematics/TRExTools.h"

#include <algorithm>
#include <cmath>
#include <forward_list>

#include <TAxis.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include "TDirectoryFile.h"
#include "TFile.h"

// **********************************************************************
// Smoothing definitions used in WSMaker
// **********************************************************************

using namespace std;

float statError(TH1 *hnom, int beg, int end) { // end is excluded
  double err = 0;
  float nomInt = hnom->IntegralAndError(beg, end - 1, err);
  return fabs(err / nomInt);
}

float computeChi2(TH1 *hnom, TH1 *hsys, int beg, int end) {
  float ratio = hsys->Integral(beg, end) / hnom->Integral(beg, end);
  float chi2 = 0;
  for (int i = beg; i < end + 1; i++) {
    if (hnom->GetBinContent(i) != 0) {
      float iratio = hsys->GetBinContent(i) / hnom->GetBinContent(i);
      float err = hnom->GetBinError(i) / hnom->GetBinContent(i);
      chi2 += ((iratio - ratio) / err) * ((iratio - ratio) / err);
    }
  }
  return chi2;
}

int findSmallerChi2(TH1 *hnom, TH1 *hsys, const std::vector<int> &extrema) {
  int pos = 0;
  float minval = 99999;
  for (unsigned int i = 0; i < extrema.size() - 1; i++) {
    float chi2 = computeChi2(hnom, hsys, extrema[i], extrema[i + 1]);
    if (chi2 < minval) {
      pos = i;
      minval = chi2;
    }
  }
  return pos;
}

void getRatioHist(TH1 *hnom, TH1 *hsys, const std::vector<int> &bins,
                  TH1 *res) {
  for (unsigned int iRefBin = 0; iRefBin < bins.size() - 1; iRefBin++) {
    float nomInt = hnom->Integral(bins.at(iRefBin), bins.at(iRefBin + 1) - 1);
    float varInt = hsys->Integral(bins.at(iRefBin), bins.at(iRefBin + 1) - 1);
    for (int b = bins.at(iRefBin); b < bins.at(iRefBin + 1); b++) {
      if (nomInt != 0) {
        res->SetBinContent(b, varInt / nomInt);
      } else {
        res->SetBinContent(b, 0);
      }
    }
  }
}

TH1 *getRatioHist(TH1 *hnom, TH1 *hsys, const std::vector<int> &bins) {
  TH1 *res = (TH1 *)hsys->Clone();
  getRatioHist(hnom, hsys, bins, res);
  return res;
}

// inclusive in lo and hi
void mergeBins(int lo, int hi, std::vector<int> &bins) {
  std::vector<int>::iterator beg =
      std::upper_bound(bins.begin(), bins.end(), lo);
  // +1 because inclusive merge
  std::vector<int>::iterator last =
      std::lower_bound(bins.begin(), bins.end(), hi + 1);
  bins.erase(beg, last);
}

std::vector<int> findExtrema(TH1 *h) {
  std::vector<int> res;
  res.push_back(1);
  int status = 0; // 1: potential max, -1: potential min
  int k = 1;
  for (int i = 2; i < h->GetNbinsX() + 1; i++) {
    // special rule for bins with 0 stat. Keep going on, until one finds another
    // bin to compare to
    if (h->GetBinContent(i) < 1.e-6) {
      continue;
    }
    if (status == 1 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      res.push_back(i - 1);
      status = -1;
    }
    if (status == -1 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      res.push_back(i - 1);
      status = 1;
    }
    if (status == 0 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      status = -1;
    }
    if (status == 0 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      status = 1;
    }
    k = i;
  }
  res.push_back(h->GetNbinsX());

  return res;
}

std::vector<int> SmoothHist::getLocalExtremaBinning(TH1 *hnom, TH1 *hsys,
                                                    unsigned int nmax) {
  // This implementation is iterative.
  // A faster (say analytic) implementation is possible if this one proves to be
  // too slow. This one is however easier to write and read
  std::vector<int> res;
  double err = 0;
  float sum = hnom->IntegralAndError(0, hnom->GetNbinsX()+1, err);
  //float sum = hnom->IntegralAndError(1, hnom->GetNbinsX(), err);
  // too large stat unc: no shape
  if (sum>0. && fabs(err / sum) > StatErr) {
    res.push_back(1);
    res.push_back(hnom->GetNbinsX() + 1);
    return res;
  }

  // normal case. Then, beginning with no rebinning
  for (int i = 1; i < hnom->GetNbinsX() + 2; i++) {
    res.push_back(i);
  }

  TH1 *ratio = getRatioHist(hnom, hsys, res);
  std::vector<int> extrema = findExtrema(ratio);

  while (extrema.size() > nmax + 2) {
    int pos = findSmallerChi2(hnom, hsys, extrema);
    mergeBins(extrema[pos], extrema[pos + 1], res);
    getRatioHist(hnom, hsys, res, ratio);
    extrema = findExtrema(ratio);
  }

  // second pass to avoid bins with too large stat uncertainty
  std::vector<int>::iterator fst = res.end();
  std::vector<int>::iterator lst = res.end();
  std::vector<int> to_remove;
  --lst;
  --fst;
  while (fst != res.begin()) {
    if (fst == lst) {
      --fst;
    } else {
      float statE = statError(hnom, *fst, *lst);
      if (statE > StatErr || statE != statE) {
        to_remove.push_back(fst - res.begin());
        --fst;
      } else {
        lst = fst;
      }
    }
  }
  for (int i : to_remove) {
    res.erase(res.begin() + i);
  }

  delete ratio;
  return res;
}


TH1 *SmoothHist::smoothHistogram(TH1 *hnom, TH1 *hsys, bool smooth) {
  
  if (rebinFactor > 0){
    hnom->Rebin(rebinFactor);
    hsys->Rebin(rebinFactor);
  }

  const std::vector<int> bins = getLocalExtremaBinning(hnom, hsys, Nmax);
  float norm_init = hsys->Integral();
  TH1 *ratio = getRatioHist(hnom, hsys, bins);

  if (smooth && ratio->GetNbinsX() > 2) {
    std::vector<float> vals(ratio->GetNbinsX() - 2);
    for (int i = 2; i < ratio->GetNbinsX(); i++) {
      vals[i - 2] =
          (2. * ratio->GetBinContent(i) + ratio->GetBinContent(i - 1) +
           ratio->GetBinContent(i + 1)) /
          4.;
    }
    for (int i = 2; i < ratio->GetNbinsX(); i++) {
      ratio->SetBinContent(i, vals[i - 2]);
    }
  }

  for (int i = 1; i < hsys->GetNbinsX() + 1; i++) {
    if (hnom->GetBinContent(i) != 0) {
      hsys->SetBinContent(i, ratio->GetBinContent(i) * hnom->GetBinContent(i));
    } else {
      hsys->SetBinContent(i, 0);
    }
  }
  hsys->Scale(norm_init / hsys->Integral());
  delete ratio;
  if (smoothedHistErrors == SmoothedHistErrors::NoErrors) {
    // set bin errors to 0 for systematics. Easier later when doing chi2 tests
    for (int i = 0; i < hsys->GetNbinsX() + 2; i++) {
      hsys->SetBinError(i, 0);
    }
  }
  return hsys;
}

// **********************************************************************
// Function for caclulating chi2 from TH1F
// **********************************************************************

double chi2(TH1 *const hnom, TH1 *const hsmooth) {

  double chi2 = 0.0;
  for (int i = 1; i < hnom->GetNbinsX() + 1; ++i) {
    chi2 += std::pow(hsmooth->GetBinContent(i) - hnom->GetBinContent(i), 2);
    // if ( h_nom->GetBinError(i) != 0.0 ) chi2 /=
    // std::pow(h_nom->GetBinError(i),2);
  }

  // chi2 += std::pow(maxDelta(h_smooth),2);

  return fabs(chi2);
}

double SmoothHist::getBestSpan(TH1 *const hnom, TH1 *const hsys) {
  // Loop over all points to be removed
  double bestSpan(0.);
  int npoints = hnom->GetNbinsX();
  if (npoints < 2)
    return bestSpan; // don't do anything if bins <2.
  for (int i = 0; i < npoints; ++i) {

    double bestChi2Point(-1.);
    double bestSpanPoint(-1.);

    // loop over spans
    for (auto span : spans) {

      double chi2_up;
      // Smooth
      smoothRemovePoint(hnom, hsys, span, chi2_up, i);

      // Store best chi2 and span
      if (bestSpanPoint < 0. || chi2_up < bestChi2Point) {
        bestChi2Point = chi2_up;
        bestSpanPoint = span;
      }
      // std::cout << "span = " << span << " , chi2_up = " << chi2_up <<
      // std::endl;

    } // span loop

    // printf("---> Point %2d : Best span up = %4.3f\n", i, bestSpanPoint);
    bestSpan += bestSpanPoint;

  } // point loop

  // Calculate best span
  bestSpan /= npoints;
  // printf("===> %-60s :   Best span = %4.3f\n", sysName.Data(), bestSpan_up);
  return bestSpan;
}

// **********************************************************************
// Function for smoothing deltas
// **********************************************************************
TH1 *SmoothHist::smoothRemovePoint(TH1 *const hnom, TH1 *const hsys,
                                   float span_up, double &chi2_up,
                                   const int pointRemoved) {
  // if nbins < 2, do nothing.
  if (hsys->GetNbinsX() < 2)
    return hsys;
  // otherwise perfome smoothing.
  //
  TH1F *h_cl = nullptr;
  TH1F *h_up_cl = nullptr;

  h_cl = (TH1F *)hnom->Clone(TString::Format("%d", pointRemoved));
  h_up_cl = (TH1F *)hsys->Clone(TString::Format("%d", pointRemoved));
  //No rebinning here. Already done.
  /*if (rebinFactor > -1) {
    h_cl->Rebin(rebinFactor);
    h_up_cl->Rebin(rebinFactor);
  } else if (rebinFactor == -2) {
    std::vector<Double_t> binning = {
        -1.000, -0.867, -0.733, -0.600, -0.467, -0.333, -0.200, -0.067,
        0.067,  0.200,  0.333,  0.467,  0.600,  0.733,  0.867,  1.000};
    Double_t *binarray = binning.data();

    h_cl = (TH1F *)hnom->Rebin(15, "_rebinned", binarray);
    h_up_cl = (TH1F *)hsys->Rebin(15, "_rebinned", binarray);
  } else {
  }*/ 

  if (relative) {
    span_up = span_up * h_cl->GetBinWidth(1);
  }

  const int npoints = h_cl->GetNbinsX();

  // Get Delta or ratio
  TH1F *h_delta_up = (TH1F *)h_up_cl->Clone();
  if (DeltaOrRatio == SmoothType::SmoothDelta) {
    h_delta_up->Add(h_cl, -1.); // h_delta_up = h_up - h_nom
  } else {
    const int nbins = h_delta_up->GetNbinsX();
    for (int i = 1; i < nbins + 1; ++i) {
      if (h_cl->GetBinContent(i) != 0.0)
        h_delta_up->SetBinContent(i, h_delta_up->GetBinContent(i) /
                                         h_cl->GetBinContent(i));
      else
        h_delta_up->SetBinContent(i, 1.0);
    }
  }

  // Convert TH1 to TGraph
  TGraph *g_delta_up = new TGraph(h_delta_up);

  // Get x values of graph
  Double_t *xpoints = new Double_t[npoints];
  for (int i = 0; i < npoints; ++i)
    xpoints[i] = h_delta_up->GetBinCenter(i + 1);

  if (pointRemoved != -1) { // Remove point i
    g_delta_up->RemovePoint(pointRemoved);
  }

  // Smooth without removing points
  TGraph *g_delta_up_fin = new TGraph(h_delta_up);
  TGraphSmooth gs_smooth_up_fin;
  TGraph *g_delta_up_smooth = nullptr;
  if (pointRemoved != -1) {
    g_delta_up_smooth = gs_smooth_up_fin.SmoothKern(g_delta_up, kernel.c_str(),
                                                    span_up, npoints, xpoints);
  } else {
    g_delta_up_smooth = gs_smooth_up_fin.SmoothKern(
        g_delta_up_fin, kernel.c_str(), span_up, npoints, xpoints);
  }

  // Get the y values of the smoothed histogram
  double *yval_up = new double[npoints];
  int filled = 0;
  bool extrap = false;
  for (int i = 0; i < g_delta_up_smooth->GetN(); ++i) {
    if (i == pointRemoved && !extrap) {
      double y_up_inter =
          g_delta_up_smooth->Eval(xpoints[pointRemoved], nullptr, "S");
      yval_up[pointRemoved] = y_up_inter;
      ++filled;
      --i;
      extrap = true;
    } else {
      double x(0.0), y_up(0.0);
      g_delta_up_smooth->GetPoint(i, x, y_up);
      if (fabs(x - xpoints[filled]) < 1.E-10) {
        yval_up[filled] = y_up;
        ++filled;
      }
    }
  }
  // ... and convert the graph to a histogram
  TH1F *h_delta_up_smooth = (TH1F *)h_delta_up->Clone();
  for (int i = 1; i < h_delta_up_smooth->GetNbinsX() + 1; ++i) {
    h_delta_up_smooth->SetBinContent(i, yval_up[i - 1]);
  }

  TH1F *h_delta_up_reb = nullptr;
  TH1F *h_delta_up_reb_sm = nullptr;

  // The following is to be used for custom rebinning
  // why apply rebinning agian???
  /*if (rebinFactor == -2) {
    std::vector<Double_t> binning = {
        -1.000, -0.867, -0.733, -0.600, -0.467, -0.333, -0.200, -0.067,
        0.067,  0.200,  0.333,  0.467,  0.600,  0.733,  0.867,  1.000};
    Double_t *binarray = binning.data();

    h_delta_up_reb = (TH1F *)h_delta_up->Rebin(15, "_rebinned", binarray);
    h_delta_up_reb_sm =
        (TH1F *)h_delta_up_smooth->Rebin(15, "_rebinned", binarray);
  }*/

  // The following is to be used for constant rebinning
    h_delta_up_reb = (TH1F *)h_delta_up->Clone();
    h_delta_up_reb_sm = (TH1F *)h_delta_up_smooth->Clone();


  // Calculate chi2
  chi2_up = chi2(h_delta_up_reb, h_delta_up_reb_sm);

  // Get absolute histograms from deltas
  // sysHistoFromDelta(h_cl, h_delta_up_reb_sm, smoothType);
  if (DeltaOrRatio == SmoothType::SmoothDelta)
    h_delta_up_reb_sm->Add(h_cl);
  else
    h_delta_up_reb_sm->Multiply(h_cl);

  if (smoothedHistErrors == SmoothedHistErrors::NoErrors) {    
    // Set bin errors to 0.
    for (int i = 1; i < h_delta_up->GetNbinsX() + 1; ++i)
      h_delta_up_reb_sm->SetBinError(i, 0.0);
  }

  delete [] xpoints;
  delete [] yval_up;
  return h_delta_up_reb_sm;
}
///

TH1 *SmoothHist::smoothWithKernel(TH1 *const hnom, TH1* hsys) {
  //gErrorIgnoreLevel = kFatal;
  // Check if nbins < 2. If yes, do nothing.
  if (hsys->GetNbinsX() < 2)
    return hsys;
  // Otherwise perform smoothing.
  //
  if (rebinFactor > 0){
    hnom->Rebin(rebinFactor);
    hsys->Rebin(rebinFactor);
  }
  // Get best span
  float bestSpan = getBestSpan(hnom, hsys);
  //if (rebinFactor > 0){
  //  bestSpan *= hnom->GetBinWidth(1)*rebinFactor;
  //}
  double chi2Sys;
  TH1 *hsys_smooth = smoothRemovePoint(hnom, hsys, bestSpan, chi2Sys, -1);
  // hsys is cloned with name include Smooth.
  hsys_smooth->SetName(hsys->GetName());

  hsys = (TH1*)hsys_smooth->Clone();
  delete hsys_smooth;
  return hsys;
}



std::vector<int> SmoothHist::findPeaks(TH1 *h) {
  std::vector<int> res;
  int status = 0; // 1: potential max, -1: potential min
  int k = 1;
  int extremum=-1;
  for (int i = 2; i < h->GetNbinsX()+1; i++) {
    if (h->GetBinContent(i) < 1.e-6) {
      continue;
    }
    if (status == 1 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      res.push_back(extremum);
      status = -1;
      extremum = i;
    }
    if (status == -1 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      res.push_back(extremum);
      status = 1;
      extremum = i;
    }
    if (status == 0 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      status = -1;
    }
    if (status == 0 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      status = 1;
    }
    if (h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6 || h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      extremum = i;
    }
    k = i;
  }
  return res;
}

int SmoothHist::testBins(TH1 *hnom, TH1 *hsys) {
  bool bigger=false;
  int count=0;
  int badbins=0;
  double maxdiff=0;
  for(int i=1; i < hnom->GetNbinsX()+1;i++){
    if(hnom->GetBinContent(i)>hsys->GetBinContent(i)+1e-6){
      if(!bigger) count=0;
      count++;
      if(i==1||i==hnom->GetNbinsX())count++;
      bigger=true;
    }else if (hnom->GetBinContent(i)<hsys->GetBinContent(i)-1e-6){
      if(bigger) count=0;
      count--;
      if(i==1||i==hnom->GetNbinsX())count--;
      bigger=false;
    }else{
      count=0;
    }
    if(i<hnom->GetNbinsX()&&maxdiff<abs(hnom->GetBinContent(i+1)-hnom->GetBinContent(i))*0.8){
      maxdiff=abs(hnom->GetBinContent(i+1)-hnom->GetBinContent(i))*0.8;
    }
    if(count>2) badbins++;
    if(count<-2) badbins++;
  }
  return badbins;
}



double SmoothHist::kdist(TH1 *hnom, TH1 *hsys) {
  double dist=0;
  for(int j=1; j < hnom->GetNbinsX()+1;j++){
    double sum=0;
    for(int i=1; i < j;i++){
      sum+=hnom->GetBinContent(i)-hsys->GetBinContent(i);
    }
    if(abs(sum)>dist) dist=abs(sum);
  }
  return dist;
}
//Decription of this algorithm: https://gitlab.cern.ch/atlas-phys/exot/CommonSystSmoothingTool/uploads/466f09c1bb0272f5782bcd39ad3d9f7e/smoothing.pdf
TH1 *SmoothHist::tchannelSmooth(TH1 *hnom, TH1 *hsys) {
  if (rebinFactor > 0){
    hnom->Rebin(rebinFactor);
    hsys->Rebin(rebinFactor);
  }
  
  //force up down consistency
  bool flip=false;
  if (hsys->GetBinContent(1) > hnom->GetBinContent(1)){
    TH1* diff = (TH1*)hsys->Clone();
    diff->Add(hnom,-1);
    hsys->Add(diff,-2);
    flip=true;
    delete diff;
  }
  
  TH1* hsyscopy = (TH1*)hsys->Clone();
  TH1* hnomcopy = (TH1*)hnom->Clone();

  std::vector<int> bins;
  float norm_init = hsys->Integral();
  TH1 *ratio;
  TH1 *ratioold;
  double koldist;
  std::vector<int> res;
  for (int j = 1; j < hnom->GetNbinsX() + 2; j++) {
    res.push_back(j);
  }
  std::vector<int> reso;
  for (int j = 1; j < hnom->GetNbinsX() + 2; j++) {
    reso.push_back(j);
  }
  std::vector<int> ex;
  int nvarch=0;
  ratio = getRatioHist(hnomcopy, hsyscopy, res);
  ratioold = getRatioHist(hnomcopy, hsyscopy, res);
  int numiter=0;
  do{ //iteratively merge bins
    double minkoldist=999999;
    int mbin=-1;
    ex = findPeaks(ratio); //get bins with a higher or lower value than both of their adjacent bins
    for (std::size_t i=0; i<res.size()-2; i++){
      if(std::find(ex.begin(),ex.end(),res[i])!=ex.end()||std::find(ex.begin(),ex.end(),res[i+1])!=ex.end()) {
        std::vector<int> tmpres;
        for (std::size_t j=0; j<res.size(); j++) 
          tmpres.push_back(res[j]); 
        mergeBins(res[i],res[i+1], tmpres); //merge bin pair temorarily to check resulting ratio
        getRatioHist(hnomcopy, hsyscopy, tmpres, ratio);
        getRatioHist(hnomcopy, hsyscopy, reso, ratioold);
        for (int j = 1; j < hsys->GetNbinsX() + 1; j++) {
          if (hnomcopy->GetBinContent(j) != 0) {
            hsys->SetBinContent(j, ratio->GetBinContent(j) * hnomcopy->GetBinContent(j));
          } else {
            hsys->SetBinContent(j, 0);
          }
	    }
        //prevent three consecutive bins all with a higher or all with a lower value than original ratio:
        int testbinnum=testBins(ratioold,ratio);
        //calculate "distance" of new and old ratio
        koldist = kdist(ratio,ratioold);
        if(koldist<minkoldist&&testbinnum==0){
          minkoldist=koldist;
          mbin=i;
        }
      }
    }
    if(mbin!=-1){
      int i=mbin;
      mergeBins(res[i],res[i+1], res); //merge bin with smallest distance
      getRatioHist(hnomcopy, hsyscopy, res, ratio);
      for (int j = 1; j < hsys->GetNbinsX() + 1; j++) {
	if (hnomcopy->GetBinContent(j) != 0) {
	  hsys->SetBinContent(j, ratio->GetBinContent(j) * hnomcopy->GetBinContent(j));
	} else {
	  hsys->SetBinContent(j, 0);
	}
      }
    }else{ //no more bins can be merged
      getRatioHist(hnomcopy, hsyscopy, res, ratio);
      for (int j = 1; j < hsys->GetNbinsX() + 1; j++) {
	if (hnomcopy->GetBinContent(j) != 0) {
	  hsys->SetBinContent(j, ratio->GetBinContent(j)* hnomcopy->GetBinContent(j));
	} else {
	  hsys->SetBinContent(j, 0);
	}
      }
      break;
    }
    nvarch++;
    // if(numiter>999)
    //   break;
    numiter++;
  }while(ex.size()>1);
  Int_t nbins = ratio->GetNbinsX();
  
  { //slightly modified ROOT smoothing algorithm from TH1::smooth(). Modification at ################
    Int_t firstbin = 1, lastbin = nbins;
    nbins = lastbin - firstbin + 1;
    Double_t *xx = new Double_t[nbins];
    Int_t nn = nbins;
    Int_t i;
    for (i=0;i<nbins;i++) {
      xx[i] = ratio->GetBinContent(i+firstbin);
    }
    if (ratio->GetNbinsX() > 2){
 
      Int_t ii;
      Double_t hh[5] = {0,0,0,0,0};
 
      std::vector<double> yy(nn);
      std::vector<double> zz(nn);
      std::vector<double> rr(nn);
 
      for (Int_t pass=0;pass<4;pass++) {
	// first copy original data into temp array
	std::copy(xx, xx+nn, zz.begin() );
 
	for (int noent = 0; noent < 2; ++noent) { // run algorithm two times
 
	  // //  do 353 i.e. running median 3, 5, and 3 in a single loop
	  for  (int kk = 0; kk < 3; kk++)  {
	    std::copy(zz.begin(), zz.end(), yy.begin());
	    int medianType = (kk != 1)  ?  3 : 5;
	    int ifirst      = (kk != 1 ) ?  1 : 2;
	    int ilast       = (kk != 1 ) ? nn-1 : nn -2;
	    //nn2 = nn - ik - 1;
	    // do all elements beside the first and last point for median 3
	    //  and first two and last 2 for median 5
	    for  ( ii = ifirst; ii < ilast; ii++)  {
	      //assert(ii - ifirst >= 0);
	      for  (int jj = 0; jj < medianType; jj++)   {
	        hh[jj] = yy[ii - ifirst + jj ];
	      }
	      zz[ii] = TMath::Median(medianType, hh);
	    }
 
	    if  (kk == 0)  {   // first median 3
	      // first point
	      hh[0] = zz[1];
	      hh[1] = zz[0];
	      hh[2] = 3*zz[1] - 2*zz[2];
	      if(pass>0) // ################  modified here: only change after 1. pass
		zz[0] = TMath::Median(3, hh);
	      // last point
	      hh[0] = zz[nn - 2];
	      hh[1] = zz[nn - 1];
	      hh[2] = 3*zz[nn - 2] - 2*zz[nn - 3];
	      if(pass>0) //  ################  modified here: only change after 1. pass
		zz[nn - 1] = TMath::Median(3, hh);
	    }
 
	    if  (kk == 1)  {   //  median 5
	      for  (ii = 0; ii < 3; ii++) {
	        hh[ii] = yy[ii];
	      }
	      zz[1] = TMath::Median(3, hh);
	      // last two points
	      for  (ii = 0; ii < 3; ii++) {
	        hh[ii] = yy[nn - 3 + ii];
	      }
	      zz[nn - 2] = TMath::Median(3, hh);
	    }
 
	  }  
 
	  std::copy ( zz.begin(), zz.end(), yy.begin() );
	  //quadratic interpolation for flat segments
	  for (ii = 2; ii < (nn - 2); ii++) {
	    if  (zz[ii - 1] != zz[ii]) continue;
	    if  (zz[ii] != zz[ii + 1]) continue;
	    hh[0] = zz[ii - 2] - zz[ii];
	    hh[1] = zz[ii + 2] - zz[ii];
	    if  (hh[0] * hh[1] <= 0) continue;
	    int jk = 1;
	    if  ( TMath::Abs(hh[1]) > TMath::Abs(hh[0]) ) jk = -1;
	    yy[ii] = -0.5*zz[ii - 2*jk] + zz[ii]/0.75 + zz[ii + 2*jk] /6.;
	    yy[ii + jk] = 0.5*(zz[ii + 2*jk] - zz[ii - 2*jk]) + zz[ii];
	  }
	
	  ///running means
	  std::copy(yy.begin(), yy.end(), zz.begin());
	  for  (ii = 1; ii < nn - 1; ii++) {
	    zz[ii] = 0.25*yy[ii - 1] + 0.5*yy[ii] + 0.25*yy[ii + 1];
	  }
	  zz[0] = yy[0];
	  zz[nn - 1] = yy[nn - 1];
	  if (noent == 0) {
 
	    // save computed values
	    std::copy(zz.begin(), zz.end(), rr.begin());
 
	    // COMPUTE  residuals
	    for  (ii = 0; ii < nn; ii++)  {
	      zz[ii] = xx[ii] - zz[ii];
	    }
	  }
	}  // end loop on noent
	for  (ii = 0; ii < nn; ii++) {
	  xx[ii] = TMath::Max((rr[ii] +zz[ii]),0.0 );
	}
      }
    }
    for (i=0;i<nbins;i++) {
      ratio->SetBinContent(i+firstbin,xx[i]);
    }
    delete [] xx;
  }//end of ROOT smoothing algorithm

  // calculate smoothed systematic template from smoothed ratio
  for (int j = 1; j < hsys->GetNbinsX() + 1; j++) {
    if (hnomcopy->GetBinContent(j) != 0) {
      hsys->SetBinContent(j, ratio->GetBinContent(j)* hnomcopy->GetBinContent(j));
    } else {
      hsys->SetBinContent(j, 0);
    }
  }

  hsys->Scale(norm_init / hsys->Integral()); // keep normalization
  delete ratio;
  delete ratioold;
  delete hsyscopy;
  delete hnomcopy;

  //force up down consistency
  if (flip){
    TH1* diff = (TH1*)hsys->Clone();
    diff->Add(hnom,-1);
    hsys->Add(diff,-2);
    delete diff;
  }

  if (smoothedHistErrors == SmoothedHistErrors::NoErrors) {
    // set bin errors to 0 for systematics. Easier later when doing chi2 tests
    for (int i = 0; i < hsys->GetNbinsX() + 2; i++) {
      hsys->SetBinError(i, 0);
    }
  }
  return hsys;
}




//-------------------
//For TRExFitter
//-------------------
TH1* SmoothHist::Smooth_Ttres(TH1 *hnom, TH1 *hsys){
  return TREx::Smooth_Ttres(hnom, hsys, independentVar, smoothedHistErrors); 
}

TH1* SmoothHist::Smooth_maxVariations(TH1* hnom, TH1* hsys){
  return TREx::Smooth_maxVariations(hnom, hsys, TRExNbins, TREx_tolerance, smoothedHistErrors);
}

//-------------------
//This is a common interface to all methods
//------------------
TH1* SmoothHist::Smooth(TH1* hnom, TH1* hsys, string SmoothingMethod,
                        SmoothedHistErrors smoothedErrors){

  if (!hnom ){
    std::cerr<<"ERROR: passing null pointer to nominal hists."<<std::endl;
  }
  if (!hsys ){
    std::cerr<<"ERROR: passing null pointer to systematic hists."<<std::endl;
  }

  if (smoothedErrors == SmoothedHistErrors::Propagated) {
    std::cerr<<"ERROR: error propagation is not yet implemented for smoothing";
    std::cerr << "method '" << SmoothingMethod << "'." << std::endl;
    return nullptr;
  }
    
  setSmoothedHistError(smoothedErrors);
  
  if (SmoothingMethod=="smoothRebinMonotonic"){
    setNmax(0);
    return smoothHistogram(hnom, hsys, true); 
  }else if(SmoothingMethod=="smoothRebinParabolic"){
    setNmax(1);
    return smoothHistogram(hnom, hsys, true);
  }else if (SmoothingMethod=="smoothDeltaUniformKernel"){
    setKernelOption("box");
    setSmoothType("delta");
    return smoothWithKernel(hnom, hsys);
  }else if(SmoothingMethod=="smoothRatioUniformKernel"){
    setKernelOption("box");
    setSmoothType("ratio");
    return smoothWithKernel(hnom, hsys);
  }else if (SmoothingMethod=="smoothDeltaGaussKernel"){
    setKernelOption("normal");
    setSmoothType("delta");
    return smoothWithKernel(hnom, hsys);
  }else if(SmoothingMethod=="smoothRatioGaussKernel"){
    setKernelOption("normal");
    setSmoothType("ratio");
    return smoothWithKernel(hnom, hsys);
  }else if (SmoothingMethod=="smoothTtresDependent"){
    setIndependentVar(false);
    return Smooth_Ttres(hnom, hsys); 
  }else if (SmoothingMethod=="smoothTtresIndependent"){
    setIndependentVar(true);
    return Smooth_Ttres(hnom, hsys); 
  }else if (SmoothingMethod=="smoothTRExDefault"){
    return Smooth_maxVariations(hnom, hsys);
  }else if (SmoothingMethod=="smoothTchannel"){
    return tchannelSmooth(hnom, hsys);
  }else{
    std::cerr<<"ERROR: Choose one method from smoothRebinMonotonic, smoothRebinParabolic,smoothDeltaUniformKernel,smoothRatioUniformKernel,smoothDeltaGaussKernel,smoothRatioGaussKernel, smoothTtresDependent, smoothTtresIndependent, smoothTRExDefault, smoothTchannel"<<std::endl;
    return nullptr;
  }


}

