#include "SmoothSystematics/TRExTools.h"

#include <algorithm>
#include <cmath>
#include <forward_list>
#include <TMath.h>
#include <TMathBase.h>
#include <iostream>

using namespace std;
using namespace TREx;

//-------------------------
//  TTres from TRExFitter
//--------------------------
//_________________________________________________________________________
//
int getBinWidth(TH1 *ratio){

    //
    // Returns the minimal number of consecutive bins with the same content
    //
    float prev=0;
    int count=1, mincount=99;

    for(int iBin = 1;iBin <= ratio->GetNbinsX(); ++iBin){

        if(ratio->GetBinContent(iBin)==0) continue;

        if(prev!=0){
            if( TMath::Abs(prev-ratio->GetBinContent(iBin))< 1e-05 ){
                count++;
            } else {
                if(count < mincount) mincount=count;
                count=1;
            }
        }
        prev = ratio->GetBinContent(iBin);
    }
    if(count < mincount) mincount=count;
    return mincount;
}

// check if one bin has:
// | (S(i)-N(i))/N(i) - (S(i-1)-N(i-1))/N(i-1)| <
// sqrt[ dS(i)^2/N(i)^2 + S(i)^2 dN(i)^2/N(i)^4
//       - 2 dS(i)/N(i) S(i) dN(i)/N(i)^2
//       dS(i+1)^2/N(i+1)^2 + S(i+1)^2 dN(i+1)^2/N(i+1)^4
//       - 2 dS(i+1)/N(i+1) S(i+1) dN(i+1)/N(i+1)^2 ]
bool TREx::systFluctuation(std::vector<TREx::Bin> &hist, bool independentVar) {
    auto dM_indep = [](const Bin &b) { return sqrt(b.dN2 + b.dS2); };
    auto dM_dep = [](const Bin &b) { return max(sqrt(b.dN2), sqrt(b.dS2)); };
    auto dM = independentVar?dM_indep:dM_dep;
    auto dMoverN = [dM](const Bin &b) {
        double N = b.N;
        if (N == 0) N = 1e-16;
        return dM(b)/N;
    };
    auto SoverN = [](const Bin &b) {
        double N = b.N;
        if (N == 0) N = 1e-16;
        return (b.S - N)/N;
    };
    auto dSoverN2 = [dMoverN](const Bin &b) -> double {
        double N = b.N;
        if (N == 0) N = 1e-16;
        //double r = b.dS2*pow(1.0/N, 2);
        double r = std::pow(dMoverN(b), 2);
        if (r < 0) return 0.0;
        return r;
    };
    int Nbins = hist.size();
    for (int k = 1; k < Nbins; ++k) {
        double variation_prev = fabs(SoverN(hist[k]) - SoverN(hist[k-1]));
        double sum_errors = sqrt(dSoverN2(hist[k]) + dSoverN2(hist[k-1]));
        //double sum_errors = max(sqrt(dSoverN2(hist[k])), sqrt(dSoverN2(hist[k-1])));
        if (variation_prev < sum_errors) return true;
    }
    return false;
}
TH1* TREx::Smooth_Ttres(TH1 *hnom, TH1 *hsyst, bool independentVar,
                        SmoothedHistErrors smoothedHistErrors){

    //
    // General idea: merge bins with large relative stat. error until systematic variation is larger than stat. error in all bins or only one bin is left
    //
    float systIntegral = hsyst->Integral();
    std::vector<Bin> hist;
    int Nbins = hnom->GetNbinsX();
    for (int k = 1; k <= Nbins; ++k) {
        hist.push_back(Bin(hnom->GetBinContent(k), hsyst->GetBinContent(k), pow(hnom->GetBinError(k), 2), pow(hsyst->GetBinError(k), 2), hnom->GetXaxis()->GetBinLowEdge(k)));
    }
    std::string temp1 = hnom->GetName();
    std::string temp2 = hsyst->GetName();

    auto dM_indep = [](const Bin &b) { return sqrt(b.dN2 + b.dS2); };
    auto dM_dep = [](const Bin &b) { return max(sqrt(b.dN2), sqrt(b.dS2)); };
    auto dM = independentVar?dM_indep:dM_dep;
    auto dMoverN = [dM](const Bin &b) {
        double N = b.N;
        if (N == 0) N = 1e-16;
        return dM(b)/N;
    };




    // merge until all bins satisfy:
    // | (S(i)-N(i))/N(i) - (S(i-1)-N(i-1))/N(i-1)| >
    // sqrt[ dS(i)^2/N(i)^2
    // remove this:      + S(i)^2 dN(i)^2/N(i)^4
    // remove this:      - 2 dS(i)/N(i) S(i) dN(i)/N(i)^2
    //       dS(i+1)^2/N(i+1)^2
    // remove this:      + S(i+1)^2 dN(i+1)^2/N(i+1)^4
    // remove this:      - 2 dS(i+1)/N(i+1) S(i+1) dN(i+1)/N(i+1)^2 ]
    // here merge first bins with largest relative difference to the nominal between two bins
    auto SoverN = [](const Bin &b) {
        double N = b.N;
        if (N == 0) N = 1e-16;
        return (b.S - N)/N;
    };
    auto dSoverN2 = [dMoverN](const Bin &b) -> double {
        double N = b.N;
        if (N == 0) N = 1e-16;
        //double r = b.dS2*pow(1.0/N, 2);
        double r = pow(dMoverN(b), 2);  // DANILO
        if (r < 0) return 0.0;
        return r;
    };
    while (systFluctuation(hist, independentVar) && (hist.size() > 1)) {
        // first check if a bin is larger than 100%
        bool mergedLarge100 = false;
        for (unsigned int k = 0; k < hist.size(); ++k) {
            if (fabs(hist[k].S - hist[k].N) >= hist[k].N) {
                std::vector<Bin>::iterator toMergeItr = hist.begin() + k;
                std::vector<Bin>::iterator toMergeSecond = std::prev(toMergeItr);
                if (k == hist.size()-1) toMergeSecond = std::prev(toMergeItr);
                else if (k == 0) toMergeSecond = std::next(toMergeItr);
            else {
                if (sqrt(dSoverN2(*std::prev(toMergeItr))) > sqrt(dSoverN2(*std::next(toMergeItr)))) {
                    toMergeSecond = std::prev(toMergeItr);
                } else {
                    toMergeSecond = std::next(toMergeItr);
                }
            }
            toMergeItr->N += toMergeSecond->N;
            toMergeItr->S += toMergeSecond->S;
            toMergeItr->dN2 += toMergeSecond->dN2;
            toMergeItr->dS2 += toMergeSecond->dS2;
            if (toMergeItr->edge > toMergeSecond->edge)
                toMergeItr->edge = toMergeSecond->edge;
            hist.erase(toMergeSecond);
            mergedLarge100 = true;
            break;
        }
    }
    if (mergedLarge100) continue;

        // at least a pair of bins have a difference of relative differences to the nominal larger than the stat error
        // find the pair of bins with largest error/relative difference
        std::vector<double> relDiff(hist.size()-1);
        std::vector<Bin>::iterator it = std::next(hist.begin());
        std::generate(relDiff.begin(), relDiff.end(), [&it,SoverN,dSoverN2]() -> double {
             double r = sqrt(dSoverN2(*it) + dSoverN2(*std::prev(it)))/fabs(SoverN(*it) - SoverN(*std::prev(it)));
             //double r = max(sqrt(dSoverN2(*it)), sqrt(dSoverN2(*std::prev(it))))/fabs(SoverN(*it) - SoverN(*std::prev(it)));
             it++;
             return r;
             });

        std::vector<double>::iterator toMergeRelItr = std::max_element(relDiff.begin(), relDiff.end());
        int binToMerge = (int) (toMergeRelItr - relDiff.begin()) + 1; // difference always taken between current and previous, so index 0, means merging 0 and 1
        std::vector<Bin>::iterator toMergeItr = hist.begin() + binToMerge; // get iterator
        std::vector<Bin>::iterator toMergeSecond = std::prev(toMergeItr); // always merge with previous
        toMergeItr->N += toMergeSecond->N;
        toMergeItr->S += toMergeSecond->S;
        toMergeItr->dN2 += toMergeSecond->dN2;
        toMergeItr->dS2 += toMergeSecond->dS2;
        if (toMergeItr->edge > toMergeSecond->edge)
            toMergeItr->edge = toMergeSecond->edge;
        hist.erase(toMergeSecond);
    }

    //
    // Define the array for new bin ranges and a template histogram with this binning
    //
    double* Bins;
    Bins=new double[hist.size()+1];
    for ( unsigned int i=0; i < hist.size(); i++) {
        Bins[i]=hist[i].edge;
    }
    Bins[hist.size()] = hnom->GetXaxis()->GetBinUpEdge(hnom->GetNbinsX());

    TH1F* hnomBinned = new TH1F("hnomBinned", "", hist.size(), Bins);
    hnomBinned->Sumw2();
    TH1F* hsystBinned = new TH1F("hsystBinned", "", hist.size(), Bins);
    hsystBinned->Sumw2();
    for (int V = 1; V <= hnomBinned->GetNbinsX(); V++) {
        hnomBinned->SetBinContent(V, hist[V-1].N);
        hnomBinned->SetBinError(V, sqrt(hist[V-1].dN2));
        hsystBinned->SetBinContent(V, hist[V-1].S);
        hsystBinned->SetBinError(V, sqrt(hist[V-1].dS2));
    }
    hsystBinned->Divide(hnomBinned);//now, hsystBinned is the relative systematic uncertainty

    //
    // Modify the systematic uncertainty histogram based on the "stat reliable" ratio to avoid statistical fluctuations
    //
    std::vector<double> err;
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        // systematic_new = nominal(with normal binning) * relative_systematic_uncertainty(with new binning)
        hsyst->SetBinContent(i,hnom->GetBinContent(i)*hsystBinned->GetBinContent(hsystBinned->FindBin( hsyst->GetBinCenter(i))));
        err.push_back(hsyst->GetBinError(i));
    }

    delete hnomBinned;
    delete hsystBinned;
    delete [] Bins;

    // use TH1::Smooth() to apply 353QH method on |S-N|/N histogram
    TH1F *ratio = (TH1F *) hsyst->Clone();
    ratio->Divide(hnom);
    hsyst->Add(hnom,-1);
    hsyst->Divide(hnom);

    //First non-empty bin
    int minbin = 0;
    for(int i=1; i < ratio->GetNbinsX()+1;i++){
        if(ratio->GetBinContent(i)!=0){
            minbin = i;
            break;
        }
    }

    //Last non-empty bin
    int maxbin = ratio->GetNbinsX();
    for(int i=maxbin; i>= 1;i--){
        if(ratio->GetBinContent(i)!=0){
            maxbin = i;
            break;
        }
    }

    // Smooth only works well for positive entries: shifts all entries by an offset of 100
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 1000 );
    }

    // Due to the rebinning, some bins can have the same content. Call the ROOT smooth function to avoid this.
    int binwidth = getBinWidth(ratio);
    if(binwidth>1) binwidth=1;
    hsyst->GetXaxis()->SetRange(minbin,maxbin);
    if(binwidth<maxbin-minbin){
        if (independentVar) {
            hsyst->Smooth(4, "R");
        } else {
            hsyst->Smooth(binwidth, "R");
        }
    }

    // Removes the 100 offset
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 1000 );
    }
    hsyst->Multiply(hnom);
    hsyst->Add(hnom);

    if(hsyst->Integral()!=0){
        hsyst->Scale(systIntegral/hsyst->Integral());
    }
    // Checks if any bin with < 0 content exists
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        double content = hsyst->GetBinContent( i );
        hsyst->SetBinError(i, err[i-1]);
        if(content < 0){
            hsyst -> SetBinContent(i, 0.);
        }
    }
    delete ratio;

    if (smoothedHistErrors == SmoothedHistErrors::NoErrors) {
      // set bin errors to 0 for systematics
      for (int i = 0; i < hsyst->GetNbinsX() + 2; i++) {
        hsyst->SetBinError(i, 0);
      }
    }
    
    return hsyst;
}

//-----------------
//This is a default smoothing in TRExFitter
//----------------
int TREx::get_nVar(TH1* hratio){

    //
    // Counts the number of slope changes
    //

    int nVar=0;
    bool thisUp = true;
    bool goingUp = true;
    double prevBin = 0, thisBin = 0;
    int usedBins = 0;
    bool hasChange = false;

    for (int bin=1; bin <=hratio->GetNbinsX(); bin++) {

        if(hratio->GetBinContent(bin)==0) continue;

        prevBin = thisBin;
        thisBin = hratio->GetBinContent(bin);

        if(fabs(thisBin-prevBin)<1e-7) continue;//skip bins with same content

        usedBins++;

        thisUp = (thisBin > prevBin);
        if (usedBins > 2 && ((thisUp && !goingUp) || (!thisUp && goingUp))){
            nVar++;
            if(hasChange)
                nVar++;
            hasChange = true;
        }
        else
            hasChange = false;


        goingUp = thisUp;
    }

    return nVar;
}


//_________________________________________________________________________
//
int TREx::rebin_getMaxVar(TH1* hnom, TH1* hsyst, double tolerance){

  //
  // Recompute the systematic histogram based on a new binning (built based on the relative stat uncertainty to suppress
  // statistical fluctuations)
  //

  //for(int i=1;i<=hsyst->GetNbinsX();i++){
  //  WriteVerboseStatus("HistoTools::rebin_getMaxVar", "In: " + std::to_string(hnom->GetBinContent(i)) + " " + std::to_string(hsyst->GetBinContent(i)));
  //}

  std::vector<double> binLimit;
  binLimit.push_back( hnom->GetXaxis()->GetXmin() );
  double relErr=20000;
  double cumulIntSyst=0;
  double cumulInt=0;
  double cumulErr=0;
  int thisBin=0;

  do {

    do { //while (relErr > tolerance && thisBin!=hnom->GetNbinsX() );

          //
          // Compute the relative statistical uncertainty of a group of bins. Performs this operation until
          // the relative statistical uncertainty is lower than the tolerance (or the number of bins)
          //

          thisBin++;

          cumulInt+=fabs(hnom->GetBinContent(thisBin));
          cumulErr+=hnom->GetBinError(thisBin)*hnom->GetBinError(thisBin);
          cumulIntSyst += hsyst->GetBinContent(thisBin);
          if ( cumulInt!=0 && cumulIntSyst!=0 ){
              relErr= sqrt(cumulErr)/cumulInt;
          }
          if (relErr==0) relErr=20000;


        } while (relErr > tolerance && thisBin!=hnom->GetNbinsX() );

        if((relErr < tolerance) || binLimit.size() == 1){//a group of bins with a sufficiently low stat error has been found, let's add it
            binLimit.push_back(hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2);
        } else {//no such group of bins has been found: merge with the last found bin
            binLimit.back() = hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2;
        }


        cumulInt=0;
        cumulErr=0;
        relErr=20000;

    } while ( thisBin!=hnom->GetNbinsX() );

    //
    // Define the array for new bin ranges and a template histogram with this binning
    //
    double* Bins;
    Bins=new double[binLimit.size()];
    for ( unsigned int i=0; i< binLimit.size(); i++) {
        Bins[i]=binLimit[i];
    }

    TH1F* rebinTemplate=new TH1F( "binRef", "binRef", binLimit.size()-1, Bins);
    rebinTemplate->Sumw2();
    for (int V=1; V<=rebinTemplate->GetNbinsX(); V++) {
        rebinTemplate->SetBinContent(V,V);
    }

    //
    // Performs a rebin "by hand" of the nominal and systematic histograms
    //
    //nominal
    TH1F* hnomBinned=new TH1F(*rebinTemplate);
    hnomBinned->Reset();
    for (int bin=0; bin <=hnom->GetNbinsX(); bin++) {
        int bigBin=hnomBinned->FindBin( hnom->GetBinCenter(bin));
        hnomBinned->SetBinContent( bigBin, hnomBinned->GetBinContent(bigBin)+ hnom->GetBinContent(bin) );
        hnomBinned->SetBinError( bigBin, sqrt( pow(hnomBinned->GetBinError(bigBin),2) + pow(hnom->GetBinError(bin),2) ) );
    }
    //systematics
    TH1F* hsystBinned=new TH1F(*rebinTemplate);
    hsystBinned->Reset();
    for (int bin=0; bin <=hsyst->GetNbinsX(); bin++) {
        int bigBin=hsystBinned->FindBin( hsyst->GetBinCenter(bin));
        hsystBinned->SetBinContent( bigBin, hsystBinned->GetBinContent(bigBin)+ hsyst->GetBinContent(bin) );
        hsystBinned->SetBinError( bigBin, sqrt( pow(hsystBinned->GetBinError(bigBin),2) + pow(hsyst->GetBinError(bin),2) ) );
    }
    hsystBinned->Divide(hnomBinned);//now, hsystBinned is the relative systematic uncertainty

    //
    // Modify the systematic uncertainty histogram based on the "stat reliable" ratio to avoid statistical fluctuations
    //
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        // systematic_new = nominal(with normal binning) * relative_systematic_uncertainty(with new binning)
        hsyst->SetBinContent(i,hnom->GetBinContent(i)*hsystBinned->GetBinContent(hsystBinned->FindBin( hsyst->GetBinCenter(i))));
    }

    //
    // Computes the number of slope variations in the new systematic histogram
    //
    int nVar = get_nVar(hsystBinned);

    delete rebinTemplate;
    delete hnomBinned;
    delete hsystBinned;
    delete [] Bins;

    return nVar;
}


//_________________________________________________________________________
//
TH1* TREx::Smooth_maxVariations(TH1* hnom, TH1* hsyst, int nbins, double tolerance,
                                SmoothedHistErrors smoothedHistErrors){

    //
    // General idea: avoid having more than "nbins" slope variations in the systematic histogram
    //
    float systIntegral = hsyst->Integral();

    //double tolerance = 0.08;

    int nVar = rebin_getMaxVar(hnom,hsyst,tolerance);


    //
    // Iterates the smoothing of the systematic histogram until the number a slope changes is lower than "nbins"
    //
    while (nVar > nbins){
        tolerance = tolerance/2.;
        nVar = rebin_getMaxVar(hnom,hsyst,tolerance);
        if(tolerance==0){
            std::string temp1 = hnom->GetName();
            std::string temp2 = hnom->GetTitle();
            std::cerr<<"HistoTools::Smooth_maxVariations---"<< "Buuuuuuuuuuug: infinite while"<<std::endl;
            std::cerr<<"HistoTools::Smooth_maxVariations---"<<temp1 <<" " << temp2 <<" nbins: " << std::to_string(nbins)<<std::endl;
            break;
        }
    }

    TH1F *ratio = (TH1F *) hsyst->Clone();
    ratio->Divide(hnom);
    hsyst->Add(hnom,-1);
    hsyst->Divide(hnom);

    //First non-empty bin
    int minbin = 0;
    for(int i=1; i < ratio->GetNbinsX()+1;i++){
        if(ratio->GetBinContent(i)!=0){
            minbin = i;
            break;
        }
    }

    //Last non-empty bin
    int maxbin = ratio->GetNbinsX();
    for(int i=maxbin; i>= 1;i--){
        if(ratio->GetBinContent(i)!=0){
            maxbin = i;
            break;
        }
    }

    // Smooth only works well for positive entries: shifts all entries by an offset of 100
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 1000 );
    }

    // Due to the rebinning, some bins can have the same content. Call the ROOT smooth function to avoid this.
    int binwidth = getBinWidth(ratio);
    if(binwidth>4) binwidth=4;
    hsyst->GetXaxis()->SetRange(minbin,maxbin);
    if(binwidth*2<maxbin-minbin){
        hsyst->Smooth(binwidth*2,"R");
    }

    // Removes the 100 offset
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 1000 );
    }
    hsyst->Multiply(hnom);
    hsyst->Add(hnom);
    if(hsyst->Integral()!=0){
        hsyst->Scale(systIntegral/hsyst->Integral());
    }
    // Checks if any bin with < 0 content exists
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        double content = hsyst->GetBinContent( i );
        if(content < 0){
            hsyst -> SetBinContent(i, 0.);
        }
    }
    delete ratio;

    if (smoothedHistErrors == SmoothedHistErrors::NoErrors) {
      // set bin errors to 0 for systematics
      for (int i = 0; i < hsyst->GetNbinsX() + 2; i++) {
        hsyst->SetBinError(i, 0);
      }
    }
    
    return hsyst;
}
