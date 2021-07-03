#ifndef __TREXTOOLS_H_
#define __TREXTOOLS_H_

#include "TH1.h"
#include "SmoothSystematics/SmoothHist.h"

#include <vector>
namespace TREx{


//-------------------------
//  TTres from TRExFitter
//--------------------------

    struct Bin {
        double N;
        double S;
        double dN2;
        double dS2;
        double edge;
        Bin(double _N, double _S, double _dN2, double _dS2, double _edge) { N = _N; S = _S; dN2 = _dN2; dS2 = _dS2; edge = _edge; }
    };
// check if one bin has:
// | (S(i)-N(i))/N(i) - (S(i-1)-N(i-1))/N(i-1)| <
// sqrt[ dS(i)^2/N(i)^2 + S(i)^2 dN(i)^2/N(i)^//       - 2 dS(i)/N(i) S(i) dN(i)/N(i)^2
//       - 2 dS(i)/N(i) S(i) dN(i)/N(i)^2
//       dS(i+1)^2/N(i+1)^2 + S(i+1)^2 dN(i+1)^2/N(i+1)^4
//       - 2 dS(i+1)/N(i+1) S(i+1) dN(i+1)/N(i+1)^2 ]
bool systFluctuation(std::vector<Bin> &hist, bool independentVar);

/**
 * This is an alternative smoothing in TRExFitter\n
 * */
  TH1* Smooth_Ttres(TH1 *hnom, TH1 *hsys, bool independentVar,
                    SmoothedHistErrors smoothedHistErrors = SmoothedHistErrors::Original);

int get_nVar(TH1* hratio);
int rebin_getMaxVar(TH1* hnom,TH1* hsys, double tolerance);

/**
 * This a default smoothing in TRExFitter\n
 * */
  TH1* Smooth_maxVariations(TH1* hnom, TH1* hsys, int nbins, double tolerance=0.08,
                            SmoothedHistErrors smoothedHistErrors = SmoothedHistErrors::Original);

}

#endif
