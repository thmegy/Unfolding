#include "SmoothSystematics/TRExTools.h"
#include "SmoothSystematics/SmoothHist.h"
#include "SmoothSystematics/PlotHist.h"
#include "SmoothSystematics/SmoothingTool.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ function TREx::Smooth_Ttres(TH1*, TH1*, bool)+;
#pragma link C++ function TREx::Smooth_maxVariations(TH1*, TH1*, int , double)+;
#pragma link C++ class SmoothHist+;
#pragma link C++ class PlotHist+;
#pragma link C++ class SmoothingTool+;

#endif
