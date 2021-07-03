#ifndef SAMPLEHIST_H
#define SAMPLEHIST_H

/// Framework includes
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/PruningUtil.h"
 
/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class TFile;
class TH1;
class TPad;
class Sample;
class NormFactor;
class ShapeFactor;
class SystematicHist;

class SampleHist {
public:
    explicit SampleHist();
    explicit SampleHist(Sample *sample,TH1 *hist);
    explicit SampleHist(Sample *sample, const std::string& histoName, const std::string& fileName);
    ~SampleHist();
    //SampleHist(const SampleHist& s) = delete;
    SampleHist(SampleHist&& s) = delete;
    //SampleHist& operator=(const SampleHist& s) = delete;
    SampleHist& operator=(SampleHist&& s) = delete;
    
    SampleHist(const SampleHist& s) = default;  // copy constructor      
    //SampleHist(SampleHist&& s) = default;
    SampleHist& operator=(const SampleHist& s) = default;
    //SampleHist& operator=(SampleHist&& s) = default;

    TH1* GetHist() const;
    const Sample* GetSample() const {return fSample;}
    std::shared_ptr<SystematicHist> AddOverallSyst(const std::string& name,const std::string& storedName,double up,double down);
    std::shared_ptr<SystematicHist> AddStatSyst(const std::string& name,const std::string& storedName,int i_bin);
    std::shared_ptr<SystematicHist> AddHistoSyst(const std::string& name,const std::string& storedName,TH1* h_up,TH1* h_down);
    std::shared_ptr<SystematicHist> AddHistoSyst(const std::string& name,const std::string& storedName, const std::string& histoName_up,
                                                 const std::string& fileName_up, const std::string& histoName_down,
                                                 const std::string& fileName_down, int pruned=0);
    std::shared_ptr<SystematicHist> GetSystematic(const std::string& systName) const;
    std::shared_ptr<SystematicHist> GetSystFromNP(const std::string& NuisParName) const;
    std::shared_ptr<NormFactor> AddNormFactor(const std::string& name,double nominal, double min, double max);
    std::shared_ptr<NormFactor> AddNormFactor(std::shared_ptr<NormFactor> normFactor);
    std::shared_ptr<NormFactor> GetNormFactor(const std::string& name) const;
    std::shared_ptr<ShapeFactor> AddShapeFactor(const std::string& name,double nominal, double min, double max);
    std::shared_ptr<ShapeFactor> AddShapeFactor(std::shared_ptr<ShapeFactor> shapeFactor);
    std::shared_ptr<ShapeFactor> GetShapeFactor(const std::string& name) const;

    bool HasSyst(const std::string& name) const;
    bool HasNorm(const std::string& name) const;
    bool HasShapeFactor(const std::string& name) const;

    void WriteToFile(const std::vector<int>& blindedBins, const std::vector<double>& scales, const double gammaThreshold, std::shared_ptr<TFile> f=nullptr,bool reWriteOrig=true);
    void ReadFromFile();

    void FixEmptyBins(const bool suppress);
    void NegativeTotalYieldWarning(TH1* hist, double yield) const;

    void Print() const;

    void Rebin(int ngroup = 2, const Double_t* xbins = 0);
    void DrawSystPlot( const std::string &syst="all", TH1* h_data=0x0,
                       bool SumAndData=false, bool bothPanels=false ) const;
    void SmoothSyst(const HistoTools::SmoothOption &opt, const bool useAlternativeShapeHistFactory, std::string syst="all", bool force=false);

    void Divide(  SampleHist* sh);
    void Multiply(SampleHist* sh);
    void Add(     SampleHist* sh,double scale=1.);
    void Scale(double scale);

    void SampleHistAdd(SampleHist* h, double scale = 1.);
    void SampleHistAddNominal(SampleHist* h, double scale);
    void CloneSampleHist(SampleHist* h, const std::set<std::string>& names, double scale = 1.);
    void SystPruning(PruningUtil *pu,TH1* hTot=nullptr);

    void DrawSystPlotUpper(TPad* pad0,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData,
                           bool bothPanels) const;

    void DrawSystPlotRatio(TPad* pad1,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData) const;

    std::vector<double> GetDataScales() const;

    std::string fName;
    Sample *fSample;
    std::unique_ptr<TH1> fHist;
    std::unique_ptr<TH1> fHist_orig;
    std::unique_ptr<TH1> fHist_regBin;
    std::unique_ptr<TH1> fHist_preSmooth; // new - to use only for syst plots
    std::shared_ptr<TH1> fHist_postFit;
    std::string fFileName;
    std::string fHistoName;
    bool fIsData;
    bool fIsSig;
    std::map<std::string,bool> fIsMorph;

    std::vector < std::shared_ptr<SystematicHist> > fSyst;

    std::vector <  std::string > fNormFactorNames;
    std::vector <  std::string > fShapeFactorNames;
    
    // other useful info
    std::string fFitName;
    std::string fRegionName;
    std::string fRegionLabel;
    std::string fVariableTitle;
    bool fSystSmoothed;
};

#endif

