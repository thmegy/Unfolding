#ifndef TRExFIT_H
#define TRExFIT_H

/// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TRExPlot.h"

// Unfolding includes
#include "UnfoldingCode/UnfoldingCode/FoldingManager.h"

// RooFit
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"
#include "RooMultiVarGaussian.h"

// RooStats
#include "RooStats/HistFactory/Measurement.h"

/// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forwards class declaration
class ConfigParser;
class FitResults;
class FittingTool;
class NormFactor;
class RooDataSet;
class RooWorkspace;
class Region;
class Sample;
class SampleHist;
class ShapeFactor;
class TGraphAsymmErrors;
class TruthSample;
class TFile;
class UnfoldingSample;
class UnfoldingSystematic;

class TRExFit {
public:

    enum FitType {
        UNDEFINED = 0,
        SPLUSB = 1,
        BONLY = 2 ,
        UNFOLDING = 3,
        EFT = 4
    };

    enum FitRegion {
        CRONLY = 1,
        CRSR = 2
    };

    enum InputType {
        HIST = 0,
        NTUP = 1
    };

    enum LimitType {
        ASYMPTOTIC = 0,
        TOYS = 1
    };

    enum TemplateInterpolationOption{
        LINEAR = 0,
        SMOOTHLINEAR = 1,
        SQUAREROOT = 2
    };

    enum PruningType{
        SEPARATESAMPLE = 0,
        BACKGROUNDREFERENCE = 1,
        COMBINEDREFERENCE = 2
    };

    struct TemplateWeight{
        std::string function;
        std::string range;
        std::string name;
        // cppcheck-suppress unusedStructMember
        double value;
    };

    TRExFit(std::string name="MyMeasurement");
    ~TRExFit();

    TRExFit(const TRExFit& t) = delete;
    TRExFit(TRExFit&& t) = delete;
    TRExFit& operator=(const TRExFit& t) = delete;
    TRExFit& operator=(TRExFit&& t) = delete;

    void SetPOI(std::string name="SigXsecOverSM",std::string unit="");
    void AddPOI(std::string name="SigXsecOverSM",std::string unit="");
    void SetStatErrorConfig(bool useIt=true, double thres=0.05, std::string cons="Poisson");
    void SetLumiErr(const double err);
    void SetLumi(const double lumi);
    void SetFitType(FitType type);
    void SetLimitType( LimitType type );
    void SetFitRegion(FitRegion region);

    std::shared_ptr<Sample> NewSample(const std::string& name,int type=0);
    std::shared_ptr<Systematic> NewSystematic(const std::string& name);
    Region* NewRegion(const std::string& name);

    // ntuple stuff
    void AddNtuplePath(const std::string& path);
    void SetMCweight(const std::string& weight);
    void SetSelection(const std::string& selection);
    void SetNtupleName(const std::string& name);
    void SetNtupleFile(const std::string& name);
    void ComputeBinning(int regIter);
    void DefineVariable(int regIter);

    // histogram stuff
    void AddHistoPath(const std::string& path);
    void SetHistoName(const std::string& name);

    void SmoothSystematics(std::string syst="all");

    // create new root file with all the histograms
    void CreateRootFiles();
    void WriteHistos(bool reWriteOrig=true) const;

    void DrawSystPlots() const;
    void DrawSystPlotsSumSamples() const;

    // read from ..
    void CloseInputFiles();
    void CorrectHistograms();

    void DrawAndSaveAll(std::string opt="");

    // separation plots
    void DrawAndSaveSeparationPlots() const;

    std::shared_ptr<TRExPlot> DrawSummary(std::string opt="", std::shared_ptr<TRExPlot> = nullptr);
    void DrawMergedPlot(std::string opt="",std::string group="") const;
    void BuildYieldTable(std::string opt="",std::string group="") const;

    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0) const;
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions) const;
    void DrawPieChartPlot(const std::string &opt="", int nCols=0,int nRows=0) const;
    void DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions) const;

    void CreateCustomAsimov() const;

    /**
      * Runs code that replaces asimov data with custom asimov for unfolding
      */
    void UnfoldingAlternativeAsimov();

    // turn to RooStats::HistFactory
    void ToRooStats(bool createWorkspace=true, bool exportOnly=true) const;

    RooStats::HistFactory::Channel OneChannelToRooStats(RooStats::HistFactory::Measurement* meas, const int ichan) const;

    RooStats::HistFactory::Sample OneSampleToRooStats(RooStats::HistFactory::Measurement* meas,
                                                      const SampleHist* h,
                                                      const int i_ch,
                                                      const int i_smp) const;

    void SystPruning() const;
    void DrawPruningPlot() const;

    /**
      * A helper function to draw the plots with normalisation for each systematic
      */
    void DrawSystematicNormalisationSummary() const;

    // Extras for EFT
    void ProcessEFTInputs(bool overwrite=false);

    // fit etc...
    void Fit(bool isLHscanOnly);
    RooDataSet* DumpData( RooWorkspace *ws, std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, std::map < std::string, double > &poiValues);
    std::map < std::string, double > PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType=SPLUSB, bool save=false);
    std::unique_ptr<RooWorkspace> PerformWorkspaceCombination( std::vector < std::string > &regionsToFit ) const;

    void PlotFittedNP();
    void PlotCorrelationMatrix();
    void PlotUnfoldedData() const;
    void GetLimit();
    void GetSignificance();
    void GetLikelihoodScan(RooWorkspace *ws, const std::string& varName, RooDataSet* data) const;
    void Get2DLikelihoodScan(RooWorkspace *ws, const std::vector<std::string>& varName, RooDataSet* data) const;

    // get fit results from txt file
    void ReadFitResults(const std::string& fileName);

    void PrintConfigSummary() const;

    Region* GetRegion(const std::string& name) const;
    std::shared_ptr<Sample> GetSample(const std::string& name) const;
    std::size_t GetSampleIndex(const std::string& name) const;

    void ProduceNPRanking(const std::string& NPnames);
    void PlotNPRanking(const bool flagSysts, const bool flagGammas) const;
    void PlotNPRankingManager() const;

    void PrintSystTables(std::string opt="") const;

    void MergeSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter
    void CombineSpecialSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter

    // for template fitting
    void AddTemplateWeight(const std::string& name, double);

    std::vector<TemplateWeight> GetTemplateWeightVec(const TemplateInterpolationOption& opt);

    std::string GetWeightFunction(std::vector<std::pair<double,std::string> > templatePair, unsigned int itemp, const TemplateInterpolationOption& opt) const;

    /**
     * Function that returns string that represents smoothed abs value function
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSmoothLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to calualte numerical correction to the smoothed linear function
     * @param parameter in the argument of the hyperbolic tangent function
     * @param size of the x axis interval
     * @param central position of the function
     * @param left position of the function on x axis
     * @param right position of the function on x axis
     * @param parameter of iteration, set to 1 for the first iteration
     * @return correction
     */
    double GetCorrection(double k, double width, double x_mean, double x_left, double init = 1.) const;

    /**
     * Helper function to approximate absolute value by sqrt(x^2+e)
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSquareRootLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to apply correction to square root aproximation
     * @param will return value for a from -a*sqrt(x^2+epsilon) +b
     * @param will return value for b from -a*sqrt(x^2+epsilon) +b
     * @param central position of the function
     * @param left position of the function on x axis
     * @param epsilon = precision of the approximation
     */
    void GetSquareCorrection(double *a, double *b, double x_i, double x_left, double epsilon) const;

    /**
     * Helper function to smooth morphing templates according to a given functional form (bin-by-bin)
     * @param name of the morphing parameter
     * @param functional form
     * @param pointer to an array of parameter values
     */
    void SmoothMorphTemplates(const std::string& name,const std::string& formula="pol1",double *p=0x0) const;

    /**
     * Helper function that draws plots with morphijng templates
     * @param name of the morphing parameter
     */
    void DrawMorphingPlots(const std::string& name) const;

    bool MorphIsAlreadyPresent(const std::string& name, const double value) const;

    // for grouped impact evaluation
    void ProduceSystSubCategoryMap();

    void BuildGroupedImpactTable() const;

    /**
     * Helper function that runs toys experiments
     */
    void RunToys();

    /**
     * Helper function to compute the variable string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string Variable(Region *reg,Sample *smp);

    /**
     * Helper function to compute the selection string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string FullSelection(Region *reg,Sample *smp);

    /**
     * Helper function to compute the weight string to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::string FullWeight(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true);

    /**
     * Helper function to compute the full paths to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::vector<std::string> FullNtuplePaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true);

    /**
     * Helper function to compute the full paths to be used when reading histograms, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::vector<std::string> FullHistogramPaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true, const bool isFolded = false);

    /**
     * A helper function to compute the fgull paths for a response matrix
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @return the full path
     */
    std::vector<std::string> FullResponseMatrixPaths(const Region* reg,
                                                     const UnfoldingSample* smp,
                                                     const UnfoldingSystematic* syst = nullptr,
                                                     const bool isUp = true) const;

    /**
     * A helper function to compute the fgull paths for a migration matrix
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @return the full path
     */
    std::vector<std::string> FullMigrationMatrixPaths(const Region* reg,
                                                      const UnfoldingSample* smp,
                                                      const UnfoldingSystematic* syst = nullptr,
                                                      const bool isUp = true) const;

    /**
     * A helper function to compute the fgull paths for acceptance
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @return the full path
     */
    std::vector<std::string> FullAcceptancePaths(const Region* reg,
                                                 const UnfoldingSample* smp,
                                                 const UnfoldingSystematic* syst = nullptr,
                                                 const bool isUp = true) const;

    /**
     * A helper function to compute the fgull paths for selection efficiency
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @return the full path
     */
    std::vector<std::string> FullSelectionEffPaths(const Region* reg,
                                                   const UnfoldingSample* smp,
                                                   const UnfoldingSystematic* syst = nullptr,
                                                   const bool isUp = true) const;

    /**
      * A helper function to combine paths for the truth distributions
      * @return a vector of the paths
      */
    std::vector<std::string> FullTruthPaths() const;

    /**
    * A helper function to get SampleHisto from a region that matches a name of the sample
    * @param Region
    * @@param name
    * @return SampleHist
    */
    std::shared_ptr<SampleHist> GetSampleHistFromName(const Region* const reg, const std::string& name) const;

    /**
     * A helper function to Copy a smoothing from a reference histogram to other histograms bin by bin
     * @param SampleHist
     * @param nominal histogram
     * @param up variation of histogram
     * @param down variation
     * @param flag if we want up or down variation returned
     * @return Up/Down variation for the new histogram
     */
    TH1* CopySmoothedHisto(const SampleHist* const sh, const TH1* const nominal, const TH1* const up, const TH1* const down, const bool isUp) const;

    /**
     * A helper function to get an index of systemaati variation that matches some name
     * @param SampleHist
     * @param name of the systematics
     * @return index
     */
    int GetSystIndex(const SampleHist* const sh, const std::string& name) const;

    std::shared_ptr<SystematicHist> CombineSpecialHistos(std::shared_ptr<SystematicHist> orig,
                                                         const std::vector<std::shared_ptr<SystematicHist> >& vec,
                                                         Systematic::COMBINATIONTYPE type,
                                                         const SampleHist* sh) const;

    /**
      *  A helper function to get the list of unique names of non-gamma systematics
      *  @return the list of unique non-gamma systematics
      */
    std::vector<std::string> GetUniqueSystNamesWithoutGamma() const;

    /**
      * A helper function to get the vector of non-validation regions
      * @return the vector on non-validation regions
      */
    std::vector<Region*> GetNonValidationRegions() const;

    /**
      * A helper function to get the vector of non-data, non-ghost samples
      * @return the vector of non-data, non-ghost samples
      */
    std::vector<std::shared_ptr<Sample> > GetNonDataNonGhostSamples() const;

    /**
      * A helper function that does the dropping of the bins
      */
    void DropBins();

    /**
      * A function that prepares signal inputs for unfolding.
      * Folded distributions are created
      */
    void PrepareUnfolding();

    /**
      * A helper function to fold systematic distributions needed for unfolding
      * @param Folding manager
      * @param output file
      * @param Region
      * @param UnfoldingSample
      * @param Current UnfoldingSystystematics
      * @param Nominal migration matrix
      */
    void ProcessUnfoldingSystematics(FoldingManager* manager,
                                     TFile* file,
                                     const Region* reg,
                                     const UnfoldingSample* sample,
                                     const UnfoldingSystematic* syst,
                                     const TH2* nominal) const;

    /** A helper function that does the actual plotting of unfolded data
      * @param unfoded data
      * @param error band
      */
    void PlotUnfold(TH1D* data,
                    TGraphAsymmErrors* band) const;

    /**
      * A helper function to plot migration or reposne matrix
      * @param matrix
      * @param flag if sample is migration
      * @param name of the region
      * @param name of the systematic
      */
    void PlotMigrationResponse(const TH2* matrix,
                               const bool isMigration,
                               const std::string& regionName,
                               const std::string& systematicName) const;

    /**
      * A helper function to force shape on some systematics
      */ 
    void RunForceShape();
    
    /**
     * A helper function to check if the fit/limit/significance is done on a mixture of real-data and Asimov-data
     */
    bool DoingMixedFitting() const;
    
    /**
     * A helper function to create a list of all regions to fit
     * @param Flag to consider FitRegion setting
     * @param If set, only regions with this DataType are included
     */
    std::vector < std:: string > ListRegionsToFit(const bool useFitRegion, int dataType=-1) const;
    
    /**
     * A helper function to create a map with specified regions and corresponding DataType
     * @param List of regions to consider
     * @param Set to true to force to use ASIMOV everywhere
     */
    std::map < std::string, int > MapRegionDataTypes(const std::vector<std::string>& regionList,bool isBlind=false) const;
        
    /**
     * A helper function to read a fit result txt file and store the NP values into a map
     * @param Text file containing fit results
     */
    std::map < std::string, double > NPValuesFromFitResults(const std::string& fitResultsFile);

    /**
     * A helper function to plot pull plots from toys when unfolding is used
     * @param Results in each bin
     * @param output file to save the result
     */
    void DrawToyPullPlot(std::vector<TH1D>& hist, TFile* out) const;
    
    /**
     * A helper function to replace the norm-factor expression of the last bin of a truth distribution when performing normalized-cross-section unfolding
     */
    void FixUnfoldingExpressions();

    /**
      * Function to run Limit estimation using toys
      * @param Data
      * @param Workspace
      */  
    void RunLimitToys(RooAbsData* data, RooWorkspace* ws) const;
    // -------------------------

    std::string fName;
    std::string fDir;
    std::string fLabel;
    std::string fInputFolder;
    std::string fInputName;
    std::string fFitResultsFile;

    std::vector < std::shared_ptr<TFile> > fFiles;

    std::vector < Region* > fRegions;
    std::vector < std::shared_ptr<Sample> > fSamples;
    std::vector < std::shared_ptr<Systematic> > fSystematics;
    std::vector < std::shared_ptr<NormFactor> >fNormFactors;
    std::vector < std::shared_ptr<ShapeFactor> > fShapeFactors;
    std::vector < std::string > fSystematicNames;
    std::vector < std::string > fNormFactorNames;
    std::vector < std::string > fShapeFactorNames;

    std::vector<std::string> fPOIs;
    std::map<std::string,std::string> fPOIunit;
    std::string fPOIforLimit;
    std::string fPOIforSig;
    bool fUseStatErr;
    double fStatErrThres;
    std::string fStatErrCons;
    bool fUseGammaPulls;

    double fLumi;
    double fLumiScale;
    double fLumiErr;

    double fThresholdSystPruning_Normalisation;
    double fThresholdSystPruning_Shape;
    double fThresholdSystLarge;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleNames;
    std::string fMCweight;
    std::string fSelection;

    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixNamesNominal;
    std::vector<std::string> fAcceptanceNames;
    std::vector<std::string> fAcceptanceFiles;
    std::vector<std::string> fAcceptancePaths;
    std::vector<std::string> fAcceptanceNamesNominal;
    std::vector<std::string> fSelectionEffNames;
    std::vector<std::string> fSelectionEffFiles;
    std::vector<std::string> fSelectionEffPaths;
    std::vector<std::string> fSelectionEffNamesNominal;
    std::vector<std::string> fMigrationNames;
    std::vector<std::string> fMigrationFiles;
    std::vector<std::string> fMigrationPaths;
    std::vector<std::string> fMigrationNamesNominal;

    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNamesNominal;

    FitResults *fFitResults;

    bool fWithPullTables;

    int fIntCode_overall;
    int fIntCode_shape;

    int fInputType; // 0: histo, 1: ntup

    bool fSystDataPlot_upFrame;
    bool fStatOnly;
    bool fGammasInStatOnly;
    bool fStatOnlyFit;
    bool fFixNPforStatOnlyFit;

    std::vector<std::string> fRegionsToPlot;
    std::vector<std::string> fSummaryPlotRegions;
    std::vector<std::string> fSummaryPlotLabels;
    std::vector<std::string> fSummaryPlotValidationRegions;
    std::vector<std::string> fSummaryPlotValidationLabels;

    double fYmin;
    double fYmax;
    double fRatioYmin;
    double fRatioYmax;
    double fRatioYminPostFit;
    double fRatioYmaxPostFit;
    std::string fRatioYtitle;
    TRExPlot::RATIOTYPE fRatioType;

    std::string fLumiLabel;
    std::string fCmeLabel;

    std::string fSuffix;
    std::string fSaveSuffix;

    bool fUpdate;
    bool fKeepPruning;

    double fBlindingThreshold;
    Common::BlindingType fBlindingType;

    int fRankingMaxNP;
    std::string fRankingOnly;
    std::string fRankingPlot;
    std::string fImageFormat;
    std::string fAtlasLabel;

    bool fDoSummaryPlot;
    bool fDoMergedPlot;
    bool fDoTables;
    bool fDoSignalRegionsPlot;
    bool fDoPieChartPlot;

    std::string fGroupedImpactCategory;

    std::string fSummaryPrefix;

    //
    // Fit caracteristics
    //
    FitType fFitType;
    FitRegion fFitRegion;
    std::map< std::string, double > fFitNPValues;
    std::map< std::string, double > fFitFixedNPs;
    std::string fFitNPValuesFromFitResults;
    bool fInjectGlobalObservables;
    std::map< std::string, double > fFitPOIAsimov;
    bool fFitIsBlind;
    bool fUseRnd;
    double fRndRange;
    long int fRndSeed;
    std::vector<std::string> fVarNameLH;
    std::vector<std::vector<std::string> > fVarName2DLH;
    double fLHscanMin;
    double fLHscanMax;
    int fLHscanSteps;
    double fLHscanMinY;
    double fLHscanMaxY;
    int fLHscanStepsY;
    bool fParal2D;
    int fParal2Dstep;
    std::vector<std::string> fVarNameMinos;
    std::vector<std::string> fVarNameHide;
    std::string fWorkspaceFileName;
    bool fDoGroupedSystImpactTable;
    std::map<std::string, std::string> fSubCategoryImpactMap;

    //
    // Limit parameters
    //
    LimitType fLimitType;
    bool fLimitIsBlind;
    bool fSignalInjection;
    double fSignalInjectionValue;
    std::string fLimitParamName;
    double fLimitParamValue;
    std::string fLimitOutputPrefixName;
    double fLimitsConfidence;

    //
    // Significance parameters
    //
    bool fSignificanceIsBlind;
    bool fSignificanceDoInjection;
    double fSignificancePOIAsimov;
    std::string fSignificanceParamName;
    double fSignificanceParamValue;
    std::string fSignificanceOutputPrefixName;

    bool fCleanTables;
    bool fSystCategoryTables;

    std::vector< std::string > fRegionGroups;

    bool fKeepPrefitBlindedBins;
    std::unique_ptr<TH1D> fBlindedBins;

    std::string fCustomAsimov;

    std::string fTableOptions;

    bool fGetGoodnessOfFit;
    int fGetChi2;

    HistoTools::SmoothOption fSmoothOption;

    bool fSuppressNegativeBinWarnings;

    std::vector<std::string> fAddAliases;

    std::vector<std::string> fCustomFunctions;
    std::vector<std::string> fCustomIncludePaths;
    std::vector<std::string> fCustomFunctionsExecutes;

    std::vector<std::string> fMorphParams;
    std::vector<std::pair<double,std::string> > fTemplatePair;
    std::vector<TRExFit::TemplateWeight> fTemplateWeightVec;
    TemplateInterpolationOption fTemplateInterpolationOption;

    std::string fBootstrap;
    std::string fBootstrapSyst;
    std::string fBootstrapSample;
    int fBootstrapIdx;

    std::vector<std::string> fDecorrSysts;
    std::string fDecorrSuff;

    bool fDoNonProfileFit;
    double fNonProfileFitSystThreshold;
    int fFitToys;
    int fToysHistoNbins;
    std::string fToysPseudodataNP;
    double fToysPseudodataNPShift;
    std::string fSmoothMorphingTemplates;
    int fPOIPrecision;

    std::string fRankingPOIName;
    bool fUseATLASRounding;
    bool fUseATLASRoundingTxt;
    bool fUseATLASRoundingTex;
    bool fuseGammasForCorr;
    bool fPropagateSystsForMorphing;
    PruningType fPruningType;

    std::vector<int> fPrePostFitCanvasSize;
    std::vector<int> fSummaryCanvasSize;
    std::vector<int> fMergeCanvasSize;
    std::vector<int> fPieChartCanvasSize;
    std::vector<int> fNPRankingCanvasSize;

    std::vector<std::string> fBlindedParameters;

    double fLabelX;
    double fLabelY;
    double fLegendX1;
    double fLegendX2;
    double fLegendY;

    double fLabelXSummary;
    double fLabelYSummary;
    double fLegendX1Summary;
    double fLegendX2Summary;
    double fLegendYSummary;

    double fLabelXMerge;
    double fLabelYMerge;
    double fLegendX1Merge;
    double fLegendX2Merge;
    double fLegendYMerge;

    int fLegendNColumns;
    int fLegendNColumnsSummary;
    int fLegendNColumnsMerge;

    bool fShowRatioPad;
    bool fShowRatioPadSummary;
    bool fShowRatioPadMerge;

    std::string fExcludeFromMorphing;

    std::vector<std::string> fScaleSamplesToData;

    bool fSaturatedModel;

    bool fDoSystNormalizationPlots;

    int fDebugNev;
    
    int fCPU;

    std::vector< std::string > fSeparationPlot;
    FoldingManager::MATRIXORIENTATION fMatrixOrientation;

    std::string fTruthDistributionPath;
    std::string fTruthDistributionFile;
    std::string fTruthDistributionName;
    int fNumberUnfoldingTruthBins;
    int fNumberUnfoldingRecoBins;
    std::vector<std::unique_ptr<UnfoldingSample> > fUnfoldingSamples;
    std::vector<std::unique_ptr<UnfoldingSystematic> > fUnfoldingSystematics;
    double fUnfoldingResultMin;
    double fUnfoldingResultMax;
    bool fHasAcceptance;
    std::string fUnfoldingTitleX;
    std::string fUnfoldingTitleY;
    double fUnfoldingRatioYmax;
    double fUnfoldingRatioYmin;
    double fUnfoldingScaleRangeY;
    bool fUnfoldingLogX;
    bool fUnfoldingLogY;
    double fUnfoldingTitleOffsetX;
    double fUnfoldingTitleOffsetY;
    std::vector<std::unique_ptr<TruthSample> > fTruthSamples;
    std::string fNominalTruthSample;
    std::string fAlternativeAsimovTruthSample;
    std::string fMigrationTitleX;
    std::string fMigrationTitleY;
    bool fMigrationLogX;
    bool fMigrationLogY;
    double fMigrationTitleOffsetX;
    double fMigrationTitleOffsetY;
    bool fPlotSystematicMigrations;
    double fMigrationZmin;
    double fMigrationZmax;
    double fResponseZmin;
    double fResponseZmax;
    bool fMigrationText;
    bool fUnfoldingDivideByBinWidth;
    double fUnfoldingDivideByLumi;
    int fRegularizationType;
    PruningUtil::SHAPEOPTION fPruningShapeOption;
    bool fSummaryLogY;
    /// This variable is needed only for multifit
    bool fUseInFit;
    bool fUseInComparison;
    bool fReorderNPs;
    bool fBlindSRs;
    bool fHEPDataFormat;
    bool fAlternativeShapeHistFactory;
    int fFitStrategy;
    bool fBinnedLikelihood;
    bool fRemoveLargeSyst;
    bool fRemoveSystOnEmptySample;
    bool fValidationPruning;
    bool fUnfoldNormXSec;
    int fUnfoldNormXSecBinN;
    bool fUsePOISinRanking;
    bool fUseHesseBeforeMigrad;
    bool fUseNllInLHscan;
    int fLimitToysStepsSplusB;
    int fLimitToysStepsB;
    int fLimitToysScanSteps;
    double fLimitToysScanMin;
    double fLimitToysScanMax;
    int fToysSeed;
    int fLimitToysSeed;
    bool fLimitPlot;
    bool fLimitFile;
    bool fDataWeighted;
};

#endif
