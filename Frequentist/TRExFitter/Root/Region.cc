// Class include
#include "TRExFitter/Region.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/YamlConverter.h"

// ROOT includes
#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TMatrixD.h"
#include "TPaveText.h"
#include "TSystem.h"

// c++ includes
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

// -------------------------------------------------------------------------------------------------
// class Region

//__________________________________________________________________________________
//
Region::Region(const string& name) :
    fName(name),
    fVariableTitle(""),
    fYTitle(""),
    fLabel(name),
    fShortLabel(name),
    fTexLabel(""),
    fFitName(""),
    fRegionType(CONTROL),
    fRegionDataType(REALDATA),
    fHasData(false),
    fData(nullptr),
    fHasSig(false),
    fYmaxScale(0),
    fYmin(0),
    fYmax(0),
    fRatioYmin(0),
    fRatioYmax(2.),
    fRatioYminPostFit(0.5),
    fRatioYmaxPostFit(1.5),
    fRatioYtitle(""),
    fRatioType(TRExPlot::RATIOTYPE::DATAOVERMC),
    fTot(nullptr),
    fErr(nullptr),
    fTot_postFit(nullptr),
    fErr_postFit(nullptr),
    fBinTransfo(""),
    fTransfoDzBkg(0.),
    fTransfoDzSig(0.),
    fTransfoFzBkg(0.),
    fTransfoFzSig(0.),
    fTransfoJpar1(0.),
    fTransfoJpar2(0.),
    fTransfoJpar3(0.),
    fVariable(""),
    fCorrVar1(""),
    fCorrVar2(""),
    fNbins(0),
    fXmin(0),
    fXmax(0),
    fSelection("1"),
    fMCweight("1"),
    fHistoNBinsRebin(-1),
    fHistoNBinsRebinPost(-1),
    fPlotPreFit(nullptr),
    fPlotPostFit(nullptr),
    fUseStatErr(false),
    fIntCode_overall(4),
    fIntCode_shape(0),
    fFitType(TRExFit::SPLUSB),
    fFitLabel(""),
    fLumiLabel(""),
    fCmeLabel(""),
    fLumiScale(1.),
    fLogScale(false),
    fBinWidth(0),
    fBlindingThreshold(-1),
    fBlindingType(Common::SOVERB),
    fSkipSmoothing(false),
    fATLASlabel("Internal"),
    fSuffix(""),
    fGroup(""),
    fKeepPrefitBlindedBins(false),
    fGetChi2(0),
    fChi2val(-1),
    fNDF(-1),
    fChi2prob(-1),
    fUseGammaPulls(false),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fLegendNColumns(2),
    fNumberUnfoldingRecoBins(0),
    fNormalizeMigrationMatrix(true),
    fHasAcceptance(false),
    fFolder(""),
    fHEPDataFormat(false),
    fAutomaticDropBins(false) {


    int canvasWidth = 600;
    int canvasHeight = 700;
    std::string cName = "c_"+fName;
    if(TRExFitter::OPTION["CanvasWidth"]!=0)  canvasWidth  = TRExFitter::OPTION["CanvasWidth"];
    if(TRExFitter::OPTION["CanvasHeight"]!=0) canvasHeight = TRExFitter::OPTION["CanvasHeight"];
    fPlotPreFit = std::make_shared<TRExPlot>(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);
    fPlotPreFit->fShowYields = TRExFitter::SHOWYIELDS;
    cName = "c_"+fName+"_postFit";
    fPlotPostFit = std::make_shared<TRExPlot>(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);
    fPlotPostFit->fShowYields = TRExFitter::SHOWYIELDS;
}

//__________________________________________________________________________________
//
Region::~Region(){
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::SetSampleHist(Sample *sample, string histoName, string fileName){
    fSampleHists.emplace_back(new SampleHist( sample, histoName, fileName ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists.back();
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else if(sample->fType==Sample::EFT){
        WriteDebugStatus("Region::SetSampleHist", "Adding EFT sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists.back()->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists.back()->fRegionName = fName;
    fSampleHists.back()->fRegionLabel = fLabel;
    fSampleHists.back()->fFitName = fFitName;
    fSampleHists.back()->fVariableTitle = fVariableTitle;
    return fSampleHists.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::SetSampleHist(Sample *sample, TH1* hist ){
    fSampleHists.emplace_back(new SampleHist( sample, hist ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists.back();
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else if(sample->fType==Sample::EFT){
        WriteDebugStatus("Region::SetSampleHist", "Adding EFT sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists.back()->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists.back()->fRegionName = fName;
    fSampleHists.back()->fRegionLabel = fLabel;
    fSampleHists.back()->fFitName = fFitName;
    fSampleHists.back()->fVariableTitle = fVariableTitle;
    return fSampleHists.back();
}

//__________________________________________________________________________________
//
void Region::AddSample(Sample* sample){
    fSamples.emplace_back(std::move(sample));
}

//__________________________________________________________________________________
//
void Region::SetBinning(int N, double *bins){
    fNbins = fHistoNBinsRebin = N;
    fHistoBins.resize(N+1);
    for(int i=0; i<=N; ++i) fHistoBins[i] = bins[i];
}

//__________________________________________________________________________________
//
void Region::Rebin(int N){
    fNbins = fHistoNBinsRebin = N;
}

//__________________________________________________________________________________
//
void Region::SetRebinning(int N, double *bins){
    fHistoNBinsRebinPost = N;
    fHistoBinsPost.resize(N+1);
    for(int i=0; i<=N; ++i) fHistoBinsPost[i] = bins[i];
}

//__________________________________________________________________________________
//
void Region::SetRegionType( RegionType type ){
    fRegionType = type;
}

//__________________________________________________________________________________
//
void Region::SetRegionDataType( DataType type ){
    fRegionDataType = type;
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::GetSampleHist(const std::string &sampleName) const{
    for(const auto& isample : fSampleHists) {
        if(isample->fName == sampleName) return isample;
    }
    return nullptr;
}

//__________________________________________________________________________________
//
void Region::BuildPreFitErrorHist(){
    WriteInfoStatus("Region::BuildPreFitErrorHist", "Building pre-fit plot for region " + fName + " ...");
    //
    fSystNames.clear();
    fNpNames.clear();
    std::map<std::string,std::string> origShapeSystName;
    std::map<std::string,bool> systIsThere;

    //
    // Collect all the systematics on all the samples
    // (in parallel also colect all the nuisance parameters)
    //
    for(const auto& isample : fSampleHists) {
        if(isample->fSample->fType == Sample::DATA) continue;
        if(isample->fSample->fType == Sample::GHOST) continue;
        if(isample->fSample->fType == Sample::EFT) continue;

        //
        // non-SHAPE Systematics
        //
        for(std::size_t i_syst = 0; i_syst < isample->fSyst.size(); ++i_syst){
            const std::string systName = isample->fSyst[i_syst]->fName;
            if(isample->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
            std::string systNuisPar;
            if(isample->fSyst[i_syst]->fSystematic!=nullptr) {
                systNuisPar = isample->fSyst[i_syst]->fSystematic->fNuisanceParameter;
            }
            if(Common::FindInStringVector(fNpNames,systNuisPar)<0){
                fNpNames.push_back(systNuisPar);
            }
        }

        //
        // SHAPE Systematics
        // (but remove the MC-stat ones (for samples with fSeparateGammas)!!
        //
        for(std::size_t i_syst = 0; i_syst < isample->fSyst.size(); ++i_syst) {
            const std::string systName = isample->fSyst[i_syst]->fName;
            if(isample->fSyst[i_syst]->fSystematic->fType!=Systematic::SHAPE) continue;
            if(systName.find("stat_")!=std::string::npos) continue;
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("shape_%s_%s_bin_%d",systName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                    origShapeSystName[gammaName] = systName;
                    TRExFitter::NPMAP[gammaName] = gammaName;
                }
                std::string systNuisPar;
                if(isample->fSyst[i_syst]->fSystematic!=nullptr) {
                    systNuisPar = gammaName;
                }
                if(Common::FindInStringVector(fNpNames,systNuisPar)<0){
                    fNpNames.push_back(systNuisPar);
                }
            }
        }
    }

    //
    // Build pre-fit error hists (for each sample and each systematic)
    //

    // - loop on samples
    //
    for(auto& isample : fSampleHists) {

        // skip data
        if(isample->fSample->fType==Sample::DATA) continue;
        if(isample->fSample->fType==Sample::GHOST) continue;
        if(isample->fSample->fType==Sample::EFT) continue;
        if(isample->fSample->fType==Sample::SIGNAL && (!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) || fRatioType==TRExPlot::RATIOTYPE::DATAOVERB)) continue;

        WriteDebugStatus("Region::BuildPreFitErrorHist", "  Sample: " + isample->fName);

        TH1* hNom = isample->fHist.get();

        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            WriteDebugStatus("Region::BuildPreFitErrorHist", "    Systematic: " + fSystNames[i_syst]);
            const std::string systName = fSystNames[i_syst];

            // get SystematicHist
            std::shared_ptr<SystematicHist> sh = isample->GetSystematic(systName);

            // hack: add a systematic hist if not there... FIXME
            if(sh==nullptr){
                isample->AddHistoSyst(systName,systName,hNom,hNom);
                sh = isample->GetSystematic(systName);
                // initialize the up and down variation histograms
                // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
                sh->fHistUp.reset(static_cast<TH1*>(hNom->Clone(Form("%s_%s_Up",  hNom->GetName(),systName.c_str()))));
                sh->fHistDown.reset(static_cast<TH1*>(hNom->Clone(Form("%s_%s_Down",hNom->GetName(),systName.c_str()))));
                sh->fHistUp->SetDirectory(nullptr);
                sh->fHistDown->SetDirectory(nullptr);
            }

            std::shared_ptr<Systematic> syst = sh->fSystematic;

            // store hist up and down
            TH1* hUp   = sh->fHistUp.get();
            TH1* hDown = sh->fHistDown.get();

            // modify them dropping shape or norm (due to pruning or shape/acc decorrelation)
            if(syst) {
                if(syst->fIsNormOnly || sh->fShapePruned) {
                    Common::DropShape(hUp,hDown,hNom);
                } 
                if(syst->fIsShapeOnly || sh->fNormPruned) {
                    Common::DropNorm(hUp,hDown,hNom);
                }
            }

            //
            // - loop on bins
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                double diffUp(0.);
                double diffDown(0.);
                double yieldUp(0.);
                double yieldDown(0.);
                const double yieldNominal = hNom->GetBinContent(i_bin);  // store nominal yield for this bin
                // if it's a systematic (NB: skip Norm-Factors!!)
                if(isample->HasSyst(fSystNames[i_syst])){
                    if(hUp!=nullptr)    yieldUp     = hUp  ->GetBinContent(i_bin);
                    else                yieldUp     = yieldNominal;
                    if(hDown!=nullptr)  yieldDown   = hDown->GetBinContent(i_bin);
                    else                yieldDown   = yieldNominal;
                    diffUp   += yieldUp   - yieldNominal;
                    diffDown += yieldDown - yieldNominal;
                }
                // for shape systematics
                if(origShapeSystName[systName]!="" && systName.find(Form("%s_bin_%d",fName.c_str(),i_bin-1))!=std::string::npos){
                    std::shared_ptr<SystematicHist> shOrig = isample->GetSystematic(origShapeSystName[systName]);
                    if(!shOrig) continue;
                    yieldUp   = shOrig->fHistUp->GetBinContent(i_bin);
                    yieldDown = 2*yieldNominal - yieldUp;
                    sh->fHistUp  ->SetBinContent(i_bin,yieldUp);
                    sh->fHistDown->SetBinContent(i_bin,yieldDown);
                    diffUp   += yieldUp   - yieldNominal;
                    diffDown += yieldDown - yieldNominal;
                }
                WriteDebugStatus("Region::BuildPreFitErrorHist", "        Bin " + std::to_string(i_bin) + ":  " + " \t +" + std::to_string(100*diffUp/yieldNominal)
                 + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");
            }
        }
    }

    // at this point all the sample-by-sample pre-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    // - loop on systematics
    fTotUp.resize(fSystNames.size());
    fTotDown.resize(fSystNames.size());
    for(std::size_t i_syst=0; i_syst<fSystNames.size(); ++i_syst){
        const std::string systName = fSystNames[i_syst];

        // initialize the tot variation hists
        fTotUp[i_syst].reset(static_cast<TH1*>(fTot->Clone(Form("h_%s_tot_%s_Up",  fName.c_str(), systName.c_str()))));
        fTotDown[i_syst].reset(static_cast<TH1*>(fTot->Clone(Form("h_%s_tot_%s_Down",fName.c_str(), systName.c_str()))));
        fTotUp[i_syst]->SetDirectory(nullptr);
        fTotDown[i_syst]->SetDirectory(nullptr);
        // - loop on bins
        for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
            // - loop on samples
            double diffUp(0.);
            double diffDown(0.);
            for(const auto& isample : fSampleHists) {
                TH1* hNom = isample->fHist.get();
                //
                // scale according to NormFactors
                double scale = 1.;
                for(unsigned int i_nf=0;i_nf<isample->fSample->fNormFactors.size();i_nf++){
                    const NormFactor *nf = isample->fSample->fNormFactors[i_nf].get();
                    // if this norm factor is a morphing one
                    if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                        const std::vector<double> scales = Common::CalculateExpression(nf, nf->fName, false, nullptr, nullptr);
                        if (scales.empty()) {
                            WriteErrorStatus("Region::BuildPreFitErrorHist", "Scale vector is empty");
                            exit(EXIT_FAILURE);
                        } 
                        scale *= scales.at(0);
                    }
                    else {
                        if (std::find(nf->fRegions.begin(), nf->fRegions.end(), fName) != nf->fRegions.end()) {
                            scale *= isample->fSample->fNormFactors[i_nf]->fNominal;
                        }
                    }
                }
                //
                // skip data
                if(isample->fSample->fType==Sample::DATA) continue;
                if(isample->fSample->fType==Sample::GHOST) continue;
                if(isample->fSample->fType==Sample::EFT) continue;
                if(isample->fSample->fType==Sample::SIGNAL && (!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) || fRatioType==TRExPlot::RATIOTYPE::DATAOVERB)) continue;
                // get SystematicHist
                std::shared_ptr<SystematicHist> sh = isample->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                const double yieldNominal = hNom->GetBinContent(i_bin);
                diffUp   += (sh->fHistUp  ->GetBinContent(i_bin) - yieldNominal)*scale;
                diffDown += (sh->fHistDown->GetBinContent(i_bin) - yieldNominal)*scale;
            }
            // add the proper bin content to the variation hists
            fTotUp[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }

    WriteDebugStatus("Region::BuildPreFitErrorHist", "----");

    //
    // build the vectors of variations (sum histograms for systematics with the same NP)
    std::vector< std::shared_ptr<TH1> > h_up;
    std::vector< std::shared_ptr<TH1> > h_down;
    for(size_t i_syst=0;i_syst< fSystNames.size(); ++i_syst){
        // look for all the already stored systematics, to find if one had the same NP
        bool found = false;
        for(size_t j_syst=0;j_syst<i_syst;++j_syst){
            if(TRExFitter::NPMAP[fSystNames[i_syst]]==TRExFitter::NPMAP[fSystNames[j_syst]]){
                found = true;
                const int whichsyst = Common::FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]]);
                if (whichsyst < 0) {
                    WriteErrorStatus("Region::BuildPreFitErrorHist", "Systematic not found in the list of systematics. This should not happen...");
                    exit(EXIT_FAILURE);
                }
                auto h_diff_up   = std::unique_ptr<TH1> (static_cast<TH1*> (h_up[whichsyst]  ->Clone(Form("%s_%s","clone_",h_up[whichsyst]  ->GetName()))));
                auto h_diff_down = std::unique_ptr<TH1> (static_cast<TH1*> (h_down[whichsyst]->Clone(Form("%s_%s","clone_",h_down[whichsyst]->GetName()))));
                h_diff_up  ->Add(fTotUp[  i_syst].get(),fTot.get(),1,-1);
                h_diff_down->Add(fTotDown[i_syst].get(),fTot.get(),1,-1);
                h_up[   Common::FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_up.get());
                h_down[ Common::FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_down.get());
                break;
            }
        }
        if(found) continue;
        //
        // if shape-only systematics to show, normalize each syst variation
        if(TRExFitter::OPTION["ShapeOnlySystBand"]>0){
            fTotUp[i_syst]  ->Scale(fTot->Integral()/fTotUp[i_syst]  ->Integral());
            fTotDown[i_syst]->Scale(fTot->Integral()/fTotDown[i_syst]->Integral());
        }
        h_up.  emplace_back(fTotUp[i_syst] );
        h_down.emplace_back(fTotDown[i_syst]);
    }
    fErr = BuildTotError( fTot.get(), h_up, h_down, fNpNames );
    fErr->SetName("g_totErr");
    // at this point fTot and fErr should be ready

    //
    // Goodness of pre-fit
    if(fGetChi2){
        // remove blinded bins
        if (!fHasData){
            WriteWarningStatus("Region::BuildPreFitErrorHist", "Data histogram is nullptr, cannot calculate Chi2 agreement.");
            WriteWarningStatus("Region::BuildPreFitErrorHist", "Maybe you do not have data sample defined?");
            return;
        }
        std::unique_ptr<TH1> h_data(static_cast<TH1*>(fData->fHist->Clone()));
        if(fGetChi2==1) fNpNames.clear();
        std::pair<double,int> res = GetChi2Test(false, h_data.get(), fTot.get(), h_up, fNpNames );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- PRE-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }
}

//__________________________________________________________________________________
//
std::shared_ptr<TRExPlot> Region::DrawPreFit(const std::vector<int>& canvasSize, string opt){

    std::shared_ptr<TRExPlot> p(nullptr);
    if (canvasSize.size() == 0){
        int canvasWidth = 600;
        int canvasHeight = 700;
        const std::string cName = "c_"+fName;
        if(TRExFitter::OPTION["CanvasWidth"]!=0)  canvasWidth  = TRExFitter::OPTION["CanvasWidth"];
        if(TRExFitter::OPTION["CanvasHeight"]!=0) canvasHeight = TRExFitter::OPTION["CanvasHeight"];
        p = std::make_shared<TRExPlot>(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);
        p->fShowYields = TRExFitter::SHOWYIELDS;
    } else {
        p = std::make_shared<TRExPlot>(("c_"+fName).c_str(), canvasSize.at(0), canvasSize.at(1),TRExFitter::NORATIO);
        p->fShowYields = TRExFitter::SHOWYIELDS;
    }
    p->SetXaxisRange(fXaxisRange);
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmax;
    p->fRatioYmin = fRatioYmin;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TRExFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        if (fLabel != "none") p->AddLabel(fLabel);
        p->AddLabel("#font[52]{Pre-fit}");
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        if (fLabel != "none") p->AddLabel(fLabel);
        if(TRExFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Pre-Fit");
    }
    YamlConverter::PlotContainer container;
    container.region = fName;
    container.xAxis = fVariableTitle;
    container.yAxis = fYTitle;

    //
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLegendNColumns = fLegendNColumns;
    if(fBlindingThreshold>=0) {
        container.blindedBins = fBlindedBins;
        p->SetBinBlinding(fBlindedBins);
    }
    p->SetDropBins(fDropBins);

    if (!fBinLabels.empty()) {
        p->ResizeBinLabel(fBinLabels.size() + 1);
        for(std::size_t i_bin=0; i_bin<fBinLabels.size(); i_bin++) {
            p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
        }
    }

    //
    // build h_tot
    //
    fTot.reset(nullptr);
    if(fHasData && opt.find("blind")==string::npos) {
        for (int ibin = 1; ibin <= fData->fHist->GetNbinsX(); ++ibin) {
            container.data.emplace_back(fData->fHist->GetBinContent(ibin));
        }
        p->SetData(fData->fHist.get(),fData->fSample->fTitle);
    }
    for(const auto& isig : fSig) {
        container.samples.emplace_back(isig->fSample->fTitle);
        std::string title = isig->fSample->fTitle;
        if(isig->fSample->fGroup != "") title = isig->fSample->fGroup;
        std::unique_ptr<TH1> h(static_cast<TH1*>(isig->fHist->Clone()));
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!isig->fSample->fUseMCStat && !isig->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf < isig->fSample->fNormFactors.size();i_nf++){
            const NormFactor *nf = isig->fSample->fNormFactors[i_nf].get();
            if (!nf->fRegions.empty() && nf->fRegions.at(0) != "none" && Common::FindInStringVector(nf->fRegions, fName) < 0) continue;
            if (Common::FindInStringVector(nf->fExclude, fName) >= 0) continue;
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                const std::vector<double> scales = Common::CalculateExpression(nf, nf->fName, false, nullptr, nullptr);
                if (scales.empty()) {
                    WriteErrorStatus("Region::DrawPreFit", "Scales are empty");
                    exit(EXIT_FAILURE);
                }
                const double scale = scales.at(0);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + isig->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                if (std::find(nf->fRegions.begin(), nf->fRegions.end(), fName) != nf->fRegions.end()) {
                    h->Scale(nf->fNominal);
                    WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + isig->fSample->fName + " by " + std::to_string(isig->fSample->fNormFactors[i_nf]->fNominal));
                }
            }
        }
        if(TRExFitter::SHOWSTACKSIG) p->AddSignal(h.get(),title);
        if(TRExFitter::SHOWNORMSIG) {
            if( (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
             || !TRExFitter::OPTION["NormSigSRonly"] )
                p->AddNormSignal(h.get(),title);
        }
        else{
            if(TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(h.get(),title);
        }
        std::vector<double> tmp;
        for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
            tmp.emplace_back(h->GetBinContent(ibin));
        }

        if(TRExFitter::SHOWOVERLAYSIG) p->AddOverSignal(h.get(),title);
        if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG){
            if(fTot==nullptr) fTot.reset(static_cast<TH1*>(h->Clone("h_tot")));
            else              fTot->Add(h.get());
            fTot->SetDirectory(nullptr);
        }
        container.signalYields.emplace_back(std::move(tmp));
    }
    for(const auto& ibkg : fBkg) {
        container.samples.emplace_back(ibkg->fSample->fTitle);
        std::string title = ibkg->fSample->fTitle;
        if(ibkg->fSample->fGroup != "") title = ibkg->fSample->fGroup;
        std::unique_ptr<TH1> h(static_cast<TH1*>(ibkg->fHist->Clone()));
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!ibkg->fSample->fUseMCStat && !ibkg->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf < ibkg->fSample->fNormFactors.size();i_nf++){
            const NormFactor *nf = ibkg->fSample->fNormFactors[i_nf].get();
            if (!nf->fRegions.empty() && nf->fRegions.at(0) != "none" && Common::FindInStringVector(nf->fRegions, fName) < 0) continue;
            if (Common::FindInStringVector(nf->fExclude, fName) >= 0) continue;
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                const std::vector<double> scales = Common::CalculateExpression(nf, nf->fName, false, nullptr, nullptr);
                if (scales.empty()) {
                    WriteErrorStatus("Region::DrawPreFit", "Scales are empty");
                    exit(EXIT_FAILURE);
                }
                const double scale = scales.at(0);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + ibkg->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                if (std::find(nf->fRegions.begin(), nf->fRegions.end(), fName) != nf->fRegions.end()) {
                    h->Scale(nf->fNominal);
                    WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + ibkg->fSample->fName + " by " + std::to_string(ibkg->fSample->fNormFactors[i_nf]->fNominal));
                }
            }
        }
        std::vector<double> tmp;
        for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
            tmp.emplace_back(h->GetBinContent(ibin));
        }
        p->AddBackground(h.get(),title);
        if(fTot==nullptr) fTot.reset(static_cast<TH1*>(h->Clone("h_tot")));
        else          fTot->Add(h.get());
        fTot->SetDirectory(nullptr);
        container.backgroundYields.emplace_back(std::move(tmp));
    }

    //
    // set error to 0 if no MCstat
    //
    if(!fUseStatErr){
        for(int i_bin=1;i_bin<=fTot->GetNbinsX();i_bin++){
            fTot->SetBinError(i_bin,0);
        }
    }
    else{
        for(int i_bin=0;i_bin<fTot->GetNbinsX()+2;i_bin++){
            if(fTot->GetBinError(i_bin)>10*fTot->GetBinContent(i_bin)) fTot->SetBinError(i_bin,fTot->GetBinContent(i_bin));
        }
    }

    p->SetTotBkg(fTot.get());
    p->BlindData();
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);

    //
    // Computes the uncertainty bands arround the h_tot histogram
    //
    BuildPreFitErrorHist();

    //
    // Save the uncertainty band and h_tot to {region}{suffix}_preFit.root file.
    //
    SavePreFitUncertaintyAndTotalMCObjects();

    //
    // Print chi2 info
    //
    if(fGetChi2 && TRExFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);

    //
    // Sets the last ingredients in the TRExPlot object
    //
    container.errors = fErr.get();
    p->SetTotBkgAsym(fErr.get());
    p->fATLASlabel = fATLASlabel;
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType==TRExPlot::RATIOTYPE::DATAOVERMC){
        p->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
    }
    p->fLabelX = fLabelX;
    p->fLabelY = fLabelY;
    p->fLegendX1 = fLegendX1;
    p->fLegendX2 = fLegendX2;
    p->fLegendY = fLegendY;
    if(fLogScale) opt += " log";
    p->Draw(opt);

    // Yaml
    YamlConverter converter{};
    converter.WritePlot(container, fFolder, false);
    if (fHEPDataFormat) {
        converter.WritePlotHEPData(container, fFolder, false);
    }

    if(TRExFitter::PREFITONPOSTFIT) {
        fPlotPreFit = p;
    }

    return p;
}

//__________________________________________________________________________________
//
double Region::GetMultFactors( FitResults* fitRes,
                               std::ofstream& pullTex,
                               const int i /*sample*/,
                               const int i_bin /*bin number*/,
                               const double binContent0,
                               const std::string &var_syst_name,
                               const bool isUp ) const{
    double multNorm(1.);
    double multShape(0.);
    double systValue(0.);
    const SampleHist *sh = fSampleHists[i].get();
    for(std::size_t i_syst = 0; i_syst < sh->fSyst.size(); ++i_syst){
        std::shared_ptr<SystematicHist> syh = sh->fSyst[i_syst];
        std::string systName = syh->fName;
        TString systNameNew(systName); // used in pull tables
        std::shared_ptr<Systematic> syst = syh->fSystematic;
        bool isOverall = syh->fIsOverall && !syh->fNormPruned;
        bool isShape   = syh->fIsShape && !syh->fShapePruned;
        if(syst){
            systName = syst->fNuisanceParameter;
            if(syst->fIsShapeOnly) isOverall = false;
            if(syst->fIsNormOnly)  isShape   = false;
            if(syst->fType==Systematic::SHAPE){
                continue;
            }
            else {
                if(var_syst_name!="" && systName==var_syst_name){
                    if(isUp) systValue = fitRes->GetNuisParValue(systName) + fitRes->GetNuisParErrUp(systName);
                    else     systValue = fitRes->GetNuisParValue(systName) + fitRes->GetNuisParErrDown(systName);
                }
                else {
                    systValue = fitRes->GetNuisParValue(systName);
                }
            }
        }
        else{
            systValue = fitRes->GetNuisParValue(systName);
        }

        //
        // Normalisation component: use the exponential interpolation and the multiplicative combination
        //
        if(isOverall){
            const double binContentUp   = (syh->fNormUp+1) * binContent0;
            const double binContentDown = (syh->fNormDown+1) * binContent0;
            const double factor = GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall);
            multNorm *= factor;
            if (fSampleHists[i]->fSample->fBuildPullTable>0){
                if (((factor > 1.01) || (factor < 0.99)) && (i_bin==1)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has norm effect " + std::to_string( factor*100 ));
                    pullTex << setprecision(2) << "norm " << systNameNew.ReplaceAll("_","-") << "&" << (factor-1)*100 << " \\% \\\\\n" << std::endl;
                }
            }
        }

        //
        // Shape component: use the linear interpolation and the additive combination
        //
        if(isShape){
            const double binContentUp   = syh->fHistShapeUp->GetBinContent(i_bin);
            const double binContentDown = syh->fHistShapeDown->GetBinContent(i_bin);
            const double factor = GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape);
            multShape += factor - 1;
            if (fSampleHists[i]->fSample->fBuildPullTable==2){
                if (((factor-1) > 0.03) || ((factor-1) < - 0.03)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has shape effect " + std::to_string( (factor-1)*100 ));
                    pullTex << setprecision(2) << "shape " << systNameNew.ReplaceAll("_","-") << " bin " << i_bin << "&" << (factor-1)*100  << " \\% \\\\\n" << std::endl;
                }
            }
        }
    }
    return multNorm*(multShape+1.);
}

//__________________________________________________________________________________
//
void Region::BuildPostFitErrorHist(FitResults *fitRes, const std::vector<std::string>& morph_names){

    WriteInfoStatus("Region::BuildPostFitErrorHist", "Building post-fit plot for region " + fName + " ...");

    //
    // 0) Collect all the systematics on all the samples
    //
    fSystNames.clear();
    std::map<string,bool> systIsThere;


    for(const auto& isample : fSampleHists) {
        if(isample->fSample->fType == Sample::DATA) continue;
        if(isample->fSample->fType == Sample::GHOST) continue;
        if(isample->fSample->fType == Sample::EFT) continue;

        //
        // Norm factors
        //
        for(const auto& inorm : isample->fSample->fNormFactors) {
            const std::string systName = inorm->fName;
            // if this norm factor is a morphing one => save the nuis.par
            if(fFitType==TRExFit::BONLY && Common::FindInStringVector(fPOIs,systName)>0) continue;
            if(inorm->fConst) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }

        //
        // Shape factors
        //
        // extract number of bins
        // loop over shape factors
        for(const auto& ishape : isample->fSample->fShapeFactors) {
            const std::string systName = ishape->fName;
            // add syst name for each bin
            for(int i_bin = 0; i_bin < isample->fHist->GetNbinsX(); i_bin++){
                const std::string systNameSF = systName + "_bin_" + std::to_string(i_bin);
                // the shape factor naming used i_bin - 1 for the first bin
                // add it as one syst per bin
                if(!systIsThere[systNameSF]){
                    fSystNames.push_back(systNameSF);
                    systIsThere[systNameSF] = true;
                }
            }
        }

        //
        // Systematics
        //
        for(std::size_t i_syst = 0; i_syst < isample->fSyst.size(); ++i_syst) {
            if(!isample->fSyst[i_syst]->fSystematic) continue;
            const std::string systName = isample->fSyst[i_syst]->fName;
            if(isample->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }

        //
        // SHAPE Systematics
        //
        for(std::size_t i_syst = 0; i_syst < isample->fSyst.size(); ++i_syst) {
            if(!isample->fSyst[i_syst]->fSystematic) continue;
            const std::string systName = isample->fSyst[i_syst]->fName;
            if(systName.find("stat_")!=std::string::npos) continue; // fSeparateGammas already added later
            if(isample->fSyst[i_syst]->fSystematic->fType!=Systematic::SHAPE) continue;
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                const std::string gammaName = Form("shape_%s_%s_bin_%d",systName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                }
            }
        }

        //
        // Gammas
        //
        if(fUseGammaPulls && (isample->fSample->fUseMCStat || isample->fSample->fSeparateGammas)){
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(isample->fSample->fSeparateGammas) {
                    gammaName = Form("shape_stat_%s_%s_bin_%d",isample->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                }
                if(!systIsThere[gammaName] && (fitRes->GetNuisParValue(gammaName)>0)){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                }
            }
        }
    }

    //
    // 1) Build post-fit error hists (for each sample and each systematic):
    //
    std::vector<double> morph_scale;
    std::vector<double> morph_scale_nominal;
    PrepareMorphScales(fitRes, &morph_scale, &morph_scale_nominal);

    // - loop on systematics
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        WriteVerboseStatus("Region::BuildPostFitErrorHist", "    Systematic: " + fSystNames[i_syst]);

        int i_morph_sample = 0;
        //
        // Get fit result
        //
        const std::string systName = fSystNames[i_syst];
        if(systName.find("saturated_model_")!=std::string::npos) continue;
        if(TRExFitter::NPMAP[systName]=="") TRExFitter::NPMAP[systName] = systName;

        // Before checking if a systematic is there in the fit results, needs first to identify which name to look for:
        // - use NuisanceParameer for systematics (and normal norm factors)
        // - use the name and NOT the NPMAP for morphing factors (NPMAP contains the morphing parameter)
        std::string systToCheck = systName;
        if(systName.find("morph_")==std::string::npos){
            systToCheck = TRExFitter::NPMAP[systName];
        }
        const double systValue   = fitRes->GetNuisParValue(systToCheck);
        const double systErrUp   = fitRes->GetNuisParErrUp(systToCheck);
        const double systErrDown = fitRes->GetNuisParErrDown(systToCheck);

        WriteVerboseStatus("Region::BuildPostFitErrorHist", "      alpha = " + std::to_string(systValue) + " +" + std::to_string(systErrUp) + " " + std::to_string(systErrDown));

        // needed for morph samples
        const int nbins = fTot_postFit->GetNbinsX();

        std::vector<double> morph_nominal(nbins, 0.0);
        std::vector<double> morph_nominal_postfit(nbins, 0.0);
        std::vector<double> morph_up_postfit(nbins, 0.0);
        std::vector<double> morph_down_postfit(nbins, 0.0);

        std::vector<double> morph_syst_up(nbins, 0.0);
        std::vector<double> morph_syst_down(nbins, 0.0);

        bool isMorph = false;
        int morph_index = -1;

        // - loop on samples
        for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
            // skip data
            if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
            if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
            if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
            WriteVerboseStatus("Region::BuildPostFitErrorHist", "  Sample: " + fSampleHists[i]->fName);

            //
            // Get SystematicHist
            //
            std::shared_ptr<SystematicHist> sh = fSampleHists[i]->GetSystematic(systName);

            // hack: add a systematic hist if not there
            if(!sh){
                fSampleHists[i]->AddHistoSyst(systName,systName,fSampleHists[i]->fHist.get(),fSampleHists[i]->fHist.get());
                sh = fSampleHists[i]->GetSystematic(systName);
            }

            //
            // initialize the up and down variation histograms
            // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
            //
            sh->fHistUp_postFit.reset(static_cast<TH1*>(fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Up_postFit",  fSampleHists[i]->fHist->GetName(),systName.c_str()))));
            sh->fHistDown_postFit.reset(static_cast<TH1*>(fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Down_postFit",fSampleHists[i]->fHist->GetName(),systName.c_str()))));
            sh->fHistUp_postFit->SetDirectory(nullptr);
            sh->fHistDown_postFit->SetDirectory(nullptr);

            // check if the sample is morph sample
            for (const auto& i_morph : morph_names){
                if (fSampleHists[i]->fIsMorph[i_morph]){
                    isMorph = true;
                    morph_index = i;
                    break;
                }
            }

            // - loop on bins
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                double diffUp(0.);
                double diffDown(0.);
                const double yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                const double yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post fit

                if (isMorph){
                    const double scaleNom  = morph_scale_nominal.at(i_morph_sample);
                    const double scaleNom_postfit  = morph_scale.at(i_morph_sample);
                    morph_nominal_postfit.at(i_bin-1)+= yieldNominal*scaleNom_postfit;
                    morph_nominal.at(i_bin-1)+= yieldNominal*scaleNom;
                }

                const size_t posTmp = systName.find("_bin_");
                const std::string gammaName      = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                const std::string gammaNameShape = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                //
                // if it's a gamma
                if(gammaName==fSystNames[i_syst] && fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                    if (isMorph){
                        morph_up_postfit.at(i_bin-1)  += yieldNominal_postFit*(systErrUp+1);
                        morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                    }
                }
                //
                // if it's a specific-sample gamma
                else if(gammaNameShape==fSystNames[i_syst] && fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                    if (isMorph){
                        morph_up_postfit.at(i_bin-1)  += yieldNominal_postFit*(systErrUp+1);
                        morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                    }
                }
                //
                // if it's a shape-systematic gamma
                else if(fSystNames[i_syst].find("shape_")!=std::string::npos && fSystNames[i_syst].find(Form("%s_bin_%d",fName.c_str(),i_bin-1))!=std::string::npos){
                    for(const auto& syh : fSampleHists[i]->fSyst){
                        std::shared_ptr<Systematic> syst = syh->fSystematic;
                        if(!syst) continue;
                        if(syst->fType==Systematic::SHAPE){
                            const std::string gammaNameShapeSyst = Form("shape_%s_%s_bin_%d",syst->fName.c_str(),fName.c_str(),i_bin-1);
                            if(gammaNameShapeSyst==fSystNames[i_syst]){
                                diffUp   += yieldNominal_postFit*systErrUp;
                                diffDown += yieldNominal_postFit*systErrDown;
                                if (isMorph){
                                    morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                                    morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                                }
                            }
                        }
                    }
                }
                //
                // if it's a norm factor
                else if(fSampleHists[i]->HasNorm(fSystNames[i_syst])){

                    // if this norm factor is a morphing one
                    if(fSystNames[i_syst].find("morph_")!=string::npos || fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fExpression.first!=""){

                        const std::vector<double> scales = Common::CalculateExpression(nullptr, fSystNames[i_syst], true, fSampleHists[i].get(), fitRes);
                        if (scales.size() != 3) {
                            WriteErrorStatus("Region::BuildPostFitErrorHist", "Scales size is not 3");
                            exit(EXIT_FAILURE);
                        }
                        morph_syst_up.at(i_bin-1)   += yieldNominal*scales.at(1);
                        morph_syst_down.at(i_bin-1) += yieldNominal*scales.at(2);
                    }
                    else{
                        diffUp   += yieldNominal_postFit*systErrUp/systValue;
                        diffDown += yieldNominal_postFit*systErrDown/systValue;
                        if (isMorph){
                            morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                            morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                        }
                    }
                }
                //
                // ShapeFactor have to get NP per bin
                else if(posTmp != std::string::npos){
                    // get the shape factor name without bin index
                    const std::string systNameSF = systName.substr(0, posTmp);
                    // get the shape factor bin as integer
                    const int iBinSF = std::atoi(systName.substr(posTmp + 5).c_str()) + 1;
                    // FIXME could still be a problem with pruning?
                    if(iBinSF == i_bin && fSampleHists[i]->HasShapeFactor(systNameSF)){
                        diffUp   += yieldNominal_postFit*systErrUp/systValue;
                        diffDown += yieldNominal_postFit*systErrDown/systValue;
                        if (isMorph){
                            morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                            morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                        }
                    }
                }
                //
                // Systematics treatment
                //
                else if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    std::ofstream dummy;
                    const double multNom  = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal );
                    const double multUp   = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, TRExFitter::NPMAP[systName], true);
                    const double multDown = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, TRExFitter::NPMAP[systName], false);
                    
                    if (isMorph){
                        if (std::abs(multNom) > 1e-9) {
                            morph_up_postfit.at(i_bin-1) += (multUp/multNom)*yieldNominal_postFit;
                            morph_down_postfit.at(i_bin-1) += (multDown/multNom)*yieldNominal_postFit;
                        }
                    } else {
                        if (std::abs(multNom) > 1e-9 ) {
                            diffUp   += (multUp/multNom   - 1.)*yieldNominal_postFit;
                            diffDown += (multDown/multNom - 1.)*yieldNominal_postFit;
                        }
                    }

                }

                WriteVerboseStatus("Region::BuildPostFitErrorHist", "        Bin " + std::to_string(i_bin) + ":   " + "\t +" + std::to_string(100*diffUp/yieldNominal)
                    + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");

                //
                // Add the proper bin content to the variation hists (coming from post-fit total histogram)
                //
                if (!isMorph){
                    sh->fHistUp_postFit  ->AddBinContent( i_bin, diffUp   );
                    sh->fHistDown_postFit->AddBinContent( i_bin, diffDown );
                }

            } // loop over bins
        } // loop over samples

        if (isMorph) {
            i_morph_sample++;
            // now apply the corrections from morph
            // Use it only ONCE for all combines morph templates
            std::shared_ptr<SystematicHist> sh = fSampleHists[morph_index]->GetSystematic(systName);

            // hack: add a systematic hist if not there
            if(sh==nullptr){
                fSampleHists[morph_index]->AddHistoSyst(systName,systName,fSampleHists[morph_index]->fHist.get(),fSampleHists[morph_index]->fHist.get());
                sh = fSampleHists[morph_index]->GetSystematic(systName);
            }

            //
            // initialize the up and down variation histograms
            // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
            //
            sh->fHistUp_postFit.reset(static_cast<TH1*>(fSampleHists[morph_index]->fHist_postFit->Clone(Form("%s_%s_Up_postFit",  fSampleHists[morph_index]->fHist->GetName(),systName.c_str()))));
            sh->fHistDown_postFit.reset(static_cast<TH1*>(fSampleHists[morph_index]->fHist_postFit->Clone(Form("%s_%s_Down_postFit",fSampleHists[morph_index]->fHist->GetName(),systName.c_str()))));
            sh->fHistUp_postFit->SetDirectory(nullptr);
            sh->fHistDown_postFit->SetDirectory(nullptr);

            // loop over bins
            for (int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                sh->fHistUp_postFit    ->AddBinContent( i_bin, (morph_up_postfit.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                sh->fHistDown_postFit  ->AddBinContent( i_bin, (morph_down_postfit.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
            }

            // uncertainty from the morph uncertainty
            if (fSystNames[i_syst].find("morph_")!=string::npos){
                for (int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                    sh->fHistUp_postFit    ->AddBinContent( i_bin, (morph_syst_up.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                    sh->fHistDown_postFit  ->AddBinContent( i_bin, (morph_syst_down.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                }
            }
        }
    }// loop over systs

    // at this point all the sample-by-sample post-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    //

    // - loop on systematics
    fTotUp_postFit.resize(fSystNames.size());
    fTotDown_postFit.resize(fSystNames.size());
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        const std::string systName = fSystNames[i_syst];
        //
        // Initialize the tot variation hists
        //
        fTotUp_postFit[i_syst].reset(static_cast<TH1*>(fTot_postFit->Clone(Form("h_tot_%s_Up_postFit",  systName.c_str()))));
        fTotUp_postFit[i_syst] -> Scale(0); //initialising the content to 0
        fTotDown_postFit[i_syst].reset(static_cast<TH1*>(fTot_postFit->Clone(Form("h_tot_%s_Down_postFit",systName.c_str()))));
        fTotDown_postFit[i_syst] -> Scale(0); //initialising the content to 0
        fTotUp_postFit[i_syst]->SetDirectory(nullptr);
        fTotDown_postFit[i_syst]->SetDirectory(nullptr);

        // - loop on bins
        for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
            double diffUp(0.);
            double diffDown(0.);
            // - loop on samples
            for(const auto& isample : fSampleHists) {
                // skip data
                if(isample->fSample->fType==Sample::DATA) continue;
                if(isample->fSample->fType==Sample::GHOST) continue;
                if(isample->fSample->fType==Sample::EFT) continue;
                if(isample->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
                // skip signal if Bkg only
                if(fFitType==TRExFit::BONLY && isample->fSample->fType==Sample::SIGNAL) continue;
                // get SystematicHist
                std::shared_ptr<SystematicHist> sh = isample->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                if(!sh) continue;
                diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin);
                diffDown += sh->fHistDown_postFit->GetBinContent(i_bin);
            }
            // add the proper bin content to the variation hists
            fTotUp_postFit[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown_postFit[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }

    // at this point all the total expectation post-fit variation histograms should be filled
    //
    // Build the vectors of variations (necessary to call BuildTotError)
    //
    std::vector< std::shared_ptr<TH1> > h_up;
    std::vector< std::shared_ptr<TH1> > h_down;
    std::vector<std::string> systNuisPars;
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        h_up.  emplace_back(fTotUp_postFit[i_syst]  );
        h_down.emplace_back(fTotDown_postFit[i_syst]);
        systNuisPars.push_back(TRExFitter::NPMAP[fSystNames[i_syst]]);
    }
    fErr_postFit = BuildTotError( fTot_postFit.get(), h_up, h_down, systNuisPars, fitRes->fCorrMatrix.get());
    fErr_postFit->SetName("g_totErr_postFit");
    // at this point fTot and fErr _postFit should be ready

    //
    // Goodness of post-fit
    if(fGetChi2){
        if (!fHasData){
            WriteWarningStatus("Region::BuildPostFitErrorHist", "Data histogram is nullptr, cannot calculate Chi2 agreement.");
            WriteWarningStatus("Region::BuildPostFitErrorHist", "Maybe you do not have data sample defined?");
            return;
        }
        // remove blinded bins
        std::unique_ptr<TH1> h_data(static_cast<TH1*>(fData->fHist->Clone()));
        if(fGetChi2==1) fSystNames.clear();
        std::pair<double,int> res = GetChi2Test(true, h_data.get(), fTot_postFit.get(), h_up, fSystNames, fitRes->fCorrMatrix.get() );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- POST-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }

}

//__________________________________________________________________________________
//
std::shared_ptr<TRExPlot> Region::DrawPostFit(FitResults* fitRes,
                                              ofstream& pullTex,
                                              const std::vector<std::string> &morph_names,
                                              const std::vector<int>& canvasSize,
                                              string opt) {

    if(TRExFitter::PREFITONPOSTFIT){
        fPlotPostFit->h_tot_bkg_prefit = static_cast<TH1*>(fPlotPreFit->GetTotBkg()->Clone("h_tot_bkg_prefit"));
    }

    std::shared_ptr<TRExPlot> p(nullptr);
    if (canvasSize.size() == 0){
        p = fPlotPostFit;
        p->fShowYields = TRExFitter::SHOWYIELDS;
    } else {
        p = std::make_shared<TRExPlot>(("c_"+fName).c_str(), canvasSize.at(0), canvasSize.at(1),TRExFitter::NORATIO);
    }

    p->SetXaxisRange(fXaxisRange);
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmaxPostFit;
    p->fRatioYmin = fRatioYminPostFit;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TRExFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        if (fLabel != "none") p->AddLabel(fLabel);
        p->AddLabel("#font[52]{Post-fit}");
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        if (fLabel != "none") p->AddLabel(fLabel);
        if(TRExFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Post-Fit");
    }
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLegendNColumns = fLegendNColumns;

    if (!fBinLabels.empty()) {
        p->ResizeBinLabel(fBinLabels.size() + 1);
        for(std::size_t i_bin=0; i_bin<fBinLabels.size(); i_bin++) {
            p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
        }
    }

    YamlConverter::PlotContainer container;
    container.region = fName;
    container.xAxis = fVariableTitle;
    container.yAxis = fYTitle;

    //
    // 0) Create a new hist for each sample
    //
    std::vector<std::shared_ptr<TH1> > hSmpNew(fSampleHists.size());
    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        hSmpNew[i].reset(static_cast<TH1*>(fSampleHists[i]->fHist->Clone()));
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if((!fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas) || fUseGammaPulls){
            for(int i_bin=0;i_bin<hSmpNew[i]->GetNbinsX()+2;i_bin++) hSmpNew[i]->SetBinError(i_bin,0.);
        }
    }

    //
    // 1) Propagates the post-fit NP values to the central value (pulls)
    //
    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
        //
        if (fSampleHists[i]->fSample->fBuildPullTable>0){
            WriteDebugStatus("Region::DrawPostFit", "Propagating post-fit to Sample " + fSampleHists[i]->fSample->fTitle);
            TString sampleTex= fSampleHists[i]->fSample->fTexTitle;
            pullTex << "\\hline\n" << endl;
            pullTex << "{\\color{blue}{$\\rightarrow \\,$ "<< sampleTex << "}} & \\\\\n"<< endl;
        }
        std::unique_ptr<TH1> hNew(static_cast<TH1*>(hSmpNew[i]->Clone()));
        for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
            const double binContent0 = hSmpNew[i]->GetBinContent(i_bin);
            const double mult_factor = GetMultFactors(fitRes, pullTex, i, i_bin, binContent0);

            //
            // Final computation
            //
            double binContentNew = binContent0*mult_factor;

            //
            // stat gammas
            if(fUseGammaPulls && (fSampleHists[i]->fSample->fUseMCStat || fSampleHists[i]->fSample->fSeparateGammas)){
                // find the gamma for this bin of this distribution in the fit results
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(fSampleHists[i]->fSample->fSeparateGammas) {
                    gammaName = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                }
                WriteDebugStatus("Region::DrawPostFit", "Looking for gamma " + gammaName);
                const double gammaValue = fitRes->GetNuisParValue(gammaName);
                WriteDebugStatus("Region::DrawPostFit", "  -->  pull = " + std::to_string(gammaValue));
                // linear effect
                if(gammaValue>0) binContentNew *= gammaValue;
            }
            // gammas from SHAPE systematics
            for(const auto& syh : fSampleHists[i]->fSyst){
                std::shared_ptr<Systematic> syst = syh->fSystematic;
                if(!syst) continue;
                if(syst->fType==Systematic::SHAPE){
                    const std::string gammaName = Form("shape_%s_%s_bin_%d",syst->fName.c_str(),fName.c_str(),i_bin-1);
                    WriteDebugStatus("Region::DrawPostFit", "Looking for gamma " + gammaName);
                    const double gammaValue = fitRes->GetNuisParValue(gammaName);
                    WriteDebugStatus("Region::DrawPostFit", "  -->  pull = " + std::to_string(gammaValue));
                    // linear effect
                    if(gammaValue>0) binContentNew *= gammaValue;
                }
            }

            //
            // Setting to new values
            //
            hNew->SetBinContent(i_bin,binContentNew);
        }
        hSmpNew[i].reset(static_cast<TH1*>(hNew->Clone()));
        fSampleHists[i]->fHist_postFit = hSmpNew[i];
        if (fSampleHists[i]->fHist_postFit) fSampleHists[i]->fHist_postFit->SetDirectory(nullptr);
    }
    //
    // 2) Scale all samples by norm factors
    //    Done after the propagation of the NP (avoids nans due to "0" value of some NormFactors)
    //    Seems consistent with application in Roostats
    //
    double nfValue;
    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
        for(const auto& inorm : fSampleHists[i]->fSample->fNormFactors) {
            const std::string nfName = inorm->fName;
            if (!inorm->fRegions.empty() && inorm->fRegions.at(0) != "none" && Common::FindInStringVector(inorm->fRegions, fName) < 0) continue;
            if (Common::FindInStringVector(inorm->fExclude, fName) >= 0) continue;
            if(inorm->fConst) nfValue = inorm->fNominal;
            else              nfValue = fitRes->GetNuisParValue(TRExFitter::NPMAP[nfName]);
            //
            // if this norm factor is a morphing one
            if(inorm->fName.find("morph_")!=string::npos || inorm->fExpression.first!=""){
                const std::vector<double> scales = Common::CalculateExpression(inorm.get(), nfName, true, nullptr, fitRes);
                if (scales.empty()) {
                    WriteErrorStatus("Region::DrawPostFit", "Scales are empty");
                    exit(EXIT_FAILURE);
                }
                double scale = scales.at(0);
                hSmpNew[i]->Scale(scale);
            }
            else{
                if (std::find(inorm->fRegions.begin(), inorm->fRegions.end(), fName) != inorm->fRegions.end()) {
                    hSmpNew[i]->Scale(nfValue);
                }
            }
        }
    }


    //
    // 2)b) Scale all samples by shape factors factors
    //

    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
        for(const auto& ishape : fSampleHists[i]->fSample->fShapeFactors) {
            const std::string sfName = ishape->fName;
            if (!ishape->fRegions.empty() && ishape->fRegions.at(0) != "none" && Common::FindInStringVector(ishape->fRegions, fName) < 0) continue;
            if (Common::FindInStringVector(ishape->fExclude, fName) >= 0) continue;
            if(sfName.find("saturated_model_sf")!=std::string::npos) continue;
            // loop over bins
            // there should be a NP per bin in the fit file
            // already checked by GetNuisParValue()
            // the shape factor naming used i_bin - 1 for the first bin
            for(int i_bin = 1; i_bin <= hSmpNew[i]->GetNbinsX(); i_bin++){
                const int iBinSF = i_bin - 1;
                const std::string sfNameBin = sfName + "_bin_" + std::to_string(iBinSF);
                double sfValue;
                if(ishape->fConst) sfValue = ishape->fNominal;
                else               sfValue = fitRes->GetNuisParValue(sfNameBin);
                // scale bin content by shape factor
                const double binContentSFNew = hSmpNew[i]->GetBinContent(i_bin) * sfValue;
                // Setting to new value
                hSmpNew[i]->SetBinContent(i_bin, binContentSFNew);
            }
        }
    }

    // Scale samples acording to fScaleSamplesToData:
    if(fScaleSamplesToData.size()>0){
        std::unique_ptr<TH1> hTot(nullptr);
        for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
            if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
            if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
            if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
            if(hTot==nullptr) hTot.reset(static_cast<TH1*>(hSmpNew[i]->Clone("hTotPostFit")));
            else              hTot->Add(hSmpNew[i].get());
            hTot->SetDirectory(nullptr);
        }
        const double totPred = hTot->Integral();
        const double totData = fData->fHist->Integral();
        double totToScale = 0;
        std::vector<std::size_t> shIdxToScale;
        for(std::size_t i = 0; i < fSampleHists.size(); ++i){
            const SampleHist *sh = fSampleHists[i].get();
            if(sh->fHist==nullptr) continue;
            if(sh->fSample->fType==Sample::GHOST){
                WriteWarningStatus("Region::DrawPostFit","Requested to scale to data a GHOST sample, " + sh->fSample->fName + ". Skipping this sample.");
                continue;
            }
            if(sh->fSample->fType==Sample::EFT){
                WriteWarningStatus("Region::DrawPostFit","Requested to scale to data a EFT sample, " + sh->fSample->fName + ". Skipping this sample.");
                continue;
            }
            if(Common::FindInStringVector(fScaleSamplesToData,sh->fSample->fName)>=0){
                shIdxToScale.emplace_back(i);
                totToScale += hSmpNew[i]->Integral();
            }
        }
        if(totToScale>0 && shIdxToScale.size()>0){
            const double scale = (totData-(totPred-totToScale))/totToScale;
            for(const auto& idx : shIdxToScale){
                WriteInfoStatus("Region::DrawPostFit","Scaling sample " + fSampleHists[idx]->fSample->fName + " by " + std::to_string(scale) + " in region " + fName);
                hSmpNew[idx]->Scale(scale);
            }
        }
    }

    //
    // 3) Add the new Sig and Bkg to plot
    //
    std::vector<std::shared_ptr<TH1> > hBkgNew;
    std::vector<std::shared_ptr<TH1> > hSigNew;
    {
        for(std::size_t i = 0; i < fSampleHists.size(); ++i){
            if(fSampleHists[i]->fSample->fType==Sample::SIGNAL){
                std::vector<double> tmp;
                for (int ibin = 1; ibin <= hSmpNew[i]->GetNbinsX(); ++ibin) {
                    tmp.emplace_back(hSmpNew[i]->GetBinContent(ibin));
                }
                container.signalYields.emplace_back(std::move(tmp));
                hSigNew.emplace_back(hSmpNew[i]);
            }
            if(fSampleHists[i]->fSample->fType==Sample::BACKGROUND){
                std::vector<double> tmp;
                for (int ibin = 1; ibin <= hSmpNew[i]->GetNbinsX(); ++ibin) {
                    tmp.emplace_back(hSmpNew[i]->GetBinContent(ibin));
                }
                container.backgroundYields.emplace_back(std::move(tmp));
                hBkgNew.emplace_back(hSmpNew[i]);
            }
        }
        if(fHasData && opt.find("blind")==string::npos) {
            for (int ibin = 1; ibin <= fData->fHist->GetNbinsX(); ++ibin) {
                container.data.emplace_back(fData->fHist->GetBinContent(ibin));
            }
            p->SetData(fData->fHist.get(),fData->fSample->fTitle);
        }
        for(std::size_t i = 0; i < fSig.size(); ++i) {
            container.samples.emplace_back(fSig[i]->fSample->fTitle);
            std::string title = fSig[i]->fSample->fTitle;
            if(fSig[i]->fSample->fGroup != "") title = fSig[i]->fSample->fGroup;
            if(TRExFitter::SHOWSTACKSIG)    p->AddSignal(    hSigNew[i].get(),title);
            if(TRExFitter::SHOWNORMSIG){
                if( (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
                 || !TRExFitter::OPTION["NormSigSRonly"] )
                    p->AddNormSignal(hSigNew[i].get(),title);
            }
            else{
                if(TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(hSigNew[i].get(),title);
            }
            if(TRExFitter::SHOWOVERLAYSIG){
                Common::ScaleNominal(fSig[i].get(),hSigNew[i].get());
                p->AddOverSignal(hSigNew[i].get(),title);
            }
        }
        for(std::size_t i = 0; i < fBkg.size(); ++i) {
            container.samples.emplace_back(fBkg[i]->fSample->fTitle);
            std::string title = fBkg[i]->fSample->fTitle;
            if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
            p->AddBackground(hBkgNew[i].get(),title);
        }
    }

    //
    // 4) Build post-fit error band
    //    Build total histogram (fTot_postFit)
    //
    int j = 0;
    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::EFT) continue;
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
        if(j==0) fTot_postFit.reset(static_cast<TH1*>(hSmpNew[i]->Clone("h_tot_postFit")));
        else fTot_postFit->Add(hSmpNew[i].get());
        j++;
    }
    fTot_postFit->SetDirectory(nullptr);

    //
    // MCstat set to 0 if disabled
    //
    if(!fUseStatErr || fUseGammaPulls){
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
            fTot_postFit->SetBinError(i_bin,0);
        }
    }

    p->SetTotBkg(fTot_postFit.get());
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);

    //
    // blinding bins
    //
    if(fBlindingThreshold>=0){

        if (fKeepPrefitBlindedBins) {
            container.blindedBins = fBlindedBins;
            p->SetBinBlinding(fBlindedBins);
        } else {
            auto CombineHistos = [](const std::vector<std::shared_ptr<TH1> >& vec) {
                std::unique_ptr<TH1D> result(nullptr);
                for (const auto& ihist : vec) {
                    if (!ihist) continue;
                    if (result) {
                        result->Add(ihist.get());
                    } else {
                        result.reset(static_cast<TH1D*>(ihist->Clone()));
                        result->SetDirectory(nullptr);
                    }
                }
                return result;
            };

            std::unique_ptr<TH1D> signal = CombineHistos(hSigNew);
            std::unique_ptr<TH1D> bkg = CombineHistos(hBkgNew);

            fBlindedBinsPostFit = Common::ComputeBlindedBins(signal.get(),
                                                             bkg.get(),
                                                             fBlindingType,
                                                             fBlindingThreshold);
            container.blindedBins = fBlindedBinsPostFit;
            p->SetBinBlinding(fBlindedBinsPostFit);
        }
    }
    p->SetDropBins(fDropBins);
    p->BlindData();

    //
    // Build error band
    //
    BuildPostFitErrorHist(fitRes, morph_names);

    //
    // Print chi2 info
    //
    if(fGetChi2 && TRExFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);

    //
    // 5) Finishes configuration of TRExPlot objects
    //
    container.errors = fErr_postFit.get();
    p->SetTotBkgAsym(fErr_postFit.get());
    p->fATLASlabel = fATLASlabel;
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType==TRExPlot::RATIOTYPE::DATAOVERMC){
        p->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
    }
    p->fLabelX = fLabelX;
    p->fLabelY = fLabelY;
    p->fLegendX1 = fLegendX1;
    p->fLegendX2 = fLegendX2;
    p->fLegendY = fLegendY;
    if(fLogScale) opt += " log";

    p->Draw(opt);

    // Yaml
    YamlConverter converter{};
    converter.WritePlot(container, fFolder, true);
    if (fHEPDataFormat) {
        converter.WritePlotHEPData(container, fFolder, true);
    }

    //
    // Print bin content and errors
    //
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
    WriteDebugStatus("Region::DrawPostFit", "Final bin contents");
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
    for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++) {
        WriteDebugStatus("Region::DrawPostFit", std::to_string(i_bin) + ":\t" + std::to_string(fTot_postFit->GetBinContent(i_bin)) + " +" +
        std::to_string( fErr_postFit->GetErrorYhigh(i_bin-1)) + " -" + std::to_string(fErr_postFit->GetErrorYlow(i_bin-1)));
    }

    //
    // Save in a root file...
    //
    gSystem->mkdir((fFitName+"/Histograms").c_str());
    WriteInfoStatus("Region::DrawPostFit", "Writing file " + fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root");
    std::unique_ptr<TFile> f(TFile::Open((fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root").c_str(),"RECREATE"));
    f->cd();
    fErr_postFit->Write("",TObject::kOverwrite);
    fTot_postFit->Write("",TObject::kOverwrite);
    if( fPlotPostFit->h_tot_bkg_prefit) {
      fPlotPostFit->h_tot_bkg_prefit->Write("",TObject::kOverwrite);
    }
    for(std::size_t i_syst=0; i_syst<fSystNames.size(); ++i_syst){
        if(fTotUp_postFit[i_syst])   fTotUp_postFit[i_syst]  ->Write("",TObject::kOverwrite);
        if(fTotDown_postFit[i_syst]) fTotDown_postFit[i_syst]->Write("",TObject::kOverwrite);
    }
    for(const auto& isample : fSampleHists) {
        if(isample->fSample->fType == Sample::DATA){
            isample->fHist->Write(Form("h_%s",isample->fName.c_str()),TObject::kOverwrite);
            continue;
        }
        if(isample->fHist_postFit){
            isample->fHist_postFit->Write(Form("h_%s_postFit",isample->fName.c_str()),TObject::kOverwrite);
            for(std::size_t i_syst=0;i_syst<isample->fSyst.size();i_syst++) {
                if(isample->fSyst[i_syst]) {
                    if(isample->fSyst[i_syst]->fHistUp_postFit) {
                        isample->fSyst[i_syst]->fHistUp_postFit  ->Write(
                          Form("h_%s_%s_Up_postFit",isample->fName.c_str(),isample->fSyst[i_syst]->fName.c_str()), TObject::kOverwrite);
                    }
                    if(isample->fSyst[i_syst]->fHistDown_postFit) {
                        isample->fSyst[i_syst]->fHistDown_postFit->Write(
                          Form("h_%s_%s_Down_postFit",isample->fName.c_str(),isample->fSyst[i_syst]->fName.c_str()),TObject::kOverwrite);
                    }
                }
            }
        }
    }
    f->Close();
    //
    return p;
}

//__________________________________________________________________________________
//
void Region::AddSelection(const std::string& selection){
    if(selection=="") return;
    if(fSelection=="1" || fSelection=="") fSelection = selection;
    else fSelection += " && "+selection;
}

//__________________________________________________________________________________
//
void Region::AddMCweight(const std::string& weight){
    if(weight=="") return;
    if(fMCweight=="1" || fMCweight=="") fMCweight = weight;
    else fMCweight += " * "+weight;
}

//__________________________________________________________________________________
//
void Region::SetVariable(const std::string& variable,int nbin,double xmin,double xmax,string corrVar1,string corrVar2){
    fVariable = variable;
    fCorrVar1 = corrVar1;
    fCorrVar2 = corrVar2;
    fNbins = nbin;
    fXmin = xmin;
    fXmax = xmax;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeVariable(const std::string& variable, const std::string& sample){
    fAlternativeVariables[sample] = variable;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeSelection(const std::string& selection, const std::string& sample){
    fAlternativeSelections[sample] = selection;
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeVariable(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(const auto& tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
    }
    if (Common::FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    return true;
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeSelection(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(const auto& tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
    }
    if (Common::FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    return true;
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeVariable(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(const auto& tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    const int idx = Common::FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    return tmpVec2[idx];
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeSelection(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(const auto& tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    const int idx = Common::FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    return tmpVec2[idx];
}

//__________________________________________________________________________________
//
void Region::SetVariableTitle(const std::string& name){
    fVariableTitle = name;
}

//__________________________________________________________________________________
//
void Region::SetLabel(const std::string& label, std::string shortLabel){
    fLabel = label;
    if(shortLabel=="") fShortLabel = label;
    else fShortLabel = shortLabel;
}

//__________________________________________________________________________________
//
void Region::Print() const{
    WriteInfoStatus("Region::Print", "    Region: " + fName);
    for(const auto& isample : fSampleHists) {
        isample->Print();
    }
}

//__________________________________________________________________________________
//
void Region::PrintSystTable(FitResults *fitRes, string opt) const{
    bool isPostFit  = false; if(opt.find("post")!=string::npos)     isPostFit  = true;
    bool doClean    = false; if(opt.find("clean")!=string::npos)    doClean    = true;
    bool doCategory = false; if(opt.find("category")!=string::npos) doCategory = true;
    bool standalone = false; if(opt.find("standalone")!=string::npos)standalone= true;
    bool landscape  = false; if(opt.find("landscape")!=string::npos) landscape = true;
    bool footnotesize=false; if(opt.find("footnotesize")!=string::npos)footnotesize=true;;
    //
    ofstream out;
    ofstream texout;
    ofstream out_cat;
    ofstream texout_cat;
    gSystem->mkdir(fFitName.c_str());
    gSystem->mkdir((fFitName+"/Tables").c_str());
    if(isPostFit){
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.tex").c_str());
        }
    }
    else{
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.tex").c_str());
        }
    }

    out << " | ";
    if (standalone) {
        texout << "\\documentclass[10pt]{article}" << endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
        texout << "\\begin{document}" << endl;
    }
    if (landscape) texout << "\\begin{landscape}" << endl;
    texout << "\\begin{table}[htbp]" << endl;
    texout << "\\begin{center}" << endl;
    if (footnotesize) texout << "\\footnotesize" << endl;
    texout << "\\begin{tabular}{|c" ;

    if(doCategory){
        out_cat << " | ";
        if (standalone) {
            texout_cat << "\\documentclass[10pt]{article}" << endl;
            texout_cat << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
            texout_cat << "\\begin{document}" << endl;
        }
        texout_cat << "\\begin{table}[htbp]" << endl;
        texout_cat << "\\begin{center}" << endl;
        texout_cat << "\\begin{tabular}{|c" ;
    }


    double Ncol = 2.;
    for(std::size_t i_smp=0;i_smp<fSampleHists.size();i_smp++){
        const SampleHist* sh = fSampleHists[i_smp].get();
        const Sample* s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        if(s->fType==Sample::EFT) continue;
        texout << "|c";
        if(doCategory) texout_cat << "|c";
        Ncol+=1;
    }
    texout << "|}" << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        texout_cat << "|}" << endl;
        texout_cat << "\\hline " << endl;
    }
    //
    double i_col = 1.;
    for(std::size_t i_smp=0;i_smp<fSampleHists.size();i_smp++){
        const SampleHist* sh = fSampleHists[i_smp].get();
        const Sample* s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        if(s->fType==Sample::EFT) continue;
        std::string title = s->fTitle;
        if(s->fTexTitle!="") title = s->fTexTitle;
        out << "      | " << s->fTitle;
        texout << "      & " << title;
        if(doCategory){
            out_cat << "      | " << s->fTitle;
            texout_cat << "      & " << title;
        }
        i_col+=1;
    }
    out << " |" << endl;
    texout << " \\\\ " << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        out_cat << " |" << endl;
        texout_cat << " \\\\ " << endl;
        texout_cat << "\\hline " << endl;
    }

    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        if(TRExFitter::SYSTMAP[fSystNames[i_syst]]!="") out << " | " << TRExFitter::SYSTMAP[fSystNames[i_syst]];
        else                                           out << " | " << fSystNames[i_syst];
        if(TRExFitter::SYSTTEX[fSystNames[i_syst]]!="")      texout << "  " << TRExFitter::SYSTTEX[fSystNames[i_syst]];
        else if(TRExFitter::SYSTMAP[fSystNames[i_syst]]!=""){
            std::string fixedTitle = TRExFitter::SYSTMAP[fSystNames[i_syst]];
            fixedTitle = Common::ReplaceString(fixedTitle,"#geq","$\\geq$");
            texout << "  " << fixedTitle;
        }
        else {
            texout << " " << fSystNames[i_syst];
        }
        for(std::size_t i_smp=0;i_smp<fSampleHists.size();i_smp++){
            const SampleHist* sh = fSampleHists[i_smp].get();
            const Sample* s = sh->fSample;
            if(s->fType==Sample::DATA) continue;
            if(s->fType==Sample::GHOST) continue;
            if(s->fType==Sample::EFT) continue;
            std::shared_ptr<SystematicHist> syh = sh->GetSystematic(fSystNames[i_syst]);
            if(syh==nullptr){
                out << " |    nan   ";
                texout << " &    nan   ";
            }
            else{
              double normUp(0);
              double normDown(0);
              if(isPostFit){
                normUp   = (syh->fHistUp_postFit->Integral()   - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
                normDown = (syh->fHistDown_postFit->Integral() - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
              }
              else{
                normUp = syh->fNormUp;
                normDown = syh->fNormDown;
              }

              out << " | " << normUp;
              texout << setprecision(3) << " & " << normUp;
              out << " / " << normDown;
              texout << setprecision(3) << " / " << normDown;
              //                 pt[i_smp]->AddText(Form("%.2f / %.2f",syh->fNormUp,syh->fNormDown) );
            }
        }

        out << " |" << endl;
        texout << " \\\\ " << endl;
    }


    if(doCategory){
        //--- Systematic tables per category:
        //--- Get systematic names per category and sample
        std::set<std::string> category_names;
        std::map<std::string, std::map<std::string, std::vector<std::string> > > category_syst_names;
        for(std::size_t i_smp=0;i_smp<fSampleHists.size();i_smp++){
            const SampleHist* sh = fSampleHists[i_smp].get();
            const Sample* s = sh->fSample;
            for(int i_samplesyst=0; i_samplesyst<(int)s->fSystematics.size();i_samplesyst++){
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;
                if(s->fType==Sample::EFT) continue;
                const std::string category = s->fSystematics[i_samplesyst]->fCategory;
                if (category!=""){
                    category_names.insert(category);
                    // Check for systematics existing not in all regions...
                    if (std::find(fSystNames.begin(), fSystNames.end(), s->fSystematics.at(i_samplesyst)->fName)!=fSystNames.end()){
                        const std::vector<std::string> sample_syste = category_syst_names[s->fName][category];
                        if (std::find(sample_syste.begin(), sample_syste.end(), s->fSystematics.at(i_samplesyst)->fName) == sample_syste.end()){
                            category_syst_names[s->fName][category].push_back(s->fSystematics.at(i_samplesyst)->fName);
                        }
                    }
                }
            }
        }
        //--- Loop over categories and samples, to get systematic effects using BuildTotError function
        for (auto category : category_names){
            out_cat << " | " << category;
            texout_cat << " " << category;
            for(std::size_t i_smp=0;i_smp<fSampleHists.size();i_smp++){
                const SampleHist* sh = fSampleHists[i_smp].get();
                const Sample* s = sh->fSample;
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;
                if(s->fType==Sample::EFT) continue;

                std::vector<std::shared_ptr<TH1> > category_histo_up;
                std::vector<std::shared_ptr<TH1> > category_histo_down;
                const std::vector<std::string> sample_syste = category_syst_names[s->fName][category];
                for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
                    if(!sh->HasSyst(fSystNames[i_syst])) continue;

                    if (std::find(sample_syste.begin(), sample_syste.end(), fSystNames.at(i_syst)) == sample_syste.end()) continue;

                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(fSystNames[i_syst]);

                    if(isPostFit){
                        TH1 *h_up = syh->fHistUp_postFit.get();
                        h_up->SetDirectory(nullptr);
                        TH1 *h_down = syh->fHistDown_postFit.get();
                        h_down->SetDirectory(nullptr);
                        category_histo_up.emplace_back(static_cast<TH1*>(h_up->Clone()));
                        category_histo_down.emplace_back(static_cast<TH1*>(h_down->Clone()));
                    }
                    else{
                        TH1 *h_up = syh->fHistUp.get();
                        h_up->SetDirectory(0);
                        TH1 *h_down = syh->fHistDown.get();
                        h_down->SetDirectory(0);
                        category_histo_up.emplace_back(static_cast<TH1*>(h_up->Clone()));
                        category_histo_down.emplace_back(static_cast<TH1*>(h_down->Clone()));
                    }
                }

                std::unique_ptr<TGraphAsymmErrors> g_err(nullptr);
                double err(0.);
                if (isPostFit){
                    g_err = BuildTotError(sh->fHist_postFit.get(), category_histo_up, category_histo_down, category_syst_names[s->fName][category], fitRes->fCorrMatrix.get());
                    if (category_histo_up.size()>0 && sh->fHist_postFit->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist_postFit->GetNbinsX()+1; ibin++){
                            const double tmp = g_err->GetErrorYhigh(ibin-1) * g_err->GetErrorYhigh(ibin-1) - sh->fHist_postFit->GetBinError(ibin) * sh->fHist_postFit->GetBinError(ibin);
                            if (tmp >=0.){ // dummy check
                                //Need to substract the statistical unc.
                                err += (std::sqrt(tmp))/sh->fHist_postFit->Integral();
                            }
                        }
                    }
                }
                else{
                    g_err = BuildTotError(sh->fHist.get(), category_histo_up, category_histo_down, category_syst_names[s->fName][category]);
                    if (category_histo_up.size()>0 && sh->fHist->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist->GetNbinsX()+1; ibin++){
                            const double tmp = g_err->GetErrorYhigh(ibin-1) * g_err->GetErrorYhigh(ibin-1) - sh->fHist->GetBinError(ibin) * sh->fHist->GetBinError(ibin);
                            if (tmp>=0.){
                                //Need to substract the statistical unc.
                                err += (std::sqrt(tmp))/sh->fHist->Integral();
                            }
                        }
                    }
                }
                out_cat << setprecision(3) << " | " << err;
                texout_cat << setprecision(3) << " & " << err;
            }
            out_cat << " |" << endl;
            texout_cat << " \\\\ " << endl;
        }
    }

    texout << "\\hline " << endl;
    texout << "\\end{tabular} " << endl;
    texout << "\\caption{Relative effect of each systematic on the yields.} " << endl;
    texout << "\\end{center} " << endl;
    texout << "\\end{table} " << endl;

    if(doCategory){
        texout_cat << "\\hline " << endl;
        texout_cat << "\\end{tabular} " << endl;
        texout_cat << "\\caption{Realtive effect of each group of systematics on the yields.} " << endl;
        texout_cat << "\\end{center} " << endl;
        texout_cat << "\\end{table} " << endl;
    }
    if (landscape) texout << "\\end{landscape}" << endl;
    if (standalone) {
        texout << "\\end{document}" << endl;
    }

    if(doClean){
        std::string shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
        shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst_clean.tex";
        gSystem->Exec(shellcommand.c_str());
        if(doCategory){
            shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
            shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += "_clean.tex";
            gSystem->Exec(shellcommand.c_str());
        }
    }
}

// --------------- Functions --------------- //

double GetDeltaN(double alpha, double Iz, double Ip, double Imi, int intCode){
    // protection against negative values
    if(Ip<=0)  Ip  = 0.00000001*Iz;
    if(Imi<=0) Imi = 0.00000001*Iz;

    double deltaN = 0.;
    if(alpha>0)      deltaN = Ip;
    else if(alpha<0) deltaN = Imi;
    else             return 1.;

    if(intCode==4){

        //////////////////////////////////////////////////////////////
        // NORMALISATION
        // =============
        // Exponential-polynomial
        // Equation solved with Mathematica
        //////////////////////////////////////////////////////////////

        if(std::fabs(alpha)>1){
            deltaN /= Iz; // divide h_tmp by the nominal
            deltaN = std::pow( deltaN, std::fabs(alpha) );  // d -> d^(|a|)
        } else {
            const double logImiIz = std::log(Imi/Iz);
            const double logImiIzSqr = logImiIz*logImiIz;
            const double logIpIz = std::log(Ip/Iz);
            const double logIpIzSqr = logIpIz*logIpIz;
            // polinomial: equations solved with Mathematica
            const double a1 = -(15*Imi - 15*Ip - 7*Imi*logImiIz + Imi*logImiIzSqr + 7*Ip*logIpIz - Ip*logIpIzSqr)/(16.*Iz);
            const double a2 = -3 + (3*Imi)/(2.*Iz) + (3*Ip)/(2.*Iz) - (9*Imi*logImiIz)/(16.*Iz) + (Imi*logImiIzSqr)/(16.*Iz) -
                (9*Ip*logIpIz)/(16.*Iz) + (Ip*logIpIzSqr)/(16.*Iz);
            const double a3 = (5*Imi)/(8.*Iz) - (5*Ip)/(8.*Iz) - (5*Imi*logImiIz)/(8.*Iz) + (Imi*logImiIzSqr)/(8.*Iz) + (5*Ip*logIpIz)/(8.*Iz) -
                (Ip*logIpIzSqr)/(8.*Iz);
            const double a4 = 3 - (3*Imi)/(2.*Iz) - (3*Ip)/(2.*Iz) + (7*Imi*logImiIz)/(8.*Iz) -
                (Imi*logImiIzSqr)/(8.*Iz) + (7*Ip*logIpIz)/(8.*Iz) - (Ip*logIpIzSqr)/(8.*Iz);
            const double a5 = (-3*Imi)/(16.*Iz) + (3*Ip)/(16.*Iz) + (3*Imi*logImiIz)/(16.*Iz) - (Imi*logImiIzSqr)/(16.*Iz) -
                (3*Ip*logIpIz)/(16.*Iz) + (Ip*logIpIzSqr)/(16.*Iz);
            const double a6 = -1 + Imi/(2.*Iz) + Ip/(2.*Iz) - (5*Imi*logImiIz)/(16.*Iz) + (Imi*logImiIzSqr)/(16.*Iz) - (5*Ip*logIpIz)/(16.*Iz) +
                (Ip*logIpIzSqr)/(16.*Iz);
            const double a = alpha;
            deltaN = 1 + a1*a + a2*a*a + a3*a*a*a + a4*a*a*a*a + a5*a*a*a*a*a + a6*a*a*a*a*a*a;
        }

    }
    else if (intCode==0) {

        //////////////////////////////////////////////////////////////
        // SHAPES
        // ======
        // Linear-polynomial
        // Extracted and adapted from RooFit
        //////////////////////////////////////////////////////////////

        if(std::fabs(alpha)>1){
            deltaN = 1. + std::fabs(alpha)*(deltaN - Iz)/Iz;
        } else {
            const double eps_plus = Ip - Iz;
            const double eps_minus = Iz - Imi;
            const double S = 0.5 * (eps_plus + eps_minus);
            const double A = 0.0625 * (eps_plus - eps_minus);
            double val = Iz + alpha * (S + alpha * A * ( 15 + alpha * alpha * (-10 + alpha * alpha * 3  ) ) );
            if (val < 0) val = 0.;
            if(Iz == Iz && Iz>0) deltaN = val/Iz;
            else deltaN = 0.;
        }
    }

    if(deltaN!=deltaN) deltaN = 1;  // to avoid nan
    if(deltaN<=0) deltaN = 0; //protection against negative values (can happen for linear extrapolation)

    return deltaN;
}

//___________________________________________________________
// function to get pre/post-fit agreement
std::pair<double,int> Region::GetChi2Test(const bool isPostFit,
                                          const TH1* h_data,
                                          const TH1* h_nominal,
                                          const std::vector< std::shared_ptr<TH1> >& h_up,
                                          const std::vector< string >& systNames,
                                          CorrelationMatrix *matrix ){
    const unsigned int nbins = h_nominal->GetNbinsX();
    int ndf = 0;
    for(unsigned int i=0;i<nbins;++i){
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos (e.g. because X-axis range contains non-physical values)
        ndf++;
    }
    //
    //Speed Up: remove irrelevant systematics (which would give in any case 0 correlation)
    std::vector< string > EffectiveSystNames;
    std::vector< unsigned int > EffectiveSystIndex;
    for(unsigned int n=0;n<systNames.size();++n){
      if(matrix!=nullptr){
        if (matrix->fNuisParIsThere[systNames[n]]) {
           EffectiveSystNames.push_back(systNames[n]);
           EffectiveSystIndex.push_back(n);
        }
        else WriteDebugStatus("GetChi2Test"," will skip syst. "+systNames[n]);
      }
      else {
        EffectiveSystNames.push_back(systNames[n]);
        EffectiveSystIndex.push_back(n);
      }
    }
    const unsigned int nsyst=EffectiveSystNames.size();
    //
    TMatrixD C(ndf,ndf);
    //
    int ibin = 0;
    int jbin = 0;
    for(unsigned int i=0;i<nbins;++i){
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
        const double ynom_i = h_nominal->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;++j){
            double sum = 0.;
            if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), j+1) != fBlindedBins.end()) continue;
            if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), j+1) != fBlindedBinsPostFit.end()) continue;
            if (std::find(fDropBins.begin(), fDropBins.end(), j+1) != fDropBins.end()) continue;
            if (h_data->GetBinContent(j+1) == 0. && h_nominal->GetBinContent(j+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
            const double ynom_j = h_nominal->GetBinContent(j+1);
            for(unsigned int n=0;n<nsyst;++n){ //n!=m, run only across correlated systs
                if (EffectiveSystNames[n].find("saturated_model") != std::string::npos) continue;
                const double ysyst_i_n = h_up[EffectiveSystIndex[n]]->GetBinContent(i+1);
                for(unsigned int m=0;m<nsyst;++m){
                    if (n==m) continue;
                    const double ysyst_j_m = h_up[EffectiveSystIndex[m]]->GetBinContent(j+1);
                    // more than Bill's suggestion: add correlation between systematics!!
                    double corr(0.);
                    if(matrix) corr = matrix->GetCorrelation(EffectiveSystNames[n],EffectiveSystNames[m]);
                    else continue;
                    sum += (ysyst_i_n-ynom_i) * corr * (ysyst_j_m-ynom_j);
                }
            }
            for(unsigned int n=0;n<systNames.size();++n){ // n==m, all systs, corr=1
                if (systNames[n].find("saturated_model") != std::string::npos) continue;
                const double ysyst_i_n = h_up[n]->GetBinContent(i+1), ysyst_j_n = h_up[n]->GetBinContent(j+1);
                sum += (ysyst_i_n-ynom_i) /* * 1.0 */ * (ysyst_j_n-ynom_j);
            }
            if(i==j && ynom_i>0) {
                sum += ynom_i;  // add stat uncertainty to diagonal
                sum += h_nominal->GetBinError(i+1) * h_nominal->GetBinError(i+1); // add MC stat as well
            }
            C[ibin][jbin] = sum;
            ++jbin;
        }
        ++ibin;
    }
    //
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    // Invert the matrix
    C.Invert();
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    double chi2 = 0.;
    ibin = 0;
    jbin = 0;
    for(unsigned int i=0;i<nbins;++i){
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
        const double ynom_i = h_nominal->GetBinContent(i+1);
        const double ydata_i = h_data->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;j++){
            if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), j+1) != fBlindedBins.end()) continue;
            if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), j+1) != fBlindedBinsPostFit.end()) continue;
            if (std::find(fDropBins.begin(), fDropBins.end(), j+1) != fDropBins.end()) continue;
            if (h_data->GetBinContent(j+1) == 0. && h_nominal->GetBinContent(j+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
            const double ynom_j = h_nominal->GetBinContent(j+1);
            const double ydata_j = h_data->GetBinContent(j+1);
            chi2 += (ydata_i - ynom_i)*C[ibin][jbin]*(ydata_j-ynom_j);
            ++jbin;
        }
        ++ibin;
    }
    //
    return std::make_pair(chi2,ndf);
}

//___________________________________________________________
// function to build total error band from:
// - a nominal histo (tot exepcted)
// - syst variation histos (eventually already scaled by post-fit pulls)
// - correlation matrix
// Note: if matrix = nullptr => no correlation, i.e. matrix = 1 (used for pre-fit, or to neglect correlation)
std::unique_ptr<TGraphAsymmErrors> BuildTotError( const TH1* const h_nominal,
                                                  const std::vector< std::shared_ptr<TH1> >& h_up,
                                                  const std::vector< std::shared_ptr<TH1> >& h_down,
                                                  const std::vector< string >& fSystNames,
                                                  CorrelationMatrix *matrix ){
    if(!h_nominal){
        WriteErrorStatus("BuildTotError","h_nominal not defined.");
        exit(EXIT_FAILURE);
    }
    if(h_up.size()!=h_down.size()){
        WriteErrorStatus("BuildTotError","h_up and h_down have different size.");
        exit(EXIT_FAILURE);
    }
    if(h_up.size()!=fSystNames.size()){
        WriteErrorStatus("BuildTotError","h_up and fSystNames have different size.");
        exit(EXIT_FAILURE);
    }
    //
    //Speed Up: remove irrelevant systematics (which would give in any case 0 correlation)
    std::vector< std::string > EffectiveSystNames;
    std::vector< std::size_t > EffectiveSystIndex;
    for(std::size_t n=0; n < fSystNames.size(); ++n){
        if(matrix!=nullptr){
            if (matrix->fNuisParIsThere[fSystNames[n]]) {
                EffectiveSystNames.push_back(fSystNames[n]);
                EffectiveSystIndex.push_back(n);
            }
        }
        else {
            EffectiveSystNames.push_back(fSystNames[n]);
            EffectiveSystIndex.push_back(n);
        }
    }
    //
    auto g_totErr = std::make_unique<TGraphAsymmErrors>(h_nominal);
    //
    // - loop on bins
    for(int i_bin=1; i_bin < h_nominal->GetNbinsX()+1; ++i_bin){
        double finalErrPlus(0.);
        double finalErrMinus(0.);
        // yieldNominal = h_nominal->GetBinContent(i_bin);
        // - loop on the syst, two by two, to include the correlations
        if (matrix!=nullptr) {
            double corr(0.);
            // if the correlation matrix is not provided, there are no off-diagonal contributions
            // to include, so these loops can be skipped
            for(std::size_t i_syst=0; i_syst < EffectiveSystNames.size(); ++i_syst){
                for(std::size_t j_syst=0; j_syst < i_syst; ++j_syst){
                    // the inner loop only runs up to i_syst-1, effectively going over half
                    // the correlation matrix (the other half is symmetric)
                    corr = matrix->GetCorrelation(EffectiveSystNames[i_syst],EffectiveSystNames[j_syst]);

                    const double errUp_i   = h_up[EffectiveSystIndex[i_syst]]  ->GetBinContent(i_bin);
                    const double errDown_i = h_down[EffectiveSystIndex[i_syst]]->GetBinContent(i_bin);
                    const double errUp_j   = h_up[EffectiveSystIndex[j_syst]]  ->GetBinContent(i_bin);
                    const double errDown_j = h_down[EffectiveSystIndex[j_syst]]->GetBinContent(i_bin);

                    //
                    // Symmetrize (seems to be done in Roostats ??)
                    //
                    const double err_i = (errUp_i - errDown_i)/2.;
                    const double err_j = (errUp_j - errDown_j)/2.;

                    //
                    // Compute the + and - variations
                    //
                    // needs a factor 2 at the end since we are looping over half the correlation matrix
                    finalErrPlus  += err_i * err_j * corr * 2;
                    finalErrMinus += err_i * err_j * corr * 2;
                }
            }
        }
        // now all diagonal el. of all systematics, corr = 1;
        for(std::size_t i_syst=0; i_syst<fSystNames.size(); ++i_syst){
            const double errUp_i   = h_up[i_syst]  ->GetBinContent(i_bin);
            const double errDown_i = h_down[i_syst]->GetBinContent(i_bin);

            //
            // Symmetrize (seems to be done in Roostats ??)
            //
            const double err_i = (errUp_i - errDown_i)/2.;

            //
            // Compute the + and - variations
            //
            finalErrPlus  += err_i * err_i ;
            finalErrMinus += err_i * err_i ;
        }
        // add stat uncertainty, which should have been stored as orignal bin errors in the h_nominal (if fUseStatErr is true)
        finalErrPlus  += h_nominal->GetBinError(i_bin) * h_nominal->GetBinError(i_bin);
        finalErrMinus += h_nominal->GetBinError(i_bin) * h_nominal->GetBinError(i_bin);


        g_totErr->SetPointEYhigh(i_bin-1,std::sqrt(finalErrPlus ));
        g_totErr->SetPointEYlow( i_bin-1,std::sqrt(finalErrMinus));
    }

    return g_totErr;
}

//___________________________________________________________
//
void Region::PrepareMorphScales(FitResults *fitRes, std::vector<double> *morph_scale, std::vector<double> *morph_scale_nominal) const{
    for(const auto& isample : fSampleHists) {
        // skip data
        if(isample->fSample->fType==Sample::DATA) continue;
        if(isample->fSample->fType==Sample::GHOST) continue;
        if(isample->fSample->fType==Sample::EFT) continue;
        for(std::size_t i_syst=0; i_syst < fSystNames.size(); ++i_syst){
            std::string systName    = fSystNames[i_syst];
            if(TRExFitter::NPMAP[systName]=="") TRExFitter::NPMAP[systName] = systName;

            if(isample->HasNorm(fSystNames[i_syst])){
                // if this norm factor is a morphing one
                if(fSystNames[i_syst].find("morph_")!=string::npos || isample->GetNormFactor(fSystNames[i_syst])->fExpression.first!=""){
                    std::vector<double> scales = Common::CalculateExpression(nullptr, fSystNames[i_syst], true, isample.get(), fitRes);
                    if (scales.empty()) {
                        WriteErrorStatus("Region::PrepareMorphScales", "Scales are empty");
                        exit(EXIT_FAILURE);
                    }
                    const double scaleNom = scales.at(0);
                    morph_scale->emplace_back(scaleNom);
                }
            }
            for(unsigned int i_nf=0;i_nf<isample->fSample->fNormFactors.size();i_nf++){
                NormFactor *nf = isample->fSample->fNormFactors[i_nf].get();
                // if this norm factor is a morphing one
                if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                    std::vector<double> scales = Common::CalculateExpression(nf, nf->fName, false, nullptr, nullptr);
                    if (scales.empty()) {
                        WriteErrorStatus("Region::PrepareMorphScales", "Scales are empty");
                        exit(EXIT_FAILURE);
                    }
                    const double scaleNom = scales.at(0);
                    morph_scale_nominal->emplace_back(scaleNom);
                }
            }
        }
    }
}

//___________________________________________________________
//
void Region::SystPruning(PruningUtil *pu){
    std::unique_ptr<TH1> hTot(nullptr);
    if(pu->GetStrategy() == 1){
        hTot = GetTotHist(false); // don't include signal
    }
    else if(pu->GetStrategy() == 2){
        hTot = GetTotHist(true); // include signal
    }
    for(auto& sh : fSampleHists){
        sh->SystPruning(pu,hTot.get());
        //
        // flag overall systematics as no shape also for pruning purposes
        for(auto& syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            if(syh->fSystematic->fNoPruning) {
                syh->fShapePruned = false;
                syh->fNormPruned  = false;
                continue;
            }
            if(syh->fSystematic->fType==Systematic::OVERALL){
                syh->fShapePruned = true;
            }
        }
        //
        // add by-hand dropping of shape or norm
        for(auto& syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            if(syh->fSystematic->fNoPruning) continue;
            if( Common::FindInStringVector(syh->fSystematic->fDropShapeIn,fName)>=0  ||
                Common::FindInStringVector(syh->fSystematic->fDropShapeIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropShapeIn, "all")>=0
            ){
                syh->fShapePruned = true;
            }
            if( Common::FindInStringVector(syh->fSystematic->fDropNormIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropNormIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropNormIn, "all")>=0
            ){
                syh->fNormPruned = true;
            }
        }
    }
    //
    // reference pruning
    for(auto& sh : fSampleHists){
        for(auto& syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            std::shared_ptr<Systematic> syst = syh->fSystematic;
            if(syst->fReferencePruning == "") continue;
            std::shared_ptr<SampleHist> refSmpH = GetSampleHist(syst->fReferencePruning);
            if(refSmpH==nullptr){
                WriteWarningStatus("Region::SystPruning", "Cannot find reference pruning sample: " + syst->fReferencePruning + " in region: " + fName);
                continue;
            }
            std::shared_ptr<SystematicHist> refSysH = refSmpH->GetSystematic(syst->fName);
            if(refSysH==nullptr){
                WriteWarningStatus("Region::SystPruning", "Cannot find systematic " + syst->fName + " for reference pruning sample " + syst->fReferencePruning);
                continue;
            }
            syh->fNormPruned = refSysH->fNormPruned;
            syh->fShapePruned = refSysH->fShapePruned;
            syh->fBadShape = refSysH->fBadShape;
            syh->fBadNorm = refSysH->fBadNorm;
        }
    }
}

//___________________________________________________________
//
std::unique_ptr<TH1> Region::GetTotHist(bool includeSignal) {
    std::unique_ptr<TH1> hTot(nullptr);
    for(auto& sh : fSampleHists){
        if(!sh->fSample) continue;
        if(sh->fSample->fType==Sample::GHOST) continue;
        if(sh->fSample->fType==Sample::DATA) continue;
        if(sh->fSample->fType==Sample::EFT) continue;
        if(!includeSignal && sh->fSample->fType==Sample::SIGNAL) continue;
        if(!sh->fHist) continue;
        TH1* hTmp = static_cast<TH1*>(sh->fHist->Clone(("hTot_"+fName).c_str()));
        // scale accoring to nominal SF (considering morphing as well)
        const double& scale = Common::GetNominalMorphScale(sh.get());
        hTmp->Scale(scale);
        if(!hTot) {
            hTot.reset(hTmp);
        } else {
            hTot->Add(hTmp);
        }
    }
    hTot->SetDirectory(nullptr);
    return hTot;
}

//___________________________________________________________
//
void Region::SavePreFitUncertaintyAndTotalMCObjects() {
    auto originalDir = gDirectory;
    TString uncBandFileName = TString::Format("%s/Histograms/%s%s_preFit.root", fFitName.c_str(), fName.c_str(), fSuffix.c_str());
    gSystem->mkdir((fFitName + "/Histograms").c_str());
    WriteDebugStatus("Region::DrawPreFit", TString::Format("Writing file %s", uncBandFileName.Data()).Data());
    std::unique_ptr<TFile> f(TFile::Open(uncBandFileName, "RECREATE"));
    f->cd();
    fErr->Write("", TObject::kOverwrite);
    fTot->Write("", TObject::kOverwrite);
    originalDir->cd();
    f->Close();
}

