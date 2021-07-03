#include "TRExFitter/LikelihoodScanManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/StatusLogbook.h"

#include "TRandom3.h"

#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

#include <algorithm>

LikelihoodScanManager::LikelihoodScanManager() :
    fScanMinX(999999),
    fScanMinY(-999999),
    fStepsX(30),
    fScanMaxX(999999),
    fScanMaxY(-999999),
    fStepsY(30),
    fUseOffset(true),
    fCPU(1),
    fParal2D(false),
    fParal2Dstep(-1),
    fUseNllInLHscan(true),
    fMinValX(99999),
    fMinValY(99999),
    fMaxValX(-99999),
    fMaxValY(-99999)
{
    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Generation);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::LinkStateMgmt);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Optimization);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Tracing);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Contents);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
        RooMsgService::instance().setStreamStatus(1,false);
    }

}

//__________________________________________________________________________________
//
LikelihoodScanManager::~LikelihoodScanManager() {
}
    
//__________________________________________________________________________________
//
LikelihoodScanManager::scanResult1D LikelihoodScanManager::Run1DScan(const RooWorkspace* ws,
                                                                     const std::string& varName,
                                                                     RooDataSet* data) const {

    if (!ws || !data) {
        WriteErrorStatus("LikelihoodScanManager::Run1DScan", "Passed nullptr for ws");
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("LikelihoodScanManager::Run1DScan", "Passed nullptr for mc");
        exit(EXIT_FAILURE);
    }
    FitUtils::FloatAllParameters(mc, true);
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    double min(-3);
    double max(3);

    bool found(false);
    RooRealVar* var(nullptr);
    {
        for (auto var_tmp : *mc->GetParametersOfInterest()) {
            var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varName) {
                WriteInfoStatus("LikelihoodScanManager::Run1DScan", "GetLikelihoodScan for NP = " + vname);
                found=true;
                min = var->getMin();
                max = var->getMax();
                break;
            }
        }
    }

    if (!found) {
        if (mc->GetNuisanceParameters()) {
            for (auto var_tmp : *mc->GetNuisanceParameters()) {
                var = static_cast<RooRealVar*>(var_tmp);
                const std::string vname = var->GetName();
                if (vname == varName || vname == "alpha_"+varName) {
                    WriteInfoStatus("LikelihoodScanManager::Run1DScan", "GetLikelihoodScan for NP = " + vname);
                    found=true;

                    if ((vname.find("gamma_") != std::string::npos)) {
                        min = 0.5;
                        max = 1.5;
                    }
                    break;
                }
            }
        }
    }

    if (fScanMinX < 99999) { // is actually set
        min = fScanMinX;
    }

    if (fScanMaxX > -99999) { // is actually set
        max = fScanMaxX;
    }

    LikelihoodScanManager::scanResult1D result;
    if (!found) {
        WriteWarningStatus("LikelihoodScanManager::Run1DScan", "Cannot find NP: " + varName);
        return result;
    }

    result.first.resize(fStepsX);
    result.second.resize(fStepsX);

    const auto offset = fUseOffset ? kTRUE : kFALSE;
    const RooArgSet* glbObs = mc->GetGlobalObservables();

    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (mc->GetNuisanceParameters()) {
        nll.reset(simPdf->createNLL(*data,
                                    RooFit::Constrain(*mc->GetNuisanceParameters()),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(offset),
                                    NumCPU(fCPU, RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE)));
    } else {
        nll.reset(simPdf->createNLL(*data,
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(offset),
                                    NumCPU(fCPU, RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE)));
    }
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
    
    double mnll = 9999999;
    var->setConstant(kTRUE); // make POI constant in the fit

    const std::size_t nFree = FitUtils::NumberOfFreeParameters(mc);

    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {

        RooMinimizer m(*nll); // get MINUIT interface of fit
        m.setPrintLevel(-1);
        m.setStrategy(1);
        m.optimizeConst(2);
        m.setEps(tol);

        WriteInfoStatus("LikelihoodScanManager::Run1DScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fStepsX) + " points");
        result.first[ipoint] = min+ipoint*(max-min)/(fStepsX - 1);
        *var = result.first[ipoint]; // set POI
        std::unique_ptr<RooFitResult> r(nullptr);
        if (nFree != 0) {
            m.migrad(); // minimize again with new posSigXsecOverSM value
            r.reset(m.save()); // save fit result
        }
        const double nllval = nll->getVal();
        if (fUseNllInLHscan || !r) {
            result.second[ipoint] = nllval;
        } else {
            result.second[ipoint] = r->minNll();
        }
        if (result.second[ipoint] < mnll) mnll = result.second[ipoint];
        WriteDebugStatus("LikelihoodScanManager::Run1DScan", "Point: " + std::to_string(result.first[ipoint]) + ", nll: " + std::to_string(result.second[ipoint]));
    }
    var->setConstant(kFALSE);

    for (auto & iY : result.second) {
        iY = iY - mnll;
    }

    TRandom3 rand(1234567);
    if (std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varName) != fBlindedParameters.end()) {
        const double rndNumber = rand.Uniform(5);
        for (auto& ix : result.first) {
            ix += rndNumber;
        }
    }

    return result;
}

//__________________________________________________________________________________
//
LikelihoodScanManager::Result2D LikelihoodScanManager::Run2DScan(const RooWorkspace* ws,
                                                                 const std::pair<std::string, std::string>& varNames,
                                                                 RooDataSet* data) {

    if (!ws || !data) {
        WriteErrorStatus("LikelihoodScanManager::Run2DScan", "Passed nullptr for ws");
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("LikelihoodScanManager::Run2DScan", "Passed nullptr for mc");
        exit(EXIT_FAILURE);
    }
    FitUtils::FloatAllParameters(mc, true);
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    int count = 0;
    RooRealVar* varX = nullptr;
    RooRealVar* varY = nullptr;
    if (mc->GetNuisanceParameters()) {
        for (auto var_tmp : *mc->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var;
                ++count;
            }
            if (count == 2) break;
        }
    }

    if (count != 2) {
        for (auto var_tmp : *mc->GetParametersOfInterest()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var;
                ++count;
            }
            if (count == 2) break;
        }
    }

    LikelihoodScanManager::Result2D result;

    if (count != 2) {
        WriteWarningStatus("LikelihoodScanManager::Run2DScan","Did not find the two parameters you want to use in the 2D likelihood scan");
        return result;
    }

    double minValX = varX->getMin();
    double maxValX = varX->getMax();
    double minValY = varY->getMin();
    double maxValY = varY->getMax();
    
    if (fScanMinX < 99999) { // is actually set
        minValX = fScanMinX;
    }
    if (fScanMinY < 99999) { // is actually set
        minValY = fScanMinY;
    }
    if (fScanMaxX > -99999) { // is actually set
        maxValX = fScanMaxX;
    }
    if (fScanMaxY > -99999) { // is actually set
        maxValY = fScanMaxY;
    }
    
    const auto offset = fUseOffset ? kTRUE : kFALSE;
    const RooArgSet* glbObs = mc->GetGlobalObservables();
    
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (mc->GetNuisanceParameters()) {
        nll.reset(simPdf->createNLL(*data,
                                    RooFit::Constrain(*mc->GetNuisanceParameters()),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(offset),
                                    NumCPU(fCPU, RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE)));
    } else {
        nll.reset(simPdf->createNLL(*data,
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(offset),
                                    NumCPU(fCPU, RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE)));
    }
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance();

    varX->setConstant(kTRUE); // make POI constant in the fit
    varY->setConstant(kTRUE); // make POI constant in the fit
    
    const std::size_t nFree = FitUtils::NumberOfFreeParameters(mc);

    result.x.resize(fStepsX);
    result.y.resize(fStepsY);
    result.z.resize(fStepsX);
    for (auto& iz : result.z) {
        iz.resize(fStepsY);
    }

    //values for parameter1, parameter2 and the NLL value
    double zmin = 9999999;
    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {
        //Set both POIs to constant
        if (fParal2D && ipoint != fParal2Dstep) continue;
        WriteInfoStatus("LikelihoodScanManager::Run2DScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fStepsX) + " points");
        result.x[ipoint] = minValX + ipoint * (maxValX - minValX) / (fStepsX - 1);
        *varX = result.x[ipoint]; // set POI
        for (int jpoint = 0; jpoint < fStepsY; ++jpoint) {
            RooMinimizer m(*nll); // get MINUIT interface of fit
            m.optimizeConst(2);
            m.setErrorLevel(-1);
            m.setPrintLevel(-1);
            m.setStrategy(1); // set precision to high
            m.setEps(tol);
            WriteInfoStatus("LikelihoodScanManager::Run2DScan","Running LHscan for subpoint " + std::to_string(jpoint+1) + " out of " + std::to_string(fStepsY) + " points");
            result.y[jpoint] = minValY + jpoint * (maxValY - minValY) / (fStepsY - 1);
            *varY = result.y[jpoint]; // set POI
            std::unique_ptr<RooFitResult> r(nullptr);
            if (nFree != 0) {
                m.migrad(); // minimize again with new posSigXsecOverSM value
                r.reset(m.save()); // save fit result
            }
            const double nllval = nll->getVal();
            double z_tmp(0);
            if (fUseNllInLHscan || !r) {
                z_tmp = nllval;
            } else {
                z_tmp = r->minNll();
            }
            result.z[ipoint][jpoint] = z_tmp;
            WriteDebugStatus("LikelihoodScanManager::Run2DScan", "Point x: " + std::to_string(result.x[ipoint]) + ", y: " + std::to_string(result.y[jpoint]) + ", nll: " + std::to_string(result.z[ipoint][jpoint]));

            // save the best values
            if (z_tmp < zmin) {
                zmin = z_tmp;
            }
        }
    }
    varX->setConstant(kFALSE);
    varY->setConstant(kFALSE);

    TRandom3 rand(1234567);
    const bool blindX = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.first) != fBlindedParameters.end();
    const bool blindY = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.second) != fBlindedParameters.end();
    const double rndNumber = rand.Uniform(5);
    if (blindX) {
        minValX += rndNumber;
        maxValX += rndNumber;
        for (auto& ix : result.x) {
            ix += rndNumber;
        }
    }

    if (blindY) {
        minValY += rndNumber;
        maxValY += rndNumber;
        for (auto& iy : result.y) {
            iy += rndNumber;
        }
    }

    fMinValX = minValX;
    fMinValY = minValY;
    fMaxValX = maxValX;
    fMaxValY = maxValY;

    return result;
}
