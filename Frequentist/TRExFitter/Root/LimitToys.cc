#include "TRExFitter/LimitToys.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "RooRandom.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"

#include "RooStats/SamplingDistPlot.h"

#include <iostream>
#include <memory>

LimitToys::LimitToys() :
    fNtoysSplusB(1000),
    fNtoysB(1000),
    fLimit(0.95),
    fScanSteps(21),
    fScanMin(0),
    fScanMax(10),
    fToysSeed(1234),
    fPlot(true),
    fFile(true),
    fOutputPath("")
{    
}


void LimitToys::RunToys(RooAbsData* data,
                        RooStats::ModelConfig* mcSplusb,
                        RooStats::ModelConfig* mcB) const {

    if (!data || !mcSplusb || !mcB) {
        WriteErrorStatus("LimitToys::RunToys", "Data or one of the ModelConfigs is nullptr!");
        return;
    }

    RooRandom::randomGenerator()->SetSeed(fToysSeed);

    RooRealVar* poi = static_cast<RooRealVar*>(mcB->GetParametersOfInterest()->first());
    mcB->SetSnapshot(*poi);

    RooStats::FrequentistCalculator freqCalc(*data, *mcB, *mcSplusb);
    
    std::unique_ptr<RooStats::ProfileLikelihoodTestStat> plr = std::make_unique<RooStats::ProfileLikelihoodTestStat>(*mcSplusb->GetPdf());
    plr->SetOneSided(true);

    const RooArgSet* glbObs = mcSplusb->GetGlobalObservables();
    RooStats::ToyMCSampler *toymcs = static_cast<RooStats::ToyMCSampler*>(freqCalc.GetTestStatSampler());
    if (glbObs) {
        toymcs->SetGlobalObservables(*glbObs);
    }
    toymcs->SetTestStatistic(plr.get());

    if (!mcSplusb->GetPdf()->canBeExtended()) {
        toymcs->SetNEventsPerToy(1);
    }

    freqCalc.SetToys(fNtoysSplusB, fNtoysB);

    RooStats::HypoTestInverter inverter(freqCalc);
    inverter.SetConfidenceLevel(fLimit);
    inverter.UseCLs(true);
    inverter.SetVerbose(false);
    inverter.SetFixedScan(fScanSteps,fScanMin,fScanMax);
    
    WriteInfoStatus("LimitToys::RunToys","Running " + std::to_string(fNtoysSplusB) + "/" + std::to_string(fNtoysB) + " toys for limits.");
    RooStats::HypoTestInverterResult* result = inverter.GetInterval();

    WriteInfoStatus("LimitToys::RunToys", "----------------------------------------------------");
    WriteInfoStatus("LimitToys::RunToys", "Results of the toys for limits");
    WriteInfoStatus("LimitToys::RunToys", "----------------------------------------------------");
    WriteInfoStatus("LimitToys::RunToys",std::to_string(100*inverter.ConfidenceLevel()) + "%  upper limit : " + std::to_string(result->UpperLimit()));
    WriteInfoStatus("LimitToys::RunToys","Expected upper limits, using the B (alternate) model : ");
    WriteInfoStatus("LimitToys::RunToys"," expected limit (median) " + std::to_string(result->GetExpectedUpperLimit(0)));
    WriteInfoStatus("LimitToys::RunToys"," expected limit (-1 sig) " + std::to_string(result->GetExpectedUpperLimit(-1)));
    WriteInfoStatus("LimitToys::RunToys"," expected limit (+1 sig) " + std::to_string(result->GetExpectedUpperLimit(1)));
    WriteInfoStatus("LimitToys::RunToys"," expected limit (-2 sig) " + std::to_string(result->GetExpectedUpperLimit(-2)));
    WriteInfoStatus("LimitToys::RunToys"," expected limit (+2 sig) " + std::to_string(result->GetExpectedUpperLimit(2)));
    WriteInfoStatus("LimitToys::RunToys", "----------------------------------------------------");

    if (fPlot) {
        WriteInfoStatus("LimitToys::RunToys", "Printing output plot to " + fOutputPath);
        TCanvas c{};
        auto plot = std::make_unique<RooStats::HypoTestInverterPlot>("HTI_Result_Plot","Limits from toys",result);
        plot->Draw("CLb 2CL");
        c.Draw();
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            c.Print((fOutputPath+"/LimitToysResult."+format).c_str());
        }

        WriteInfoStatus("LimitToys::RunToys","Printing test-statistic output plot to "+fOutputPath);
        TCanvas cteststat{};
        const int nEntries = result->ArraySize();
        if (nEntries > 1) {
            int ny = (int) std::ceil( std::sqrt(nEntries) );
            int nx = (int) std::ceil( double(nEntries)/ny );
            cteststat.Divide( nx,ny);
        }
        for (int i=0; i<nEntries; i++) {
            if (nEntries > 1) cteststat.cd(i+1);           
            auto pl = plot->MakeTestStatPlot(i);
            pl->SetLogYaxis(true);
            pl->Draw();
        }
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            cteststat.Print((fOutputPath+"/LimitToysSampleDistr."+format).c_str());
        }
        WriteInfoStatus("LimitToys::RunToys","Printing test-statistic output plot to "+fOutputPath);
    }

    if (fFile) {
        RootOutput(result);
    }
}

void LimitToys::RootOutput(RooStats::HypoTestInverterResult* result) const {
    std::unique_ptr<TFile> file(TFile::Open((fOutputPath+"/Toys.root").c_str(), "RECREATE")); 
    if (!file) {
        WriteErrorStatus("LimitToys::RootOutput", "Cannot open file at: " + fOutputPath+"/Toys.root");
        return;
    }
    result->Write();
    
    TTree *tree = new TTree("stats","LimitToys");
    Float_t mu;
    Float_t CLb;
    Float_t CLs;
    Float_t CLsplusb;
    Float_t CLsplusbError;
    Float_t upperLimit;
    Float_t expectedLimit;
    Float_t expectedLimit_plus1;
    Float_t expectedLimit_plus2;
    Float_t expectedLimit_minus1;
    Float_t expectedLimit_minus2;

    tree->Branch("mu", &mu);
    tree->Branch("CLb", &CLb);
    tree->Branch("CLs", &CLs);
    tree->Branch("CLsplusb", &CLsplusb);
    tree->Branch("CLsplusbError", &CLsplusbError);
    tree->Branch("observedLimit", &upperLimit);
    tree->Branch("expectedLimit", &expectedLimit);
    tree->Branch("expectedLimit_plus1", &expectedLimit_plus1);
    tree->Branch("expectedLimit_plus2", &expectedLimit_plus2);
    tree->Branch("expectedLimit_minus1", &expectedLimit_minus1);
    tree->Branch("expectedLimit_minus2", &expectedLimit_minus2);

    for (int i = 0; i < result->ArraySize(); ++i) {
        mu                   = result->GetXValue(i);
        CLb                  = result->CLb(i);
        CLs                  = result->CLs(i);
        CLsplusb             = result->CLsplusb(i);
        CLsplusbError        = result->CLsplusbError(i);
        upperLimit           = result->UpperLimit();
        expectedLimit        = result->GetExpectedUpperLimit(0);
        expectedLimit_plus1  = result->GetExpectedUpperLimit(1);
        expectedLimit_plus2  = result->GetExpectedUpperLimit(2);
        expectedLimit_minus1 = result->GetExpectedUpperLimit(-1);
        expectedLimit_minus2 = result->GetExpectedUpperLimit(-2);
        tree->Fill();
    }

    tree->SetDirectory(file.get());
    file->Write();
    file->Close();
}
