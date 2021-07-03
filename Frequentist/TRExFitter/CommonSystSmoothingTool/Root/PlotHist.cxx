
//#include <io.h>   // For access().
#include <sys/stat.h>  // For stat().
#include <sys/types.h> // For stat().
#include <algorithm>

#include "SmoothSystematics/PlotHist.h"

//#include "TPad.h"
#include "TLegend.h"
#include "TObject.h"

#include "TCanvas.h"

using namespace std;
PlotHist::PlotHist(string pFile) {
  if (pFile.empty()) {
    pFile = "plots.root";
  }

  struct stat st;
  if (stat("plots/", &st) != 0) {
    mkdir("plots/", S_IRUSR | S_IWUSR | S_IXUSR);
  }
  pFile = "plots/" + pFile;

  OutputPlotFile = TFile::Open(pFile.data(), "recreate");
}

PlotHist::~PlotHist() {
  if (OutputPlotFile->IsOpen())
    OutputPlotFile->Close();
  // delete OutputPlotFile;
}

void PlotHist::plot(TH1 *phnom, TH1 *hsys, TH1 *hsysSmooth) {
  OutputPlotFile->cd();

  // open canvas for plots without ratio
  TCanvas *Canv = new TCanvas("Canv", "", 800, 600);
  Canv->cd();
  if (RebinFactor == 0) {
    RebinFactor = phnom->GetNbinsX() / 10;
  }
  else if(RebinFactor < 0){
    RebinFactor = 1;
  }
  phnom->Rebin(RebinFactor);
  hsys->Rebin(RebinFactor);
  hsysSmooth->Rebin(RebinFactor);

  if (Normalize) {
    phnom->Scale(1. / phnom->Integral());
    hsys->Scale(1. / hsys->Integral());
    hsysSmooth->Scale(1. / hsysSmooth->Integral());
  }
  float max = phnom->GetMaximum();
  phnom->SetMaximum(max * 1.5);
  phnom->SetMinimum(0);
  // set colors
  hsys->SetLineColor(kRed);
  hsysSmooth->SetLineColor(30);
  hsysSmooth->SetLineStyle(10);

  phnom->Draw("Hist");
  hsys->Draw("Hist same");
  hsysSmooth->Draw("hist same");

//    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
//    leg->AddEntry(phnom, "nominal","l");
//    leg->AddEntry(hsys,"original", "l");
//    leg->AddEntry(hsysSmooth, "smoothed", "l");
//    leg->Draw();
    

  string outputName = hsys->GetName();
  outputName = "plots/" + outputName;
  string forPDF = outputName + ".pdf";
  Canv->Print(forPDF.data());
  OutputPlotFile->cd();
  Canv->Write(outputName.data());
  Canv->Clear();
  delete Canv;
}

TH1* PlotHist::getPull(TH1* hNOM, TH1* hSYS, bool fillErr){
  TH1F* hPull =(TH1F*)hNOM->Clone("hPull"); //clone empty hist.
  hPull->Reset("ICESM"); 
  int nBins = hNOM->GetNbinsX();
  for (int iBin=0; iBin<=nBins+1; iBin++){
    float bConNom = hNOM->GetBinContent(iBin);
    float bConSys = hSYS->GetBinContent(iBin);
    if(bConNom == 0 ){
      hPull->SetBinContent(iBin,0.);
    }else{
      hPull->SetBinContent(iBin,(bConSys-bConNom)*100./bConNom);
    }
    if(fillErr){
      float bErrNom = hNOM->GetBinError(iBin);
      if(bConNom == 0 ){
        hPull->SetBinError(iBin,0.);
      }else{
        hPull->SetBinError(iBin,bErrNom*100./bConNom);
      }
    }else{
      hPull->SetBinError(iBin,0);
    }
  }
  return hPull;
}

void PlotHist::plotWithRatio(TH1 *phnom, TH1 *hsys, TH1 *hsysSmooth) {
  OutputPlotFile->cd();

  // open canvas for plots with ratio
  TCanvas *CanvRatio = new TCanvas("CanvRatio", "", 800, 800);
  CanvRatio->cd();

  if(hsysSmooth->GetNbinsX() < hsys->GetNbinsX()){
    int tmpRebinFactor = hsys->GetNbinsX()/hsysSmooth->GetNbinsX();
    hsys->Rebin(tmpRebinFactor);
  }
  if(hsysSmooth->GetNbinsX() < phnom->GetNbinsX()){
    int tmpRebinFactor = phnom->GetNbinsX()/hsysSmooth->GetNbinsX();
    phnom->Rebin(tmpRebinFactor);
  }
  
  if (RebinFactor == 0) {
    RebinFactor = phnom->GetNbinsX() / 10;
    hsysSmooth->Rebin(RebinFactor);
    phnom->Rebin(RebinFactor);
    hsys->Rebin(RebinFactor);
  }
  else if(RebinFactor > 0){ 
    hsysSmooth->Rebin(RebinFactor);
    phnom->Rebin(RebinFactor);
    hsys->Rebin(RebinFactor);
  }
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetTicky(2);
  //pad1->SetTickx(2);
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();

  float max = phnom->GetMaximum();
  phnom->GetYaxis()->SetRangeUser(0.0, max * 1.5);

  // set colors
  hsys->SetLineColor(kRed);
  string hTitle = hsys->GetTitle();
  string toErase = phnom->GetTitle();
  size_t pos = hTitle.find(toErase);
  hTitle.erase(pos,toErase.length()+1 );
  hsys->SetTitle(hTitle.data());
  hsysSmooth->SetLineColor(30);
  hsysSmooth->SetLineStyle(10);
  hsysSmooth->SetTitle("Smoothed");
  
  if(Normalize){
    phnom->DrawNormalized("Hist");
    hsys->DrawNormalized("Hist same");
    hsysSmooth->DrawNormalized("hist same");
  }else {
    phnom->Draw("Hist");
    hsys->Draw("Hist same");
    hsysSmooth->Draw("hist same");
  }
  TLegend* tleg = pad1->BuildLegend(0.6,0.6,0.9,0.9);
    //leg.AddEntry(phnom, "nominal","l");
    //leg.AddEntry(hsys,"original", "l");
    //leg.AddEntry(hsysSmooth, "smoothed", "l");
  tleg->Draw();
  
  CanvRatio->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.298);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->SetTicky(2);
  pad2->SetTickx(2);
  pad2->Draw();
  pad2->cd();
  TH1* hnomR = getPull(phnom,phnom, true);
  hnomR->SetFillStyle(3001);
  hnomR->SetFillColor(kGray);
  hnomR->GetXaxis()->SetTitle(phnom->GetXaxis()->GetTitle());
  TH1 *hsysR = getPull(phnom,hsys); //(TH1 *)hsys->Clone("Sys");
  //hsysR->Scale(100.); // trun to %
  //hnomR->Scale(100.);
  float minR = 1.3*abs(hsysR->GetMinimum());
  float maxR = 1.3*abs(hsysR->GetMaximum());
  maxR = std::max(minR, maxR);
  if(maxR>1.){
    hnomR->GetYaxis()->SetRangeUser(-maxR, maxR);
    //hnomR->SetMinimum(-maxR);  // Define Y ..
    //hnomR->SetMaximum(maxR); // .. range
  }
  hnomR->SetBit(TH1::kNoTitle, true);
  hnomR->SetLineColor(kBlack);
  hnomR->GetYaxis()->SetLabelSize(15);
  hnomR->GetYaxis()->SetLabelFont(43);
  hnomR->GetYaxis()->SetTitleFont(43);
  hnomR->GetYaxis()->SetTitleSize(20);
  hnomR->GetYaxis()->SetTitleOffset(1.55);
  hnomR->GetYaxis()->SetTitle("#frac{Sys-Nom}{Nom} [%]");
  hnomR->GetXaxis()->SetLabelSize(15);
  hnomR->GetXaxis()->SetLabelFont(43);
  //hnomR->GetXaxis()->SetTitle();
  hnomR->GetXaxis()->SetTitleFont(43);
  hnomR->GetXaxis()->SetTitleSize(25);
  hnomR->GetXaxis()->SetTitleOffset(4.0);

  hnomR->Draw("hist e2");
  
  hsysR->SetLineColor(kRed);
  hsysR->SetLineStyle(2);
  hsysR->SetLineWidth(2);
  hsysR->Draw("hist same");
  TH1 *hsmoothR = getPull(phnom,hsysSmooth); //(TH1 *)hsysSmooth->Clone("smoothedR");
  //hsmoothR->Scale(100.);
  hsmoothR->SetLineColor(30);
  hsmoothR->SetLineStyle(7);
  hsmoothR->SetLineWidth(2);
  hsmoothR->Draw("hist same");

  CanvRatio->cd();
  string outputName = hsys->GetName();
  outputName = "plots/" + outputName + "Ratio";
  string forPDF = outputName + ".pdf";
  CanvRatio->Print(forPDF.data());
  OutputPlotFile->cd();
  CanvRatio->Write(outputName.data());
  CanvRatio->Clear();
  delete CanvRatio;
}
