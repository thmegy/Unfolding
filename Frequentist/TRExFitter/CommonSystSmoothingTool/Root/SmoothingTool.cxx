#include "SmoothSystematics/SmoothingTool.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

#include <algorithm>
#include <fstream>

#include "SmoothSystematics/PlotHist.h"
//#include "TLegend.h"

SmoothingTool::~SmoothingTool() {
}

void SmoothingTool::addInputFileNom(const string NomFileName,
                                    const string TdirName) {
  TdirectoryNom = TdirName;
  TFile *testOpen = TFile::Open(NomFileName.data());
  if (!testOpen->IsOpen()) {
    cerr << "ERROR: can not open file: " << NomFileName << endl;
  } else {
    testOpen->Close();
  }
  delete testOpen;
  InputFileNomName = NomFileName;
  if (InputFileSysName.empty())
    InputFileSysName = NomFileName;
  if (OutputFileSysName.empty())
    OutputFileSysName = NomFileName;
}

void SmoothingTool::addInputFileSys(const string SysFileName,
                                    const string TdirName) {
  TdirectorySys = TdirName;
  TFile *testOpen = TFile::Open(SysFileName.data());

  if (!testOpen->IsOpen()) {
    cerr << "ERROR: can not open file: " << SysFileName << endl;
  } else {
    testOpen->Close();
  }
  delete testOpen;
  // will write smoothed hist to original file
  InputFileSysName = SysFileName;
  if (InputFileNomName.empty())
    InputFileNomName = SysFileName;
  if (OutputFileSysName.empty())
    OutputFileSysName = SysFileName;
}

void SmoothingTool::addOutputFileSys(const string SysFileName) {
  // TdirectorySys = TdirName;
  OutputFileSysName = SysFileName;
}

void SmoothingTool::addHist(string hnomName, vector<string> hsysNameList) {
  ListOfHistToSmooth[hnomName] = hsysNameList;
}
bool SmoothingTool::skip(string checkSkip) {
  if (checkSkip.empty())
    return false;
  for (auto &iStr : SkipList) {
    if (iStr == checkSkip)
      return true;
    std::size_t found = checkSkip.find(iStr);
    if (found != std::string::npos)
      return true;
  }
  return false;
}
//MainStr contain at least one of SubStr.
bool findSubStr(string MainStr, vector<string> SubStr){
  std::size_t found = std::string::npos;
  bool result = false;
  for(auto &iStr : SubStr){
    found = MainStr.find(iStr);
    if (found != std::string::npos){
      result = true;
      return result;
    }
  }
  return result;
}

void SmoothingTool::generateList(int FirstN) {

  TFile *InputFileNom = TFile::Open(InputFileNomName.data(), "read");

  // Check if input file for systematic hists are same as nom hists.
  TFile *InputFileSys = InputFileNom;
  if (InputFileSysName != InputFileNomName)
    InputFileSys = TFile::Open(InputFileSysName.data(), "read");
  TObject *objNom;
  TObject *objSys;
  TKey *keyNom;
  TKey *keySys;
  TDirectory *tdNom = InputFileNom;
  TDirectory *tdSys = InputFileSys;

  TIter nextNom(InputFileNom->GetListOfKeys());
  if (!TdirectoryNom.empty()) {
    tdNom = (TDirectory *)InputFileNom->Get(TdirectoryNom.data());
    if (tdNom != nullptr)
      nextNom = tdNom->GetListOfKeys();
    else {
      cerr << "Can not get TDirectory: " << TdirectoryNom << endl;
      exit(-1);
    }
  }

  TIter nextSys(InputFileSys->GetListOfKeys());
  if (!TdirectorySys.empty()) {
    tdSys = (TDirectory *)InputFileSys->Get(TdirectorySys.data());
    if (tdSys != nullptr)
      nextSys = tdSys->GetListOfKeys();
    else {
      cerr << "Can not get TDirectory: " << TdirectorySys << endl;
      exit(-1);
    }
  }
  // loop over Nom.
  int counter = 0;
  tdNom->cd();
  while ((keyNom = (TKey *)nextNom())) {
    objNom = keyNom->ReadObj();
    // TClass *clNom = gROOT->GetClass(objName.data());
    if (!objNom->InheritsFrom("TH1"))
      continue;
    // check if this is in skip list.
    string NomName = objNom->GetName();
    if (skip(NomName)) continue;

    //if not a Nominal hist.
    if (findSubStr(NomName,{"Sys", "sys", "up", "down", "Up", "Down"})) continue;
    // loop over Sys
    tdSys->cd(); 
    while ((keySys = (TKey *)nextSys())) {
      if ((FirstN > 0) && (counter >= FirstN))
        return;
      objSys = keySys->ReadObj();
      if (!objSys->InheritsFrom("TH1"))
        continue;
      // check skip list.
      string SysName = objSys->GetName();
      if (skip(SysName))
        continue;
      
      //if it is not a Systematics hist.
      if (!findSubStr(SysName,{"Sys", "sys", "up", "down", "Up", "Down"})) continue;
      // check if this is a sys of this Nominal sample.
      std::size_t found = SysName.find_first_of(NomName + "_");
      if (found == std::string::npos)
        continue;
      // if every thing is Ok, store the list of names of nominal ans
      // systimatice histograms
      ListOfHistToSmooth[NomName].push_back(
          SysName); // here we have name of Nominal and corresponding systimatic
                    // hists names.
      // count number of histograms.
      if (FirstN > 0)
        counter++;
    } // end of Sys loop.
  }   // end of Nom loop
  if (InputFileNom->IsOpen())
    InputFileNom->Close();
  if (InputFileSys->IsOpen())
    InputFileSys->Close();
} // end of generateList.

// Run Smoothing
bool SmoothingTool::runSmoothing() {
  
  TFile *InputFileNom = TFile::Open(InputFileNomName.data(), "update");

  // Check if input file for systematic hists are same as nom hists.
  TFile *InputFileSys = InputFileNom;
  if (InputFileSysName != InputFileNomName)
    InputFileSys = TFile::Open(InputFileSysName.data(), "update");
  // Check if output file different.
  TFile *OutputFileSys = InputFileSys;
  if (OutputFileSysName != InputFileSysName) {
    OutputFileSys = TFile::Open(OutputFileSysName.data(), "recreate");
    OutputFileSys->cd(); // this is become a current dir.
    // copy original hist
    copyHist(InputFileNom);
    OutputFileSys->cd();
    // OutputFileSys->Write();
  }

  // Loop over list of Hists for smooth.
  TDirectory *tmpTdNom = InputFileNom;
  if (!TdirectoryNom.empty()) {
    tmpTdNom = (TDirectory *)InputFileNom->Get(TdirectoryNom.data());
    if (tmpTdNom == nullptr) { //  hnom = (TH1*)tmpTdNom->Get(NomName.data());
      cerr << "Can not get TDirectory: " << TdirectorySys << endl;
      exit(-1);
    }
  }

  TDirectory *tmpTdSys = InputFileSys;
  if (!TdirectorySys.empty()) {
    tmpTdSys = (TDirectory *)InputFileSys->Get(TdirectorySys.data());
    if (tmpTdSys == nullptr) { // hsys = (TH1*)tmpTD->Get(isysName.data());
      cerr << "Can not get TDirectory: " << TdirectorySys << endl;
      exit(-1);
    }
  }

  TDirectory *outDir = OutputFileSys;
  if (!TdirectorySys.empty()) {
    outDir = (TDirectory *)OutputFileSys->Get(TdirectorySys.data());
    if (outDir == nullptr) {
      cerr << "Can not get TDirectory: " << TdirectorySys << endl;
      exit(-1);
    }
  }

  for (auto const &ipair : ListOfHistToSmooth) {
    string NomName = ipair.first;
    tmpTdNom->cd();
    TH1 *hnom = (TH1 *)tmpTdNom->Get(NomName.data());
    if (hnom == nullptr) {
      cerr << "ERROR: can not get Nominal Hist: " << NomName << endl;
      continue;
    }
    
    if(rebinFactor > 0){
      hnom->Rebin(rebinFactor);
    }

    // loop over sys names.
    for (auto const &isysName : ipair.second) {
      // TH1* hsys = nullptr;
      tmpTdSys->cd();
      TH1 *hsys = (TH1 *)tmpTdSys->Get(isysName.data());

      if (hsys == nullptr) {
        cerr << "ERROR: can not get Systematic Hist: " << isysName << endl;
        continue;
      }
      if(rebinFactor){
        hsys->Rebin(rebinFactor);
      }
      string hsysNameSmooth = string(hsys->GetName()) + SmoothedHistName;
      // Here call smothHist alg.
      TH1 *hsmoothed = SmoothAlg.Smooth(hnom, (TH1 *)hsys->Clone(hsysNameSmooth.c_str()), SmoothingOption, smoothedHistErrors);
      /*
      if (SmoothingOption == "merge") {
        hsmoothed = SmoothAlg.smoothHistogram(
            hnom, (TH1 *)hsys->Clone(hsysNameSmooth.c_str()), true);
      } else if (SmoothingOption == "kernel") {
        hsmoothed = SmoothAlg.smoothWithKernel(
            hnom, (TH1 *)hsys->Clone(hsysNameSmooth.c_str()));
      } else if (SmoothingOption == "Ttres"){
        hsmoothed = SmoothAlg.Smooth_Ttres(hnom, (TH1 *)hsys->Clone(hsysNameSmooth.c_str()), independentVar);
      
      } else if (SmoothingOption == "TRExDefault"){
        hsmoothed = SmoothAlg.Smooth_maxVariations(hnom, (TH1 *)hsys->Clone(hsysNameSmooth.c_str()), TREx_nbins);
      
      }else {
        cerr << "Wrong Smoothing option, choose merge, kernel, Ttres or TRExDefault!" << endl;
        exit(-1);
      }*/
      if (!hsmoothed){
        exit(-1);
      }

      // here write out the smoothed hist to file.
      outDir->cd();
      hsmoothed->Write();

      cout << "Smoothed sys: " << hsmoothed->GetName() << endl;
      // plot hist befor delete them.
      if (MakePlot)
        PlotAlg->plotWithRatio(hnom, hsys, hsmoothed);
      delete hsmoothed;
      delete hsys;
      // cout<<"Smoothed sys: "<<isysName<<endl;
    } // over sys.
    // if(tmpTdSys) delete tmpTdSys;
    if(hnom) delete hnom;
  } // over list of Hist.
  if (OutputFileSys->IsOpen()) {
    OutputFileSys->Write();
    OutputFileSys->Close();
  }
  if (InputFileNom->IsOpen())
    InputFileNom->Close();
  if (InputFileSys->IsOpen())
    InputFileSys->Close();
  // end of plot.
  if (PlotAlg) {
    PlotAlg->~PlotHist();
    //delete PlotAlg;
  }

  return true;
}

// Copy source directory to target.
void SmoothingTool::copyHist(TDirectory *source) {
  // copy things to current directory.
  TDirectory *savdir = gDirectory;
  // TDirectory *adir = (TDirectory*) savdir->Get(source->GetName());
  // if(!adir) adir = savdir->mkdir(source->GetName());
  // adir->cd();
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey *)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl)
      continue;
    if (cl->InheritsFrom("TDirectory")) {
      TDirectory *tmpDir = (TDirectory *)key->ReadObj();
      TDirectory *adir = savdir->mkdir(tmpDir->GetName());
      adir->cd();
      copyHist(tmpDir);
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      savdir->cd();
      obj->Write();
      delete obj;
    }

  } // end while
  savdir->SaveSelf(kTRUE);
  savdir->cd(); // return to target dir.
}
