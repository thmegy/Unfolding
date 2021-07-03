import ROOT
import sys

def openRootFile(fName, option="read"):
  oFile = ROOT.TFile.Open(fName, option)
  if oFile == 0:
    if option == "read":
      print("ERROR: Cannot open file: %s"%fName)
    else: 
      print("ERROR: Cannot create file: %s"%fName)
    sys.exit(0)
  return oFile


def copyHists(source, target, hList):
  for iH in hList:
    source.cd()
    tmpH = source.Get(iH)
    if tmpH == 0:
      continue
    #copy hist
    target.cd()
    tmpH.Write()


def main():

  iFileName = "/afs/cern.ch/user/y/yabulait/workarea/multibAnalysis/stats/testHist/MVA/13TeV_TwoLepton_2tag2pjet_150ptv_SR_mva.root"
  outFileName = "/afs/cern.ch/user/y/yabulait/workarea/multibAnalysis/stats/testHist/MVA/13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root"

  iFile = openRootFile(iFileName, "read")
  outFile = openRootFile(outFileName, "recreate")

  #Copy nominal
  copyHists(iFile, outFile, ["Zbb"])

  #Copy Sys
  sysTDir = iFile.Get("Systematics")
  if sysTDir == 0:
    print("ERROR: can not get TDir: Systematics")
    iFile.Close()
    outFile.Close()
    sys.exit(0)

  sysSaveTDir = outFile.mkdir("Systematics")
  sysSaveTDir.cd()

  sysList = ["Zbb_SysJET_21NP_JET_EffectiveNP_2__1up", "Zbb_SysJET_JvtEfficiency__1up", "Zbb_SysJET_JvtEfficiency__1down"]
  copyHists(sysTDir, sysSaveTDir, sysList)
  
  iFile.cd()
  iFile.Close()
  outFile.cd()
  outFile.Close()

if __name__ == "__main__":
  main()

