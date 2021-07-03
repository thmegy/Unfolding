
import sys,os,argparse
import ROOT
ROOT.gROOT.SetBatch()
import copy 

ENV_WORK_DIR=os.environ["WorkDir_DIR"]

ROOT.gSystem.Load("%s/lib/libCommonSystSmoothingToolLib.so"%ENV_WORK_DIR)

ROOT.gStyle.SetOptStat(0)

parse = argparse.ArgumentParser(description='Test Smoothing methods.',prog='testSmooth.py', add_help=True)
parse.add_argument('--inputFile', dest='inputFile', default='13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root',type=str, help='Input root file.')
parse.add_argument('--hnom', dest='hnom', default='Nominal',type=str, help='Name of nominal histogram')
parse.add_argument('--hsys', dest='hsys',default='', type=str, help='Name of systematic histogram')
parse.add_argument('--subdir', dest='subdir', default='',type=str, help='subdirectory in root file')
parse.add_argument('--methods', dest='methods', default='smoothRebinParabolic', type=str, choices=['smoothRebinMonotonic', 'smoothRebinParabolic', 'smoothDeltaUniformKernel','smoothTtresDependent','smoothTRExDefault'])

parse.add_argument('--outputFile', dest='outputFile', type=str, default='output.root', help="Name of plot <name.root>")
parse.add_argument('--doPlot',action='store_true', default=False)

parser=parse.parse_args()

smoothTool=ROOT.SmoothHist()
plotTool = ROOT.PlotHist()

smoothTool.setStatErrThreshold(0.05)
smoothTool.setTRExTolerance(0.08)
smoothTool.setTRExNbins(1)

class Method:
  def __init__(self, Name, LineStyle, LineColor, RatioLineStyle):
    self.Name=Name
    self.LineStyle=LineStyle
    self.LineColor=LineColor
    self.RatioLineStyle=RatioLineStyle

  def hStyle(self,hist):
    hist.SetLineColor(self.LineColor)
    hist.SetLineStyle(self.LineStyle)
    hist.SetTitle("h%s"%self.Name)
  def hRatioStyle(self,hist):
    hist.SetLineStyle(self.RatioLineStyle)
    hist.SetLineColor(self.LineColor)
    hist.SetLineWidth(3)


MethodCollection={"smoothRebinParabolic": Method("smoothRebinParabolic",10,30,7),
    "smoothRebinMonotonic":Method("smoothRebinMonotonic",10,31,7),
    "smoothRatioUniformKernel":Method("smoothRebinMonotoni",10,ROOT.kGreen+2,7),
    "smoothDeltaUniformKernel":Method("smoothDeltaUniformKernel",10,ROOT.kGreen+4,7),
    "smoothRatioGaussKernel":Method("smoothRatioGaussKernel",10,ROOT.kGreen+6,7),
    "smoothDeltaGaussKernel": Method("smoothDeltaGaussKernel", 10,ROOT.kGreen,7),
    "smoothTtresDependent": Method("smoothTtresDependent",10,ROOT.kMagenta+2,7),
    "smoothTtresIndependent": Method("smoothTtresIndependent",10,ROOT.kMagenta+4,7),
    "smoothTRExDefault": Method("smoothTRExDefault",10,ROOT.kYellow+2,7)}

#plot
def plot_hist(hnom, hsys, hsmooth):
    canv = ROOT.TCanvas("canvTest","Nom/Sys",800,800)
    pad1 = ROOT.TPad("pad1Test", "pad1", 0, 0.5, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetTicky(2)
    pad1.Draw()
    pad1.cd()
    hnom.GetXaxis().SetRangeUser(730, 4000.)
    hnom.Draw("Hist")
    
    hsys.SetLineColor(ROOT.kRed)
    hsys.Draw("Hist same")
    
    MethodCollection[parser.methods].hRatioStyle(hsmooth)
    hsmooth.Draw("Hist same")
    
    leg = ROOT.TLegend(0.6,0.3,0.9,0.9)
    leg.AddEntry(hnom, "nominal","l")
    leg.AddEntry(hsys,hsys.GetName(), "l")
    leg.AddEntry(hsmooth, hsmooth.GetName(),"l")

    leg.Draw()

    canv.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.498)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx() 
    pad2.SetGridy()
    pad2.SetTicky(2)
    pad2.SetTickx(2)
    pad2.Draw()
    pad2.cd()

    hnomR = plotTool.getPull(hnom,hnom, True)
    hnomR.SetFillStyle(3001)
    hnomR.SetFillColor(ROOT.kGray)
    hnomR.GetXaxis().SetTitle(hnom.GetXaxis().GetTitle())
    hsysR = plotTool.getPull(hnom,hsys)
    minR = 1.5*abs(hsysR.GetMinimum())
    maxR = 1.5*abs(hsysR.GetMaximum())
    maxR = max(minR, maxR)
    if maxR>10.: maxR=25
    if maxR>1.: 
      hnomR.SetMinimum(-maxR)
      hnomR.SetMaximum(maxR)
    hnomR.SetBit(ROOT.TH1.kNoTitle, True)
    hnomR.SetLineColor(ROOT.kBlack)
    hnomR.GetYaxis().SetLabelSize(15)
    hnomR.GetYaxis().SetLabelFont(43)
    hnomR.GetYaxis().SetTitleFont(43)
    hnomR.GetYaxis().SetTitleSize(20)
    hnomR.GetYaxis().SetTitleOffset(1.55)
    hnomR.GetYaxis().SetTitle("#frac{Sys-Nom}{Nom} [%]")
    hnomR.GetXaxis().SetLabelSize(15)
    hnomR.GetXaxis().SetLabelFont(43)
    hnomR.GetXaxis().SetTitleFont(43)
    hnomR.GetXaxis().SetTitleSize(25)
    hnomR.GetXaxis().SetTitleOffset(4.0)
    hnomR.GetXaxis().SetRangeUser(730,4000.)
    hnomR.Draw("hist e2")

    hsysR.SetLineColor(ROOT.kRed)
    hsysR.SetLineStyle(2)
    hsysR.SetLineWidth(3)
    hsysR.Draw("hist same")


    #for ihR in hRlist:
    #    ihR.Draw("hist same")
    hsmoothR=plotTool.getPull(hnom, hsmooth)
    MethodCollection[parser.methods].hRatioStyle(hsmoothR)
    hsmoothR.Draw("hist same")

    canv.Print("%s.pdf"%hsys.GetName())


def smooth(hnom, hSyslist):
    """
    h_list=smooth(hnom, hSyslist)
    #hnom --> nominal histogram
    #hSyslist --> systematic histogram
    #return smoothed histogram list
    """
    #ih.Clone("h%s"%(ih.GetName()+parser.methods))
    return [smoothTool.Smooth(hnom, ih,parser.methods)
            for ih in hSyslist]

def save_hist(hnom, hSyslist, subdir=None):
    """
    save_hist(hnom, hRlist)
    #hnom --> nominal histogram
    #hRlist --> smoothed by 
                several smoothing method 
    #save them to a root file. 
    """
    output_root_file = ROOT.TFile.Open(parser.outputFile,
                                       "recreate")
    output_root_file.cd()
    if subdir is not None:
        saveDir = ROOT.gDirectory.mkdir(subdir)
        saveDir.cd()
        
    hnom.Write()
    for ih in hSyslist:
         ih.Write()
            
    output_root_file.Close()
 

#start from here

subdir_ = parser.subdir #"OneTagBOneProbeB"
inFile = ROOT.TFile.Open(parser.inputFile)

hnom = inFile.Get(subdir_+"/"+parser.hnom)
#hnom.SetLineColor(ROOT.kBlue)

tdSys = inFile.Get(subdir_)
dkeys = tdSys.GetListOfKeys()
hSyslist_org_name=[dkeys.At(i).GetName() for i in range(dkeys.GetEntries())]
hSyslist_org=list()
for ihname_ in hSyslist_org_name:
    if ihname_ is "Nominal": continue
    hSyslist_org.append(tdSys.Get(ihname_))

# we got all hists

hSyslist_sm = smooth(hnom, copy.deepcopy(hSyslist_org))

if parser.doPlot:
    for org_,sm_ in list(zip(hSyslist_org, hSyslist_sm)):
        plot_hist(hnom, org_,sm_)

save_hist(hnom, hSyslist_sm, subdir_)