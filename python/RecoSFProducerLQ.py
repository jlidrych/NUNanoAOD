import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class RecoSFProducerLQ(Module):
    def __init__(self,targetfile,name="MuonRecoSF",norm=True,verbose=False,doSysVar=True):
        self.targetfile = targetfile
        self.name = name
        self.norm = norm
        self.verbose = verbose
        self.doSysVar = doSysVar
    def loadHisto(self,filename,hname):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(None)
        tf.Close()
        return hist
    def beginJob(self):
	    pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        self.out.branch(self.name+"Up", "F")
        self.out.branch(self.name+"Down", "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        #        print "here we are.."
        """process event, return True (go to next module) or False (fail, go to next event)"""
        lep_cat = 0
        if hasattr(event,"lep_category"):
            lep_cat = int(getattr(event,"lep_category"))
        if lep_cat < 1 :
            l1_pt = 0
            l1_eta = 0
            l2_pt = 0
            l2_eta = 0
            l1_flavor = 0
        else:
            l1_pt = float(getattr(event,"lead_lep_p"))
            l1_eta = abs(float(getattr(event,"lead_lep_eta")))
            l2_pt = float(getattr(event,"trail_lep_p"))
            l2_eta = abs(float(getattr(event,"trail_lep_eta")))
        weight = 1
        weightError = 0

        if lep_cat==1:# or lep_cat==5 or lep_cat==7 : #these are MM. MML and MMLL lepton categories
            hist = self.loadHisto(self.targetfile,"RecoEffPvsEta")
            nxBins = 8
            nyBins = 2
            searchbin1x = -1
            searchbin1y = -1
            searchbin2x = -1
            searchbin2y = -1

            for xbin in range(1,nxBins+1):
                if l1_pt < hist.GetXaxis().GetBinLowEdge(1):
                    searchbin1x = 1
                    break
                if l1_pt > hist.GetXaxis().GetBinLowEdge(xbin) and l1_pt < (hist.GetXaxis().GetBinLowEdge(xbin) + hist.GetXaxis().GetBinWidth(xbin) ):
                    searchbin1x = xbin
                    break
            if l1_pt > hist.GetXaxis().GetBinLowEdge(nxBins):
                searchbin1x = nxBins

            for ybin in range(1,nyBins+1):
                if l1_eta < hist.GetYaxis().GetBinLowEdge(1):
                    searchbin1y = 1
                    break
                if l1_eta > hist.GetYaxis().GetBinLowEdge(ybin) and l1_eta < (hist.GetYaxis().GetBinLowEdge(ybin) + hist.GetYaxis().GetBinWidth(ybin) ):
                    searchbin1y = ybin
                    break
            if l1_eta > hist.GetYaxis().GetBinLowEdge(nyBins):
                searchbin1y = nyBins

            for xbin in range(1,nxBins+1):
                if l2_pt < hist.GetXaxis().GetBinLowEdge(1):
                    searchbin2x = 1
                    break
                if l2_pt > hist.GetXaxis().GetBinLowEdge(xbin) and l2_pt < (hist.GetXaxis().GetBinLowEdge(xbin) + hist.GetXaxis().GetBinWidth(xbin) ):
                    searchbin2x = xbin
                    break
            if l2_pt > hist.GetXaxis().GetBinLowEdge(nxBins):
                searchbin2x = nxBins

            for ybin in range(1,nyBins+1):
                if l2_eta < hist.GetYaxis().GetBinLowEdge(1):
                    searchbin2y = 1
                    break
                if l2_eta > hist.GetYaxis().GetBinLowEdge(ybin) and l2_eta < (hist.GetYaxis().GetBinLowEdge(ybin) + hist.GetYaxis().GetBinWidth(ybin) ):
                    searchbin2y = ybin
                    break
            if l2_eta > hist.GetYaxis().GetBinLowEdge(nyBins):
                searchbin2y = nyBins

            weight1 = hist.GetBinContent(searchbin1x,searchbin1y)
            weight1Err = hist.GetBinError(searchbin1x,searchbin1y)

            weight2 = hist.GetBinContent(searchbin2x,searchbin2y)
            weight2Err = hist.GetBinError(searchbin2x,searchbin2y)

            weight = weight1*weight2
            
            self.out.fillBranch(self.name, weight)
            if self.doSysVar:
                weightError = (weight1+weight1Err)*(weight2+weight2Err)
                self.out.fillBranch(self.name+"Up", weightError)
                weightError = (weight1-weight1Err)*(weight2-weight2Err)
                self.out.fillBranch(self.name+"Down",weightError)
            return True

        self.out.fillBranch(self.name,1.)
        if self.doSysVar:
            self.out.fillBranch(self.name+"Up", 1.)
            self.out.fillBranch(self.name+"Down", 1.)
        return True       

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
RecoSF_UL2017 = "%s/src/PhysicsTools/MonoZ/data/MuonRecoSF/MuonRecoSF-UL2017.root" % os.environ['CMSSW_BASE']
MuonRecoSF_UL2017 = lambda : RecoSFProducerLQ(RecoSF_UL2017,verbose=False, doSysVar=True)
