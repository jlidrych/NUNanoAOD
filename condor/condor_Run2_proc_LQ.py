import os, sys, re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# import PSet
import yaml
#Import the NanoAOD-tools that we will need
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme import jetRecalib
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
#Import the MonoZ analysis tools
from PhysicsTools.MonoZ.LeptoquarkProducer import *
from PhysicsTools.MonoZ.GenWeightProducer import *
from PhysicsTools.MonoZ.RecoSFProducerLQ import *
from PhysicsTools.MonoZ.TriggerSFProducerLQ import *
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-isMC'   , '--isMC'   , type=int, default=1     , help="")
parser.add_argument('-jobNum' , '--jobNum' , type=int, default=1     , help="")
parser.add_argument('-era'    , '--era'    , type=str, default="2018", help="")
parser.add_argument('-doSyst' , '--doSyst' , type=int, default=0     , help="")
parser.add_argument('-infile' , '--infile' , type=str, default=None  , help="")
parser.add_argument('-dataset', '--dataset', type=str, default="X"   , help="")
parser.add_argument('-nevt'   , '--nevt'   , type=str, default=-1    , help="")
parser.add_argument('-json'   , '--json'   , type=str, default=None  , help="")

options  = parser.parse_args()

def inputfile(nanofile):
   tested   = False
   forceaaa = False
   pfn=os.popen("edmFileUtil -d %s"%(nanofile)).read()
   pfn=re.sub("\n","",pfn)
   print nanofile," -> ",pfn
   if (os.getenv("GLIDECLIENT_Group","") != "overflow" and
       os.getenv("GLIDECLIENT_Group","") != "overflow_conservative" and not
       forceaaa ):
      if not tested:
         print "Testing file open"
         testfile=ROOT.TFile.Open(pfn)
         if testfile and testfile.IsOpen() :
            print "Test OK"
            nanofile=pfn
            testfile.Close()
         else:
            if "root://cms-xrd-global.cern.ch/" not in nanofile:
               nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
            forceaaa=True
      else:
         nanofile = pfn
   else:
       if "root://cms-xrd-global.cern.ch/" not in nanofile:
          nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
   return nanofile

options.infile = inputfile(options.infile)


print "---------------------------"
print " -- options  = ", options
print " -- is MC    = ", options.isMC
print " -- jobNum   = ", options.jobNum
print " -- era      = ", options.era
print " -- in file  = ", options.infile
print " -- dataset  = ", options.dataset
print "---------------------------"

if options.era=="UL2016_preVFP":
   era1 = "2016"
if options.era=="UL2016":
   era1 = "2016"
if options.era=="UL2017":
   era1= "2017"
if options.era=="UL2018":
   era1 = "2018"


xsection = 1.0
#nevents = 1
if options.isMC:
   condtag_ = "NANOAODSIM"
   if options.dataset == "X":
      options.dataset = options.infile
      options.dataset = options.dataset.split('/store')[1].split("/")
      condtag_ = options.dataset[5]
      options.dataset = options.dataset[3]
      runperiod_ = "None"
   print "[check] condtag_ == ", condtag_
   print "[check] dataset  == ", options.dataset
else:
   if options.dataset == "X":
      print "[check] ",options.dataset
      options.dataset = options.infile
      options.dataset = options.dataset.split('/store')[1].split("/")
      condtag_ = options.dataset[2]
      print "[check] condtag_ == ", condtag_
      print "[check] era == ", options.era
      print "[check] dataset: ",options.dataset[3]
      options.dataset = options.dataset[3]
      runperiod_ = condtag_.split(era1)[1][:1]
   else:
      options.dataset = options.dataset.split("/")
      condtag_ = options.dataset[2]
      options.dataset = options.dataset[1]
      runperiod_ = condtag_.split(era1)[1][:1]
   print "[check] condtag_ == ", condtag_
   print "[check] dataset  == ", options.dataset

dataset = options.dataset 

if options.isMC:
#   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:
       MC_xsecs = yaml.safe_load(file)
   try:
       xsection *= MC_xsecs[dataset]["xsec"]
       xsection *= MC_xsecs[dataset]["kr"]
       xsection *= MC_xsecs[dataset]["br"]
       print("The dataset name: ", dataset)
       print("The xsection: ", xsection)
   except:
       print("WARNING: I did not find the xsection for that MC sample. Check the dataset name and the relevant yaml file")

pre_selection = "Sum$(Muon_pt>40 && abs(Muon_eta)<2.5) > 1 "

if float(options.nevt) > 0:
   print " passing this cut and : ", options.nevt
   pre_selection += ' && (Entry$ < {})'.format(options.nevt)


if "QCD" in dataset:
    modules_era = [ GenWeightProducer( isMC = options.isMC, xsec = xsection, dopdf =  False, do_xsecscale = True) ]
else:
    modules_era = [ GenWeightProducer( isMC = options.isMC, xsec = xsection, dopdf =  True, do_xsecscale = True) ]


pro_syst = [ "jer"]
ext_syst = [ "puWeight"]#, "PDF", "MuonSF", "ElecronSF", "EWK", "nvtxWeight","TriggerSFWeight","btagEventWeight", "QCDScale0w", "QCDScale1w", "QCDScale2w"]

if options.era=="UL2017":
   for ix in ['jesAbsolute', 'jesAbsolute_2017', 'jesBBEC1', 'jesBBEC1_2017', 'jesEC2', 'jesEC2_2017', 'jesFlavorQCD', 'jesHF', 'jesHF_2017', 'jesRelativeBal', 'jesRelativeSample_2017']:
      pro_syst.append(ix)

jetmetCorrector = createJMECorrector(isMC = options.isMC, dataYear = options.era, runPeriod= runperiod_ , jesUncert="Merged", applyHEMfix = (True if options.era=="UL2018" else False)) # add METfixEE2017
if options.isMC:
   print "sample : ", options.dataset, " candtag : ", condtag_

   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2_LQ.yaml"))
   except yaml.YAMLError as exc:
      print(exc)
   if options.era=="UL2017":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2017.SingleMuon", "") + ")"

   if options.era=="UL2017":
      modules_era.append(puAutoWeight_UL2017())
      modules_era.append(jetmetCorrector())
#      modules_era.append(btagSFProducer("UL2017", "deepcsv",['shape_corr']))
      modules_era.append(muonScaleResUL2017())
      modules_era.append(lepSF_UL2017_HighPt())

   modules_era.append(LeptoquarkProducer(isMC=options.isMC, era=str(options.era), do_syst=0, syst_var=''))


   if options.era=="UL2017":
      modules_era.append(MuonRecoSF_UL2017())
      modules_era.append(TriggerSF_UL2017());

   # for shift-based systematics
   if options.doSyst:
      for sys in pro_syst:
         for var in ["Up", "Down"]:
	    # if "jesTotal" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
	    # if "jer" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
            modules_era.append(LeptoquarkProducer(options.isMC, str(options.era), do_syst=options.doSyst, syst_var=sys+var))
else:
   if options.era=="UL2016_preVFP":
      modules_era.append(muonScaleResUL2016_preVPF())
   if options.era=="UL2016":
      modules_era.append(muonScaleResUL2016())
   if options.era=="UL2017":
      modules_era.append(muonScaleResUL2017())
   if options.era=="UL2018":
      modules_era.append(muonScaleResUL2018())

   print "sample : ", options.dataset, " candtag : ", condtag_

   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2_LQ.yaml"))
   except yaml.YAMLError as exc:
      print(exc)
   if options.era=="UL2017":
      if 'Run2017B' in condtag_:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017.RunB", "") + ")"
      else:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017.SingleMuon", "") + ")"

   print " -- era : ",

   runPeriod = condtag_.split(era1)[1][:1]
   jmeCorrections = createJMECorrector(isMC=False, 
                                       dataYear=options.era, 
                                       runPeriod=runPeriod, 
                                       jesUncert="Total")   
   modules_era.append(jetmetCorrector())

   modules_era.append(LeptoquarkProducer(isMC=options.isMC, era=str(options.era), do_syst=0, syst_var=''))
   if options.era=="UL2017":
       options.json = "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
   print "---- JSON used is : ", options.json


for i in modules_era:
   print "modules : ", i

print "Selection : ", pre_selection

p = PostProcessor(
   ".", [options.infile],
   cut=pre_selection,
   branchsel="keep_and_drop.txt",
   outputbranchsel="keep_and_drop_post.txt",
   haddFileName="tree_%s.root" % str(options.jobNum),
   modules=modules_era,
   provenance=True,
   noOut=False,
   fwkJobReport=True,
   jsonInput=options.json
)
p.run()
