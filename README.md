# LQ NanoAOD
Instructions for LQ analysis using NanoAOD based on HH analysis


## Quick start

```bash
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv

mkdir PhysicsTools/
git clone https://github.com/jlidrych/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone --branch LQanalysis-UL https://github.com/jlidrych/NUNanoAOD.git PhysicsTools/MonoZ

scram b -j 10

## Example: running over a signal MC file
To run locally, edit condor_Run2_proc_LQ.py so that the xsections yaml is from the data directory.
Uncomment:
   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
Comment out: 
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:

cd $CMSSW_BASE/src/PhysicsTools/MonoZ/condor/
voms-proxy-init -voms cms --valid 168:00

python condor_Run2_proc_LQ.py --isMC=1 --doSyst=0 --era=UL2017 --nevt=1000 --infile=root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/130000/5E245D45-B8B0-AE44-B783-FB3C1B0624C5.root

## Example: submitting jobs to condor - UL 2017
To run on condor, edit condor_Run2_proc_LQ.py 
Comment out:
   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
Uncomment out: 
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:
Can also turn off systematics by changing --doSyst default value

python runLQ.py --isMC=1 --era=UL2017 --tag=[tagname] --input=../data/list_2017_MC_LQ.txt
python runLQ.py --isMC=0 --era=UL2017 --tag=[tagname] --input=../data/list_2017_Data_Mu_NANOv9.txt

```
