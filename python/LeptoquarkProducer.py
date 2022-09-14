import ROOT
import sys, os
import numpy as np
import math
from importlib import import_module
import itertools
from copy import deepcopy
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import PhysicsTools.NanoAODTools.postprocessing.tools as tk

ROOT.PyConfig.IgnoreCommandLineOptions = True


class LeptoquarkProducer(Module):
    def __init__(self, isMC, era, do_syst=False, syst_var=''):
        self.isMC = isMC
        self.era = era
        self.do_syst = do_syst
        self.syst_var = syst_var
        self.zmass = 91.1873
        self.Hmass = 125.10
	if self.syst_var !='':
          self.syst_suffix = '_sys_' + self.syst_var if self.do_syst else ''
	else:
          self.syst_suffix = syst_var

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("met_pt{}".format(self.syst_suffix), "F")
        self.out.branch("met_phi{}".format(self.syst_suffix), "F")
        self.out.branch("lep_category{}".format(self.syst_suffix), "I") 
        # 1 = dimuon channel, 2 = dielectron channel
        self.out.branch("event_category{}".format(self.syst_suffix), "I") 
        # for ABCD method
        # 1 = A OS+ISO
        # 2 = B OS+non-ISO
        # 3 = C SS+ISO
        # 4 = D SS+non-ISO

        self.out.branch("ngood_muons{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_electrons{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_taus{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_jets{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjets{}".format(self.syst_suffix), "I")

        self.out.branch("Z_pt{}".format(self.syst_suffix), "F")
        self.out.branch("Z_eta{}".format(self.syst_suffix), "F")
        self.out.branch("Z_phi{}".format(self.syst_suffix), "F")
        self.out.branch("Z_mass{}".format(self.syst_suffix), "F")

        self.out.branch("lead_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("lead_lep_p{}".format(self.syst_suffix), "F")
        self.out.branch("lead_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("lead_lep_phi{}".format(self.syst_suffix), "F")

        self.out.branch("trail_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trail_lep_p{}".format(self.syst_suffix), "F")
        self.out.branch("trail_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trail_lep_phi{}".format(self.syst_suffix), "F")

        self.out.branch("lead_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("lead_jet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("lead_jet_phi{}".format(self.syst_suffix), "F")

        self.out.branch("lead_bjet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("lead_bjet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("lead_bjet_phi{}".format(self.syst_suffix), "F")

        self.out.branch("met_filter{}".format(self.syst_suffix), "I")

        self.out.branch("jetHT{}".format(self.syst_suffix), "F")
        self.out.branch("ST{}".format(self.syst_suffix), "F")

        if self.isMC and len(self.syst_suffix)==0:
            self.out.branch("w_muon_SF{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFUp{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFDown{}".format(self.syst_suffix), "F")

        if self.isMC:
            self.out.branch("w_btag_SF{}".format(self.syst_suffix),"F")     

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def electron_id(self, electron, wp):
        pass_id = 0
        if (wp == "80"):
            try:
                pass_id = electron.mvaFall17V2Iso_WP80
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WP80
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WP80
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "
                        
            return pass_id
        elif (wp == "90"):
            try:
                pass_id = electron.mvaFall17V2Iso_WP90
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WP90
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WP90
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id
        elif (wp == "WPL"):
            try:
                pass_id = electron.mvaFall17V2Iso_WPL
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WPL
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WPL
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "
            return pass_id

    def btagDeepJet_id(self, wp):
        # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        if (self.era == "2016" and wp == "loose"):
            return 0.0614
        elif (self.era == "2016" and wp == "medium"):
            return 0.3093
        elif (self.era == "2016" and wp == "tight"):
            return 0.7221
        elif (self.era == "2017" and wp == "loose"):
            return 0.0521
        elif (self.era == "2017" and wp == "medium"):
            return 0.3033
        elif (self.era == "2017" and wp == "tight"):
            return 0.7489
        elif (self.era == "UL2017" and wp == "loose"):
            return 0.0532
        elif (self.era == "UL2017" and wp == "medium"):
            return 0.3040
        elif (self.era == "UL2017" and wp == "tight"):
            return 0.7476
        elif (self.era == "UL2018" and wp == "loose"):
            return 0.0490
        elif (self.era == "UL2018" and wp == "medium"):
            return 0.2783
        elif (self.era == "UL2018" and wp == "tight"):
            return 0.7100

    def btagDeepCSV_id(self, wp):
        # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        if (self.era == "UL2016" and wp == "medium"):
            return 0.5847
        elif (self.era == "UL2016_preVFP" and wp == "medium"):
            return 0.6001
        elif (self.era == "UL2017" and wp == "medium"):
            return 0.4506
        elif (self.era == "UL2018" and wp == "medium"):
            return 0.4168


    def met_filter(self, flag, filter_mask=True):
        if self.era =="UL2018" or self.era =="UL2017":
            filter_mask = flag.ecalBadCalibFilter
        return filter_mask and (
              (flag.HBHENoiseFilter)
           and (flag.HBHENoiseIsoFilter)
           and (flag.EcalDeadCellTriggerPrimitiveFilter)
           and (flag.goodVertices)
           and (flag.eeBadScFilter)
           and (flag.globalSuperTightHalo2016Filter)
           and (flag.BadPFMuonFilter)
           and (flag.BadPFMuonDzFilter)
        )

    def duplicate_removal(self):
        """
        For data, same event could come from different datasets
        FIXME: need to be implemented check the source from
        the old MonoZ code
        https://github.com/NEUAnalyses/monoZ_Analysis/blob/master/src/MonoZSelector.cc#L463
        """
        pass

    def analyze(self, event):
        """
        process event, return True (go to next module)
        or False (fail, go to next event)
        """
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        jets = list(Collection(event, "Jet"))
        taus = list(Collection(event, "Tau"))
        flag = Object(event, "Flag")
        met = Object(event, "MET")

        # in case of systematic take the shifted values are default
        # For the central values, need to include jetMetTool all the time
        # Jet systematics
        if self.syst_var == "":
            syst_var = "nom"
        else:
            syst_var = self.syst_var
        # checking something
        try:
            var_jet_pts = getattr(event,  "Jet_pt_{}".format(syst_var), None)
            if var_jet_pts:
                for i,jet in enumerate(jets):
                    jet.pt = var_jet_pts[i]
            else:
                print 'WARNING: jet pts with variation {}'
                'not available, using the nominal value'.format(syst_var)
        except:
            var_jet_pts = getattr(event,  "Jet_pt_nom", None)
            for i,jet in enumerate(jets):
                jet.pt = var_jet_pts[i]
        try:
            var_met_pt  = getattr(event,  "MET_pt_{}".format(syst_var), None)
            var_met_phi = getattr(event, "MET_phi_{}".format(syst_var), None)
            if var_met_pt:
                met.pt = var_met_pt
            else:
                print 'WARNING: MET pt with variation '
                '{} not available, using the nominal value'.format(syst_var)
            if var_met_phi:
                met.phi = var_met_phi
            else:
                print 'WARNING: MET phi with variation {}'
                'not available, using the nominal value'.format(syst_var)
        except:
            var_met_pt  = getattr(event,  "MET_T1_pt", None)
            var_met_phi = getattr(event, "MET_T1_phi", None)
            if var_met_pt:
                met.pt = var_met_pt
            if var_met_phi:
                met.phi = var_met_phi

        met_p4 = ROOT.TLorentzVector()
        met_p4.SetPtEtaPhiM(met.pt,0.0,met.phi, 0.0)


        # Electrons Energy
        if "ElectronEn" in self.syst_var:
            (met_px, met_py) = ( met.pt*np.cos(met.phi), met.pt*np.sin(met.phi) )
            if "Up" in self.syst_var:
                for i, elec in enumerate(electrons):
                    met_px = met_px + (elec.energyErr)*np.cos(elec.phi)/math.cosh(elec.eta)
                    met_py = met_py + (elec.energyErr)*np.sin(elec.phi)/math.cosh(elec.eta)
                    elec.pt = elec.pt + elec.energyErr/math.cosh(elec.eta)
            else:
                for i, elec in enumerate(electrons):
                    met_px = met_px - (elec.energyErr)*np.cos(elec.phi)/math.cosh(elec.eta)
                    met_py = met_py - (elec.energyErr)*np.sin(elec.phi)/math.cosh(elec.eta)
                    elec.pt = elec.pt - elec.energyErr/math.cosh(elec.eta)
            met.pt  = math.sqrt(met_px**2 + met_py**2)
            met.phi = math.atan2(met_py, met_px)

        # Muons Energy
        #if self.isMC:
        muons_pts = getattr(event, "Muon_corrected_pt")
        for i, muon in enumerate(muons):
            muon.pt = muons_pts[i]

        if "MuonEn" in self.syst_var:
            (met_px, met_py) = ( met.pt*np.cos(met.phi), met.pt*np.sin(met.phi) )
            if "Up" in self.syst_var:
                muons_pts = getattr(event, "Muon_correctedUp_pt")
                for i, muon in enumerate(muons):
                    met_px = met_px - (muons_pts[i] - muon.pt)*np.cos(muon.phi)
                    met_py = met_py - (muons_pts[i] - muon.pt)*np.sin(muon.phi)
                    muon.pt = muons_pts[i]
            else:
                muons_pts = getattr(event, "Muon_correctedDown_pt")
                for i, muon in enumerate(muons):
                    met_px =met_px - (muons_pts[i] - muon.pt)*np.cos(muon.phi)
                    met_py =met_py - (muons_pts[i] - muon.pt)*np.sin(muon.phi)
                    muon.pt = muons_pts[i]
            met.pt  = math.sqrt(met_px**2 + met_py**2)
            met.phi = math.atan2(met_py, met_px)            


        # filling and contructing the event categorisation
        self.out.fillBranch("met_pt{}".format(self.syst_suffix), met.pt)
        self.out.fillBranch("met_phi{}".format(self.syst_suffix), met.phi)

        pass_met_filter = self.met_filter(flag, True)
        self.out.fillBranch("met_filter{}".format(self.syst_suffix), pass_met_filter)

        # count electrons and muons
        good_leptons = []
        good_muons = []
        good_electrons = []

        muons.sort(key=lambda muon: muon.pt, reverse=True)
        electrons.sort(key=lambda el: el.pt, reverse=True)
        # Choose loose/medium-quality e/mu for event categorization

        for idx,mu in enumerate(muons):
            pass_ips = abs(mu.dxy) < 0.2 and abs(mu.dz) < 0.5
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= 52
            pass_ids = (mu.highPtId > 0) and mu.isGlobal #and mu.triggerIdLoose
            if pass_fid and pass_ids and pass_ips:
                good_muons.append(mu)

        for idy,el in enumerate(electrons):
            id_CB = el.cutBased
            # changing to MVA based ID :
            if el.pt >= 30 and abs(el.eta) <= 2.5 and self.electron_id(el, "WP90"):
                good_electrons.append(el)

        # let sort the muons in pt
        good_muons.sort(key=lambda x: x.pt, reverse=True)
        good_electrons.sort(key=lambda x: x.pt, reverse=True)

        good_leptons = good_muons + good_electrons
        ngood_muons = len(good_muons)
        ngood_electrons = len(good_electrons)
        ngood_leptons = len(good_leptons)

        self.out.fillBranch("ngood_muons{}".format(self.syst_suffix), ngood_muons)
        self.out.fillBranch("ngood_electrons{}".format(self.syst_suffix), ngood_electrons)

        if ngood_muons >=2:
            self.out.fillBranch("lep_category{}".format(self.syst_suffix), 1)
            if (good_muons[0].pdgId * good_muons[1].pdgId == -13*13):
                if(good_muons[0].triggerIdLoose and good_muons[1].triggerIdLoose):
                    self.out.fillBranch("event_category{}".format(self.syst_suffix),1)
                else:
                    self.out.fillBranch("event_category{}".format(self.syst_suffix),2)
            if (good_muons[0].pdgId * good_muons[1].pdgId == 13*13):
                if(good_muons[0].triggerIdLoose and good_muons[1].triggerIdLoose):
                    self.out.fillBranch("event_category{}".format(self.syst_suffix),3)
                else:
                    self.out.fillBranch("event_category{}".format(self.syst_suffix),4)

        if self.isMC:
            w_muon_SF     = 1.0
            w_muon_SFUp   = 1.0
            w_muon_SFDown = 1.0
            if ngood_leptons >= 2:
                if abs(good_leptons[0].pdgId) == 13:
                    w_muon_SF     *=  good_leptons[0].SF
                    w_muon_SFUp   *= (good_leptons[0].SF + good_leptons[0].SFErr)
                    w_muon_SFDown *= (good_leptons[0].SF - good_leptons[0].SFErr)
                if abs(good_leptons[1].pdgId) == 13:
                    w_muon_SF     *=  good_leptons[1].SF
                    w_muon_SFUp   *= (good_leptons[1].SF + good_leptons[1].SFErr)
                    w_muon_SFDown *= (good_leptons[1].SF - good_leptons[1].SFErr)
            self.out.fillBranch("w_muon_SF"        , w_muon_SF        )
            self.out.fillBranch("w_muon_SFUp"      , w_muon_SFUp      )
            self.out.fillBranch("w_muon_SFDown"    , w_muon_SFDown    )

        # process jet
        good_jets = []
        good_bjets = []
        btagSF = 1
        jetHT = 0

        for jet in jets:
            if jet.pt < 50.0 or abs(jet.eta) > 2.4:
                continue
            if jet.jetId < 2:
                continue
            if tk.closest(jet, good_leptons)[1] < 0.4:
                continue
            good_jets.append(jet)
            jetHT += jet.pt

            if self.isMC:
                try:
                    btagSF *=jet.btagSF_deepcsv_shape
                except:
                    btagSF *=1
            # Count b-tag with loose WP DeepCSV
            # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            if jet.btagDeepB > self.btagDeepCSV_id("medium"):
                good_bjets.append(jet)

        good_jets.sort(key=lambda jet: jet.pt, reverse=True)
        good_bjets.sort(key=lambda jet: jet.btagDeepB, reverse=True)

       #We will remove jets later so better count them now
        num_jets = len(good_jets)
        num_bjets = len(good_bjets)

        self.out.fillBranch("jetHT{}".format(self.syst_suffix), jetHT)
        self.out.fillBranch("ngood_jets{}".format(self.syst_suffix), num_jets)
        self.out.fillBranch("ngood_bjets{}".format(self.syst_suffix),num_bjets)
        if self.isMC:
            self.out.fillBranch("w_btag_SF{}".format(self.syst_suffix),btagSF)

        if num_bjets > 0:
            self.out.fillBranch("lead_bjet_pt{}".format(self.syst_suffix),good_bjets[0].pt)
            self.out.fillBranch("lead_bjet_eta{}".format(self.syst_suffix),good_bjets[0].eta)
            self.out.fillBranch("lead_bjet_phi{}".format(self.syst_suffix),good_bjets[0].phi)
        else:
            self.out.fillBranch("lead_bjet_pt{}".format(self.syst_suffix),-99)
            self.out.fillBranch("lead_bjet_eta{}".format(self.syst_suffix),-99)
            self.out.fillBranch("lead_bjet_phi{}".format(self.syst_suffix),-99)            

        if num_jets > 0:
            self.out.fillBranch("lead_jet_pt{}".format(self.syst_suffix),good_jets[0].pt)
            self.out.fillBranch("lead_jet_eta{}".format(self.syst_suffix),good_jets[0].eta)
            self.out.fillBranch("lead_jet_phi{}".format(self.syst_suffix),good_jets[0].phi)
        else:
            self.out.fillBranch("lead_jet_pt{}".format(self.syst_suffix),-99)
            self.out.fillBranch("lead_jet_eta{}".format(self.syst_suffix),-99)
            self.out.fillBranch("lead_jet_phi{}".format(self.syst_suffix),-99)

        Z_cand = ROOT.TLorentzVector()
        STcalc = 0
        STcalc += jetHT
        
        if ngood_muons >= 2:
            Z_cand = good_muons[0].p4() + good_muons[1].p4()
            STcalc += good_muons[0].pt
            STcalc += good_muons[1].pt
            self.out.fillBranch("Z_pt{}".format(self.syst_suffix),Z_cand.Pt())
            self.out.fillBranch("Z_eta{}".format(self.syst_suffix),Z_cand.Eta())
            self.out.fillBranch("Z_phi{}".format(self.syst_suffix),Z_cand.Phi())
            self.out.fillBranch("Z_mass{}".format(self.syst_suffix),Z_cand.M())
        
            self.out.fillBranch("lead_lep_pt{}".format(self.syst_suffix), good_muons[0].pt)
            self.out.fillBranch("lead_lep_p{}".format(self.syst_suffix), good_muons[0].p4().P())
            self.out.fillBranch("lead_lep_eta{}".format(self.syst_suffix), good_muons[0].eta)
            self.out.fillBranch("lead_lep_phi{}".format(self.syst_suffix), good_muons[0].phi)
            
            self.out.fillBranch("trail_lep_pt{}".format(self.syst_suffix), good_muons[1].pt)
            self.out.fillBranch("trail_lep_p{}".format(self.syst_suffix), good_muons[1].p4().P())
            self.out.fillBranch("trail_lep_eta{}".format(self.syst_suffix), good_muons[1].eta)
            self.out.fillBranch("trail_lep_phi{}".format(self.syst_suffix), good_muons[1].phi)
            self.out.fillBranch("ST{}".format(self.syst_suffix), STcalc)


        # process taus
        had_taus = []
        for tau in taus:
            if tk.closest(tau, good_leptons)[1] < 0.4:
                continue
            # only hadronic tau decay
            if tau.decayMode != 5:
                continue
            if tau.pt > 18 and abs(tau.eta) <= 2.3:
                had_taus.append(tau)

        self.out.fillBranch("ngood_taus{}".format(self.syst_suffix), len(had_taus))

        if( (ngood_muons >= 2) and (num_jets > 0) ):
            return True
        else:
            return False
