'''
Plotting stage of the Stau analysis
'''
import os
import ROOT
# use this : source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
intLumi        = 2.05e8
###If scaleSig=0 or scaleBack=0, we don't apply any additional scaling, on top of the normalization to cross section and integrated luminosity, as defined in finalSel.py
###If scaleSig or scaleBack is not defined, plots will be normalized to 1
scaleSig       = 1
scaleBkg       = 1.

ana_tex        = ''
delphesVersion = ''
energy         = 240.
collider       = 'FCC-ee'

# Input/output directories
inputDir       = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final/"
outdir         = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final/plots"
os.makedirs(outdir, exist_ok=True)

formats  = ['pdf']
yaxis    = ['log', 'lin']
stacksig = ['nostack']

splitLeg = True

# Variables to plots
variables = [
            # Gen-level stau
            # "n_GenStau",
            # "GenStau_vx",
            # "GenStau_vy",
            # "GenStau_vz",
            # "GenStau_Lxy",
            # "GenStau_Lxyz",

            # # Tau kinematics
            # "GenTau_e",
            # "GenTau_pt",
            # "GenTau_eta",
            # "GenTau_phi",
            # "GenTau_px",
            # "GenTau_py",
            # "GenTau_pz",

            # # Tau vertex
            # "GenTau_vx",
            # "GenTau_vy",
            # "GenTau_vz",

            # # Final-state muons
            # "n_FSGenMuon",
            # "FSGenMuon_e",
            # "FSGenMuon_pt",
            # "FSGenMuon_px",
            # "FSGenMuon_py",
            # "FSGenMuon_pz",
            # "FSGenMuon_eta",
            # "FSGenMuon_phi",
            # "FSGenMuon_charge",
            # "FSGenMuon_vx",
            # "FSGenMuon_vy",
            # "FSGenMuon_vz",

            # # Final-state electrons
            # "n_FSGenElectron",
            # "FSGenElectron_e",
            # "FSGenElectron_p",
            # "FSGenElectron_pt",
            # "FSGenElectron_px",
            # "FSGenElectron_py",
            # "FSGenElectron_pz",
            # "FSGenElectron_eta",
            # "FSGenElectron_phi",
            # "FSGenElectron_charge",
            # "FSGenElectron_vx",
            # "FSGenElectron_vy",
            # "FSGenElectron_vz",

            # Reco Jets
            # "n_RecoJets",
            # "RecoJet_e",
            # "RecoJet_pt",
            # "RecoJet_px",
            # "RecoJet_py",
            # "RecoJet_pz",
            # "RecoJet_eta",
            # "RecoJet_phi",
            # "RecoJet_charge",

            # # Reco Electrons
            # "n_RecoElectrons",
            # "RecoElectrons_e",
            # "RecoElectrons_pt",
            # "RecoElectrons_px",
            # "RecoElectrons_py",
            # "RecoElectrons_pz",
            # "RecoElectrons_eta",
            # "RecoElectrons_phi",
            # "RecoElectrons_charge",

            # # Reco Muons
            # "n_RecoMuons",
            # "RecoMuons_e",
            # "RecoMuons_pt",
            # "RecoMuons_px",
            # "RecoMuons_py",
            # "RecoMuons_pz",
            # "RecoMuons_eta",
            # "RecoMuons_phi",
            # "RecoMuons_charge",

            # # MET
            # "RecoMissingEnergy_e",
            # "RecoMissingEnergy_pt",
            # "RecoMissingEnergy_eta",
            # "RecoMissingEnergy_phi",

            "n_seltracks_DVs",
            "n_trks_seltracks_DVs",
            "invMass_seltracks_DVs",
            "DV_evt_seltracks_chi2",
            "DV_evt_seltracks_normchi2",

            "Reco_DVs_merged_Lxy",
            "Reco_DVs_merged_Lxyz",
            "Reco_ee_energy",
            "Reco_mumu_energy",
            "Reco_ee_invMass",
            "Reco_mumu_invMass",
            "Reco_ee_px",
            "Reco_mumu_px",
            "Reco_ee_py",
            "Reco_mumu_py",
            "Reco_ee_pz",
            "Reco_mumu_pz",
            "muon_electron_overlap",
]
# Define selections, labels, colors, plots, legends
selections = {}
selections[''] = [
    "selNone",
    "sel",
]

extralabel = {}
extralabel['selNone'] = "Before selection"
extralabel["sel"] = "Missing Energy > 18 GeV"

color_wheel = [
    # colors from DESY color guide
    "#EB5E2D",  # Red
    "#8CBE23",  # Light green
    "#00B1AA",  # Turquoise
    "#D2006E",  # Magenta
    "#917DB9",  # Violet
    "#C3B700",  # Olive
    "#FAC800",  # Yellow
    "#B92D41",  # Dark red
    "#00A64B",  # Green
    "#006987",  # Petrol
    "#8C3C5B",  # Aubergine
    "#504F8F",  # Purple
    "#828F2B",  # Dark olive
    "#004A6F"   # Dark blue
]

colors = {}

# Signal
colors['FCCee_100_stau_10mm_ctau_ecm_240'] = ROOT.TColor.GetColor(color_wheel[0]) 
colors['FCCee_110_stau_10mm_ctau_ecm_240'] = ROOT.TColor.GetColor(color_wheel[1])

# Background
colors['p8_ee_WW_ecm240'] = ROOT.TColor.GetColor(color_wheel[1]) 
colors['p8_ee_ZZ_ecm240'] = ROOT.TColor.GetColor(color_wheel[2]) 
colors['mgp8_ee_zh_ecm240_hbb'] = ROOT.TColor.GetColor(color_wheel[3])
colors['wzp6_ee_nuenueH_Htautau_ecm240'] = ROOT.TColor.GetColor(color_wheel[4])
colors['wzp6_ee_bbH_Htautau_ecm240'] = ROOT.TColor.GetColor(color_wheel[5])

# Plot and legend structure
plots = {}
plots[''] = {
    'signal': {
        'FCCee_100_stau_10mm_ctau_ecm_240': ['FCCee_100_stau_10mm_ctau_ecm_240'],
        'FCCee_110_stau_10mm_ctau_ecm_240': ['FCCee_110_stau_10mm_ctau_ecm_240'],
    },
    'backgrounds': {
        # 'p8_ee_WW_mumu_ecm240': ['p8_ee_WW_mumu_ecm240'],
        # 'p8_ee_WW_ee_ecm240': ['p8_ee_WW_ee_ecm240'],
        'p8_ee_WW_ecm240': ['p8_ee_WW_ecm240'],
        'p8_ee_ZZ_ecm240': ['p8_ee_ZZ_ecm240'],
        'mgp8_ee_zh_ecm240_hbb': ['mgp8_ee_zh_ecm240_ßhbb'],
        'wzp6_ee_nuenueH_Htautau_ecm240': ['wzp6_ee_nuenueH_Htautau_ecm240'],
        'wzp6_ee_bbH_Htautau_ecm240': ['wzp6_ee_bbH_Htautau_ecm240'],
    }
}

legend = {}
legend['FCCee_100_stau_10mm_ctau_ecm_240'] = 'm_{stau} = 100 GeV, ctau_0 = 10 mm'
legend['FCCee_110_stau_10mm_ctau_ecm_240'] = 'm_{stau} = 110 GeV, ctau_0 = 10 mm'
# BACKRGOUNDS
legend['p8_ee_WW_ecm240'] = 'e^{+}e^{-} #rightarrow #mu^{+}#mu^{-}'
legend['p8_ee_ZZ_ecm240'] = 'e^{+}e^{-} #rightarrow #mu^{+}#mu^{-}'
legend['mgp8_ee_zh_ecm240_hbb'] = 'e^{+}e^{-} #rightarrow ZH, H #rightarrow b#bar{b}'
legend['wzp6_ee_nuenueH_Htautau_ecm240'] = 'e^{+}e^{-} #rightarrow H #rightarrow #tau^{+}#tau^{-}'
legend['wzp6_ee_bbH_Htautau_ecm240'] = 'e^{+}e^{-} #rightarrow H #rightarrow #tau^{+}#tau^{-}'
