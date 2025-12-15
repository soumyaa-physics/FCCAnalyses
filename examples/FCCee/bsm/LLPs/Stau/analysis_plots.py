'''
Plotting stage of the Stau analysis
'''
import ROOT

intLumi        = 2.05e8
scaleSig       = -1.
scaleBkg       = -1.

ana_tex        = ''
delphesVersion = ''
energy         = 240.
collider       = 'FCC-ee'

inputDir       = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final/"
outdir         = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final/plots"


formats  = ['pdf']
yaxis    = ['log']
stacksig = ['nostack']

splitLeg = True

variables = [
                # Gen-level stau
            "GenStau",
            "n_GenStau",
            "GenStau_observed_lifetime_xyz",
            "GenStau_vx",
            "GenStau_vy",
            "GenStau_vz",
            "GenStau_Lxy",
            "GenStau_Lxyz",

            # Stau daughters
            "StauDaughters",

            # Tau kinematics
            "GenTau_e",
            "GenTau_pt",
            "GenTau_eta",
            "GenTau_phi",
            "GenTau_px",
            "GenTau_py",
            "GenTau_pz",

            # Tau vertex
            "GenTau_vx",
            "GenTau_vy",
            "GenTau_vz",

            # Final-state pions
            "FSGenPion",
            "n_FSGenPion",
            "FSGenPion_e",
            "FSGenPion_p",
            "FSGenPion_pt",
            "FSGenPion_px",
            "FSGenPion_py",
            "FSGenPion_pz",
            "FSGenPion_eta",
            "FSGenPion_theta",
            "FSGenPion_phi",
            "FSGenPion_charge",
            "FSGenPion_vertex_x",
            "FSGenPion_vertex_y",
            "FSGenPion_vertex_z",

            # Final-state muons
            "FSGenMuon",
            "n_FSGenMuon",
            "FSGenMuon_e",
            "FSGenMuon_p",
            "FSGenMuon_pt",
            "FSGenMuon_px",
            "FSGenMuon_py",
            "FSGenMuon_pz",
            "FSGenMuon_eta",
            "FSGenMuon_theta",
            "FSGenMuon_phi",
            "FSGenMuon_charge",
            "FSGenMuon_vertex_x",
            "FSGenMuon_vertex_y",
            "FSGenMuon_vertex_z",

            # Final-state electrons
            "FSGenElectron",
            "n_FSGenElectron",
            "FSGenElectron_e",
            "FSGenElectron_p",
            "FSGenElectron_pt",
            "FSGenElectron_px",
            "FSGenElectron_py",
            "FSGenElectron_pz",
            "FSGenElectron_eta",
            "FSGenElectron_theta",
            "FSGenElectron_phi",
            "FSGenElectron_charge",
            "FSGenElectron_vertex_x",
            "FSGenElectron_vertex_y",
            "FSGenElectron_vertex_z",

            # Reco Jets
            "n_RecoJets",
            "RecoJet_e",
            "RecoJet_p",
            "RecoJet_pt",
            "RecoJet_px",
            "RecoJet_py",
            "RecoJet_pz",
            "RecoJet_eta",
            "RecoJet_theta",
            "RecoJet_phi",
            "RecoJet_charge",

            # Reco Electrons
            "RecoElectrons",
            "n_RecoElectrons",
            "RecoElectron_e",
            "RecoElectron_p",
            "RecoElectron_pt",
            "RecoElectron_px",
            "RecoElectron_py",
            "RecoElectron_pz",
            "RecoElectron_eta",
            "RecoElectron_theta",
            "RecoElectron_phi",
            "RecoElectron_charge",

            # Reco Muons
            "RecoMuons",
            "n_RecoMuons",
            "RecoMuons_e",
            "RecoMuons_p",
            "RecoMuons_pt",
            "RecoMuons_px",
            "RecoMuons_py",
            "RecoMuons_pz",
            "RecoMuons_eta",
            "RecoMuons_theta",
            "RecoMuons_phi",
            "RecoMuons_charge",

            # MET
            "RecoMissingEnergy_e",
            "RecoMissingEnergy_p",
            "RecoMissingEnergy_pt",
            "RecoMissingEnergy_eta",
            "RecoMissingEnergy_theta",
            "RecoMissingEnergy_phi",
]

selections = {}
selections['Stau'] = [
    "selNone",
]

extralabel = {}
extralabel['selNone'] = "Before selection"

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
colors['FCC_100stau_240com_000'] = ROOT.TColor.GetColor(color_wheel[0]) 

# Background
# colors['background_3photons_cme_91p188'] = ROOT.TColor.GetColor(color_wheel[13])

plots = {}
plots['stau'] = {
    'signal': {
        'FCC_100stau_240com_000': ['FCC_100stau_240com_000'],
    },
}

legend = {}
legend['FCC_100stau_240com_000'] = 'm_{stau} = 100 GeV, ctau_0 = 10 mm'

