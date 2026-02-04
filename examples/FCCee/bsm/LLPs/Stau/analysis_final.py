'''
Final stage of the stau analysis
'''

# Input/output directories
inputDir  = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/output_stage1"
outputDir = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final/"

# List of datasets used in the analysis
processList = {
        #######################################################
        #             CME: 240 GeV (ZH)- 10mm                 #
        #######################################################
        'FCCee_100_stau_10mm_ctau_ecm_240'  : {'fraction': 1.0},
        'FCCee_110_stau_10mm_ctau_ecm_240' : {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH)- 20cm                 #
        #######################################################
        'FCCee_110_stau_20cm_ctau_ecm_240'  : {'fraction': 1.0},

        #######################################################
        #               WINTER 2023                           #
        #######################################################
        # 'p8_ee_WW_ecm240': {'fraction': 1.0,'chunks':100},
        # 'p8_ee_ZZ_ecm240': {'fraction': 1.0,'chunks':100},
        # 'mgp8_ee_zh_ecm240_hbb': {'fraction': 1.0,'chunks':100},
        # 'wzp6_ee_nuenueH_Htautau_ecm240': {'fraction': 1.0,'chunks':100},
        # 'wzp6_ee_bbH_Htautau_ecm240': {'fraction': 1.0,'chunks':100},

        # old samples
        # 'p8_ee_WW_mumu_ecm240': {'fraction': 0.5,'chunks':100},
        # 'p8_ee_WW_ee_ecm240': {'fraction': 0.5,'chunks':100},
        # 'p8_ee_Zbb_ecm91_EvtGen_Bd2KstarTauTau': {'fraction': 1.0,'chunks':100},
        } 

prodTag = "FCCee/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
# procDict = "FCCee_procDict_spring2021_IDEA.json"


# Add samples which are not part of the offical process
procDictAdd = {
        #######################################################
        #             CME: 240 GeV (ZH)- 10mm                 #
        #######################################################
        'FCCee_100_stau_10mm_ctau_ecm_240': {"numberOfEvents": 10000, "sumOfWeights": 10000, "crossSection": 0.07733 ,     "kfactor": 1.0, "matchingEfficiency": 1.0},
        # # 'FCCee_105_stau_10mm_ctau_ecm_240': {"numberOfEvents": 10000, "sumOfWeights": 10000, "crossSection": 0.02048,     "kfactor": 1.0, "matchingEfficiency": 1.0},
        'FCCee_110_stau_10mm_ctau_ecm_240': {"numberOfEvents": 10000, "sumOfWeights": 10000, "crossSection":  0.02923,     "kfactor": 1.0, "matchingEfficiency": 1.0},
        # # 'FCCee_115_stau_10mm_ctau_ecm_240': {"numberOfEvents": 10000, "sumOfWeights": 10000, "crossSection": 0.02048,     "kfactor": 1.0, "matchingEfficiency": 1.0},

        #######################################################
        #             CME: 240 GeV (ZH)- 20cm                 #
        #######################################################
        'FCCee_110_stau_20cm_ctau_ecm_240': {"numberOfEvents": 10000, "sumOfWeights": 10000, "crossSection":  0.02923,     "kfactor": 1.0, "matchingEfficiency": 1.0},

        
        }

intLumi = 2.05e8
#intLumi = 1.92e7
#intLumi = 1.08e7
#intLumi = 2.7e6
doScale = True
# saveMetaData = True

# # Required for plots
# writeOutScales = True
# writeMetaDataToFile = True

# Number of threads to use
# nCPUS = 2

# Whether to produce ROOT TTrees, default is False
doTree = False

# Save cut yields and efficiencies in LaTeX table
saveTabular = True

# Save cut yields and efficiencies in JSON file
saveJSON = True

# Dictionary with the list of cuts. The key is the name of the selection that will be added to the output file
cutList = {

    "selNone": "n_RecoTracks > -1",

    "sel_MET": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 18."
    ),

    # "sel_nEle": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && n_RecoElectrons < 3"
    # ),

    # "sel_nEle_nMu": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && n_RecoElectrons < 3"
    #     " && n_RecoMuons < 3"
    # ),

    # "sel_nEle_nMu_ee": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && n_RecoElectrons < 3"
    #     " && n_RecoMuons < 3"
    #     " && Reco_ee_invMass == -1."
    # ),

    # "sel_final": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && n_RecoElectrons < 3"
    #     " && n_RecoMuons < 3"
    #     " && Reco_ee_invMass == -1."
    #     " && Reco_mumu_invMass == -1."
    # ),
    
    # "sel_semiLepHad": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     "&& ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1) || (n_RecoElectrons == 0 && n_RecoMuons == 0))"
    # ),

    "sel_semiLeptonic": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 18."
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
    ),

    "sel_semiLep_jetcut": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 18."
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
        " && n_RecoJets < 4"
    ),

    # "sel_semiLep_DVcut": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
    #     " && n_RecoJets < 4"
    #     # exactly one DV from selected tracks
    #     " && n_seltracks_DVs == 1"
    #     # with <= 3 tracks
    #     # " && n_trks_seltracks_DVs.size() == 1" # ensure that there is only one DV
        # " && n_trks_seltracks_DVs[0] <= 3" # cut on number of tracks in that one DV
    # )

    "sel_semiLep_ntrks_cut": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 18."
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
        " && n_RecoJets < 4"
        # exactly one DV from selected tracks
        # " && n_seltracks_DVs == 1"
        # with <= 3 tracks
        " && n_trks_seltracks_DVs.size() == 1" # ensure that there is only one DV
        " && n_trks_seltracks_DVs[0] <= 3" # cut on number of tracks in that one DV
    ),


    "sel_semiLep_jet_cut": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 18."
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
        " && n_RecoJets < 4"
        # exactly one DV from selected tracks
        # " && n_seltracks_DVs == 1"
        # with <= 3 tracks
        " && n_trks_seltracks_DVs.size() == 1" # ensure that there is only one DV
        " && n_trks_seltracks_DVs[0] <= 3" # cut on number of tracks in that one DV
        " && RecoJet_e[0] < 120" # cut on leading jet energy to further suppress hadronic background
    ),


    # "fullyHadronic": (
    #     "n_RecoTracks > -1"
    #     " && RecoMissingEnergy_e > 18."
    #     " && n_RecoElectrons == 0"
    #     " && n_RecoMuons == 0"
    # ),

    "hadronic_cuts": (
        "n_RecoTracks > -1"
        " && RecoMissingEnergy_e > 40."
        " && n_RecoElectrons == 0"
        " && n_RecoMuons == 0"
        " && n_RecoJets < 4"
    ),

}
cutLabels = {

    "selNone": "Before selection",
    "sel_MET": "MET > 18 GeV",
    # "sel_nEle": "Number of reconstructed electrons < 3",
    # "sel_nEle_nMu": "Number of reconstructed muons < 3",
    # "sel_nEle_nMu_ee": "Invariant mass of ee < 85 GeV",
    # "sel_final": "Invariant mass of mm < 85 GeV",
    # "sel_semiLepHad": "MET > 18 GeV + (0 or 1 reco lepton",
    "sel_semiLeptonic": "MET > 18 GeV + 1 reco lepton ",
    "sel_semiLep_jetcut": "MET > 18 GeV +  1 reco lepton + < 4 jets",
    # "sel_semiLep_DVcut": "MET > 18 GeV +  1 reco lepton + < 4 jets + DV cut",
    "sel_semiLep_ntrks_cut": "MET > 18 GeV +  1 reco lepton + < 4 jets + DV with <= 3 tracks",
    "sel_semiLep_jet_cut": "MET > 18 GeV +  1 reco lepton + < 4 jets + DV with <= 3 tracks + leading jet_e < 120 GeV",
    # "fullyHadronic": "MET > 18 GeV + 0 reco leptons",
    "hadronic_cuts": "MET > 40 GeV + 0 leptons + < 4 jets",
}

'''
Dictionary for the output variables/histograms. 
The key is the name of the variable in the output files. 
"name" is the name of the variable in the input file, 
"title" is the x-axis label of the histogram, 
"bin" the number of bins of the histogram, 
"xmin" the minimum x-axis value and 
"xmax" the maximum x-axis value
'''

histoList = {
    # Gen-level stau
    "n_GenStau":           {"name":"n_GenStau",          "title":"Number of gen staus",             "bin":5,   "xmin":-0.5,  "xmax":4.5},
    "GenStau_vx":          {"name":"GenStau_vx",         "title":"Stau vertex x",                   "bin":50,  "xmin":-10,   "xmax":10},
    "GenStau_vy":          {"name":"GenStau_vy",         "title":"Stau vertex y",                   "bin":50,  "xmin":-10,   "xmax":10},
    "GenStau_vz":          {"name":"GenStau_vz",         "title":"Stau vertex z",                   "bin":50,  "xmin":-50,   "xmax":50},
    "GenStau_Lxy":         {"name":"GenStau_Lxy",        "title":"Transverse decay length Lxy",     "bin":50,  "xmin":0,     "xmax":10},
    "GenStau_Lxyz":        {"name":"GenStau_Lxyz",       "title":"3D decay length Lxyz",            "bin":50,  "xmin":0,     "xmax":10},
    # "GenStau_observed_lifetime_xyz":   {"name":"GenStau_observed_lifetime_xyz",   "title":"Observed stau lifetime (xyz)",    "bin":50,  "xmin":0,     "xmax":100},
    # "n_StauDaughters":       {"name":"n_StauDaughters",        "title":"Number of Stau Daughters",        "bin":5,   "xmin":-0.5, "xmax":9.5},
    "GenGravitino_e":        {"name":"GenGravitino_e",        "title":"Gen Gravitino energy",            "bin":100,  "xmin":0,    "xmax":200},
    # Tau
    "GenTau_e":              {"name":"GenTau_e",              "title":"Gen Tau energy",                  "bin":100,  "xmin":0,    "xmax":200},
    "GenTau_px":             {"name":"GenTau_px",             "title":"Gen Tau px",                      "bin":100,  "xmin":-200, "xmax":200},
    "GenTau_py":             {"name":"GenTau_py",             "title":"Gen Tau py",                      "bin":100,  "xmin":-200, "xmax":200},
    "GenTau_pz":             {"name":"GenTau_pz",             "title":"Gen Tau pz",                      "bin":100,  "xmin":-500, "xmax":500},
    "GenTau_pt":             {"name":"GenTau_pt",             "title":"Gen Tau pt",                      "bin":100,  "xmin":0,    "xmax":200},
    "GenTau_eta":            {"name":"GenTau_eta",            "title":"Gen Tau eta",                     "bin":100,  "xmin":-5,   "xmax":5},
    "GenTau_phi":            {"name":"GenTau_phi",            "title":"Gen Tau phi",                     "bin":100,  "xmin":-3.2, "xmax":3.2},

    "GenTau_vx":             {"name":"GenTau_vx",             "title":"Gen Tau vx",                      "bin":100,  "xmin":-10,  "xmax":10},
    "GenTau_vy":             {"name":"GenTau_vy",             "title":"Gen Tau vy",                      "bin":100,  "xmin":-10,  "xmax":10},
    "GenTau_vz":             {"name":"GenTau_vz",             "title":"Gen Tau vz",                      "bin":100,  "xmin":-50,  "xmax":50},
    "GenTau_cTau":           {"name":"GenTau_cTau",           "title":"cTau",                    "bin":120,   "xmin":0,    "xmax":60},
    # "decayLengthTau":        {"name": "decayLengthTau",    "title":"Gen Tau decay length",           "bin":50,   "xmin":0,    "xmax":50},

    # Final-state muons
    "n_FSGenMuon":       {"name":"n_FSGenMuon",     "title":"Number of FS muons", "bin":11,   "xmin":-0.5,  "xmax":10.5},
    "FSGenMuon_e":      {"name":"FSGenMuon_e",      "title":"FS muon energy",     "bin":50, "xmin":0, "xmax":200},
    "FSGenMuon_px":     {"name":"FSGenMuon_px",     "title":"FS muon px",         "bin":50, "xmin":-100, "xmax":100},
    "FSGenMuon_py":     {"name":"FSGenMuon_py",     "title":"FS muon py",         "bin":50, "xmin":-100, "xmax":100},
    "FSGenMuon_pz":     {"name":"FSGenMuon_pz",     "title":"FS muon pz",         "bin":50, "xmin":-200, "xmax":200},
    "FSGenMuon_pt":     {"name":"FSGenMuon_pt",     "title":"FS muon pt",         "bin":50, "xmin":0, "xmax":100},
    "FSGenMuon_eta":    {"name":"FSGenMuon_eta",    "title":"FS muon eta",        "bin":50, "xmin":-5, "xmax":5},
    "FSGenMuon_phi":    {"name":"FSGenMuon_phi",    "title":"FS muon phi",        "bin":64, "xmin":-3.2, "xmax":3.2},
    "FSGenMuon_vx":{"name":"FSGenMuon_vx","title":"FS muon vertex x","bin":50,"xmin":-10,"xmax":10},
    "FSGenMuon_vy":{"name":"FSGenMuon_vy","title":"FS muon vertex y","bin":50,"xmin":-10,"xmax":10},
    "FSGenMuon_vz":{"name":"FSGenMuon_vz","title":"FS muon vertex z","bin":50,"xmin":-50,"xmax":50},
    "FSGenMuon_charge": {"name":"FSGenMuon_charge", "title":"FS muon charge", "bin":3, "xmin":-1.5, "xmax":1.5},

    # Final-state electrons
    "n_FSGenElectron": {"name":"n_FSGenElectron","title":"Number of FS electrons","bin":11,"xmin":-0.5,"xmax":10.5},
    "FSGenElectron_e":    {"name":"FSGenElectron_e",    "title":"FS electron energy", "bin":50, "xmin":0, "xmax":200},
    "FSGenElectron_px":   {"name":"FSGenElectron_px",   "title":"FS electron px",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenElectron_py":   {"name":"FSGenElectron_py",   "title":"FS electron py",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenElectron_pz":   {"name":"FSGenElectron_pz",   "title":"FS electron pz",     "bin":50, "xmin":-200, "xmax":200},
    "FSGenElectron_pt":   {"name":"FSGenElectron_pt",   "title":"FS electron pt",     "bin":50, "xmin":0, "xmax":100},
    "FSGenElectron_eta":  {"name":"FSGenElectron_eta",  "title":"FS electron eta",    "bin":50, "xmin":-5, "xmax":5},
    "FSGenElectron_phi":  {"name":"FSGenElectron_phi",  "title":"FS electron phi",    "bin":64, "xmin":-3.2, "xmax":3.2},
    "FSGenElectron_vx": {"name":"FSGenElectron_vx","title":"FS electron vertex x","bin":50,"xmin":-10,"xmax":10},
    "FSGenElectron_vy": {"name":"FSGenElectron_vy","title":"FS electron vertex y","bin":50,"xmin":-10,"xmax":10},
    "FSGenElectron_vz": {"name":"FSGenElectron_vz","title":"FS electron vertex z","bin":50,"xmin":-50,"xmax":50},
    "FSGenElectron_charge": {"name":"FSGenElectron_charge", "title":"FS electron charge", "bin":3, "xmin":-1.5, "xmax":1.5},

    # Final-state neutrinos
    "n_FSGenNeutrino": {"name":"n_FSGenNeutrino","title":"Number of FS neutrinos","bin":11,"xmin":-0.5,"xmax":10.5},
    "FSGenNeutrino_e":    {"name":"FSGenNeutrino_e",    "title":"FS neutrino energy", "bin":160, "xmin":0, "xmax":160},
    "FSGenNeutrino_px":     {"name":"FSGenNeutrino_px",   "title":"FS neutrino px",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenNeutrino_py":     {"name":"FSGenNeutrino_py",   "title":"FS neutrino py",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenNeutrino_pz":     {"name":"FSGenNeutrino_pz",   "title":"FS neutrino pz",     "bin":50, "xmin":-200, "xmax":200},
    "FSGenNeutrino_pt":     {"name":"FSGenNeutrino_pt",   "title":"FS neutrino pt",    "bin":50, "xmin":0, "xmax":100},
    "FS_GenNeutrino_eta":    {"name":"FSGenNeutrino_eta",  "title":"FS neutrino eta",    "bin":50, "xmin":-5, "xmax":5},
    "FSGenNeutrino_phi":    {"name":"FSGenNeutrino_phi",  "title":"FS neutrino phi",    "bin":64, "xmin":-3.2, "xmax":3.2},

    # Reco Jets
    "n_RecoJets":     {"name":"n_RecoJets","title":"Number of reconstructed jets","bin":10,"xmin":-0.5,"xmax":9.5},
    "RecoJet_e":      {"name":"RecoJet_e",      "title":"Jet energy", "bin":50,"xmin":0,"xmax":200},
    "RecoJet_px":     {"name":"RecoJet_px",     "title":"Jet px",     "bin":50,"xmin":-100,"xmax":100},
    "RecoJet_py":     {"name":"RecoJet_py",     "title":"Jet py",     "bin":50,"xmin":-100,"xmax":100},
    "RecoJet_pz":     {"name":"RecoJet_pz",     "title":"Jet pz",     "bin":50,"xmin":-200,"xmax":200},
    "RecoJet_pt":     {"name":"RecoJet_pt",     "title":"Jet pt",     "bin":50,"xmin":0,"xmax":100},
    "RecoJet_eta":    {"name":"RecoJet_eta",    "title":"Jet eta",    "bin":50,"xmin":-5,"xmax":5},
    "RecoJet_phi":    {"name":"RecoJet_phi",    "title":"Jet phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
    "RecoJet_charge": {"name":"RecoJet_charge", "title":"Jet charge", "bin":3,"xmin":-1.5,"xmax":1.5},

    # Reco Electrons
    "n_RecoElectrons":     {"name":"n_RecoElectrons",   "title":"Number of reconstructed electrons","bin":5,"xmin":-0.5,"xmax":4.5},
    "RecoElectrons_e":      {"name":"RecoElectrons_e",      "title":"Reco electron energy", "bin":50,"xmin":0,"xmax":200},
    "RecoElectrons_px":     {"name":"RecoElectrons_px",     "title":"Reco electron px",     "bin":50,"xmin":-100,"xmax":100},
    "RecoElectrons_py":     {"name":"RecoElectrons_py",     "title":"Reco electron py",     "bin":50,"xmin":-100,"xmax":100},
    "RecoElectrons_pz":     {"name":"RecoElectrons_pz",     "title":"Reco electron pz",     "bin":50,"xmin":-200,"xmax":200},
    "RecoElectrons_pt":     {"name":"RecoElectrons_pt",     "title":"Reco electron pt",     "bin":50,"xmin":0,"xmax":100},
    "RecoElectrons_eta":    {"name":"RecoElectrons_eta",    "title":"Reco electron eta",    "bin":50,"xmin":-5,"xmax":5},
    "RecoElectrons_phi":    {"name":"RecoElectrons_phi",    "title":"Reco electron phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
    "RecoElectrons_charge": {"name":"RecoElectrons_charge", "title":"Reco electron charge","bin":3,"xmin":-1.5,"xmax":1.5},

    # # Reco Muons
    "n_RecoMuons":      {"name":"n_RecoMuons","title":"Number of reconstructed muons","bin":5,"xmin":-0.5,"xmax":4.5},
    "RecoMuons_e":      {"name":"RecoMuons_e",      "title":"Reco muon energy", "bin":50,"xmin":0,"xmax":200},
    "RecoMuons_px":     {"name":"RecoMuons_px",     "title":"Reco muon px",     "bin":50,"xmin":-100,"xmax":100},
    "RecoMuons_py":     {"name":"RecoMuons_py",     "title":"Reco muon py",     "bin":50,"xmin":-100,"xmax":100},
    "RecoMuons_pz":     {"name":"RecoMuons_pz",     "title":"Reco muon pz",     "bin":50,"xmin":-200,"xmax":200},
    "RecoMuons_pt":     {"name":"RecoMuons_pt",     "title":"Reco muon pt",     "bin":50,"xmin":0,"xmax":100},
    "RecoMuons_eta":    {"name":"RecoMuons_eta",    "title":"Reco muon eta",    "bin":50,"xmin":-5,"xmax":5},
    "RecoMuons_phi":    {"name":"RecoMuons_phi",    "title":"Reco muon phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
    "RecoMuons_charge": {"name":"RecoMuons_charge", "title":"Reco muon charge", "bin":3,"xmin":-1.5,"xmax":1.5},

    # MET
    "RecoMissingEnergy_e":   {"name":"RecoMissingEnergy_e",   "title":"Missing ET energy", "bin":50,"xmin":0,"xmax":200},
    "RecoMissingEnergy_pt":  {"name":"RecoMissingEnergy_pt",  "title":"Missing ET pt",     "bin":50,"xmin":0,"xmax":200},
    "RecoMissingEnergy_eta": {"name":"RecoMissingEnergy_eta", "title":"Missing ET eta",    "bin":50,"xmin":-5,"xmax":5},
    "RecoMissingEnergy_phi": {"name":"RecoMissingEnergy_phi", "title":"Missing ET phi",    "bin":64,"xmin":-3.2,"xmax":3.2},

    # Reco DVs from selected tracks
    "n_seltracks_DVs": {"name":"n_seltracks_DVs",   "title":"Number of reconstructed DVs",  "bin":11 , "xmin":-0.5, "xmax":10.5},
    "n_nonprimary_tracks": {"name":"n_nonprimary_tracks",   "title":"Number of non-primary tracks DVs",  "bin":20, "xmin":-0.5, "xmax":19.5},
    "n_RecoedPrimaryTracks": {"name":"n_RecoedPrimaryTracks",   "title":"Number of primary tracks DVs",  "bin":25, "xmin":-0.5, "xmax":50.5},
    "sel_tracks_pt_DV": {"name": "sel_tracks_pt_DV", "title": "pt of non-primary tracks", "bin":40, "xmin":0.0, "xmax": 80},
    "sel_tracks_D0_DV": {"name": "sel_tracks_D0_DV", "title": "D0 of non-primary tracks", "bin":200, "xmin":0.0, "xmax": 200.0},
    "sel_tracks_Z0_DV": {"name": "sel_tracks_Z0_DV", "title": "Z0 of non-primary tracks", "bin":200, "xmin":0.0, "xmax": 200.0},

    "n_trks_seltracks_DVs": {"name":"n_trks_seltracks_DVs",  "title":"Number of tracks per DV",   "bin":10, "xmin":-0.5, "xmax":9.5},
    "invMass_seltracks_DVs": {"name":"invMass_seltracks_DVs",  "title":"DV invariant mass [GeV]",    "bin":100, "xmin":0,    "xmax":10},
    "DV_evt_seltracks_chi2": {"name":"DV_evt_seltracks_chi2",    "title":"DV fit chi2",    "bin":10, "xmin":0, "xmax":10},
    "DV_evt_seltracks_normchi2": {"name":"DV_evt_seltracks_normchi2",    "title":"DV fit normalized chi2",    "bin":50, "xmin":0, "xmax":10},

    "Reco_DVs_merged_Lxy": {"name":"Reco_DVs_merged_Lxy","title":"DV L_{xy} [mm]","bin":250,"xmin":0,"xmax":250},
    "Reco_DVs_merged_Lxyz": {"name":"Reco_DVs_merged_Lxyz","title":"DV L_{xyz} [mm]","bin":250,"xmin":0,"xmax":250},

    # Reco dielectron system
    "Reco_ee_energy":  {"name":"Reco_ee_energy",  "title":"ee energy",   "bin":50, "xmin":0, "xmax":200},
    "Reco_ee_px":      {"name":"Reco_ee_px",      "title":"ee px",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_ee_py":      {"name":"Reco_ee_py",      "title":"ee py",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_ee_pz":      {"name":"Reco_ee_pz",      "title":"ee pz",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_ee_invMass": {"name":"Reco_ee_invMass","title":"ee inv mass","bin":50,"xmin":0,"xmax":100},

    # Reco dimuon system
    "Reco_mumu_energy":  {"name":"Reco_mumu_energy",  "title":"mu mu energy",   "bin":50, "xmin":0, "xmax":200},
    "Reco_mumu_px":      {"name":"Reco_mumu_px",      "title":"mu mu px",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_mumu_py":      {"name":"Reco_mumu_py",      "title":"mu mu py",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_mumu_pz":      {"name":"Reco_mumu_pz",      "title":"mu mu pz",       "bin":50, "xmin":-100, "xmax":100},
    "Reco_mumu_invMass": {"name":"Reco_mumu_invMass", "title":"mu mu inv mass", "bin":50, "xmin":0, "xmax":100},
    "muon_electron_overlap": {
        "name":"muon_electron_overlap",
        "title":"muon-electron overlap flag",
        "bin":2,
        "xmin":0,
        "xmax":2
    },

}

