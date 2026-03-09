'''
Final stage of the stau analysis
'''
# Input/output directories
# for signal:
# inputDir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/stage1_0603"
# for background:
inputDir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/condor_0603"
outputDir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final_0603/"

# List of datasets used in the analysis
processList = {
        #######################################################
        #                      SIGNAL                         #
        #######################################################
        #######################################################
        #             CME: 240 GeV (ZH) - 20 cm               #
        #######################################################
        'FCCee_100_stau_20cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_20cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_20cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_20cm_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH) - 50 cm               #
        #######################################################
        'FCCee_100_stau_50cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_50cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_50cm_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_50cm_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH) - 1 m                 #
        #######################################################
        'FCCee_100_stau_1m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_1m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_1m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_1m_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH) - 2 m                 #
        #######################################################
        'FCCee_100_stau_2m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_2m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_2m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_2m_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH) - 3 m                 #
        #######################################################
        'FCCee_100_stau_3m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_3m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_3m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_3m_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #             CME: 240 GeV (ZH) - 4 m                 #
        #######################################################
        'FCCee_100_stau_4m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_105_stau_4m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_110_stau_4m_ctau_ecm_240': {'fraction': 1.0},
        'FCCee_115_stau_4m_ctau_ecm_240': {'fraction': 1.0},

        #######################################################
        #               WINTER 2023                           #
        #######################################################
        'p8_ee_WW_ecm240': {'fraction': 1.0,'chunks':100},
        'p8_ee_ZZ_ecm240': {'fraction': 1.0,'chunks':100},
        'mgp8_ee_zh_ecm240_hbb': {'fraction': 1.0,'chunks':100},
        'wzp6_ee_nuenueH_Htautau_ecm240': {'fraction': 1.0,'chunks':100},
        'wzp6_ee_bbH_Htautau_ecm240': {'fraction': 1.0,'chunks':100},

        } 

prodTag = "FCCee/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
# procDict = "FCCee_procDict_spring2021_IDEA.json"


# Add samples which are not part of the offical process
procDictAdd = {

    # 20 cm
    'FCCee_100_stau_20cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_20cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_20cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_20cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    # 50 cm
    'FCCee_100_stau_50cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_50cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_50cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_50cm_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    # 1 m
    'FCCee_100_stau_1m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_1m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_1m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_1m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    # 2 m
    'FCCee_100_stau_2m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_2m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_2m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_2m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    # 3 m
    'FCCee_100_stau_3m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_3m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_3m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_3m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    # 4 m
    'FCCee_100_stau_4m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.07735, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_105_stau_4m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.05196, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_110_stau_4m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.02923, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_115_stau_4m_ctau_ecm_240': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.01067, "kfactor": 1.0, "matchingEfficiency": 1.0},

    }

# from madgraph: -Cross-section :   0.02923 +- 5.14e-06 pb
intLumi = 1.08e7

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
   
    # # "semiLepHad": ( 
    # #     "n_RecoTracks > -1"
    # #     "&& ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1) || (n_RecoElectrons == 0 && n_RecoMuons == 0))"
    # # ),

    # # "semiLeptonic": (
    # #     "n_RecoTracks > -1"
    # #     " && ((n_RecoElectrons > 0 && n_RecoMuons == 0) || (n_RecoElectrons ==0 && n_RecoMuons > 0))"
    # # ),

    "semiLeptonic_KV": (
        "n_RecoTracks > -1"
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons ==0 && n_RecoMuons == 1))"
        " && nKinkCandidates_passVeto > 0 "
        # " && KinkVertex_ntracks < 3 "
    ),

    "semiLep_DV": (
        "n_RecoTracks > -1"
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
        " && nKinkCandidates_passVeto == 0"
        " && nDisplacedVertices_failInnerHitVeto > 0 &&  nDisplacedVertices_failInnerHitVeto < 3"
    ),

    # ## DV region is better for small beta gamma ctau
    "hadronic_DV": (
        "n_RecoTracks > -1"
        " && n_RecoElectrons == 0"
        " && n_RecoMuons == 0"
        " && nKinkCandidates_passVeto == 0"
        " && nDisplacedVertices_failInnerHitVeto > 0 &&  nDisplacedVertices_failInnerHitVeto < 3"
    ), 

    # ## large beta gamma ctau
    "hadronic_KV": (
        " n_RecoTracks > -1"
        # " && n_RecoElectrons == 0"
        # " && n_RecoMuons == 0"
        " && nKinkCandidates_passVeto > 0 "
    ),
        
    # "hadronic_KV_lepallowed": (
    #     " n_RecoTracks > -1"
    #     " && nKinkCandidates_passVeto > 0 "
    # ),


}
cutLabels = {
    "selNone": "selNone",
    # "semiLepHad": "(0 or 1 reco lepton)",
    # "semiLeptonic": "semiLeptonic",
    "semiLeptonic_KV": "semiLeptonic_KV",
    "semiLep_DV": "semiLep_DV",
    "hadronic_KV":  "hadronic KV",
    "hadronic_DV": "hadronic_DV ",
    # "hadronic_KV_lepallowed" : "hadronic_KV_lepallowed",
    # "hadronic_allchannels": "testing all channel versions",
}

histoList = {
    # Gen-level stau
    # "n_GenStau":           {"name":"n_GenStau",          "title":"Number of gen staus",             "bin":5,   "xmin":-0.5,  "xmax":4.5},
    # "GenStau_vx":          {"name":"GenStau_vx",         "title":"Stau vertex x",                   "bin":50,  "xmin":-10,   "xmax":10},
    # "GenStau_vy":          {"name":"GenStau_vy",         "title":"Stau vertex y",                   "bin":50,  "xmin":-10,   "xmax":10},
    # "GenStau_vz":          {"name":"GenStau_vz",         "title":"Stau vertex z",                   "bin":50,  "xmin":-50,   "xmax":50},
    # "GenStau_Lxy":         {"name":"GenStau_Lxy",        "title":"Transverse decay length Lxy",     "bin":50,  "xmin":0,     "xmax":10},
    # "GenStau_Lxyz":        {"name":"GenStau_Lxyz",       "title":"3D decay length Lxyz",            "bin":50,  "xmin":0,     "xmax":10},
    # # "GenStau_observed_lifetime_xyz":   {"name":"GenStau_observed_lifetime_xyz",   "title":"Observed stau lifetime (xyz)",    "bin":50,  "xmin":0,     "xmax":100},
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
    # "GenTau_theta":          {"name":"GenTau_theta",          "title":"Gen Tau theta",                   "bin":100,  "xmin":0,    "xmax":3.2},
    "GenTau_vx":             {"name":"GenTau_vx",             "title":"Gen Tau vx",                      "bin":100,  "xmin":-10,  "xmax":10},
    "GenTau_vy":             {"name":"GenTau_vy",             "title":"Gen Tau vy",                      "bin":100,  "xmin":-10,  "xmax":10},
    "GenTau_vz":             {"name":"GenTau_vz",             "title":"Gen Tau vz",                      "bin":100,  "xmin":-50,  "xmax":50},
    # "GenTau_cTau":           {"name":"GenTau_cTau",           "title":"cTau",                    "bin":150,   "xmin":0,    "xmax":300},
    # "decayLengthTau":        {"name": "decayLengthTau",    "title":"Gen Tau decay length",           "bin":50,   "xmin":0,    "xmax":50},
    # "GenStau_theta":          {"name":"GenStau_theta",          "title":"Gen Stau theta",                   "bin":100,  "xmin":0,    "xmax":3.2},
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
    "nDisplaced_Vertices": {"name":"nDisplaced_Vertices",   "title":"Number of reconstructed DVs",  "bin":11 , "xmin":-0.5, "xmax":10.5},
    "n_nonprimary_tracks": {"name":"n_nonprimary_tracks",   "title":"Number of non-primary tracks DVs",  "bin":20, "xmin":-0.5, "xmax":19.5},
    "n_RecoedPrimaryTracks": {"name":"n_RecoedPrimaryTracks",   "title":"Number of primary tracks DVs",  "bin":25, "xmin":-0.5, "xmax":50.5},
    
    "sel_tracks_pt_DV": {"name": "sel_tracks_pt_DV", "title": "pt of non-primary tracks", "bin":40, "xmin":0.0, "xmax": 80},
    "sel_tracks_D0_DV": {"name": "sel_tracks_D0_DV", "title": "D0 of non-primary tracks", "bin":200, "xmin":0.0, "xmax": 200.0},
    "sel_tracks_Z0_DV": {"name": "sel_tracks_Z0_DV", "title": "Z0 of non-primary tracks", "bin":200, "xmin":0.0, "xmax": 200.0},

    "nDisplacedVertices_failInnerHitVeto" :  {"name":"nDisplacedVertices_failInnerHitVeto",  "title":"#DV failing the hit veto",   "bin":10, "xmin":-0.5, "xmax":9.5},
    "nTracks_DV_failInnerHitVeto": {"name":"nTracks_DV_failInnerHitVeto",  "title":"Number of tracks per DV (fail hit veto)",   "bin":10, "xmin":-0.5, "xmax":9.5},
    "invMass_seltracks_DVs": {"name":"invMass_seltracks_DVs",  "title":"DV invariant mass [GeV]",    "bin":100, "xmin":0,    "xmax":10},
    "invMass_seltracks_DVs_zoom": {"name":"invMass_seltracks_DVs",  "title":"DV invariant mass [GeV]",    "bin":10, "xmin":0,    "xmax":2},
    "DV_evt_seltracks_chi2": {"name":"DV_evt_seltracks_chi2",    "title":"DV fit chi2",    "bin":10, "xmin":0, "xmax":10},
    "DV_evt_seltracks_normchi2": {"name":"DV_evt_seltracks_normchi2",    "title":"DV fit normalized chi2",    "bin":50, "xmin":0, "xmax":10},
    "Reco_seltracks_DVs_Lxy": {"name":"Reco_seltracks_DVs_Lxy","title":"DV L_{xy} [mm]","bin":100,"xmin":0,"xmax":250},
    "Reco_seltracks_DVs_Lxyz": {"name":"Reco_seltracks_DVs_Lxyz","title":"DV L_{xyz} [mm]","bin":100,"xmin":0,"xmax":250},

    "RecoVisibleEnergy": {"name":"RecoVisibleEnergy", "title":"Visible energy", "bin":120,"xmin":0,"xmax":240},
    "RecoMissingEnergy3D": {"name":"RecoMissingEnergy3D", "title":"Calculated Missing energy", "bin":120,"xmin":0,"xmax":240},
    
    "KinkCandidates_passInnerHitVeto": {"name":"KinkCandidates_passInnerHitVeto", "title":"Number of kink vertices", "bin":5, "xmin":-0.5, "xmax":4.5},
    "nKinkCandidates_passVeto" : {"name":"nKinkCandidates_passVeto", "title":"Number of kink vertices passing the hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
    "KinkVertex_invMass": {"name":"KinkVertex_invMass", "title":"Invariant mass of kink vertex", "bin":150, "xmin":0, "xmax":150},
    "KinkVertex_dxy": {"name":"KinkVertex_dxy", "title":"dxy of kink vertex", "bin":100, "xmin":0, "xmax":2000},
    "KinkVertex_d3d": {"name":"KinkVertex_d3d", "title":"d3d of kink vertex", "bin":100, "xmin":0, "xmax":2000},
    "KinkVertex_ntracks": {"name":"KinkVertex_ntracks", "title":"Number of tracks in kink vertex", "bin":10, "xmin":-0.5, "xmax":9.5},
    "nKinkVertices": {"name":"nKinkVertices", "title":"Number of kink vertices before hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
    "KinkAngle": {"name": "KinkAngle", "title": "angle between the r_pvkv and P_kv", "bin":180, "xmin":0, "xmax":180},
}

