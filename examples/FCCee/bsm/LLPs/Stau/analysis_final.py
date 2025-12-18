'''
Final stage of the stau analysis
'''

# Input/output directories
inputDir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/working"
outputDir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final"

# List of datasets used in the analysis
process_list = {
        #######################################################
        #               CME: 240 GeV (ZH)                     #
        #######################################################
        'FCC_100stau_240com_000'  : {'fraction': 1.0},
        } 


procDict = "FCCee_procDict_winter2023_IDEA.json"


intLumi = 2.05e8
#intLumi = 1.92e7
#intLumi = 1.08e7
#intLumi = 2.7e6
doScale = True

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
    "sel": " RecoMissingEnergy_p > 20.",
}

cutLabels = {
    "selNone": "Before selection",
    "sel": "Missing Energy > 20 GeV "
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


    # Final-state muons
    "n_FSGenMuon":       {"name":"n_FSGenMuon",     "title":"Number of FS muons", "bin":5,   "xmin":-0.5,  "xmax":4.5},
    "FSGenMuon_e":      {"name":"FSGenMuon_e",      "title":"FS muon energy",     "bin":50, "xmin":0, "xmax":200},
    "FSGenMuon_px":     {"name":"FSGenMuon_px",     "title":"FS muon px",         "bin":50, "xmin":-100, "xmax":100},
    "FSGenMuon_py":     {"name":"FSGenMuon_py",     "title":"FS muon py",         "bin":50, "xmin":-100, "xmax":100},
    "FSGenMuon_pz":     {"name":"FSGenMuon_pz",     "title":"FS muon pz",         "bin":50, "xmin":-200, "xmax":200},
    "FSGenMuon_pt":     {"name":"FSGenMuon_pt",     "title":"FS muon pt",         "bin":50, "xmin":0, "xmax":100},
    "FSGenMuon_eta":    {"name":"FSGenMuon_eta",    "title":"FS muon eta",        "bin":50, "xmin":-5, "xmax":5},
    "FSGenMuon_phi":    {"name":"FSGenMuon_phi",    "title":"FS muon phi",        "bin":64, "xmin":-3.2, "xmax":3.2},
    # "FSGenMuon_vertex_x":{"name":"FSGenMuon_vertex_x","title":"FS muon vertex x","bin":50,"xmin":-10,"xmax":10},
    # "FSGenMuon_vertex_y":{"name":"FSGenMuon_vertex_y","title":"FS muon vertex y","bin":50,"xmin":-10,"xmax":10},
    # "FSGenMuon_vertex_z":{"name":"FSGenMuon_vertex_z","title":"FS muon vertex z","bin":50,"xmin":-50,"xmax":50},
    "FSGenMuon_charge": {"name":"FSGenMuon_charge", "title":"FS muon charge", "bin":3, "xmin":-1.5, "xmax":1.5},

    # Final-state electrons
    "n_FSGenElectron": {"name":"n_FSGenElectron","title":"Number of FS electrons","bin":5,"xmin":-0.5,"xmax":4.5},
    "FSGenElectron_e":    {"name":"FSGenElectron_e",    "title":"FS electron energy", "bin":50, "xmin":0, "xmax":200},
    "FSGenElectron_px":   {"name":"FSGenElectron_px",   "title":"FS electron px",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenElectron_py":   {"name":"FSGenElectron_py",   "title":"FS electron py",     "bin":50, "xmin":-100, "xmax":100},
    "FSGenElectron_pz":   {"name":"FSGenElectron_pz",   "title":"FS electron pz",     "bin":50, "xmin":-200, "xmax":200},
    "FSGenElectron_pt":   {"name":"FSGenElectron_pt",   "title":"FS electron pt",     "bin":50, "xmin":0, "xmax":100},
    "FSGenElectron_eta":  {"name":"FSGenElectron_eta",  "title":"FS electron eta",    "bin":50, "xmin":-5, "xmax":5},
    "FSGenElectron_phi":  {"name":"FSGenElectron_phi",  "title":"FS electron phi",    "bin":64, "xmin":-3.2, "xmax":3.2},
    # "FSGenElectron_vertex_x": {"name":"FSGenElectron_vertex_x","title":"FS electron vertex x","bin":50,"xmin":-10,"xmax":10},
    # "FSGenElectron_vertex_y": {"name":"FSGenElectron_vertex_y","title":"FS electron vertex y","bin":50,"xmin":-10,"xmax":10},
    # "FSGenElectron_vertex_z": {"name":"FSGenElectron_vertex_z","title":"FS electron vertex z","bin":50,"xmin":-50,"xmax":50},
    "FSGenElectron_charge": {"name":"FSGenElectron_charge", "title":"FS electron charge", "bin":3, "xmin":-1.5, "xmax":1.5},

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
    # "RecoElectron_e":      {"name":"RecoElectron_e",      "title":"Reco electron energy", "bin":50,"xmin":0,"xmax":200},
    # "RecoElectron_px":     {"name":"RecoElectron_px",     "title":"Reco electron px",     "bin":50,"xmin":-100,"xmax":100},
    # "RecoElectron_py":     {"name":"RecoElectron_py",     "title":"Reco electron py",     "bin":50,"xmin":-100,"xmax":100},
    # "RecoElectron_pz":     {"name":"RecoElectron_pz",     "title":"Reco electron pz",     "bin":50,"xmin":-200,"xmax":200},
    # "RecoElectron_pt":     {"name":"RecoElectron_pt",     "title":"Reco electron pt",     "bin":50,"xmin":0,"xmax":100},
    # "RecoElectron_eta":    {"name":"RecoElectron_eta",    "title":"Reco electron eta",    "bin":50,"xmin":-5,"xmax":5},
    # "RecoElectron_phi":    {"name":"RecoElectron_phi",    "title":"Reco electron phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
    # "RecoElectron_charge": {"name":"RecoElectron_charge", "title":"Reco electron charge","bin":3,"xmin":-1.5,"xmax":1.5},

    # Reco Muons
    "n_RecoMuons":      {"name":"n_RecoMuons","title":"Number of reconstructed muons","bin":5,"xmin":-0.5,"xmax":4.5},
    # "RecoMuons_e":      {"name":"RecoMuons_e",      "title":"Reco muon energy", "bin":50,"xmin":0,"xmax":200},
    # "RecoMuons_px":     {"name":"RecoMuons_px",     "title":"Reco muon px",     "bin":50,"xmin":-100,"xmax":100},
    # "RecoMuons_py":     {"name":"RecoMuons_py",     "title":"Reco muon py",     "bin":50,"xmin":-100,"xmax":100},
    # "RecoMuons_pz":     {"name":"RecoMuons_pz",     "title":"Reco muon pz",     "bin":50,"xmin":-200,"xmax":200},
    # "RecoMuons_pt":     {"name":"RecoMuons_pt",     "title":"Reco muon pt",     "bin":50,"xmin":0,"xmax":100},
    # "RecoMuons_eta":    {"name":"RecoMuons_eta",    "title":"Reco muon eta",    "bin":50,"xmin":-5,"xmax":5},
    # "RecoMuons_phi":    {"name":"RecoMuons_phi",    "title":"Reco muon phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
    # "RecoMuons_charge": {"name":"RecoMuons_charge", "title":"Reco muon charge", "bin":3,"xmin":-1.5,"xmax":1.5},

    # MET
    "RecoMissingEnergy_e":   {"name":"RecoMissingEnergy_e",   "title":"Missing ET energy", "bin":50,"xmin":0,"xmax":200},
    "RecoMissingEnergy_pt":  {"name":"RecoMissingEnergy_pt",  "title":"Missing ET pt",     "bin":50,"xmin":0,"xmax":200},
    "RecoMissingEnergy_eta": {"name":"RecoMissingEnergy_eta", "title":"Missing ET eta",    "bin":50,"xmin":-5,"xmax":5},
    "RecoMissingEnergy_phi": {"name":"RecoMissingEnergy_phi", "title":"Missing ET phi",    "bin":64,"xmin":-3.2,"xmax":3.2},
}

