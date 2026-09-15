'''
Final stage of the stau analysis
'''
import argparse
import os

# Input/output directories
# for signal:
# inputDir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/stage1_0603"
# for background:
# inputDir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/condor_0603"
# combined:
parser = argparse.ArgumentParser(description="FolderNames")
parser.add_argument("-o", "--output_folder_name", type=str, default="")

args, _ = parser.parse_known_args()

inputDir  = "/afs/desy.de/user/m/mrandria/DUST/STAGE1"
outputDir = "/afs/desy.de/user/m/mrandria/DUST/FINAL"
os.makedirs(outputDir, exist_ok=True)

# List of datasets used in the analysis
processList = {
        # ###################################################
        # #                     SIGNAL                      #
        # ###################################################
        # ###################################################
        #             FCCee: 365 GeV - 0.5 m              #
        ###################################################
        "FCCee_120_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_0.5m_ctau_ecm_365": {'fraction': 1.0},

        ###################################################
        #             FCCee: 365 GeV - 1 m                #
        ###################################################
        "FCCee_120_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_1m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_1m_ctau_ecm_365": {'fraction': 1.0},

        ###################################################
        #             FCCee: 365 GeV - 2 m                #
        ###################################################
        "FCCee_120_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_2m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_2m_ctau_ecm_365": {'fraction': 1.0},

        ###################################################
        #             FCCee: 365 GeV - 5 m                #
        ###################################################
        "FCCee_120_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_5m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_5m_ctau_ecm_365": {'fraction': 1.0},

        # ###################################################
        # #             FCCee: 365 GeV - 10 m               #
        # ###################################################
        "FCCee_120_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_10m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_10m_ctau_ecm_365": {'fraction': 1.0},

        ###################################################
        #             FCCee: 365 GeV - 20 m               #
        ###################################################
        "FCCee_120_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_20m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_20m_ctau_ecm_365": {'fraction': 1.0},

        ###################################################
        #             FCCee: 365 GeV - 50 m               #
        ###################################################
        "FCCee_120_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_130_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_140_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_150_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_160_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_170_stau_50m_ctau_ecm_365": {'fraction': 1.0},
        "FCCee_180_stau_50m_ctau_ecm_365": {'fraction': 1.0},

        ##################################################
        #                BACKGROUND                    #
        ##################################################
        'p8_ee_WW_ecm365': {'fraction': 1.0,'chunks':100},
        'p8_ee_ZZ_ecm365': {'fraction': 1.0,'chunks':100},
        'wzp6_ee_nuenueH_Htautau_ecm365': {'fraction': 1.0,'chunks':100},
        'wzp6_ee_bbH_Htautau_ecm365': {'fraction': 1.0,'chunks':100},
        'p8_ee_tt_ecm365': {'fraction': 1.0,'chunks':100},
        'wzp6_ee_tautau_ecm365': {'fraction': 1.0,'chunks':100},
    }

prodTag = "FCCee/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
# procDict = "FCCee_procDict_spring2021_IDEA.json"


# Add samples which are not part of the offical process
procDictAdd = {
    'FCCee_120_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_0.5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},

    'FCCee_120_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_1m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},
    
    'FCCee_120_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_2m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},

    'FCCee_120_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_5m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},

    'FCCee_120_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_10m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},

    'FCCee_120_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_20m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},

    'FCCee_120_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.36226000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_130_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 6.75979000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_140_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 5.16167000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_150_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 3.61349000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_160_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 2.17607000e-02, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_170_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 9.40847000e-03, "kfactor": 1.0, "matchingEfficiency": 1.0},
    'FCCee_180_stau_50m_ctau_ecm_365': {"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 8.77581000e-04, "kfactor": 1.0, "matchingEfficiency": 1.0},
}

# intLumi = 0.67e6 / 4 # lumi for 1 year per IP
intLumi = 0.67e6 / 4

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

    "semiLep_KV": (
        "n_RecoTracks > -1"
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons ==0 && n_RecoMuons == 1))"
        " && (Any(RecoElectrons_pt > 10) || Any(RecoMuons_pt > 10))"
        " && Any(nKinkCandidates_passVeto > 0 && abs(KinkAngle) > 20)"
        # " && nKinkCandidates_passVeto > 0"
        # " && Any(abs(KinkAngle) > 20)"
    ),

    "semiLep_DV": (
        "n_RecoTracks > -1"
        " && ((n_RecoElectrons == 1 && n_RecoMuons == 0) || (n_RecoElectrons == 0 && n_RecoMuons == 1))"
        " && (Any(RecoElectrons_pt > 10) || Any(RecoMuons_pt > 10))"
        " && nKinkCandidates_passVeto == 0"
        " && nDisplacedVertices_failInnerHitVeto > 0 && nDisplacedVertices_failInnerHitVeto < 3"
        " && All(PV2V0Cos < 0.95)"
    ),

    "hadronic_KV": (
        " n_RecoTracks > -1"
        " && n_RecoElectrons == 0"
        " && n_RecoMuons == 0"
        " && nKinkCandidates_passVeto > 0"
        " && Any(abs(KinkAngle) > 20)"
    ),

    "hadronic_DV": (
        "n_RecoTracks > -1"
        " && n_RecoElectrons == 0"
        " && n_RecoMuons == 0"
        " && nKinkCandidates_passVeto == 0"
        " && nDisplacedVertices_failInnerHitVeto > 0 && nDisplacedVertices_failInnerHitVeto < 3"
        " && Any(PV2V0Cos < 0.95)"
    ),

    # "escaping": (
    #     "n_RecoTracks > -1"
    #     " && nKinkCandidates_passVeto == 0"
    #     " && nDisplacedVertices_failInnerHitVeto == 0"
    # ),

    "escaping_Staus": (
        "n_RecoTracks > -1"
        #" && n_RecoedPrimaryTracks < 3"
        " && nKinkCandidates_passVeto == 0"
        " && nDisplacedVertices_failInnerHitVeto == 0"
        " && Any(TOF_EscapingTracks > 12500)"
        " && All(TOF_EscapingTracks > 12500)"
        #" && Any(RecoedPrimaryTracks_pt > 50)"
    ),
}

cutLabels = {
    "selNone": "selNone",
    "semiLep_KV": "semiLep_KV",
    "semiLep_DV": "semiLep_DV",
    "hadronic_KV": "hadronic_KV",
    "hadronic_DV": "hadronic_DV",
    # "escaping": "escaping",
    "escaping_Staus": "escaping_Staus",
}


# histoList = {
#     #================ Track information ================#
#     "n_RecoedPrimaryTracks": {"name":"n_RecoedPrimaryTracks", "title":"Number of primary tracks DVs", "bin":25, "xmin":-0.5, "xmax":50.5},
#     "n_AcceptedTracks": {"name":"n_AcceptedTracks", "title":"Number of accepted tracks", "bin":50, "xmin":-0.5, "xmax":49.5},
#     "TOF_AcceptedTracks": {"name":"TOF_AcceptedTracks", "title":"Time of flight [ps]", "bin":200, "xmin":0, "xmax":50000},
#     "PrimaryVertex_ntracks": {"name":"PrimaryVertex_ntracks", "title":"Number of tracks at primary vertex", "bin":50, "xmin":-0.5, "xmax":49.5},
#     "n_RecoTracks": {"name":"n_RecoTracks", "title":"Number of reconstructed tracks", "bin":50, "xmin":-0.5, "xmax":49.5},
#     "n_nonprimary_tracks": {"name":"n_nonprimary_tracks", "title":"Number of non-primary tracks", "bin":20, "xmin":-0.5, "xmax":19.5},


#     #================ DV information ================#
#     "nDisplaced_Vertices": {"name":"nDisplaced_Vertices", "title":"Number of reconstructed DVs", "bin":11, "xmin":-0.5, "xmax":10.5},
#     "nTracks_DV": {"name":"nTracks_DV", "title":"Number of tracks per DV", "bin":10, "xmin":-0.5, "xmax":9.5},
#     "nDisplacedVertices_failInnerHitVeto": {"name":"nDisplacedVertices_failInnerHitVeto", "title":"Number of DVs failing the inner hit veto", "bin":10, "xmin":-0.5, "xmax":9.5},
#     "nTracks_DV_failInnerHitVeto": {"name":"nTracks_DV_failInnerHitVeto", "title":"Number of tracks per DV failing the inner hit veto", "bin":10, "xmin":-0.5, "xmax":9.5},
#     "invMass_seltracks_DVs": {"name":"invMass_seltracks_DVs", "title":"DV invariant mass [GeV]", "bin":100, "xmin":0, "xmax":10},
#     "DV_evt_seltracks_chi2": {"name":"DV_evt_seltracks_chi2", "title":"DV fit #chi^{2}", "bin":10, "xmin":0, "xmax":10},
#     "DV_evt_seltracks_normchi2": {"name":"DV_evt_seltracks_normchi2", "title":"DV fit normalized #chi^{2}", "bin":50, "xmin":0, "xmax":10},
#     "Reco_seltracks_DVs_Lxy": {"name":"Reco_seltracks_DVs_Lxy", "title":"DV transverse decay length L_{xy} [mm]", "bin":100, "xmin":0, "xmax":250},
#     "Reco_seltracks_DVs_Lxyz": {"name":"Reco_seltracks_DVs_Lxyz", "title":"DV 3D decay length L_{xyz} [mm]", "bin":100, "xmin":0, "xmax":250},


#     #================ Reco Electrons ================#
#     "n_RecoElectrons": {"name":"n_RecoElectrons", "title":"Number of reconstructed electrons", "bin":5, "xmin":-0.5, "xmax":4.5},
#     "RecoElectrons_e": {"name":"RecoElectrons_e", "title":"Reco electron energy [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoElectrons_p": {"name":"RecoElectrons_p", "title":"Reco electron momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoElectrons_p_NEW": {"name":"RecoElectrons_p_NEW", "title":"Reco electron momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoElectrons_pt": {"name":"RecoElectrons_pt", "title":"Reco electron p_{T} [GeV]", "bin":50, "xmin":0, "xmax":100},
#     "RecoElectrons_px": {"name":"RecoElectrons_px", "title":"Reco electron p_{x} [GeV]", "bin":50, "xmin":-100, "xmax":100},
#     "RecoElectrons_py": {"name":"RecoElectrons_py", "title":"Reco electron p_{y} [GeV]", "bin":50, "xmin":-100, "xmax":100},
#     "RecoElectrons_pz": {"name":"RecoElectrons_pz", "title":"Reco electron p_{z} [GeV]", "bin":50, "xmin":-200, "xmax":200},
#     "RecoElectrons_theta": {"name":"RecoElectrons_theta", "title":"Reco electron #theta", "bin":100, "xmin":0, "xmax":3.2},
#     "RecoElectrons_phi": {"name":"RecoElectrons_phi", "title":"Reco electron #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
#     "RecoElectrons_charge": {"name":"RecoElectrons_charge", "title":"Reco electron charge", "bin":3, "xmin":-1.5, "xmax":1.5},


#     #================ Reco Muons ================#
#     "n_RecoMuons": {"name":"n_RecoMuons", "title":"Number of reconstructed muons", "bin":5, "xmin":-0.5, "xmax":4.5},
#     "RecoMuons_e": {"name":"RecoMuons_e", "title":"Reco muon energy [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoMuons_p": {"name":"RecoMuons_p", "title":"Reco muon momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoMuons_p_NEW": {"name":"RecoMuons_p_NEW", "title":"Reco muon momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
#     "RecoMuons_pt": {"name":"RecoMuons_pt", "title":"Reco muon p_{T} [GeV]", "bin":50, "xmin":0, "xmax":100},
#     "RecoMuons_px": {"name":"RecoMuons_px", "title":"Reco muon p_{x} [GeV]", "bin":50, "xmin":-100, "xmax":100},
#     "RecoMuons_py": {"name":"RecoMuons_py", "title":"Reco muon p_{y} [GeV]", "bin":50, "xmin":-100, "xmax":100},
#     "RecoMuons_pz": {"name":"RecoMuons_pz", "title":"Reco muon p_{z} [GeV]", "bin":50, "xmin":-200, "xmax":200},
#     "RecoMuons_eta": {"name":"RecoMuons_eta", "title":"Reco muon #eta", "bin":50, "xmin":-5, "xmax":5},
#     "RecoMuons_theta": {"name":"RecoMuons_theta", "title":"Reco muon #theta", "bin":100, "xmin":0, "xmax":3.2},
#     "RecoMuons_phi": {"name":"RecoMuons_phi", "title":"Reco muon #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
#     "RecoMuons_charge": {"name":"RecoMuons_charge", "title":"Reco muon charge", "bin":3, "xmin":-1.5, "xmax":1.5},


#     # #================ Primary track variables ================#
#     "RecoedPrimaryTracks_charge": {"name":"RecoedPrimaryTracks_charge", "title":"Recoed primary tracks charge", "bin":3, "xmin":-1.5, "xmax":1.5},
#     "sel_tracks_charge": {"name":"sel_tracks_charge", "title":"Selected tracks charge", "bin":3, "xmin":-1.5, "xmax":1.5},
#     "RecoedPrimaryTracks_d0": {"name":"RecoedPrimaryTracks_d0", "title":"Recoed primary tracks d_{0} [mm]", "bin":100, "xmin":-100, "xmax":100},
#     "RecoedPrimaryTracks_phi": {"name":"RecoedPrimaryTracks_phi", "title":"Recoed primary tracks #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
#     "RecoedPrimaryTracks_theta": {"name":"RecoedPrimaryTracks_theta", "title":"Recoed primary tracks #theta", "bin":30, "xmin":-4, "xmax":4},
#     "RecoedPrimaryTracks_p": {"name":"RecoedPrimaryTracks_p", "title":"Reconstructed primary tracks: p [GeV]", "bin":100, "xmin":0, "xmax":365},
#     "RecoedPrimaryTracks_pt": {"name":"RecoedPrimaryTracks_pt", "title":"Reconstructed primary tracks: p_{T} [GeV]", "bin":100, "xmin":0, "xmax":365},


#     #================ Kinked candidates ================#
#     "KinkCandidates_passInnerHitVeto": {"name":"KinkCandidates_passInnerHitVeto", "title":"Number of kink vertices", "bin":5, "xmin":-0.5, "xmax":4.5},
#     "nKinkCandidates_passVeto": {"name":"nKinkCandidates_passVeto", "title":"Number of kink vertices passing the hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
#     "KinkVertex_invMass": {"name":"KinkVertex_invMass", "title":"Invariant mass of kink vertex [GeV]", "bin":150, "xmin":0, "xmax":150},
#     "nKinkVertices": {"name":"nKinkVertices", "title":"Number of kink vertices before hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
#     "KinkVertex_ntracks": {"name":"KinkVertex_ntracks", "title":"Number of tracks in kink vertex", "bin":10, "xmin":-0.5, "xmax":9.5},
#     "KinkAngle": {"name":"KinkAngle", "title":"Angle between the r_{PVKV} and P_{KV}", "bin":180, "xmin":0, "xmax":180},
#     "KinkVertex_dxy": {"name":"KinkVertex_dxy", "title":"d_{xy} of kink vertex [mm]", "bin":100, "xmin":0, "xmax":2000},
#     "KinkVertex_d3d": {"name":"KinkVertex_d3d", "title":"d_{3D} of kink vertex [mm]", "bin":100, "xmin":0, "xmax":2000},


#     #================ Additional track / TOF variables ================#
#     "PV2V0Cos": {"name":"PV2V0Cos", "title":"Cosine angle between displaced vertex and primary vertex", "bin":100, "xmin":-1, "xmax":1},
#     "RecoedPrimaryTrack_mass": {"name":"RecoedPrimaryTrack_mass", "title":"Mass of escaping tracks [GeV]", "bin":50, "xmin":100, "xmax":200},
#     "TOF_RecoedPrimaryTracks": {"name":"TOF_RecoedPrimaryTracks", "title":"Time of flight of recoed primary tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
#     "TOF_sel_tracks": {"name":"TOF_sel_tracks", "title":"Time of flight of selected tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
#     "TOF_EscapingTracks": {"name":"TOF_EscapingTracks", "title":"Time of flight of escaping tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
#     "EscapingTracks_mass": {"name":"EscapingTracks_mass", "title":"Mass of escaping tracks [GeV]", "bin":50, "xmin":100, "xmax":200},
#     "EscapingTracks_pt": {"name":"EscapingTracks_pt", "title":"Transverse momentum of escaping tracks [GeV]", "bin":100, "xmin":0, "xmax":365},
#     "EscapingTracks_p": {"name":"EscapingTracks_p", "title":"Momentum of escaping tracks [GeV]", "bin":100, "xmin":0, "xmax":365},
#     "EscapingTracks_theta": {"name":"EscapingTracks_theta", "title":"Theta of escaping tracks", "bin":100, "xmin":0, "xmax":3.2},
#     "EscapingTracks_phi": {"name":"EscapingTracks_phi", "title":"Phi of escaping tracks", "bin":64, "xmin":-3.2, "xmax":3.2},
#     "n_EscapingTracks": {"name":"n_EscapingTracks", "title":"Number of escaping tracks", "bin":25, "xmin":-0.5, "xmax":24.5},
#     # "RecoedPrimaryTracks_firstHitLoc": {"name":"RecoedPrimaryTracks_firstHitLoc", "title":"Reconstructed primary tracks: first hit location [mm]", "bin":100, "xmin":0, "xmax":365},
#     # "RecoedPrimaryTracks_lastHitLoc": {"name":"RecoedPrimaryTracks_lastHitLoc", "title":"Reconstructed primary tracks: last hit location [mm]", "bin":100, "xmin":0, "xmax":365},
    
# }


histoList = {
    #================ Track information ================#
    "n_RecoedPrimaryTracks": {"name":"n_RecoedPrimaryTracks", "title":"Number of primary tracks DVs", "bin":25, "xmin":-0.5, "xmax":50.5},
    "n_AcceptedTracks": {"name":"n_AcceptedTracks", "title":"Number of accepted tracks", "bin":50, "xmin":-0.5, "xmax":49.5},
    "TOF_AcceptedTracks": {"name":"TOF_AcceptedTracks", "title":"Time of flight [ps]", "bin":200, "xmin":0, "xmax":50000},
    "PrimaryVertex_ntracks": {"name":"PrimaryVertex_ntracks", "title":"Number of tracks at primary vertex", "bin":50, "xmin":-0.5, "xmax":49.5},
    "n_RecoTracks": {"name":"n_RecoTracks", "title":"Number of reconstructed tracks", "bin":50, "xmin":-0.5, "xmax":49.5},
    "n_nonprimary_tracks": {"name":"n_nonprimary_tracks", "title":"Number of non-primary tracks", "bin":20, "xmin":-0.5, "xmax":19.5},


    #================ DV information ================#
    "nDisplaced_Vertices": {"name":"nDisplaced_Vertices", "title":"Number of reconstructed DVs", "bin":11, "xmin":-0.5, "xmax":10.5},
    "nTracks_DV": {"name":"nTracks_DV", "title":"Number of tracks per DV", "bin":10, "xmin":-0.5, "xmax":9.5},
    "nDisplacedVertices_failInnerHitVeto": {"name":"nDisplacedVertices_failInnerHitVeto", "title":"Number of DVs failing the inner hit veto", "bin":10, "xmin":-0.5, "xmax":9.5},
    "nTracks_DV_failInnerHitVeto": {"name":"nTracks_DV_failInnerHitVeto", "title":"Number of tracks per DV failing the inner hit veto", "bin":10, "xmin":-0.5, "xmax":9.5},
    "invMass_seltracks_DVs": {"name":"invMass_seltracks_DVs", "title":"DV invariant mass [GeV]", "bin":100, "xmin":0, "xmax":10},
    "DV_evt_seltracks_chi2": {"name":"DV_evt_seltracks_chi2", "title":"DV fit #chi^{2}", "bin":10, "xmin":0, "xmax":10},
    "DV_evt_seltracks_normchi2": {"name":"DV_evt_seltracks_normchi2", "title":"DV fit normalized #chi^{2}", "bin":50, "xmin":0, "xmax":10},
    "Reco_seltracks_DVs_Lxy": {"name":"Reco_seltracks_DVs_Lxy", "title":"DV transverse decay length L_{xy} [mm]", "bin":100, "xmin":0, "xmax":250},
    "Reco_seltracks_DVs_Lxyz": {"name":"Reco_seltracks_DVs_Lxyz", "title":"DV 3D decay length L_{xyz} [mm]", "bin":100, "xmin":0, "xmax":250},


    #================ Hit pattern information ================#
    # "RecoParticles_firstHitLoc": {"name":"RecoParticles_firstHitLoc", "title":"First hit location [mm]", "bin":100, "xmin":0, "xmax":365},
    # "RecoParticles_lastHitLoc": {"name":"RecoParticles_lastHitLoc", "title":"Last hit location [mm]", "bin":100, "xmin":0, "xmax":365},
    # "RecoParticles_nHits": {"name":"RecoParticles_nHits", "title":"Number of hits", "bin":50, "xmin":0, "xmax":50},
    # "RecoParticles_nDriftChamberHits": {"name":"RecoParticles_nDriftChamberHits", "title":"Number of drift chamber hits", "bin":50, "xmin":0, "xmax":50},


    #================ Reco Electrons ================# 
    "n_RecoElectrons": {"name":"n_RecoElectrons", "title":"Number of reconstructed electrons", "bin":5, "xmin":-0.5, "xmax":4.5},
    "RecoElectrons_e": {"name":"RecoElectrons_e", "title":"Reco electron energy [GeV]", "bin":50, "xmin":0, "xmax":200},
    "RecoElectrons_p": {"name":"RecoElectrons_p", "title":"Reco electron momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
    "RecoElectrons_pt": {"name":"RecoElectrons_pt", "title":"Reco electron p_{T} [GeV]", "bin":50, "xmin":0, "xmax":100},
    "RecoElectrons_px": {"name":"RecoElectrons_px", "title":"Reco electron p_{x} [GeV]", "bin":50, "xmin":-100, "xmax":100},
    "RecoElectrons_py": {"name":"RecoElectrons_py", "title":"Reco electron p_{y} [GeV]", "bin":50, "xmin":-100, "xmax":100},
    "RecoElectrons_pz": {"name":"RecoElectrons_pz", "title":"Reco electron p_{z} [GeV]", "bin":50, "xmin":-200, "xmax":200},
    "RecoElectrons_eta": {"name":"RecoElectrons_eta", "title": "Reco electron #eta", "bin":50, "xmin":-5, "xmax":5},
    "RecoElectrons_theta": {"name":"RecoElectrons_theta", "title":"Reco electron #theta", "bin":100, "xmin":0, "xmax":3.2},
    "RecoElectrons_phi": {"name":"RecoElectrons_phi", "title":"Reco electron #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
    "RecoElectrons_charge": {"name":"RecoElectrons_charge", "title":"Reco electron charge", "bin":3, "xmin":-1.5, "xmax":1.5},


    #================ Reco Muons ================#
    "n_RecoMuons": {"name":"n_RecoMuons", "title":"Number of reconstructed muons", "bin":5, "xmin":-0.5, "xmax":4.5},
    "RecoMuons_e": {"name":"RecoMuons_e", "title":"Reco muon energy [GeV]", "bin":50, "xmin":0, "xmax":200},
    "RecoMuons_p": {"name":"RecoMuons_p", "title":"Reco muon momentum [GeV]", "bin":50, "xmin":0, "xmax":200},
    "RecoMuons_pt": {"name":"RecoMuons_pt", "title":"Reco muon p_{T} [GeV]", "bin":50, "xmin":0, "xmax":100},
    "RecoMuons_px": {"name":"RecoMuons_px", "title":"Reco muon p_{x} [GeV]", "bin":50, "xmin":-100, "xmax":100},
    "RecoMuons_py": {"name":"RecoMuons_py", "title":"Reco muon p_{y} [GeV]", "bin":50, "xmin":-100, "xmax":100},
    "RecoMuons_pz": {"name":"RecoMuons_pz", "title":"Reco muon p_{z} [GeV]", "bin":50, "xmin":-200, "xmax":200},
    "RecoMuons_eta": {"name":"RecoMuons_eta", "title":"Reco muon #eta", "bin":50, "xmin":-5, "xmax":5},
    "RecoMuons_theta": {"name":"RecoMuons_theta", "title":"Reco muon #theta", "bin":100, "xmin":0, "xmax":3.2},
    "RecoMuons_phi": {"name":"RecoMuons_phi", "title":"Reco muon #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
    "RecoMuons_charge": {"name":"RecoMuons_charge", "title":"Reco muon charge", "bin":3, "xmin":-1.5, "xmax":1.5},


    #================ Time variables ================#
    "GenStau_time": {"name":"GenStau_time", "title":"Stau time [ps]", "bin":50, "xmin":0, "xmax":10},
    "GenTau_time": {"name":"GenTau_time", "title":"Gen Tau time [ps]", "bin":50, "xmin":0, "xmax":10},
    "FSGenElectron_time": {"name":"FSGenElectron_time", "title":"FS electron time [ps]", "bin":50, "xmin":0, "xmax":10},
    "FSGenMuon_time": {"name":"FSGenMuon_time", "title":"FS muon time [ps]", "bin":50, "xmin":0, "xmax":10},
    "GenTau_status": {"name":"GenTau_status", "title":"Gen Tau status", "bin":20, "xmin":-0.5, "xmax":19.5},
    "GenTau_cTau": {"name":"GenTau_cTau", "title":"Gen Tau c#tau [mm]", "bin":150, "xmin":0, "xmax":4},
    # "GenStau_theta": {"name":"GenStau_theta", "title":"Gen Stau #theta", "bin":100, "xmin":0, "xmax":3.2}, # not in bg
    # "GenTau_theta": {"name":"GenTau_theta", "title":"Gen Tau #theta", "bin":100, "xmin":0, "xmax":3.2}, # not in bg


    #================ Primary track variables ================#
    "RecoedPrimaryTracks_charge": {"name":"RecoedPrimaryTracks_charge", "title":"Recoed primary tracks charge", "bin":3, "xmin":-1.5, "xmax":1.5},
    "sel_tracks_charge": {"name":"sel_tracks_charge", "title":"Selected tracks charge", "bin":3, "xmin":-1.5, "xmax":1.5},
    "RecoedPrimaryTracks_d0": {"name":"RecoedPrimaryTracks_d0", "title":"Recoed primary tracks d_{0} [mm]", "bin":100, "xmin":-100, "xmax":100},
    "RecoedPrimaryTracks_phi": {"name":"RecoedPrimaryTracks_phi", "title":"Recoed primary tracks #phi", "bin":64, "xmin":-3.2, "xmax":3.2},
    "RecoedPrimaryTracks_theta": {"name":"RecoedPrimaryTracks_theta", "title":"Recoed primary tracks #theta", "bin":30, "xmin":-4, "xmax":4},
    "RecoedPrimaryTracks_p": {"name":"RecoedPrimaryTracks_p", "title":"Reconstructed primary tracks: p [GeV]", "bin":100, "xmin":0, "xmax":365},
    "RecoedPrimaryTracks_pt": {"name":"RecoedPrimaryTracks_pt", "title":"Reconstructed primary tracks: p_{T} [GeV]", "bin":100, "xmin":0, "xmax":365},


    #================ Kinked candidates ================#
    "nKinkCandidates_passVeto": {"name":"nKinkCandidates_passVeto", "title":"Number of kink vertices passing the hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
    "KinkVertex_invMass": {"name":"KinkVertex_invMass", "title":"Invariant mass of kink vertex [GeV]", "bin":150, "xmin":0, "xmax":150},
    "nKinkVertices": {"name":"nKinkVertices", "title":"Number of kink vertices before hit veto", "bin":5, "xmin":-0.5, "xmax":4.5},
    "KinkVertex_ntracks": {"name":"KinkVertex_ntracks", "title":"Number of tracks in kink vertex", "bin":10, "xmin":-0.5, "xmax":9.5},
    "KinkAngle": {"name":"KinkAngle", "title":"Angle between the r_{PVKV} and P_{KV}", "bin":180, "xmin":0, "xmax":180},
    "KinkVertex_dxy": {"name":"KinkVertex_dxy", "title":"d_{xy} of kink vertex [mm]", "bin":100, "xmin":0, "xmax":2000},
    "KinkVertex_d3d": {"name":"KinkVertex_d3d", "title":"d_{3D} of kink vertex [mm]", "bin":100, "xmin":0, "xmax":2000},


    #================ Additional track / TOF variables ================# '
    "PV2V0Cos": {"name":"PV2V0Cos", "title":"Cosine angle between displaced vertex and primary vertex", "bin":100, "xmin":-1, "xmax":1},
    "RecoedPrimaryTrack_mass": {"name":"RecoedPrimaryTrack_mass", "title":"Mass of escaping tracks [GeV]", "bin":50, "xmin":100, "xmax":200},
    "TOF_RecoedPrimaryTracks": {"name":"TOF_RecoedPrimaryTracks", "title":"Time of flight of recoed primary tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
    # "Beta_RecoedPrimaryTracks": {"name":"Beta_RecoedPrimaryTracks", "title":"Beta of recoed primary tracks", "bin":100, "xmin":0, "xmax":1},
    "TOF_sel_tracks": {"name":"TOF_sel_tracks", "title":"Time of flight of selected tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
    "TOF_EscapingTracks": {"name":"TOF_EscapingTracks", "title":"Time of flight of escaping tracks [ps]", "bin":200, "xmin":0, "xmax":50000},
    "EscapingTracks_mass": {"name":"EscapingTracks_mass", "title":"Mass of escaping tracks [GeV]", "bin":50, "xmin":100, "xmax":200},
    "EscapingTracks_pt": {"name":"EscapingTracks_pt", "title":"Transverse momentum of escaping tracks [GeV]", "bin":100, "xmin":0, "xmax":365},
    "EscapingTracks_p": {"name":"EscapingTracks_p", "title":"Momentum of escaping tracks [GeV]", "bin":100, "xmin":0, "xmax":365},
    "EscapingTracks_theta": {"name":"EscapingTracks_theta", "title":"Theta of escaping tracks", "bin":100, "xmin":0, "xmax":3.2},
    "EscapingTracks_phi": {"name":"EscapingTracks_phi", "title":"Phi of escaping tracks", "bin":64, "xmin":-3.2, "xmax":3.2},
    "n_EscapingTracks": {"name":"n_EscapingTracks", "title":"Number of escaping tracks", "bin":25, "xmin":-0.5, "xmax":24.5},
}
