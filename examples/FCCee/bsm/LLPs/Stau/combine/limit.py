import ROOT
import subprocess
import os
import argparse

# -------------------- Signal files (full paths) --------------------
signal_files = {
    "selNone": {
        "FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits_selNone_histo.root",
        "FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits_selNone_histo.root",
        "FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits_selNone_histo.root",
    },
    "semiLeptonic_KV": {
        "FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits_semiLeptonic_KV_histo.root",
        "FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits_semiLeptonic_KV_histo.root",
        "FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits_semiLeptonic_KV_histo.root",
    },
    "hadronic_KV": {
        "FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_20cm_ctau_ecm_240_delphes_nhits_hadronic_KV_histo.root",
        "FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_1p5m_ctau_ecm_240_delphes_nhits_hadronic_KV_histo.root",
        "FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits":
           "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/FCCee_110_stau_3m_ctau_ecm_240_delphes_nhits_hadronic_KV_histo.root",
    },
}

# -------------------- Background files (full paths) --------------------
background_files = {
    "selNone": {
        "WW":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_WW_ecm240_selNone_histo.root",
        "ZZ":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_ZZ_ecm240_selNone_histo.root",
        "ZH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/mgp8_ee_zh_ecm240_hbb_selNone_histo.root",
        "nunuH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_nuenueH_Htautau_ecm240_selNone_histo.root",
        "bbH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_bbH_Htautau_ecm240_selNone_histo.root",
    },
    "semiLeptonic_KV": {
        "WW":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_WW_ecm240_semiLeptonic_KV_histo.root",
        "ZZ":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_ZZ_ecm240_semiLeptonic_KV_histo.root",
        "ZH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/mgp8_ee_zh_ecm240_hbb_semiLeptonic_KV_histo.root",
        "nunuH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_nuenueH_Htautau_ecm240_semiLeptonic_KV_histo.root",
        "bbH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_bbH_Htautau_ecm240_semiLeptonic_KV_histo.root",
    },
    "hadronic_KV": {
        "WW":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_WW_ecm240_hadronic_KV_histo.root",
        "ZZ":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/p8_ee_ZZ_ecm240_hadronic_KV_histo.root",
        "ZH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/mgp8_ee_zh_ecm240_hbb_hadronic_KV_histo.root",
        "nunuH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_nuenueH_Htautau_ecm240_hadronic_KV_histo.root",
        "bbH":"/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/NEW_FINAL/wzp6_ee_bbH_Htautau_ecm240_hadronic_KV_histo.root",
    },
}

# -------------------- Selections --------------------
selections = {
    "selNone": "Before selection",
    "semiLeptonic_KV": "1 reco lepton + at least 1 KV",
    "hadronic_KV": "0 reco lepton + at least 1 KV",
}

hist_name = "KinkVertex_d3d"
json_limits = {}

# -------------------- Argument Parsing --------------------
parser = argparse.ArgumentParser(description="Create datacards for a given selection")
parser.add_argument("--sel", required=True, choices=list(selections.keys()), help="selection")
args = parser.parse_args()

sel = args.sel

# -------------------- Sum signal histograms --------------------
sig_hist_sum = None
for proc, file_path in signal_files[sel].items():
    sig_file = ROOT.TFile.Open(file_path, "READ")
    if not sig_file or sig_file.IsZombie():
        raise FileNotFoundError(f"Cannot open signal file: {file_path}")
    h = sig_file.Get(hist_name)
    if not h:
        sig_file.ls()
        raise RuntimeError(f"Histogram '{hist_name}' not found in {file_path}")
    if sig_hist_sum is None:
        sig_hist_sum = h.Clone("signal")
    else:
        sig_hist_sum.Add(h)
    sig_file.Close()

# -------------------- Sum background histograms --------------------
bkg_hist_sum = None
for proc, file_path in background_files[sel].items():
    bkg_file = ROOT.TFile.Open(file_path, "READ")
    if not bkg_file or bkg_file.IsZombie():
        raise FileNotFoundError(f"Cannot open background file: {file_path}")
    h = bkg_file.Get(hist_name)
    if not h:
        bkg_file.ls()
        raise RuntimeError(f"Histogram '{hist_name}' not found in {file_path}")
    if bkg_hist_sum is None:
        bkg_hist_sum = h.Clone("background")
    else:
        bkg_hist_sum.Add(h)
    bkg_file.Close()

# -------------------- Set zero-bin errors --------------------
for hist in [sig_hist_sum, bkg_hist_sum]:
    for i in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i) == 0:
            hist.SetBinError(i, 1.84)

# -------------------- Write datacard ROOT file --------------------
output_root = f"./{sel}_datacard.root"
out_file = ROOT.TFile(output_root, "RECREATE")
sig_hist_sum.Write()
bkg_hist_sum.Write()
data_hist = bkg_hist_sum.Clone("data_obs")
data_hist.Write()
out_file.Close()

# -------------------- Write text datacard --------------------
datacard_path = f"./{sel}_datacard.txt"
with open(datacard_path, "w") as f:
    f.write("# Datacard for combined histograms\n")
    f.write("imax 1  number of channels\n")
    f.write("jmax 1  number of backgrounds\n")
    f.write("kmax *  number of nuisance parameters\n")
    f.write(f"shapes * * {output_root} $PROCESS $PROCESS_$SYSTEMATIC\n")
    f.write("------------\n")
    f.write("bin bin1\n")
    f.write(f"observation {int(data_hist.Integral())}\n")
    f.write("------------\n")
    f.write("bin bin1 bin1\n")
    f.write("process signal background\n")
    f.write("process 0 1\n")
    f.write(f"rate {sig_hist_sum.Integral()} {bkg_hist_sum.Integral()}\n")
    f.write("------------\n")

# -------------------- Run combine --------------------
combine_cmd = ["combine", "-M", "AsymptoticLimits", datacard_path]
result = subprocess.run(combine_cmd, capture_output=True, text=True)
print(result.stdout)