#!/usr/bin/env python3
import json
import os

input_json = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final_0303/results.json"
output_dir = "datacards"
os.makedirs(output_dir, exist_ok=True)

channels = ["hadronic_KV"]

backgrounds = [
    "p8_ee_WW_ecm240",
    "p8_ee_ZZ_ecm240",
    "mgp8_ee_zh_ecm240_hbb",
    "wzp6_ee_nuenueH_Htautau_ecm240",
    "wzp6_ee_bbH_Htautau_ecm240"
]

signals = [
    'FCCee_100_stau_20cm_ctau_ecm_240',
    'FCCee_105_stau_20cm_ctau_ecm_240',
    'FCCee_110_stau_20cm_ctau_ecm_240',
    'FCCee_115_stau_20cm_ctau_ecm_240',
    'FCCee_100_stau_50cm_ctau_ecm_240',
    'FCCee_105_stau_50cm_ctau_ecm_240',
    'FCCee_110_stau_50cm_ctau_ecm_240',
    'FCCee_115_stau_50cm_ctau_ecm_240',
    'FCCee_100_stau_1m_ctau_ecm_240',
    'FCCee_105_stau_1m_ctau_ecm_240',
    'FCCee_110_stau_1m_ctau_ecm_240',
    'FCCee_115_stau_1m_ctau_ecm_240',
    'FCCee_100_stau_2m_ctau_ecm_240',
    'FCCee_105_stau_2m_ctau_ecm_240',
    'FCCee_110_stau_2m_ctau_ecm_240',
    'FCCee_115_stau_2m_ctau_ecm_240',
    'FCCee_100_stau_3m_ctau_ecm_240',
    'FCCee_105_stau_3m_ctau_ecm_240',
    'FCCee_110_stau_3m_ctau_ecm_240',
    'FCCee_115_stau_3m_ctau_ecm_240',
    'FCCee_100_stau_4m_ctau_ecm_240',
    'FCCee_105_stau_4m_ctau_ecm_240',
    'FCCee_110_stau_4m_ctau_ecm_240',
    'FCCee_115_stau_4m_ctau_ecm_240',
]

def make_label(s):
    parts = s.split("_")
    return f"sig_{parts[1]}_{parts[3]}"

with open(input_json) as f:
    results = json.load(f)

def get_rate(proc, chan):
    val = results.get(proc, {}).get(chan, {}).get("n_events_raw", 0)
    return max(float(val), 5)  # never exactly zero

# ── One datacard per signal ───────────────────────────────────────────────────
for sig_proc in signals:
    sig_label = make_label(sig_proc)

    # Observation = sum of backgrounds (background-only hypothesis)
    obs_per_chan = {}
    for chan in channels:
        obs_per_chan[chan] = sum(get_rate(b, chan) for b in backgrounds)

    # All processes for this datacard: backgrounds + ONE signal
    # Signal gets process index 0, backgrounds get 1,2,3,...
    all_proc_names  = backgrounds + [sig_label]
    all_proc_keys   = backgrounds + [sig_proc]
    # Signal = 0, backgrounds = 1..N  (matching the template exactly)
    all_proc_nums   = list(range(1, len(backgrounds) + 1)) + [0]

    n_cols = len(channels) * len(all_proc_names)

    fname = os.path.join(output_dir, f"datacard_{sig_label}.txt")
    with open(fname, "w") as dc:

        dc.write(f"# Datacard for signal: {sig_label}\n\n")
        dc.write(f"imax {len(channels)}  number of channels\n")
        dc.write("jmax *  number of backgrounds\n")
        dc.write("kmax *  number of nuisance parameters\n")
        dc.write("------------\n")

        # --- Observation ---
        dc.write("bin".ljust(15))
        dc.write("  ".join(channels) + "\n")
        dc.write("observation".ljust(15))
        dc.write("  ".join(f"{obs_per_chan[c]:.0f}" for c in channels) + "\n")
        dc.write("------------\n")

        # --- bin / process / process / rate ---
        # Each column = one (channel, process) pair
        col_chans = [chan for chan in channels for _ in all_proc_names]
        col_names = all_proc_names * len(channels)
        col_keys  = all_proc_keys  * len(channels)
        col_nums  = all_proc_nums  * len(channels)

        dc.write("bin".ljust(15)     + "  ".join(col_chans)             + "\n")
        dc.write("process".ljust(15) + "  ".join(col_names)             + "\n")
        dc.write("process".ljust(15) + "  ".join(str(x) for x in col_nums) + "\n")

        rates = []
        for key, chan in zip(col_keys, col_chans):
            rates.append(get_rate(key, chan))

        dc.write("rate".ljust(15)
                 + "  ".join(f"{r:.3f}" for r in rates) + "\n")
        dc.write("------------\n")

        # --- Systematics ---
        dc.write("lumi  lnN  "
                 + "  ".join(["1.02"] * n_cols) + "\n")
        dc.write("# 2% lumi uncertainty on all processes\n")

    print(f"Written: {fname}")

# # ── Shell script to run all limits ───────────────────────────────────────────
# run_script = os.path.join(output_dir, "run_limits.sh")
# with open(run_script, "w") as sh:
#     sh.write("#!/bin/bash\n")
#     sh.write("mkdir -p results\n\n")
#     for sig_proc in signals:
#         label = make_label(sig_proc)
#         dc    = f"datacard_{label}.txt"
#         sh.write(f"echo 'Running limit for {label}'\n")
#         sh.write(
#             f"combine -M AsymptoticLimits {dc} "
#             f"--name _{label} "
#             f"-m 125 "           # placeholder mass; change as needed
#             f"--rMin 0 --rMax 200 "
#             f"> results/log_{label}.txt 2>&1\n\n"
#         )
# sh.write("echo 'All done. Results in higgsCombine_*.root'\n")
# os.chmod(run_script, 0o755)
# print(f"\nRun all limits with:  cd {output_dir} && ./run_limits.sh")