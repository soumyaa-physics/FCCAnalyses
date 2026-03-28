#!/usr/bin/env python3
import json
import os
import subprocess

json_limits = {}
# final_0303 has a different hadronic kV selection with the lepton veto included
input_json = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final_1day_2003/results.json"
output_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/datacards_all/datacards_2003_1day"
limits_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(limits_dir, exist_ok=True)

# lumi = 1.08e7  # pb^-1
lumi = 6480 # lumi in pb-1 for 1 day at 1 IP
# lumi = 0.9e6 # per year at IP in pb-1
# lumi = 2.7e6 # per IP at 3 years
SCALE_FACTOR = 100     # signal always divided by this in the datacard

channels = ["hadronic_KV", "semiLeptonic_KV", "semiLep_DV", "hadronic_DV"]

backgrounds = [
    "p8_ee_WW_ecm240",
    "p8_ee_ZZ_ecm240",
    "mgp8_ee_zh_ecm240_hbb",
    "wzp6_ee_nuenueH_Htautau_ecm240",
    "wzp6_ee_bbH_Htautau_ecm240"
]

mc_info = {
    "p8_ee_WW_ecm240": {"N_MC": 3800000, "sigma": 16.4385}, # 0.01 of all events
    "p8_ee_ZZ_ecm240": {"N_MC": 5700000, "sigma": 1.35899}, # 0.1 all events : 5,616,209.3
    # now processed all of these:
    "mgp8_ee_zh_ecm240_hbb": {"N_MC": 100000, "sigma": 1},
    "wzp6_ee_nuenueH_Htautau_ecm240": {"N_MC": 1200000, "sigma": 0.001219},
    "wzp6_ee_bbH_Htautau_ecm240": {"N_MC": 400000, "sigma": 0.00188},
}

# Signal MC - should also be scaled
signal_mc_info = {
    # 20 cm
    "FCCee_100_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_20cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 50 cm
    "FCCee_100_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_50cm_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 1 m
    "FCCee_100_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_1m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 2 m
    "FCCee_100_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_2m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 3 m
    "FCCee_100_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_3m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 4 m
    "FCCee_100_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_4m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 6 m
    "FCCee_100_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_6m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 10 m
    "FCCee_100_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_10m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},

    # 20 m
    "FCCee_100_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.07735},
    "FCCee_105_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.05196},
    "FCCee_110_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.02923},
    "FCCee_115_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.01067},
    "FCCee_118_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.002752},
    "FCCee_119_stau_20m_ctau_ecm_240": {"N_MC": 100000, "sigma": 0.0009792},
}

# # scaling all signals with lifetime > 2m : pythia issue
# scale_large_ctau = {
#     "20cm" : 1.0,
#     "50cm" : 1.0,
#     "1m" : 1.0,
#     "2m": 0.9996,
#     "3m": 0.9944,
#     "4m": 0.9795,
#     "6m": 0.9251,
#     "10m": 0.7888,
#     "20m": 0.5404
# }

signals = list(signal_mc_info.keys())

def make_label(s):
    parts = s.split("_")
    return f"sig_{parts[1]}_{parts[3]}"

def extract_mass_lifetime(sig_name):
    parts = sig_name.split("_")
    return parts[1], parts[3]

with open(input_json) as f:
    results = json.load(f)

def get_rate(proc, chan):
    """
    Return expected event yield for proc in chan.
    Both backgrounds and signals use: raw * (sigma * L / N_MC)
    For 0 events left: use 95% CL Poisson upper limit = 3/N_MC * sigma * L
    Signal is additionally divided by SCALE_FACTOR.
    """
    raw = float(results.get(proc, {}).get(chan, {}).get("n_events_raw", 0))

    if proc in mc_info:  # background
        info = mc_info[proc]
        if raw == 0:
            scaled = (3.0 / info["N_MC"]) * info["sigma"] * lumi
        else:
            scaled = raw * (info["sigma"] * lumi / info["N_MC"])
            # print(scaled)

    else:  # signal
        info = signal_mc_info[proc]
        # get the lifetime
        # _, lifetime = extract_mass_lifetime(proc)
        # scale_factor_ctau = scale_large_ctau.get(lifetime, 1.0)

        if raw == 0:
            scaled = 0
            # scaled = (3.0 / info["N_MC"]) * info["sigma"] * lumi / SCALE_FACTOR
        else:
            scaled_lumi = raw * (info["sigma"] * lumi / info["N_MC"])
            scaled = scaled_lumi/ SCALE_FACTOR
            
    return max(scaled, 1e-9)


for sig_proc in signals:

    sig_label = make_label(sig_proc)
    print(f"\n[INFO] Processing {sig_label}  (signal divided by {SCALE_FACTOR})")

    obs_per_chan = {ch: sum(get_rate(b, ch) for b in backgrounds) for ch in channels}

    all_proc_names = backgrounds + [sig_label]
    all_proc_keys  = backgrounds + [sig_proc]
    all_proc_nums  = list(range(1, len(backgrounds) + 1)) + [0]

    n_cols = len(channels) * len(all_proc_names)
    fname  = os.path.join(output_dir, f"datacard_{sig_label}.txt")

    with open(fname, "w") as dc:
        dc.write(f"# Datacard for signal: {sig_label}\n")
        dc.write(f"# Signal scaled down by {SCALE_FACTOR} -- multiply r by {SCALE_FACTOR} to get true signal strength\n\n")
        dc.write(f"imax {len(channels)}\n")
        dc.write("jmax *\n")
        dc.write("kmax *\n")
        dc.write("------------\n")

        dc.write("bin".ljust(15)         + "  ".join(channels) + "\n")
        dc.write("observation".ljust(15) + "  ".join(f"{obs_per_chan[c]:.0f}" for c in channels) + "\n")
        dc.write("------------\n")

        col_chans = [c for c in channels for _ in all_proc_names]
        col_keys  = all_proc_keys * len(channels)
        col_names = all_proc_names * len(channels)
        col_nums  = all_proc_nums  * len(channels)

        dc.write("bin".ljust(15)     + "  ".join(col_chans)                + "\n")
        dc.write("process".ljust(15) + "  ".join(col_names)                + "\n")
        dc.write("process".ljust(15) + "  ".join(str(x) for x in col_nums) + "\n")

        rates = [get_rate(k, c) for k, c in zip(col_keys, col_chans)]
        dc.write("rate".ljust(15) + "  ".join(f"{r:.6f}" for r in rates) + "\n")
        dc.write("------------\n")
        dc.write("lumi  lnN  " + "  ".join(["1.02"] * n_cols) + "\n") # considered 2% uncertainity on lumi- integrated
        # lumi uncertainity per day: 3\%

    # Run Combine
    result = subprocess.run(
        ["combine", "-M", "AsymptoticLimits", fname],
        capture_output=True, text=True
    )

    exp0 = exp_m1 = exp_p1 = exp_m2 = exp_p2 = obs = None
    for line in result.stdout.split("\n"):
        if   "Expected 50.0%" in line: exp0   = float(line.split()[-1])
        elif "Expected 16.0%" in line: exp_m1 = float(line.split()[-1])
        elif "Expected 84.0%" in line: exp_p1 = float(line.split()[-1])
        elif "Expected  2.5%" in line: exp_m2 = float(line.split()[-1])
        elif "Expected 97.5%" in line: exp_p2 = float(line.split()[-1])
        elif "Observed Limit:" in line: obs   = float(line.split()[-1])

    if exp0 is None:
        print(f"[WARNING] Combine failed for {sig_label}. stderr:\n{result.stderr}")

    mass, lifetime = extract_mass_lifetime(sig_proc)
    json_limits[f"{mass}_{lifetime}"] = {
        "exp-2": exp_m2 or 0,
        "exp-1": exp_m1 or 0,
        "exp0":  exp0   or 0,
        "exp+1": exp_p1 or 0,
        "exp+2": exp_p2 or 0,
        "obs":   obs    or 0,
    }

json_path = os.path.join(limits_dir, "limits_1day_2003.json") 
with open(json_path, "w") as f:
    json.dump(json_limits, f, indent=2)

print(f"\nAll limits saved to {json_path}")