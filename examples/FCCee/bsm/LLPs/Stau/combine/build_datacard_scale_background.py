#!/usr/bin/env python3
import json
import os
import subprocess

# save here to access for plotting
json_limits = {}

input_json = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/final_0303/results.json"
output_dir = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/datacards_scaledbg"
limits_dir = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(limits_dir, exist_ok=True)

channels = ["hadronic_KV","semiLeptonic_KV", "semiLep_DV", "hadronic_DV"]

backgrounds = [
    "p8_ee_WW_ecm240",
    "p8_ee_ZZ_ecm240",
    "mgp8_ee_zh_ecm240_hbb",
    "wzp6_ee_nuenueH_Htautau_ecm240",
    "wzp6_ee_bbH_Htautau_ecm240"
]

mc_info = {
    "p8_ee_WW_ecm240": {"N_MC": 20000, "sigma": 16.4385},
    "p8_ee_ZZ_ecm240": {"N_MC": 6000, "sigma": 1.35899},
    "mgp8_ee_zh_ecm240_hbb": {"N_MC": 10000, "sigma": 1},
    "wzp6_ee_nuenueH_Htautau_ecm240": {"N_MC": 4000, "sigma": 0.001219},
    "wzp6_ee_bbH_Htautau_ecm240": {"N_MC": 12000, "sigma": 0.00188},
}

signals = [
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

def extract_mass_lifetime(sig_name):
    parts = sig_name.split("_")
    return parts[1], parts[3]  # mass, lifetime

with open(input_json) as f:
    results = json.load(f)

L = 1.08e7  # pb⁻¹
# SCALE_THRESHOLD = 0.1 # can be changed
SCALE_FACTOR = 10

def get_rate(proc, chan):
    if proc in mc_info:  # background
        info = mc_info[proc]
        raw = float(results.get(proc, {}).get(chan, {}).get("n_events_raw", 0))

        if raw == 0:
            scaled = (3 / info["N_MC"]) * info["sigma"] * L # 95% confidence limit
        else:
            scaled = raw * (info["sigma"] * L / info["N_MC"])

    else:  # signal
        raw = float(results.get(proc, {}).get(chan, {}).get("n_events_raw", 0))
        scaled = raw

    return max(scaled, 1e-9)


for sig_proc in signals:

    sig_label = make_label(sig_proc)

    obs_per_chan = {}
    bkg_totals = {}

    for ch in channels:
        total_bkg = sum(get_rate(b, ch) for b in backgrounds)
        obs_per_chan[ch] = total_bkg
        bkg_totals[ch] = total_bkg

    # Decide scaling
    scale_value = 1.0
    for ch in channels:
        sig_rate = get_rate(sig_proc, ch)
        bkg_total = bkg_totals[ch]

        if bkg_total > 0:
            scale_value = SCALE_FACTOR
            print(f"[INFO] Scaling {sig_label} by 1/{SCALE_FACTOR}") # for book-keeping, in order to rescale it
            break

    all_proc_names = backgrounds + [sig_label]
    all_proc_keys = backgrounds + [sig_proc]
    all_proc_nums = list(range(1, len(backgrounds) + 1)) + [0]

    n_cols = len(channels) * len(all_proc_names)
    fname = os.path.join(output_dir, f"datacard_{sig_label}.txt")

    with open(fname, "w") as dc:

        dc.write(f"# Datacard for signal: {sig_label}\n\n")
        dc.write(f"imax {len(channels)}\n")
        dc.write("jmax *\n")
        dc.write("kmax *\n")
        dc.write("------------\n")

        dc.write("bin".ljust(15) + "  ".join(channels) + "\n")
        dc.write("observation".ljust(15) +
                 "  ".join(f"{obs_per_chan[c]:.0f}" for c in channels) + "\n")
        dc.write("------------\n")

        col_chans = [c for c in channels for _ in all_proc_names]
        col_nums = all_proc_nums * len(channels)

        dc.write("bin".ljust(15) + "  ".join(col_chans) + "\n")
        dc.write("process".ljust(15) + "  ".join(all_proc_names * len(channels)) + "\n")
        dc.write("process".ljust(15) + "  ".join(str(x) for x in col_nums) + "\n")

        rates = []
        for key, ch in zip(all_proc_keys * len(channels), col_chans):
            r = get_rate(key, ch)
            if key == sig_proc:
                r = r / scale_value
            rates.append(r)

        dc.write("rate".ljust(15) +
                 "  ".join(f"{r:.6f}" for r in rates) + "\n")
        dc.write("------------\n")

        dc.write("lumi  lnN  " + "  ".join(["1.02"] * n_cols) + "\n")

    # Run combine
    result = subprocess.run(
        ["combine", "-M", "AsymptoticLimits", fname],
        capture_output=True,
        text=True
    )

    output_text = result.stdout

    exp0 = exp_m1 = exp_p1 = exp_m2 = exp_p2 = obs = None

    for line in output_text.split("\n"):
        if "Expected 50.0%" in line: exp0 = float(line.split()[-1])
        elif "Expected 16.0%" in line: exp_m1 = float(line.split()[-1])
        elif "Expected 84.0%" in line: exp_p1 = float(line.split()[-1])
        elif "Expected  2.5%" in line: exp_m2 = float(line.split()[-1])
        elif "Expected 97.5%" in line: exp_p2 = float(line.split()[-1])
        elif "Observed Limit:" in line: obs = float(line.split()[-1])

    mass, lifetime = extract_mass_lifetime(sig_proc)
    key = f"{mass}_{lifetime}"

    json_limits[key] = {
        "exp-2": exp_m2 or 0,
        "exp-1": exp_m1 or 0,
        "exp0":  exp0  or 0,
        "exp+1": exp_p1 or 0,
        "exp+2": exp_p2 or 0,
        "obs":   obs   or 0,
    }


# Save limits
json_path = os.path.join(limits_dir, "limits.json")

with open(json_path, "w") as f:
    json.dump(json_limits, f, indent=2)

print(f"All limits saved to {json_path}")