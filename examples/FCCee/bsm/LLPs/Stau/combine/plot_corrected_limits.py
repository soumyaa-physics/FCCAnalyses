import json
import ROOT
import os
from array import array

pastel_green  = ROOT.TColor.GetColor("#A6CEE3")
pastel_yellow = ROOT.TColor.GetColor("#FFF2A8")

lifetime     = "4m"       # choose: 20cm, 50cm, 1m, 2m, 3m, 4m
SCALE_FACTOR = 0.1

cross_sections = {
    "100": 0.07735,
    "105": 0.05196,
    "110": 0.02923,
    "115": 0.01067,
}

masses = ["100", "105", "110", "115"]
xvals  = array('d', [float(m) for m in masses])

limits_dir = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"
json_path  = os.path.join(limits_dir, "limits.json")

with open(json_path) as f:
    data = json.load(f)

# For each mass, look up the entry for this lifetime
def get_entry(mass, lt):
    key = f"{mass}_{lt}"
    if key not in data:
        raise KeyError(f"Key {key} not found in limits.json")
    return data[key]

def sigma_limit(mass, key):
    entry = get_entry(mass, lifetime)
    return entry[key] * SCALE_FACTOR * cross_sections[mass]

n = len(masses)

obs_y  = array('d', [sigma_limit(m, "obs")   for m in masses])
exp_y  = array('d', [sigma_limit(m, "exp0")  for m in masses])
exp_m1 = array('d', [sigma_limit(m, "exp-1") for m in masses])
exp_p1 = array('d', [sigma_limit(m, "exp+1") for m in masses])
exp_m2 = array('d', [sigma_limit(m, "exp-2") for m in masses])
exp_p2 = array('d', [sigma_limit(m, "exp+2") for m in masses])

# Theory curve — one point per mass
theory_y = array('d', [cross_sections[m] for m in masses])

# ── Graphs ─────────────────────────────────────────────────────────────────
gr_obs    = ROOT.TGraph(n, xvals, obs_y)
gr_exp    = ROOT.TGraph(n, xvals, exp_y)
gr_theory = ROOT.TGraph(n, xvals, theory_y)

gr_1sigma = ROOT.TGraph(2*n)
for i in range(n):
    gr_1sigma.SetPoint(i,       xvals[i], exp_p1[i])
    gr_1sigma.SetPoint(2*n-i-1, xvals[i], exp_m1[i])
gr_1sigma.SetFillColor(pastel_yellow)

gr_2sigma = ROOT.TGraph(2*n)
for i in range(n):
    gr_2sigma.SetPoint(i,       xvals[i], exp_p2[i])
    gr_2sigma.SetPoint(2*n-i-1, xvals[i], exp_m2[i])
gr_2sigma.SetFillColor(pastel_green)

# ── Canvas ─────────────────────────────────────────────────────────────────
c = ROOT.TCanvas("c", "Brazil plot vs mass", 1200, 900)
c.SetLogy()

gr_2sigma.SetTitle(
    f"Stau GMSB, c#tau = {lifetime};"
    "Stau mass [GeV];"
    "95% CL Limit on #sigma [pb]"
)

gr_2sigma.Draw("AF")
gr_2sigma.GetXaxis().SetRangeUser(97, 118)
gr_2sigma.GetYaxis().SetRangeUser(1e-4, 2.0)

gr_1sigma.Draw("F same")

gr_exp.SetLineStyle(2)
gr_exp.SetLineWidth(2)
gr_exp.Draw("L same")

gr_obs.SetMarkerStyle(20)
gr_obs.SetLineWidth(2)
gr_obs.Draw("PL same")

gr_theory.SetLineColor(ROOT.kRed)
gr_theory.SetLineWidth(2)
gr_theory.SetLineStyle(1)
gr_theory.Draw("L same")

# ── Legend ─────────────────────────────────────────────────────────────────
legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(gr_obs,    "Observed",            "pl")
legend.AddEntry(gr_exp,    "Expected",            "l")
legend.AddEntry(gr_1sigma, "Expected #pm1#sigma", "f")
legend.AddEntry(gr_2sigma, "Expected #pm2#sigma", "f")
legend.AddEntry(gr_theory, "Theory #sigma",       "l")
legend.Draw()

c.SaveAs(f"brazil_plot_lifetime_all_channels{lifetime}.png")
print(f"Saved: brazil_plot_lifetime_{lifetime}.png")