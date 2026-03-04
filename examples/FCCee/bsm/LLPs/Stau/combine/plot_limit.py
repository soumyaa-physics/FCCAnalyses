import json
import ROOT
import os
from array import array

pastel_green  = ROOT.TColor.GetColor("#A6CEE3")
pastel_yellow = ROOT.TColor.GetColor("#FFF2A8")

mass_point   = "110"
SCALE_FACTOR = 10

cross_sections = {
    "100": 0.07735,
    "105": 0.05196,
    "110": 0.02923,
    "115": 0.01067,
}

limits_dir = "/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"
json_path  = os.path.join(limits_dir, "limits.json")

with open(json_path) as f:
    data = json.load(f)

filtered = {k: v for k, v in data.items() if k.startswith(mass_point + "_")}

def lifetime_to_float(s):
    s = s.strip()
    if s.endswith("cm"):
        return float(s.replace("cm", "")) / 100.0
    elif s.endswith("m"):
        return float(s.replace("m", ""))
    return float(s)

lifetimes = sorted(
    [k.split("_")[1] for k in filtered.keys()],
    key=lifetime_to_float
)

xvals        = array('d', [lifetime_to_float(x) for x in lifetimes])
sigma_theory = cross_sections[mass_point]

def sigma_limit(entry, key):
    return entry[key] * SCALE_FACTOR * sigma_theory

n = len(lifetimes)

# ── All arrays use sigma_limit ─────────────────────────────────────────────
obs_y  = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "obs")   for lt in lifetimes])
exp_y  = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "exp0")  for lt in lifetimes])
exp_m1 = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "exp-1") for lt in lifetimes])
exp_p1 = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "exp+1") for lt in lifetimes])
exp_m2 = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "exp-2") for lt in lifetimes])
exp_p2 = array('d', [sigma_limit(filtered[f"{mass_point}_{lt}"], "exp+2") for lt in lifetimes])

# ── Graphs ─────────────────────────────────────────────────────────────────
gr_obs = ROOT.TGraph(n, xvals, obs_y)
gr_exp = ROOT.TGraph(n, xvals, exp_y)

gr_1sigma = ROOT.TGraph(2*n)
for i in range(n):
    gr_1sigma.SetPoint(i,         xvals[i], exp_p1[i])
    gr_1sigma.SetPoint(2*n-i-1,   xvals[i], exp_m1[i])
gr_1sigma.SetFillColor(pastel_yellow)

gr_2sigma = ROOT.TGraph(2*n)
for i in range(n):
    gr_2sigma.SetPoint(i,         xvals[i], exp_p2[i])
    gr_2sigma.SetPoint(2*n-i-1,   xvals[i], exp_m2[i])
gr_2sigma.SetFillColor(pastel_green)

# ── Canvas ─────────────────────────────────────────────────────────────────
c = ROOT.TCanvas("c", "Brazil plot", 1200, 900)
c.SetLogy()

gr_2sigma.SetTitle(
    f"Stau GMSB, m = {mass_point} GeV;"
    "Lifetime c#tau [m];"
    "95% CL Limit on #sigma [pb]"
)

gr_2sigma.Draw("AF")
gr_2sigma.GetXaxis().SetRangeUser(xvals[0] - 0.1, xvals[-1] + 0.5)
gr_2sigma.GetYaxis().SetRangeUser(1e-3, 2.0)

gr_1sigma.Draw("F same")

gr_exp.SetLineStyle(2)
gr_exp.SetLineWidth(2)
gr_exp.Draw("L same")

gr_obs.SetMarkerStyle(20)
gr_obs.SetLineWidth(2)
gr_obs.Draw("PL same")

theory_line = ROOT.TLine(xvals[0], sigma_theory, xvals[-1], sigma_theory)
theory_line.SetLineColor(ROOT.kRed)
theory_line.SetLineWidth(2)
theory_line.Draw()

legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(gr_obs,      "Observed",            "pl")
legend.AddEntry(gr_exp,      "Expected",            "l")
legend.AddEntry(gr_1sigma,   "Expected #pm1#sigma", "f")
legend.AddEntry(gr_2sigma,   "Expected #pm2#sigma", "f")
legend.AddEntry(theory_line, "Theory #sigma",       "l")
legend.Draw()

c.SaveAs(f"brazil_plot_mass_{mass_point}.png")
print("Saved.")