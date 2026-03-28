import json
import ROOT
import os
import numpy as np

lifetime_labels = ["20cm", "50cm", "1m", "2m", "3m", "4m"]
masses = np.array([100, 105, 110, 115])
cross_sections = {
    "100": 0.07735,
    "105": 0.05196,
    "110": 0.02923,
    "115": 0.01067,
}
SCALE_FACTOR = 0.001
limits_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"
json_path = os.path.join(limits_dir, "limits.json")

with open(json_path) as f:
    data = json.load(f)

def parse_lifetime(s):
    if "cm" in s:
        return float(s.replace("cm", ""))
    if "m" in s:
        return float(s.replace("m", "")) * 100
    raise ValueError(f"Unknown lifetime format: {s}")

def get_entry(mass, lt):
    key = f"{mass}_{lt}"
    if key not in data:
        raise KeyError(f"{key} not found in limits.json")
    return data[key]

lifetimes = np.array([parse_lifetime(x) for x in lifetime_labels])

# Fixed bin edges so each (mass, lifetime) point sits at a bin centre
mass_edges = np.array([97.5, 102.5, 107.5, 112.5, 117.5], dtype=float)
lt_edges   = np.array([10., 35., 75., 150., 250., 350., 450.], dtype=float)

def make_xs_hist(name, title, key):
    h = ROOT.TH2F(name, title,
                  len(mass_edges) - 1, mass_edges,
                  len(lt_edges)   - 1, lt_edges)
    for i, lt in enumerate(lifetime_labels):
        for j, m in enumerate(masses):
            entry  = get_entry(str(m), lt)
            xs_lim = entry[key] * SCALE_FACTOR * cross_sections[str(m)]
            h.SetBinContent(j + 1, i + 1, xs_lim)
    return h

def make_ratio_hist(name, key):
    h = ROOT.TH2F(name, "",
                  len(mass_edges) - 1, mass_edges,
                  len(lt_edges)   - 1, lt_edges)
    for i, lt in enumerate(lifetime_labels):
        for j, m in enumerate(masses):
            entry = get_entry(str(m), lt)
            limit_val = entry[key] * SCALE_FACTOR
            theory_val = cross_sections[str(m)]   
            ratio = limit_val/theory_val 
            h.SetBinContent(j + 1, i + 1, ratio)
    return h

h_obs_xs   = make_xs_hist("h_obs_xs", "Observed 95% CL upper limit on #sigma [pb]", "obs")
# get the ratio with theory
h_obs_r    = make_ratio_hist("h_obs_r",    "obs")
h_exp_r    = make_ratio_hist("h_exp_r",    "exp0")
h_exp_m1_r = make_ratio_hist("h_exp_m1_r", "exp-1")
h_exp_p1_r = make_ratio_hist("h_exp_p1_r", "exp+1")
h_exp_m2_r = make_ratio_hist("h_exp_m2_r", "exp-2")
h_exp_p2_r = make_ratio_hist("h_exp_p2_r", "exp+2")

# # Keep ROOT objects alive (prevent garbage collection)
_refs = []

# Draw contour at ratio = contour_val
def draw_contour(h, contour_val, color, style, width):
    hc = h.Clone(h.GetName() + "_cont")
    hc.SetContour(1, np.array([contour_val]))
    hc.SetLineColor(color)
    hc.SetLineStyle(style)
    hc.SetLineWidth(width)
    hc.Draw("CONT3 SAME")
    _refs.append(hc)
    return hc

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBeach)
ROOT.gStyle.SetNumberContours(100)

c = ROOT.TCanvas("c", "Stau sensitivity", 950, 720)
c.SetLogy()
c.SetLogz()
c.SetRightMargin(0.16)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)

h_obs_xs.GetXaxis().SetTitle("m_{#tilde{#tau}} [GeV]")
h_obs_xs.GetYaxis().SetTitle("c#tau_{0} [cm]")
h_obs_xs.GetZaxis().SetTitle("95% CL upper limit on #sigma [pb]")
h_obs_xs.GetXaxis().SetTitleOffset(1.1)
h_obs_xs.GetYaxis().SetTitleOffset(1.3)
h_obs_xs.GetZaxis().SetTitleOffset(1.4)

# Draw colour map
h_obs_xs.Draw("COLZ")

CONTOUR = 1e-2

draw_contour(h_exp_m2_r, CONTOUR, ROOT.kRed,   style=3, width=2)
draw_contour(h_exp_p2_r, CONTOUR, ROOT.kRed,   style=3, width=2)
draw_contour(h_exp_m1_r, CONTOUR, ROOT.kRed,   style=1, width=2)
draw_contour(h_exp_p1_r, CONTOUR, ROOT.kRed,   style=1, width=2)
# draw_contour(h_exp_r,    CONTOUR, ROOT.kRed,   style=1, width=3)
draw_contour(h_obs_r,    CONTOUR, ROOT.kBlack, style=1, width=3)

leg = ROOT.TLegend(0.14, 0.67, 0.54, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.028)

d_obs = ROOT.TLine(); d_obs.SetLineColor(ROOT.kBlack); d_obs.SetLineStyle(1); d_obs.SetLineWidth(3)
# d_exp = ROOT.TLine(); d_exp.SetLineColor(ROOT.kRed);   d_exp.SetLineStyle(1); d_exp.SetLineWidth(3)
d_1s  = ROOT.TLine(); d_1s.SetLineColor(ROOT.kRed);    d_1s.SetLineStyle(1);  d_1s.SetLineWidth(2)
d_2s  = ROOT.TLine(); d_2s.SetLineColor(ROOT.kRed);    d_2s.SetLineStyle(3);  d_2s.SetLineWidth(2)

leg.AddEntry(d_obs, "Observed (#sigma^{95% CL}/#sigma_{theory} = 10^{-2})", "l")
# leg.AddEntry(d_exp, "Expected (#sigma^{95% CL}/#sigma_{theory} = 10^{-2})", "l")
leg.AddEntry(d_1s,  "Expected #pm1#sigma",                                   "l")
leg.AddEntry(d_2s,  "Expected #pm2#sigma",                                   "l")
leg.Draw()
_refs += [leg, d_obs, d_1s, d_2s]


label = ROOT.TLatex()
label.SetNDC()
label.SetTextSize(0.035)
label.DrawLatex(0.60, 0.25, "FCC-ee Simulation")
label.DrawLatex(0.60, 0.20, "#sqrt{s} = 240 GeV")

c.Update()
c.SaveAs("stau_exclusion_limit.png")
c.SaveAs("stau_exclusion_limit.pdf")
print("Saved stau_exclusion_limit.png / .pdf")