import json
import ROOT
import os
from array import array

pastel_green  = ROOT.TColor.GetColor("#A6CEE3")
pastel_yellow = ROOT.TColor.GetColor("#FFF2A8")

lifetimes = ["20cm", "50cm", "1m", "2m", "3m", "4m", "6m", "10m", "20m"]

SCALE_FACTOR = 0.01 # taken from datacards

cross_sections = {
    "100": 0.07735,
    "105": 0.05196,
    "110": 0.02923,
    "115": 0.01067,
    "118": 0.002752,
    "119": 0.0009792,  
}

masses = ["100", "105", "110", "115", "118" , "119"]
xvals  = array('d', [float(m) for m in masses])

limits_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"
json_path  = os.path.join(limits_dir, "limits_1day_2003.json")

with open(json_path) as f:
    data = json.load(f)

# For each mass, look up the entry for this lifetime
def get_entry(mass, lt):
    key = f"{mass}_{lt}"
    if key not in data:
        raise KeyError(f"Key {key} not found in limits.json")
    return data[key]

for lifetime in lifetimes:

    def sigma_limit(mass, key): # rescale the limits with the scalefactor
        entry = get_entry(mass, lifetime)
        return entry[key] * SCALE_FACTOR * cross_sections[mass]

    n = len(masses)

    # obs_y  = array('d', [sigma_limit(m, "obs")   for m in masses])
    exp_y  = array('d', [sigma_limit(m, "exp0")  for m in masses])
    exp_m1 = array('d', [sigma_limit(m, "exp-1") for m in masses])
    exp_p1 = array('d', [sigma_limit(m, "exp+1") for m in masses])
    exp_m2 = array('d', [sigma_limit(m, "exp-2") for m in masses])
    exp_p2 = array('d', [sigma_limit(m, "exp+2") for m in masses])

    # Theory curve — one point per mass
    theory_y = array('d', [cross_sections[m] for m in masses])

    # Graph
    # gr_obs    = ROOT.TGraph(n, xvals, obs_y)
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

    # canvas
    c = ROOT.TCanvas("c", "Brazil plot vs mass", 1200, 900)
    c.SetLogy()

    gr_2sigma.SetTitle(
        f"Long-lived #tilde{{#tau}} in GMSB (c#tau = {lifetime});"
        "m_{#tilde{#tau}} [GeV];"
        "95% CL limit on #sigma [pb]"
    )

    gr_2sigma.Draw("AF") # fill the polygon- this gets the 2sig band
    gr_2sigma.GetXaxis().SetRangeUser(97, 120)
    gr_2sigma.GetYaxis().SetRangeUser(1e-7, 2.0)

    gr_1sigma.Draw("F same") # on same canvas get the 1sig band

    gr_exp.SetLineStyle(2)
    gr_exp.SetLineWidth(2)
    gr_exp.Draw("L same") # draw line between points in the same canvas

    # gr_obs.SetMarkerStyle(20)
    # gr_obs.SetLineWidth(3)
    # gr_obs.Draw("PL same") # markers and line

    gr_theory.SetLineColor(ROOT.kRed)
    gr_theory.SetLineWidth(3)
    gr_theory.SetLineStyle(1)
    gr_theory.Draw("L same") # line

    # legend
    legend = ROOT.TLegend(0.70, 0.70, 0.88, 0.88) #x1, y1, x2, y2
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # legend.AddEntry(gr_obs,    "Observed",            "pl")
    legend.AddEntry(gr_exp,    "Expected",            "l")
    legend.AddEntry(gr_1sigma, "Expected #pm1#sigma", "f")
    legend.AddEntry(gr_2sigma, "Expected #pm2#sigma", "f")
    legend.AddEntry(gr_theory, "Theory #sigma",       "l")
    legend.Draw()

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.035)
    label.DrawLatex(0.18, 0.85, "FCC-ee Simulation")
    label.DrawLatex(0.18, 0.80, "#sqrt{s} = 240 GeV")

    lumi_label = ROOT.TLatex()
    lumi_label.SetNDC()
    lumi_label.SetTextSize(0.035)
    lumi_label.SetTextAlign(31)  # right-aligned
    # lumi_label.DrawLatex(0.88, 0.15, "L = 1.08e7 pb^{-1} (Integrated Lumi)")
    # lumi_label.DrawLatex(0.88, 0.15, "L = 2.7e6 pb^{-1} (3 years, 1 IP)")
    # lumi_label.DrawLatex(0.88, 0.15, "L = 0.9e6 pb^{-1} (1 year, 1 IP)")
    lumi_label.DrawLatex(0.88, 0.15, "L = 6480 pb^{-1} (1 day, 1 IP)")


    # save the file
    plots_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/plots_lumi1day_2003"
    outfile = os.path.join(plots_dir, f"brazil_plot_lifetime_all_channels_{lifetime}.png")
    c.SaveAs(outfile)
    print(f"Saved: {outfile}")