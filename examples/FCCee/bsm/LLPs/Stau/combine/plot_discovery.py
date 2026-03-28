import json
import ROOT
import os
from array import array

colors = [
    ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange,
    ROOT.kCyan+1, ROOT.kPink+1, ROOT.kViolet, ROOT.kTeal+1
]

lifetimes = ["20cm", "50cm", "1m", "6m", "10m", "20m"]

SCALE_FACTOR = 0.001

cross_sections = {
    "100": 0.07735,
    "105": 0.05196,
    "110": 0.02923,
    "115": 0.01067,
    "118": 0.002752,
    "119": 0.0009792,
}

masses = ["100", "105", "110", "115", "118", "119"]
xvals  = array('d', [float(m) for m in masses])

limits_dir = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/limits"
json_path  = os.path.join(limits_dir, "limits_3years_2003.json")
plots_dir  = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/combine/plots_lumi3years_2003"
os.makedirs(plots_dir, exist_ok=True)

with open(json_path) as f:
    data = json.load(f)

def get_entry(mass, lt):
    key = f"{mass}_{lt}"
    if key not in data:
        raise KeyError(f"Key {key} not found in limits.json")
    return data[key]

L_current = 2.7e6  # 3 years 1 IP

# Create canvas for overlay
c = ROOT.TCanvas("c_overlay", "Lumi required for Discovery", 1200, 900)
c.SetLogy()
c.SetGrid()
c.SetTickx()
c.SetTicky()
frame = ROOT.TH2F("frame", "5 #sigma discovery reach; m_{#tilde{#tau}} [GeV]; Required Lumi [pb^{-1}]", 
                  100, 95, 125, 100, 1e1, 1e10)
frame.SetStats(0)
frame.SetMinimum(1e1)
frame.SetMaximum(1e14)
frame.Draw()
c.SetLogy()  # apply log scale after drawing
graphs = []

for idx, lifetime in enumerate(lifetimes):

    exp_y = array('d', [
        get_entry(m, lifetime)["exp0"] * SCALE_FACTOR  * cross_sections[m]
        for m in masses
    ])
    theory_y = array('d', [cross_sections[m] for m in masses])


    lumi_required = []
    for i, m in enumerate(masses):
        sigma_lim = exp_y[i]   # expected limit on sigma at L_current
        sigma_th  = theory_y[i]
        max_L5    = 1e10        # cap at 1e9 pb^-1 (upper plot boundary)
        min_L5    = 1e1       # floor — points below this are "already discoverable"

        if sigma_lim <= 0:
            lumi_required.append(max_L5)
            continue

        # Z at L_current ~ sigma_th / sigma_lim  (Asimov approximation)
        Z_current = sigma_th / sigma_lim
        L5 = L_current * (5.0 / Z_current) ** 2

        # print(f"  mass={m}, ct={lifetime}: sigma_th={sigma_th:.4g} pb, "
        #       f"sigma_lim={sigma_lim:.4g} pb, Z_now={Z_current:.3f}, L5={L5:.3g} pb^-1")

        L5 = max(min_L5, min(L5, max_L5))
        lumi_required.append(L5)

    print()
    lumi_arr = array('d', lumi_required)

    gr = ROOT.TGraph(len(masses), xvals, lumi_arr)
    gr.SetLineColor(colors[idx % len(colors)])
    gr.SetLineWidth(2)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.0)
    gr.SetMarkerColor(colors[idx % len(colors)])
    gr.Draw("PL same")
    graphs.append((gr, lifetime))

# Horizontal reference line at L_current
ref_line = ROOT.TLine(95, L_current, 125, L_current)
ref_line.SetLineColor(ROOT.kBlack)
ref_line.SetLineStyle(2)
ref_line.SetLineWidth(2)
ref_line.Draw()

ref_label = ROOT.TLatex()
ref_label.SetNDC(False)
ref_label.SetTextSize(0.028)
ref_label.SetTextColor(ROOT.kBlack)
ref_label.DrawLatex(96, L_current * 0.50, "L = 2.7ab^{-1} (3 years, 1 IP)")

# 1-year line
L_1y = 0.9e6
ref_line_1y = ROOT.TLine(95, L_1y, 125, L_1y)
ref_line_1y.SetLineColor(ROOT.kBlue)
ref_line_1y.SetLineStyle(2)
ref_line_1y.SetLineWidth(2)
ref_line_1y.Draw()

ref_label_1y = ROOT.TLatex()
ref_label_1y.SetNDC(False)
ref_label_1y.SetTextSize(0.028)
ref_label_1y.SetTextColor(ROOT.kBlue)
ref_label_1y.DrawLatex(96, L_1y * 0.50, "L = 0.9ab^{-1} (1 year, 1 IP)")

# 1-day line
L_1d = 6480
ref_line_1d = ROOT.TLine(95, L_1d, 125, L_1d)
ref_line_1d.SetLineColor(ROOT.kGreen+2)
ref_line_1d.SetLineStyle(2)
ref_line_1d.SetLineWidth(2)
ref_line_1d.Draw()

ref_label_1d = ROOT.TLatex()
ref_label_1d.SetNDC(False)
ref_label_1d.SetTextSize(0.028)
ref_label_1d.SetTextColor(ROOT.kGreen+2)
ref_label_1d.DrawLatex(96, L_1d * 0.50, "L = 6480 pb^{-1} (1 day, 1 IP)")

# Legend
legend = ROOT.TLegend(0.15, 0.60, 0.40, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.028)
for gr, lt in graphs:
    legend.AddEntry(gr, f"c#tau = {lt}", "lp")
legend.Draw()

# Labels
label = ROOT.TLatex()
label.SetNDC()
label.SetTextSize(0.035)
label.DrawLatex(0.68, 0.20, "FCC-ee Simulation")
label.DrawLatex(0.68, 0.15, "#sqrt{s} = 240 GeV")

lumi_label = ROOT.TLatex()
lumi_label.SetNDC()
lumi_label.SetTextSize(0.030)
lumi_label.SetTextAlign(31)
# lumi_label.DrawLatex(0.88, 0.12, "L_{ref} = 2.7#times10^{6} pb^{-1} (3 years, 1 IP)")


c.RedrawAxis()
c.Update()

# Save
outfile = os.path.join(plots_dir, "discovery_lumi_all_lifetimes_2003.png")
c.SaveAs(outfile)
print(f"Saved overlay plot: {outfile}")