import ROOT
import math

input = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/stage1_combined_0603/FCCee_119_stau_2m_ctau_ecm_240.root"
output = "/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau//FCCee_119_stau_2m_ctau_ecm_240_REWEIGHTED_TO_20m.root"

ctau0 = 2.0   # generated
ctau1 = 20.0  # target

def lifetime_weight(ctau, ctau0, ctau1):
    return (ctau0 / ctau1) * math.exp(
        -ctau * (1/ctau1 - 1/ctau0)
    )

f_in = ROOT.TFile.Open(input)
tree_in = f_in.Get("events")

f_out = ROOT.TFile(output, "RECREATE")
tree_out = tree_in.CloneTree(0)  # clone structure only

# new branch
weight = ROOT.std.vector('float')()
tree_out.Branch("weight_ctau", weight)

for event in tree_in:
    weight.clear()

    ctau_vals = event.GenTau_cTau  # check units!

    for ctau in ctau_vals:
        w = lifetime_weight(ctau, ctau0, ctau1)
        weight.push_back(w)

    tree_out.Fill()

f_out.Write()
f_out.Close()
f_in.Close()