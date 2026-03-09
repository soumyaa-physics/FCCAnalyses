import ROOT
import glob

files = glob.glob("output/condor_0603/wzp6_ee_bbH_Htautau_ecm240/*.root")

total = 0

for f in files:
    tf = ROOT.TFile.Open(f)
    p = tf.Get("eventsProcessed")
    total += p.GetVal()

print("Total eventsProcessed =", total)