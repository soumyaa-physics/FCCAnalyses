# import ROOT

# f = ROOT.TFile.Open("/eos/user/s/svashish/MG5_aMC_v3_6_6/FCC_100stau_240com/Events/run_03/FCC_100stau_240com_000.root")
# t = f.Get("Delphes")

# f_out = ROOT.TFile.Open("FCC_100stau_240com_000_events.root", "RECREATE")
# t_clone = t.CloneTree()
# t_clone.SetName("events")
# t_clone.Write()

# df = ROOT.RDataFrame(t)

# # get column names
# cols = df.GetColumnNames()

# print("Columns in the dataframe:")
# for c in cols:
#     print(c)

    
# f_out.Close()
# f.Close()

import ROOT

# open ROOT file
f = ROOT.TFile.Open("FCC_100stau_240com_000_events.root")
t = f.Get("events")  # or 'events' if you renamed the tree

# create RDataFrame
df = ROOT.RDataFrame(t)

# get column names
cols = df.GetColumnNames()

print("Columns in the dataframe:")
for c in cols:
    print(c)
