# import ROOT
# import glob

# files = glob.glob("output/condor_0603/wzp6_ee_bbH_Htautau_ecm240/*.root")

# total = 0

# for f in files:
#     tf = ROOT.TFile.Open(f)
#     p = tf.Get("eventsProcessed")
#     total += p.GetVal()

# print("Total eventsProcessed =", total)

import ROOT
import math

m_stau = 110.0      
E_stau = 240.0 / 2  
cTau_m_list = [ 1,2,3,4, 6, 10, 20]  

p = math.sqrt(E_stau**2 - m_stau**2)
beta_gamma = p/m_stau
print(f"For stau mass: {m_stau} E = {E_stau:.2f} GeV, p = {p:.2f} GeV")
print(f"beta*gamma : {beta_gamma}")

for ctau in cTau_m_list:
    L_lab = beta_gamma * ctau  # average decay length (true) in meters
    # as pythia dumps everything at average ctau of 202cm: generator level decay length of the stau
    # L = 2.02  # meters
    # prob_decay = math.exp(-L/L_lab)
    print(f"cTau = {ctau} m -> Average lab-frame decay length: {L_lab:.2f} m")
    # print(f"Probability of decay within 202cm: {1- prob_decay:.4f}")

# cTau_m_list_new = [ 10, 20]  

# for ctau_new in cTau_m_list_new:
#     L_lab = beta_gamma * ctau  # reference
#     weight = 2/ctau_new * math.exp(-)

def average_ctau(root_file, tree_name="events", branch_name="decayLengthTau"):
    f = ROOT.TFile.Open(root_file)
    tree = f.Get(tree_name)

    sum_ctau = 0
    n_entries = 0
    
    for event in tree:
        branch = getattr(event, branch_name)
        
        # Check if branch is iterable (vector) or scalar
        try:
            n = len(branch)
            # multiple staus in event
            for val in branch:
                sum_ctau += val
                n_entries += 1
        except TypeError:
            # single value per event
            sum_ctau += branch
            n_entries += 1
    
    f.Close()
    if n_entries == 0:
        return 0
    return sum_ctau / n_entries

# file1 = "output/stage1_combined_0603/FCCee_118_stau_10m_ctau_ecm_240.root"
# file2 = "output/stage1_combined_0603/FCCee_118_stau_20m_ctau_ecm_240.root"
# file3 = "output/stage1_combined_0603/FCCee_118_stau_6m_ctau_ecm_240.root"
# file4 = "output/stage1_combined_0603/FCCee_118_stau_4m_ctau_ecm_240.root"

# file1 = "output/stage1_combined_0603/FCCee_119_stau_10m_ctau_ecm_240.root"
# file2 = "output/stage1_combined_0603/FCCee_119_stau_20m_ctau_ecm_240.root"
# file3 = "output/stage1_combined_0603/FCCee_119_stau_6m_ctau_ecm_240.root"
# file4 = "output/stage1_combined_0603/FCCee_119_stau_4m_ctau_ecm_240.root"
# file5 = "output/stage1_combined_0603/FCCee_119_stau_1m_ctau_ecm_240.root"
# file6 = "output/stage1_combined_0603/FCCee_119_stau_2m_ctau_ecm_240.root"
# file7 = "output/stage1_combined_0603/FCCee_119_stau_3m_ctau_ecm_240.root"

file1 = "output/testing/FCCee_110_stau_10m_ctau_ecm_240.root"
file2 = "output/testing/FCCee_110_stau_20m_ctau_ecm_240.root"
file3 = "output/stage1_combined_0603/FCCee_110_stau_6m_ctau_ecm_240.root"
file4 = "output/stage1_combined_0603/FCCee_110_stau_4m_ctau_ecm_240.root"
file5 = "output/stage1_combined_0603/FCCee_110_stau_1m_ctau_ecm_240.root"
file6 = "output/stage1_combined_0603/FCCee_110_stau_2m_ctau_ecm_240.root"
file7 = "output/stage1_combined_0603/FCCee_110_stau_3m_ctau_ecm_240.root"


avg1 = average_ctau(file1)
avg2 = average_ctau(file2)
avg3 = average_ctau(file3)
avg4 = average_ctau(file4)
avg5 = average_ctau(file5)
avg6 = average_ctau(file6)
avg7 = average_ctau(file7)
print(f"######### Decay position of stau of mass: 110 GeV  ###########")
print(f"Average gentau_ctau for 1m input: {avg5:.2f} cm")
print(f"Average gentau_ctau for 2m input: {avg6:.2f} cm")
print(f"Average gentau_ctau for 3m input: {avg7:.2f} cm")
print(f"Average decayLengthTau for 4m input: {avg4:.2f} cm")
print(f"Average decayLengthTau for 6m input: {avg3:.2f} cm")
print(f"Average decayLengthTau for 10m input: {avg1:.2f} cm")
print(f"Average decayLengthTau for 20m input: {avg2:.2f} cm")

