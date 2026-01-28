The EMD4HEP output is taken from: [Delphes-EDM4HEP](https://github.com/soumyaa-physics/delphes/tree/master/cards/IDEA)

Source the setup:
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

For running pre-edm4hep samples source this and switch to pre-edm4hep branch:
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
```
use setup_pre.sh now
run:
```
fccanalysis run analysis_stage_pre1.py --nevents 10000
```

analysis_stage1.py : 
```
fccanalysis run analysis_stage1.py 
```
analysis stage final:
```
fccanalysis final analysis_final.py 
```
analysis plot:
```
fccanalysis plots analysis_plots.py
```
for condor submission use: useful for backgrounds
```
fccanalysis submit analysis_stage1.py 
```
EDM4hep → analysis_stage1 → flat ntuple
flat ntuple → analysis_final → histograms
histograms → analysis_plot → plots

References:
1. [VertexUtils](https://github.com/HEP-FCC/FCCAnalyses/blob/763cb483f4c8e605b3182c8b1d076cdd920739b2/analyzers/dataframe/src/VertexingUtils.cc)
2. [vertex_github](https://hep-fcc.github.io/FCCAnalyses/doc/latest/VertexingUtils_8h.html)

condor jobs
condor_q
