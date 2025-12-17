The EMD4HEP output is taken from: [Delphes-EDM4HEP](https://github.com/soumyaa-physics/delphes/tree/master/cards/IDEA)

analysis_stage1.py : 
```
fccanalysis run analysis_stage1.py 
```
Analysis notes and current issues:

in exotic higgs decays - analysis_stage1
there is a selection criteria to select tracks for displaced vertex reco- line 119
the class called there  ` VertexingUtils::sel_pt_tracks(1)`
no longer exists
trying to find an alternate

--updates
removed displaced vertex info
and MET

References:
1. [VertexUtils](https://github.com/HEP-FCC/FCCAnalyses/blob/763cb483f4c8e605b3182c8b1d076cdd920739b2/analyzers/dataframe/src/VertexingUtils.cc)
2. [vertex_github](https://hep-fcc.github.io/FCCAnalyses/doc/latest/VertexingUtils_8h.html)
