from argparse import ArgumentParser

# Mandatory: Analysis class where the user defines the operations on the dataframe
class Analysis():
    '''
    Generator-level llp stau properties from Delphes + Pythia output.
    '''
    def __init__(self, cmdline_args):
        parser = ArgumentParser(description='Additional analysis arguments', usage='Provided after "--"')
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        # Mandatory: List of datasets used in the analysis
        self.process_list = {
            # uncomment the samples you want to run over
            #######################################################
            #               CME: 240 GeV (ZH)                     #
            #######################################################
            'FCC_100stau_240COM'  : {'fraction': 1.0},
            # Background samples: e⁺e⁻ → γγγ
            #'background_3photons_cme_91p188'    : {'fraction': 1.0},
            } 
        
        # Optional: number of threads to run on, default is 'all available'
        # self.n_threads = 4

        # Optional: batch queue name when running on HTCondor, default is 'longlunch'
        # self.batch_queue = 'longlunch'

        # Optional: computing account when running on CERN's HTCondor, default is 'group_u_FCC.local_gen'
        # self.comp_group = 'group_u_FCC.local_gen'

        # Input/output directories
        # uncomment the input directory you want to use
        self.input_dir = '/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/edm4hep_output'
        self.output_dir = '/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output'

        self.tree_name = "events"
        print(f"Using tree: {self.tree_name}")

        # Optional HTCondor settings
        # self.run_batch = True
        # self.batch_queue = 'longlunch'
        # self.comp_group = 'group_u_FCC.local_gen'

        # Optional: number of threads
        # self.n_threads = 4


    # Mandatory: analyzers function to define the analysis graph
    def analyzers(self, df):
        '''
        Generator-level particle properties.
        '''
        df = df.Alias("MCRecoAssociations0", "_RecoMCLink_from.index")
        df = df.Alias("MCRecoAssociations1", "_RecoMCLink_to.index")

        df2 = (
            df
            # --------------------------
            # Generator-level particles
            # --------------------------
            .Alias("Particle1", "_Particle_daughters.index")
            .Define("GenStau", "MCParticle::sel_pdgID(1000015,true)(Particle)")
            .Define("n_GenStau", "MCParticle::get_n(GenStau)")
            .Define("GenStau_vx", "MCParticle::get_vertex_x(GenStau)") 
            .Define("GenStau_vy", "MCParticle::get_vertex_y(GenStau)")
            .Define("GenStau_vz", "MCParticle::get_vertex_z(GenStau)")
            .Define("GenStau_Lxy", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy)")
            .Define("GenStau_Lxyz", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy + GenStau_vz*GenStau_vz)")

            # Stau daughters
            .Define("GenTau", "MCParticle::sel_pdgID(15,true)(Particle)")
            .Define("GenGravitino", "MCParticle::sel_pdgID(1000049,false)(Particle)")

            # Tau kinematics and vertex
            .Define("GenTau_px", "MCParticle::get_px(GenTau)")
            .Define("GenTau_py", "MCParticle::get_py(GenTau)")
            .Define("GenTau_pz", "MCParticle::get_pz(GenTau)")
            .Define("GenTau_pt", "MCParticle::get_pt(GenTau)")
            .Define("GenTau_eta", "MCParticle::get_eta(GenTau)")
            .Define("GenTau_phi", "MCParticle::get_phi(GenTau)")
            .Define("GenTau_e", "MCParticle::get_e(GenTau)")
            .Define("GenTau_charge", "MCParticle::get_charge(GenTau)")
            .Define("GenTau_vx", "MCParticle::get_vertex_x(GenTau)")
            .Define("GenTau_vy", "MCParticle::get_vertex_y(GenTau)")
            .Define("GenTau_vz", "MCParticle::get_vertex_z(GenTau)")

            .Define("Tau_prod", "MCParticle::get_vertex(GenTau)") 
            .Define("GenTau_mass", "MCParticle::get_mass(GenTau)")  
            .Define("decayLengthTau", 'return sqrt(GenTau_vx.at(0)*GenTau_vx.at(0) + GenTau_vy.at(0)*GenTau_vy.at(0) + GenTau_vz.at(0)*GenTau_vz.at(0))')
           
            # Final-state generator particles (status 1)
            .Define("GenElectron_PID", "MCParticle::sel_pdgID(11, true)(Particle)")
            .Define("FSGenElectron", "MCParticle::sel_genStatus(1)(GenElectron_PID)") 
            .Define("GenMuon_PID", "MCParticle::sel_pdgID(13, true)(Particle)")
            .Define("FSGenMuon", "MCParticle::sel_genStatus(1)(GenMuon_PID)")  
            .Define("GenPhoton_PID", "MCParticle::sel_pdgID(22, false)(Particle)")
            .Define("FSGenPhoton", "MCParticle::sel_genStatus(1)(GenPhoton_PID)") 


            # Kinematics helper for FSGen electrons and positrons
            .Define("n_FSGenElectron", "MCParticle::get_n(FSGenElectron)")
            .Define("FSGenElectron_e", "if (n_FSGenElectron>0) return MCParticle::get_e(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_p", "if (n_FSGenElectron>0) return MCParticle::get_p(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_pt", "if (n_FSGenElectron>0) return MCParticle::get_pt(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_px", "if (n_FSGenElectron>0) return MCParticle::get_px(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_py", "if (n_FSGenElectron>0) return MCParticle::get_py(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_pz", "if (n_FSGenElectron>0) return MCParticle::get_pz(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_eta", "if (n_FSGenElectron>0) return MCParticle::get_eta(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_theta", "if (n_FSGenElectron>0) return MCParticle::get_theta(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_phi", "if (n_FSGenElectron>0) return MCParticle::get_phi(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_charge", "if (n_FSGenElectron>0) return MCParticle::get_charge(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_vx", "if (n_FSGenElectron>0) return MCParticle::get_vertex_x( FSGenElectron ); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_vy", "if (n_FSGenElectron>0) return MCParticle::get_vertex_y( FSGenElectron ); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_vz", "if (n_FSGenElectron>0) return MCParticle::get_vertex_z( FSGenElectron ); else return MCParticle::get_genStatus(GenElectron_PID);")

            # Kinematics for FSGen muons and anti-muons
            .Define("n_FSGenMuon", "MCParticle::get_n(FSGenMuon)")
            .Define("FSGenMuon_e", "if (n_FSGenMuon>0) return MCParticle::get_e(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_p", "if (n_FSGenMuon>0) return MCParticle::get_p(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_pt", "if (n_FSGenMuon>0) return MCParticle::get_pt(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_px", "if (n_FSGenMuon>0) return MCParticle::get_px(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_py", "if (n_FSGenMuon>0) return MCParticle::get_py(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_pz", "if (n_FSGenMuon>0) return MCParticle::get_pz(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_eta", "if (n_FSGenMuon>0) return MCParticle::get_eta(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_theta", "if (n_FSGenMuon>0) return MCParticle::get_theta(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_phi", "if (n_FSGenMuon>0) return MCParticle::get_phi(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_charge", "if (n_FSGenMuon>0) return MCParticle::get_charge(FSGenMuon); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_vx", "if (n_FSGenMuon>0) return MCParticle::get_vertex_x( FSGenMuon ); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_vy", "if (n_FSGenMuon>0) return MCParticle::get_vertex_y( FSGenMuon ); else return MCParticle::get_genStatus(GenMuon_PID);")
            .Define("FSGenMuon_vz", "if (n_FSGenMuon>0) return MCParticle::get_vertex_z( FSGenMuon ); else return MCParticle::get_genStatus(GenMuon_PID);")

            #Kinematics for FSGen photons
            .Define("n_FSGenPhoton", "MCParticle::get_n(FSGenPhoton)")
            .Define("FSGenPhoton_e", "MCParticle::get_e(FSGenPhoton)")
            .Define("FSGenPhoton_p", "MCParticle::get_p(FSGenPhoton)")
            .Define("FSGenPhoton_pt", "MCParticle::get_pt(FSGenPhoton)")
            .Define("FSGenPhoton_px", "MCParticle::get_px(FSGenPhoton)")
            .Define("FSGenPhoton_py", "MCParticle::get_py(FSGenPhoton)")
            .Define("FSGenPhoton_pz", "MCParticle::get_pz(FSGenPhoton)")
            .Define("FSGenPhoton_eta", "MCParticle::get_eta(FSGenPhoton)")
            .Define("FSGenPhoton_theta", "MCParticle::get_theta(FSGenPhoton)")
            .Define("FSGenPhoton_phi", "MCParticle::get_phi(FSGenPhoton)")

            # --------------------------
            # Reconstructed particles
            # --------------------------
            # MC Primary Vertex
            .Define("MC_PrimaryVertex",  "MCParticle::get_EventPrimaryVertex(21)( Particle )" )

            # Tracks
            .Define("n_RecoTracks",  "ReconstructedParticle2Track::getTK_n(_EFlowTrack_trackStates)")

            # Vertex fitting

            # First, reconstruct a vertex from all tracks 
            # Input parameters are 1 = primary vertex, _EFlowTrack_trackStates contains all tracks, bool beamspotconstraint = true, bsc x/y/z (sigma xy, sigma_z)
            .Define("VertexObject_allTracks",  "VertexFitterSimple::VertexFitter_Tk ( 1, _EFlowTrack_trackStates, true, 4.5, 20e-3, 300)")

            # Select the tracks that are reconstructed  as primaries, so it removes all long lived decay tracks
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( _EFlowTrack_trackStates, true, 4.5, 20e-3, 300, 0., 0., 0.)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")

            # the final primary vertex : final/refined fit
            .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
            .Define("PrimaryVertex",   "VertexingUtils::get_VertexData( PrimaryVertexObject )")

            # Displaced vertex reconstruction
            # The selection functions have to be updated and taken from the Vertexing Utils.cc 
 
            # # select tracks with pT > 1 GeV
            # .Define('sel_tracks_pt', "VertexingUtils::get_p_SV(1)(_EFlowTrack_trackStates)")
            # select tracks with |d0 |> 2 mm
            # .Define('sel_tracks', 'VertexingUtils::get_dR_SV_obj(2)(sel_tracks_pt)')
            # arguments: d0sig_min, d0sig_max, z0sig_min, z0sig_max
            # .Define('sel_tracks', 'selTracks(2., 100., -1., 1e5)(ReconstructedParticle, _EFlowTrack_trackStates)')
            
            # # .Define("sel_tracks",   "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack,  RecoedPrimaryTracks )")
            # .Define("sel_tracks", "VertexFitterSimple::get_NonPrimaryTracks(_EFlowTrack_trackStates, VertexObject_allTracks)")
            # # find the DVs
            # .Define("DV_evt_seltracks", "VertexFinderLCFIPlus::get_SV_event(sel_tracks, _EFlowTrack_trackStates, PrimaryVertexObject, true, 9., 40., 5.)")
            # # number of DVs
            # .Define('n_seltracks_DVs', 'VertexingUtils::get_n_SV(DV_evt_seltracks)')
            # # number of tracks from the DVs
            # .Define('n_trks_seltracks_DVs', 'VertexingUtils::get_VertexNtrk(DV_evt_seltracks)')
            # # invariant mass at the DVs (assuming the tracks to be pions)
            # .Define('invMass_seltracks_DVs', 'VertexingUtils::get_invM(DV_evt_seltracks)')

            # # get the chi2 distributions of the DVs from selected tracks - to check if tracks originate from the same physical point
            # .Define("DV_evt_seltracks_chi2",    "VertexingUtils::get_chi2_SV(DV_evt_seltracks)")
            # .Define("DV_evt_seltracks_normchi2","VertexingUtils::get_norm_chi2_SV(DV_evt_seltracks)") # DV chi2 (normalised)

            # # get the decay radius and full 3D distance of all the DVs from selected tracks
            # .Define("Reco_seltracks_DVs_Lxy","VertexingUtils::get_dxy_SV(DV_evt_seltracks, PrimaryVertexObject)")
            # .Define("Reco_seltracks_DVs_Lxyz","VertexingUtils::get_d3d_SV(DV_evt_seltracks, PrimaryVertexObject)")
            
            # # merge vertices with position within 10*error-of-position, get the tracks from the merged vertices and refit
            # # this is hard coded in merge vertices: to only select the pair of vertices which are spatially consistent
            # .Define('merged_DVs', 'VertexingUtils::mergeVertices(DV_evt_seltracks, _EFlowTrack_trackStates)')
            # # number of merged DVs
            # .Define("merged_DVs_n", "VertexingUtils::get_n_SV(merged_DVs)")
            # # number of tracks from the merged DVs
            # .Define('n_trks_merged_DVs', 'VertexingUtils::get_VertexNtrk(merged_DVs)')
            # # invariant mass at the merged DVs
            # .Define('invMass_merged_DVs', 'VertexingUtils::get_invM(merged_DVs)')

            # # get the chi2 distributions of the merged DVs
            # .Define("merged_DVs_chi2",    "VertexingUtils::get_chi2_SV(merged_DVs)")
            # .Define("merged_DVs_normchi2","VertexingUtils::get_norm_chi2_SV(merged_DVs)") # DV chi2 (normalised)

            # get the decay radius of all the merged DVs
            # .Define("Reco_DVs_merged_Lxy","VertexingUtils::get_dxy_SV(merged_DVs, PrimaryVertexObject)")
            # .Define("Reco_DVs_merged_Lxyz","VertexingUtils::get_d3d_SV(merged_DVs, PrimaryVertexObject)")

            # .Define("RecoTauTracks", "VertexingUtils::get_tracksInJets(Jet, _EFlowTrack_trackStates, Jet_to_Track_indices, 0)")
            # .Define("RecoTauDecayVertexObject", "VertexFitterSimple::VertexFitter_Tk(2, RecoTauTracks)")
            # .Define("RecoTauDecayVertex", "VertexingUtils::get_VertexData(RecoTauDecayVertexObject)")
            # .Define("RecoTau_Lxy", "sqrt(RecoTauDecayVertex.position.x*RecoTauDecayVertex.position.x + RecoTauDecayVertex.position.y*RecoTauDecayVertex.position.y)")
            # .Define("RecoTau_Lxyz", "sqrt(RecoTauDecayVertex.position.x*RecoTauDecayVertex.position.x + RecoTauDecayVertex.position.y*RecoTauDecayVertex.position.y + RecoTauDecayVertex.position.z*RecoTauDecayVertex.position.z)")

            # JETS
            .Define("n_RecoJets", "ReconstructedParticle::get_n(Jet)")
            .Define("RecoJet_e",      "ReconstructedParticle::get_e(Jet)")
            .Define("RecoJet_p",      "ReconstructedParticle::get_p(Jet)") 
            .Define("RecoJet_pt",      "ReconstructedParticle::get_pt(Jet)")
            .Define("RecoJet_px",      "ReconstructedParticle::get_px(Jet)")
            .Define("RecoJet_py",      "ReconstructedParticle::get_py(Jet)")
            .Define("RecoJet_pz",      "ReconstructedParticle::get_pz(Jet)")
		    .Define("RecoJet_eta",     "ReconstructedParticle::get_eta(Jet)") 
            .Define("RecoJet_theta",   "ReconstructedParticle::get_theta(Jet)")
		    .Define("RecoJet_phi",     "ReconstructedParticle::get_phi(Jet)") 
            .Define("RecoJet_charge",  "ReconstructedParticle::get_charge(Jet)")
            .Define("RecoJet_mvis",     "ReconstructedParticle::get_P4vis(Jet)")
            .Define("RecoJetTrack_absD0", "return abs(ReconstructedParticle2Track::getRP2TRK_D0(Jet,_EFlowTrack_trackStates))")
            .Define("RecoJetTrack_absZ0", "return abs(ReconstructedParticle2Track::getRP2TRK_Z0(Jet,_EFlowTrack_trackStates))")
            .Define("RecoJetTrack_absD0sig", "return abs(ReconstructedParticle2Track::getRP2TRK_D0_sig(Jet,_EFlowTrack_trackStates))") 
            .Define("RecoJetTrack_absZ0sig", "return abs(ReconstructedParticle2Track::getRP2TRK_Z0_sig(Jet,_EFlowTrack_trackStates))")
            .Define("RecoJetTrack_D0cov", "ReconstructedParticle2Track::getRP2TRK_D0_cov(Jet,_EFlowTrack_trackStates)") 
            .Define("RecoJetTrack_Z0cov", "ReconstructedParticle2Track::getRP2TRK_Z0_cov(Jet,_EFlowTrack_trackStates)")

            # Electrons
            .Alias("Electron0", "Electron_objIdx.index")
            .Define('RecoElectrons',  'ReconstructedParticle::get(Electron0, ReconstructedParticles)') 
            .Define("n_RecoElectrons", "ReconstructedParticle::get_n(RecoElectrons)")
            .Define("RecoElectrons_p", "ReconstructedParticle::get_p(RecoElectrons)")
            .Define("RecoElectrons_px", "ReconstructedParticle::get_px(RecoElectrons)")
            .Define("RecoElectrons_py", "ReconstructedParticle::get_py(RecoElectrons)")
            .Define("RecoElectrons_pz", "ReconstructedParticle::get_pz(RecoElectrons)")
            .Define("RecoElectrons_pt", "ReconstructedParticle::get_pt(RecoElectrons)")
            .Define("RecoElectrons_eta", "ReconstructedParticle::get_eta(RecoElectrons)")
            .Define("RecoElectrons_phi", "ReconstructedParticle::get_phi(RecoElectrons)")
            .Define("RecoElectrons_e", "ReconstructedParticle::get_e(RecoElectrons)")
            .Define("RecoElectrons_charge", "ReconstructedParticle::get_charge(RecoElectrons)")
            .Define("RecoElectrons_theta", "ReconstructedParticle::get_theta(RecoElectrons)")
            # .Define("RecoElectrons_vx", "ReconstructedParticle::get_vertex_x(RecoElectrons)")
            # .Define("RecoElectrons_vy", "ReconstructedParticle::get_vertex_y(RecoElectrons)")
            # .Define("RecoElectrons_vz", "ReconstructedParticle::get_vertex_z(RecoElectrons)")

            # Muons
            .Alias("Muon0", "Muon_objIdx.index")
            .Define('RecoMuons',  'ReconstructedParticle::get(Muon0, ReconstructedParticles)') 
            .Define("n_RecoMuons", "ReconstructedParticle::get_n(RecoMuons)")
            .Define("RecoMuons_p", "ReconstructedParticle::get_p(RecoMuons)")
            .Define("RecoMuons_px", "ReconstructedParticle::get_px(RecoMuons)")
            .Define("RecoMuons_py", "ReconstructedParticle::get_py(RecoMuons)")
            .Define("RecoMuons_pz", "ReconstructedParticle::get_pz(RecoMuons)")
            .Define("RecoMuons_pt", "ReconstructedParticle::get_pt(RecoMuons)")
            .Define("RecoMuons_eta", "ReconstructedParticle::get_eta(RecoMuons)")
            .Define("RecoMuons_phi", "ReconstructedParticle::get_phi(RecoMuons)")
            .Define("RecoMuons_e", "ReconstructedParticle::get_e(RecoMuons)")
            .Define("RecoMuons_charge", "ReconstructedParticle::get_charge(RecoMuons)")
            .Define("RecoMuons_theta", "ReconstructedParticle::get_theta(RecoMuons)")
            # .Define("RecoMuons_vx", "ReconstructedParticle::get_vertex_x(RecoMuons)")
            # .Define("RecoMuons_vy", "ReconstructedParticle::get_vertex_y(RecoMuons)")
            # .Define("RecoMuons_vz", "ReconstructedParticle::get_vertex_z(RecoMuons)")

            # PHOTONS
            .Alias("Photon0", "Photon_objIdx.index") 
            .Define("RecoPhotons",  "ReconstructedParticle::get(Photon0, ReconstructedParticles)")
            .Define("n_RecoPhotons",  "ReconstructedParticle::get_n(RecoPhotons)") 
            .Define("RecoPhoton_e",      "ReconstructedParticle::get_e(RecoPhotons)")
            .Define("RecoPhoton_p",      "ReconstructedParticle::get_p(RecoPhotons)")
            .Define("RecoPhoton_pt",      "ReconstructedParticle::get_pt(RecoPhotons)")
            .Define("RecoPhoton_px",      "ReconstructedParticle::get_px(RecoPhotons)")
            .Define("RecoPhoton_py",      "ReconstructedParticle::get_py(RecoPhotons)")
            .Define("RecoPhoton_pz",      "ReconstructedParticle::get_pz(RecoPhotons)")
		    .Define("RecoPhoton_eta",     "ReconstructedParticle::get_eta(RecoPhotons)")
            .Define("RecoPhoton_theta",   "ReconstructedParticle::get_theta(RecoPhotons)")
		    .Define("RecoPhoton_phi",     "ReconstructedParticle::get_phi(RecoPhotons)") 
            .Define("RecoPhoton_charge",  "ReconstructedParticle::get_charge(RecoPhotons)")

            # MET
            .Define("RecoMissingEnergy_e", "ReconstructedParticle::get_e(MissingET)[0]")
		    .Define("RecoMissingEnergy_p", "ReconstructedParticle::get_p(MissingET)[0]")
		    .Define("RecoMissingEnergy_pt", "ReconstructedParticle::get_pt(MissingET)[0]")
		    .Define("RecoMissingEnergy_px", "ReconstructedParticle::get_px(MissingET)[0]") #x-component of RecoMissingEnergy
		    .Define("RecoMissingEnergy_py", "ReconstructedParticle::get_py(MissingET)[0]") #y-component of RecoMissingEnergy
		    .Define("RecoMissingEnergy_pz", "ReconstructedParticle::get_pz(MissingET)[0]") #z-component of RecoMissingEnergy
		    .Define("RecoMissingEnergy_eta", "ReconstructedParticle::get_eta(MissingET)[0]")
		    .Define("RecoMissingEnergy_theta", "ReconstructedParticle::get_theta(MissingET)[0]")
		    .Define("RecoMissingEnergy_phi", "ReconstructedParticle::get_phi(MissingET)[0]") #angle of RecoMissingEnergy

        )
        return df2

# Mandatory: output function
    def output(self):
        '''
        Output branches to be saved in the final ROOT file.
        '''
        branch_list = [
            # Gen-level stau
            "GenStau",
            "n_GenStau",
            "GenStau_vx",
            "GenStau_vy",
            "GenStau_vz",
            "GenStau_Lxy",
            "GenStau_Lxyz",

            # # Stau daughters
            "GenTau",
            "GenGravitino",

            # Tau kinematics
            "GenTau_e",
            "GenTau_pt",
            "GenTau_eta",
            "GenTau_phi",
            "GenTau_px",
            "GenTau_py",
            "GenTau_pz",
            "GenTau_charge",
            "decayLengthTau",

            # Tau vertex
            "GenTau_vx",
            "GenTau_vy",
            "GenTau_vz",

            # Final-state electrons
            "FSGenElectron",
            "n_FSGenElectron",
            "FSGenElectron_e",
            "FSGenElectron_p",
            "FSGenElectron_pt",
            "FSGenElectron_px",
            "FSGenElectron_py",
            "FSGenElectron_pz",
            "FSGenElectron_eta",
            "FSGenElectron_theta",
            "FSGenElectron_phi",
            "FSGenElectron_charge",
            "FSGenElectron_vx",
            "FSGenElectron_vy",
            "FSGenElectron_vz",

            # Final-state muons
            "FSGenMuon",
            "n_FSGenMuon",
            "FSGenMuon_e",
            "FSGenMuon_p",
            "FSGenMuon_pt",
            "FSGenMuon_px",
            "FSGenMuon_py",
            "FSGenMuon_pz",
            "FSGenMuon_eta",
            "FSGenMuon_theta",
            "FSGenMuon_phi",
            "FSGenMuon_charge",
            "FSGenMuon_vx",
            "FSGenMuon_vy",
            "FSGenMuon_vz",

            # Final-state photons
            "FSGenPhoton",
            "n_FSGenPhoton",
            "FSGenPhoton_e",
            "FSGenPhoton_p",
            "FSGenPhoton_pt",
            "FSGenPhoton_px",
            "FSGenPhoton_py",
            "FSGenPhoton_pz",
            "FSGenPhoton_eta",
            "FSGenPhoton_theta",
            "FSGenPhoton_phi",

            # # # Track information
            # "n_RecoTracks",
            # "n_seltracks_DVs",
            # "n_trks_seltracks_DVs",
            # "invMass_seltracks_DVs",
            # "DV_evt_seltracks_chi2",
            # "DV_evt_seltracks_normchi2",
            # "Reco_seltracks_DVs_Lxy",
            # "Reco_seltracks_DVs_Lxyz",

            # "merged_DVs_n",
            # "n_trks_merged_DVs",
            # "invMass_merged_DVs",
            # "merged_DVs_chi2",
            # "merged_DVs_normchi2",
            # "Reco_DVs_merged_Lxy",
            # "Reco_DVs_merged_Lxyz",

            # Reco Jets
            "n_RecoJets",
            "RecoJet_e",
            "RecoJet_p",
            "RecoJet_pt",
            "RecoJet_px",
            "RecoJet_py",
            "RecoJet_pz",
            "RecoJet_eta",
            "RecoJet_theta",
            "RecoJet_phi",
            "RecoJet_charge",
            "RecoJet_mvis",
            "RecoJetTrack_absD0",
            "RecoJetTrack_absZ0",
            "RecoJetTrack_absD0sig",
            "RecoJetTrack_absZ0sig",
            "RecoJetTrack_D0cov",
            "RecoJetTrack_Z0cov",

            # Reco Electrons
            "RecoElectrons",
            "n_RecoElectrons",
            "RecoElectrons_e",
            "RecoElectrons_p",
            "RecoElectrons_pt",
            "RecoElectrons_px",
            "RecoElectrons_py",
            "RecoElectrons_pz",
            "RecoElectrons_eta",
            "RecoElectrons_theta",
            "RecoElectrons_phi",
            "RecoElectrons_charge",
            # "RecoElectrons_vx",
            # "RecoElectrons_vy",
            # "RecoElectrons_vz",

            # Reco Muons
            "RecoMuons",
            "n_RecoMuons",
            "RecoMuons_e",
            "RecoMuons_p",
            "RecoMuons_pt",
            "RecoMuons_px",
            "RecoMuons_py",
            "RecoMuons_pz",
            "RecoMuons_eta",
            "RecoMuons_theta",
            "RecoMuons_phi",
            "RecoMuons_charge",
            # "RecoMuons_vx",
            # "RecoMuons_vy",
            # "RecoMuons_vz",

            # Reco Photons
            "RecoPhotons",
            "n_RecoPhotons",
            "RecoPhoton_e",
            "RecoPhoton_p",
            "RecoPhoton_pt",
            "RecoPhoton_px",
            "RecoPhoton_py",
            "RecoPhoton_pz",
            "RecoPhoton_eta",
            "RecoPhoton_theta",
            "RecoPhoton_phi",
            "RecoPhoton_charge",

            # MET
            "RecoMissingEnergy_e",
            "RecoMissingEnergy_p",
            "RecoMissingEnergy_pt",
            "RecoMissingEnergy_px",
            "RecoMissingEnergy_py",
            "RecoMissingEnergy_pz",
            "RecoMissingEnergy_eta",
            "RecoMissingEnergy_theta",
            "RecoMissingEnergy_phi",
        ]

        return branch_list
