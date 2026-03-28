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
            # 20 cm
            # "FCCee_100_stau_20cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_105_stau_20cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_110_stau_20cm_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_115_stau_20cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_20cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_20cm_ctau_ecm_240": {'fraction': 1.0},

            # # 50 cm
            # # "FCCee_100_stau_50cm_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_105_stau_50cm_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_110_stau_50cm_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_115_stau_50cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_50cm_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_50cm_ctau_ecm_240": {'fraction': 1.0},

            # # 1 m
            # # "FCCee_100_stau_1m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_105_stau_1m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_110_stau_1m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_115_stau_1m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_1m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_1m_ctau_ecm_240": {'fraction': 1.0},

            # # 2 m
            # # "FCCee_100_stau_2m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_105_stau_2m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_110_stau_2m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_115_stau_2m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_2m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_2m_ctau_ecm_240": {'fraction': 1.0},

            # 3 m
            # "FCCee_100_stau_3m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_105_stau_3m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_110_stau_3m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_115_stau_3m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_3m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_3m_ctau_ecm_240": {'fraction': 1.0},

            # # # 4 m
            # # # "FCCee_100_stau_4m_ctau_ecm_240": {'fraction': 1.0},
            # # # "FCCee_105_stau_4m_ctau_ecm_240": {'fraction': 1.0},
            # # # "FCCee_110_stau_4m_ctau_ecm_240": {'fraction': 1.0},
            # # # "FCCee_115_stau_4m_ctau_ecm_240": {'fraction': 1.0},
            "FCCee_118_stau_4m_ctau_ecm_240": {'fraction': 1.0},
            "FCCee_119_stau_4m_ctau_ecm_240": {'fraction': 1.0},

            # # # 6 m
            # # "FCCee_100_stau_6m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_105_stau_6m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_110_stau_6m_ctau_ecm_240": {'fraction': 1.0},
            # # "FCCee_115_stau_6m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_6m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_6m_ctau_ecm_240": {'fraction': 1.0},

            # # 10 m
            # "FCCee_100_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_105_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_110_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_115_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_10m_ctau_ecm_240": {'fraction': 1.0},
            
            # # # 20 m
            # "FCCee_100_stau_20m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_105_stau_20m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_110_stau_20m_ctau_ecm_240": {'fraction': 1.0}, ## SOME PROBLEM HERE!!
            # "FCCee_115_stau_20m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_118_stau_20m_ctau_ecm_240": {'fraction': 1.0},
            # "FCCee_119_stau_20m_ctau_ecm_240": {'fraction': 1.0},

        }

        # prodTag     = "FCCee/winter2023/IDEA/"
        # prodTag     = "FCCee/spring2021/IDEA/"

        # Input/output directories
        # uncomment the input directory you want to use

        # backgrounds are stored here:
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA'
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/'
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/'

        self.input_dir = '/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/edm4hep_output_2003'
        self.output_dir = '/eos/home-s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/stage1_combined_2003'

        # for batch
        # self.output_dir = '.output/condor_jobs/'

        # # Optional: output directory on eos, if specified files will be copied
        # # there once the batch job is done, default is empty
        # self.output_dir_eos = '/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output/condor_jobs'
        # eosType = "eosuser"

        self.tree_name = "events"
        print(f"Using tree: {self.tree_name}")

        # # self.should_transfer_files = True

        # # Optional HTCondor settings
        # # self.run_batch = True # does not exist anymore
        # self.batch_queue = 'tomorrow'
        # self.comp_group = 'group_u_FCC.local_gen'
        # # self.
        # # # Optional: number of threads
        # self.n_threads = 4

    # Mandatory: analyzers function to define the analysis graph
    def analyzers(self, df):
        '''
        Generator-level particle properties.
            '''
        cols = set(df.GetColumnNames())
        hasRecoMCLink = False

        cols = df.GetColumnNames()
        # For private samples with legacy or EDM4hep MC–Reco associations
        if "_RecoMCLink_from.index" in cols:
            print("INFO: Using legacy RecoMCLink")
            df = df.Alias("MCRecoAssociations0", "_RecoMCLink_from.index")
            df = df.Alias("MCRecoAssociations1", "_RecoMCLink_to.index")

        # For centrally-produced samples with EDM4hep MC–Reco associations
        elif "MCRecoAssociations#0.index" in cols:
            print("INFO: Using EDM4hep MCRecoAssociations")
            df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
            df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

        else:
            raise RuntimeError("No MC–Reco association collection found")
        
        hasEFlowTrackStates = (
        "_EFlowTrack_trackStates" in cols or
        "EFlowTrack.trackStates_begin" in cols
        )
        
        if "EFlowTrack_1" in cols:
            df = df.Alias("TrackStates", "EFlowTrack_1")
        elif "_EFlowTrack_trackStates" in cols:
            df = df.Alias("TrackStates", "_EFlowTrack_trackStates")
        elif "EFlowTrack.trackStates_begin" in cols:
            df = df.Alias("TrackStates", "EFlowTrack.trackStates_begin")
        else:
            raise RuntimeError("No Electron index collection found")

        if "Electron#0.index" in cols:
            df = df.Alias("ElectronIdx", "Electron#0.index")
        elif "Electron_objIdx.index" in cols:
            df = df.Alias("ElectronIdx", "Electron_objIdx.index")
        else:
            raise RuntimeError("No Electron index collection found")
        
        if "Muon#0.index" in cols:
            df = df.Alias("MuonIdx", "Muon#0.index")
        elif "Muon_objIdx.index" in cols:
            df = df.Alias("MuonIdx", "Muon_objIdx.index")
        else:
            raise RuntimeError("No Muon index collection found")
        
        if "Photon#0.index" in cols:
            df = df.Alias("PhotonIdx", "Photon#0.index")
        elif "Photon_objIdx.index" in cols:
            df = df.Alias("PhotonIdx", "Photon_objIdx.index")   
        else:
            raise RuntimeError("No Photon index collection found")

        df2 = (
            df
            # --------------------------
            # Generator-level particles
            # --------------------------
            # .Alias("Particle1", "_Particle_daughters.index")
            .Define("GenStau", "MCParticle::sel_pdgID(1000015,true)(Particle)")
            .Define("n_GenStau", "MCParticle::get_n(GenStau)")
            .Define("GenStau_status", "MCParticle::get_genStatus(GenStau)")
            .Define("GenStau_vx", "MCParticle::get_vertex_x(GenStau)") 
            .Define("GenStau_vy", "MCParticle::get_vertex_y(GenStau)")
            .Define("GenStau_e", "MCParticle::get_e(GenStau)")
            .Define("GenStau_vz", "MCParticle::get_vertex_z(GenStau)")
            .Define("GenStau_Lxy", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy)")
            .Define("GenStau_Lxyz", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy + GenStau_vz*GenStau_vz)")
            .Define("GenStau_time", "MCParticle::get_time(GenStau)")
            .Define("GenStau_theta", "MCParticle::get_theta(GenStau)")
            .Define("GenStau_phi", "MCParticle::get_phi(GenStau)")
            # Stau daughters
            # .Define("StauTauDecay", "MCParticle::get_indices_MotherByIndex(GenStau_idx[0], 15,  false,  true, false, Particle, Particle1.index)")
            # .Define("GenTau_daughters", "MCParticle::get_list_of_particles_from_decay(GenTau)")

            .Define("GenTau", "MCParticle::sel_pdgID(15,true)(Particle)")
            .Define("GenTau_status", "MCParticle::get_genStatus(GenTau)")
            .Define("GenTau_theta", "MCParticle::get_theta(GenTau)")
            .Define("GenGravitino", "MCParticle::sel_pdgID(1000049,false)(Particle)")
            .Define("n_GenGravitino", "MCParticle::get_n(GenGravitino)")
            .Define("GenGravitino_status", "MCParticle::get_genStatus(GenGravitino)")
            .Define("GenGravitino_e", "MCParticle::get_e(GenGravitino)")

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
            .Define("GenTau_time", "MCParticle::get_time(GenTau)")
            .Define("GenTau_cTau",
            "ROOT::VecOps::RVec<float> v; "
            "for (auto t : GenTau_time) v.push_back(t * 2.99792458e10); "
            "return v;")

            .Define("Tau_prod", "MCParticle::get_vertex(GenTau)") 
            .Define("GenTau_mass", "MCParticle::get_mass(GenTau)")  
            .Define("decayLengthTau", "GenTau_vx.size() > 0 ? sqrt(GenTau_vx[0]*GenTau_vx[0] + GenTau_vy[0]*GenTau_vy[0] + GenTau_vz[0]*GenTau_vz[0]) : -1.f")

            # Final-state generator particles (status 1)
            .Define("GenElectron_PID", "MCParticle::sel_pdgID(11, true)(Particle)")
            .Define("FSGenElectron", "MCParticle::sel_genStatus(1)(GenElectron_PID)") 
            .Define("GenMuon_PID", "MCParticle::sel_pdgID(13, true)(Particle)")
            .Define("FSGenMuon", "MCParticle::sel_genStatus(1)(GenMuon_PID)")  
            .Define("GenPhoton_PID", "MCParticle::sel_pdgID(22, false)(Particle)")
            .Define("FSGenPhoton", "MCParticle::sel_genStatus(1)(GenPhoton_PID)") 

            # Kinematics helper for FSGen electrons and positrons
            .Define("n_FSGenElectron", "MCParticle::get_n(FSGenElectron)")
            .Define("FSGenElectron_e", "MCParticle::get_e(FSGenElectron)")
            .Define("FSGenElectron_p", "MCParticle::get_p(FSGenElectron)")
            .Define("FSGenElectron_pt", "MCParticle::get_pt(FSGenElectron)")
            .Define("FSGenElectron_px", "MCParticle::get_px(FSGenElectron)")
            .Define("FSGenElectron_py", "MCParticle::get_py(FSGenElectron)")
            .Define("FSGenElectron_pz", "MCParticle::get_pz(FSGenElectron)")
            .Define("FSGenElectron_eta", "MCParticle::get_eta(FSGenElectron)")
            .Define("FSGenElectron_theta", "MCParticle::get_theta(FSGenElectron)")
            .Define("FSGenElectron_phi", "MCParticle::get_phi(FSGenElectron)")
            .Define("FSGenElectron_charge", "MCParticle::get_charge(FSGenElectron)")
            .Define("FSGenElectron_vx", "MCParticle::get_vertex_x(FSGenElectron)")
            .Define("FSGenElectron_vy", "MCParticle::get_vertex_y(FSGenElectron)")
            .Define("FSGenElectron_vz", "MCParticle::get_vertex_z(FSGenElectron)")
            .Define("FSGenElectron_time", "MCParticle::get_time(FSGenElectron)")

            # Kinematics for FSGen muons and anti-muons
            .Define("n_FSGenMuon", "MCParticle::get_n(FSGenMuon)")
            .Define("FSGenMuon_e", "MCParticle::get_e(FSGenMuon)")
            .Define("FSGenMuon_p", "MCParticle::get_p(FSGenMuon)")
            .Define("FSGenMuon_pt", "MCParticle::get_pt(FSGenMuon)")
            .Define("FSGenMuon_px", "MCParticle::get_px(FSGenMuon)")
            .Define("FSGenMuon_py", "MCParticle::get_py(FSGenMuon)")
            .Define("FSGenMuon_pz", "MCParticle::get_pz(FSGenMuon)")
            .Define("FSGenMuon_eta", "MCParticle::get_eta(FSGenMuon)")
            .Define("FSGenMuon_theta", "MCParticle::get_theta(FSGenMuon)")
            .Define("FSGenMuon_phi", "MCParticle::get_phi(FSGenMuon)")
            .Define("FSGenMuon_charge", "MCParticle::get_charge(FSGenMuon)")
            .Define("FSGenMuon_vx", "MCParticle::get_vertex_x(FSGenMuon)")
            .Define("FSGenMuon_vy", "MCParticle::get_vertex_y(FSGenMuon)")
            .Define("FSGenMuon_vz", "MCParticle::get_vertex_z(FSGenMuon)")
            .Define("FSGenMuon_time", "MCParticle::get_time(FSGenMuon)")
        
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

            # custon neutrino PID to include nu_e, nu_mu, nu_tau
            .Define("GenNeutrino_PID", "FCCAnalyses::MCParticle::sel_pdgID(16, true)(Particle)")
            .Define("FSGenNeutrino", "FCCAnalyses::MCParticle::sel_genStatus(1)(GenNeutrino_PID)") # gen status==1 means final state particle (FS)
            .Define("n_FSGenNeutrino", "FCCAnalyses::MCParticle::get_n(FSGenNeutrino)")
            .Define("FSGenNeutrino_e", "FCCAnalyses::MCParticle::get_e(FSGenNeutrino)")
            .Define("FSGenNeutrino_p", "FCCAnalyses::MCParticle::get_p(FSGenNeutrino)")
            .Define("FSGenNeutrino_pt", "FCCAnalyses::MCParticle::get_pt(FSGenNeutrino)")
            .Define("FSGenNeutrino_px", "FCCAnalyses::MCParticle::get_px(FSGenNeutrino)")
            .Define("FSGenNeutrino_py", "FCCAnalyses::MCParticle::get_py(FSGenNeutrino)")
            .Define("FSGenNeutrino_pz", "FCCAnalyses::MCParticle::get_pz(FSGenNeutrino)")
            .Define("FSGenNeutrino_eta", "FCCAnalyses::MCParticle::get_eta(FSGenNeutrino)")
            .Define("FSGenNeutrino_theta", "FCCAnalyses::MCParticle::get_theta(FSGenNeutrino)")
            .Define("FSGenNeutrino_phi", "FCCAnalyses::MCParticle::get_phi(FSGenNeutrino)")
            .Define("FSGenNeutrino_charge", "FCCAnalyses::MCParticle::get_charge(FSGenNeutrino)")
            
            # # JETS
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
            .Define("RecoJetTrack_absD0", "return abs(ReconstructedParticle2Track::getRP2TRK_D0(Jet,TrackStates))")
            .Define("RecoJetTrack_absZ0", "return abs(ReconstructedParticle2Track::getRP2TRK_Z0(Jet,TrackStates))")
            .Define("RecoJetTrack_absD0sig", "return abs(ReconstructedParticle2Track::getRP2TRK_D0_sig(Jet,TrackStates))") 
            .Define("RecoJetTrack_absZ0sig", "return abs(ReconstructedParticle2Track::getRP2TRK_Z0_sig(Jet,TrackStates))")
            .Define("RecoJetTrack_D0cov", "ReconstructedParticle2Track::getRP2TRK_D0_cov(Jet,TrackStates)") 
            .Define("RecoJetTrack_Z0cov", "ReconstructedParticle2Track::getRP2TRK_Z0_cov(Jet,TrackStates)")

            # Electrons
            .Alias("Electron0", "ElectronIdx")
            .Define("RecoElectrons",  "ReconstructedParticle::get(Electron0, ReconstructedParticles)") 
            .Define("n_RecoElectrons", "ReconstructedParticle::get_n(RecoElectrons)")
            .Define("RecoElectrons_p", "ReconstructedParticle::get_p(RecoElectrons)")
            # Muons
            .Alias("Muon0", "MuonIdx")
            .Define("RecoMuons",  "ReconstructedParticle::get(Muon0, ReconstructedParticles)") 
            .Define("n_RecoMuons", "ReconstructedParticle::get_n(RecoMuons)")
            .Define("RecoMuons_p", "ReconstructedParticle::get_p(RecoMuons)")

            # PHOTONS
            .Alias("Photon0", "PhotonIdx") 
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
            .Define("hasMissingET", "MissingET.size() > 0")
            .Define("RecoMissingEnergy_e", "hasMissingET ? ReconstructedParticle::get_e(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_p",  "hasMissingET ? ReconstructedParticle::get_p(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_pt",  "hasMissingET ? ReconstructedParticle::get_pt(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_px", "hasMissingET ? ReconstructedParticle::get_px(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_py", "hasMissingET ? ReconstructedParticle::get_py(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_pz", "hasMissingET ? ReconstructedParticle::get_pz(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_eta", "hasMissingET ? ReconstructedParticle::get_eta(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_theta",  "hasMissingET ? ReconstructedParticle::get_theta(MissingET)[0] : 0.f")
            .Define("RecoMissingEnergy_phi",  "hasMissingET ? ReconstructedParticle::get_phi(MissingET)[0] : 0.f")

            # Hit information for the tracks 
            .Define("RecoForTracks", "ReconstructedParticle2Track::recoParticleIndices_forTracks(TrackStates,ReconstructedParticles)")
            .Define("RecoIndices", "Utils::index_range(ReconstructedParticles)")
            .Define("RecoParticles_hitsOnTrack","VertexingUtils::getHitsOnTrack(" \
            "RecoIndices, ReconstructedParticles, Particle, MCRecoAssociations0, MCRecoAssociations1)")
            .Define("RecoParticles_hitPatterns","VertexingUtils::toHitPatterns(RecoParticles_hitsOnTrack)")
            .Define("RecoParticles_firstHitLoc","VertexingUtils::getFirstHits(RecoParticles_hitPatterns)")
            .Define("RecoParticles_lastHitLoc","VertexingUtils::getLastHits(RecoParticles_hitPatterns)")
            .Define("RecoParticles_nHits","VertexingUtils::getNHits(RecoParticles_hitPatterns)")
            .Define("RecoParticles_nDriftChamberHits","VertexingUtils::getNDCHits(RecoParticles_hitPatterns)")
            # for the kink vertices, let's ask for at least 10 hits on track, of which at least 8 in the drift chamber. 
            .Define("RecoParticles_passNhits_forKV","VertexingUtils::passHitCount(RecoParticles_hitPatterns, 10,8)")


            # --------------------------
            # Reconstructed particles
            # --------------------------
            # MC Primary Vertex
            # .Define("MC_PrimaryVertex",  "MCParticle::get_EventPrimaryVertex(21)( Particle )" )

            # Tracks
            .Define("n_RecoTracks",  "ReconstructedParticle2Track::getTK_n(TrackStates)")
            .Define("AcceptedTracks", "VertexFitterSimple::getSelectedTracks(TrackStates, 2.0)")
            .Define("n_AcceptedTracks",  "ReconstructedParticle2Track::getTK_n(AcceptedTracks)")

            # Select the tracks that are reconstructed  as primaries, so it removes all long lived decay tracks
            # this would be our stau tracks
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( AcceptedTracks, true, 4.5, 20e-3, 300, 0., 0., 0.)")
            .Define("RecoedPrimaryTracks_phi",  "ReconstructedParticle2Track::getRP2TRK_phi(ReconstructedParticles, RecoedPrimaryTracks)")
            .Define("RecoedPrimaryTracks_d0", "ReconstructedParticle2Track::getRP2TRK_D0(ReconstructedParticles, RecoedPrimaryTracks)")
            # .Define("RecoedPrimaryTracks_theta",  "ReconstructedParticle2Track::getRP2TRK_theta(ReconstructedParticles, RecoedPrimaryTracks)")
            .Define("RecoedPrimaryTracks_p",  "ReconstructedParticle2Track::getRP2TRK_mom(ReconstructedParticles, RecoedPrimaryTracks)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")

            .Define("Recoed_Primary_D0sig", "ReconstructedParticle2Track::getRP2TRK_D0_sig(ReconstructedParticles, RecoedPrimaryTracks)")
            .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300 ) ") 
            .Define("PrimaryVertex",   "VertexingUtils::get_VertexData( PrimaryVertexObject )")
            .Define("PrimaryVertex_ntracks", "FCCAnalyses::VertexingUtils::get_VertexNtrk( PrimaryVertexObject )")

            ############################################################################################
            # Additional filters:
            .Define("RecoElectronsAbove10","RecoElectrons[RecoElectrons_p > 10.]")
            .Define("n_RecoElectronsAbove10", "ReconstructedParticle::get_n(RecoElectronsAbove10)")
            .Define("RecoMuonsAbove10","RecoMuons[RecoMuons_p > 10.]")
            .Define("n_RecoMuonsAbove10", "ReconstructedParticle::get_n(RecoMuonsAbove10)")
            .Filter("n_RecoElectronsAbove10 < 2 && n_RecoMuonsAbove10 < 2")
            .Filter("n_RecoTracks < 9")
            ############################################################################################
            
            .Define("sel_tracks", "VertexFitterSimple::get_NonPrimaryTracks(AcceptedTracks, RecoedPrimaryTracks)") # 100 events/sec
            .Define("sel_tracks_D0_DV", "ReconstructedParticle2Track::getRP2TRK_D0(ReconstructedParticles, sel_tracks)")
            .Define("sel_tracks_D0sig_DV", "ReconstructedParticle2Track::getRP2TRK_D0_sig(ReconstructedParticles, sel_tracks)")
            .Define("sel_tracks_Z0_DV", "ReconstructedParticle2Track::getRP2TRK_Z0(ReconstructedParticles, sel_tracks)")
            .Define("n_nonprimary_tracks", "ReconstructedParticle2Track::getTK_n(sel_tracks)")
            .Define('sel_tracks_pt_DV', 'ReconstructedParticle2Track::getRP2TRK_mom(ReconstructedParticles ,sel_tracks)') 

            # KINK FINDER FROM SELECTED TRACKS - TARGETTING 1 PRONG TAU DECAYS 
            .Define("RecoedPrimaryTracks_charge", "ReconstructedParticle2Track::getRP2TRK_charge(ReconstructedParticles, RecoedPrimaryTracks)")
            .Define("sel_tracks_charge", "ReconstructedParticle2Track::getRP2TRK_charge(ReconstructedParticles, sel_tracks)")
            # old way:
            .Define("KinkCandidates","ReconstructedParticle2Track::findKink_candidate(ReconstructedParticles, RecoedPrimaryTracks, sel_tracks, TrackStates, RecoParticles_passNhits_forKV)")
            .Define("nKinkVertices_old", "KinkCandidates.size()")
            # this gives us the vertex object:
            .Define("KinkCandidates_VertexObject","ReconstructedParticle2Track::KinkCandidate_VertexObject(ReconstructedParticles, RecoedPrimaryTracks, sel_tracks, TrackStates, RecoParticles_passNhits_forKV )")
            .Define("nKinkVertices", "KinkCandidates_VertexObject.size()")

            .Define("KinkCandidates_passInnerHitVeto","VertexingUtils::passInnerHitVeto(KinkCandidates_VertexObject, RecoParticles_hitPatterns, RecoForTracks, true,true,false,false)")            
            .Define("nKinkCandidates_passVeto", "ROOT::VecOps::Sum(KinkCandidates_passInnerHitVeto)")
            .Define("KinkCandidates_VertexObject_passVeto", "KinkCandidates_VertexObject[KinkCandidates_passInnerHitVeto]")
            # number of tracks at the kink vertex (should be 2 for 1 prong tau decays)
            .Define("KinkVertex_ntracks", "VertexingUtils::get_VertexNtrk(KinkCandidates_VertexObject_passVeto)") 
            # invariant mass of a two track vertex   # CAUTION: m1 -> first track; m2 -> second track
            .Define("KinkVertex_invMass", "VertexingUtils::get_kink_mass(KinkCandidates_VertexObject_passVeto)")
            .Define("KinkVertex_SV" , "VertexingUtils::get_position_SV(KinkCandidates_VertexObject_passVeto)")  

            # the final primary vertex : final/refined fit

            # vector of distances of all reconstructed SV from PV (in mm in xy plane)
            .Define("KinkVertex_dxy", "VertexingUtils::get_dxy_SV(KinkCandidates_VertexObject_passVeto, PrimaryVertexObject)")
            # vector of distances of all reconstructed SV from PV (in mm in 3D)
            .Define("KinkVertex_d3d", "VertexingUtils::get_d3d_SV(KinkCandidates_VertexObject_passVeto, PrimaryVertexObject)")
            .Define("KinkCosAngle", 
                    "nKinkCandidates_passVeto > 0 ? "
                    "VertexingUtils::get_PV2V0angle(KinkCandidates_VertexObject_passVeto[0], PrimaryVertexObject) "
                    ": -2.0")  # cos of the angle between r_PVtoKV and p_KV
            .Define("KinkAngle", "TMath::ACos(KinkCosAngle)*180/3.14")            

            # DISPLACED VERTICES FROM SELECTED TRACKS- TARGETTING 3 PRONG TAU DECAY
            .Filter("sel_tracks.size()>0")
            # find the DVs from the non-primary tracks
            .Define("Displaced_Vertex", "VertexFinderLCFIPlus::get_SV_event(sel_tracks, TrackStates, PrimaryVertexObject, true, 9., 40., 5.)")
            # number of displaced vertices found
            .Define('nDisplaced_Vertices', 'VertexingUtils::get_n_SV(Displaced_Vertex)')
            # number of tracks from the DVs
            .Define('nTracks_DV', 'VertexingUtils::get_VertexNtrk(Displaced_Vertex)') 
            # momentum of the selected tracks`  `
            # Hit pattern check- only take displaced vertices which fail this criteria (if DV passes hit criteria it is actually a KV)
            # require only outgoing tracks, no incoming tracks - veto them
            .Define("DisplacedVertices_failInnerHitVeto","VertexingUtils::passInnerHitVeto(Displaced_Vertex, RecoParticles_hitPatterns, RecoForTracks, false, true, true, false)")
            # .Define("Displaced_VertexObject_failVeto", "Displaced_Vertex[!DisplacedVertices_passInnerHitVeto]")
            # fail mask = logical NOT of pass mask
            # .Define("DisplacedVertices_failInnerHitVeto", "!DisplacedVertices_passInnerHitVeto")
            # get only the vertices that failed the veto
            .Define("Displaced_VertexObject_failVeto", "Displaced_Vertex[DisplacedVertices_failInnerHitVeto]")
            # count how many failed- sum of bool
            .Define("nDisplacedVertices_failInnerHitVeto", "ROOT::VecOps::Sum(DisplacedVertices_failInnerHitVeto)")   
            .Define('nTracks_DV_failInnerHitVeto', 'VertexingUtils::get_VertexNtrk(Displaced_VertexObject_failVeto)')
           # invariant mass at the DVs (assuming the tracks to be pions)
            .Define('invMass_seltracks_DVs', 'VertexingUtils::get_invM(Displaced_VertexObject_failVeto)')
            # get the chi2 distributions of the DVs from selected tracks - to check if tracks originate from the same physical point
            .Define("DV_evt_seltracks_chi2",    "VertexingUtils::get_chi2_SV(Displaced_VertexObject_failVeto)")
            .Define("DV_evt_seltracks_normchi2","VertexingUtils::get_norm_chi2_SV(Displaced_VertexObject_failVeto)") # DV chi2 (normalised)

            # # get the decay radius and full 3D distance of all the DVs from selected tracks
            .Define("Reco_seltracks_DVs_Lxy","VertexingUtils::get_dxy_SV(Displaced_VertexObject_failVeto, PrimaryVertexObject)")
            .Define("Reco_seltracks_DVs_Lxyz","VertexingUtils::get_d3d_SV(Displaced_VertexObject_failVeto, PrimaryVertexObject)")

            # FILTERED ELECTRON COLLECTION:
            .Define("RecoElectrons_p_NEW",      "ReconstructedParticle::get_p(RecoElectronsAbove10)")
            .Define("RecoElectrons_px",     "ReconstructedParticle::get_px(RecoElectronsAbove10)")
            .Define("RecoElectrons_py",     "ReconstructedParticle::get_py(RecoElectronsAbove10)")
            .Define("RecoElectrons_pz",     "ReconstructedParticle::get_pz(RecoElectronsAbove10)")
            .Define("RecoElectrons_pt",     "ReconstructedParticle::get_pt(RecoElectronsAbove10)")
            .Define("RecoElectrons_eta",    "ReconstructedParticle::get_eta(RecoElectronsAbove10)")
            .Define("RecoElectrons_phi",    "ReconstructedParticle::get_phi(RecoElectronsAbove10)")
            .Define("RecoElectrons_e",      "ReconstructedParticle::get_e(RecoElectronsAbove10)")
            .Define("RecoElectrons_charge", "ReconstructedParticle::get_charge(RecoElectronsAbove10)")
            .Define("RecoElectrons_theta",  "ReconstructedParticle::get_theta(RecoElectronsAbove10)")

            # FILTERED MUON COLLECTION:
            .Define("RecoMuons_p_NEW",      "ReconstructedParticle::get_p(RecoMuonsAbove10)")
            .Define("RecoMuons_px",     "ReconstructedParticle::get_px(RecoMuonsAbove10)")
            .Define("RecoMuons_py",     "ReconstructedParticle::get_py(RecoMuonsAbove10)")
            .Define("RecoMuons_pz",     "ReconstructedParticle::get_pz(RecoMuonsAbove10)")
            .Define("RecoMuons_pt",     "ReconstructedParticle::get_pt(RecoMuonsAbove10)")
            .Define("RecoMuons_eta",    "ReconstructedParticle::get_eta(RecoMuonsAbove10)")
            .Define("RecoMuons_phi",    "ReconstructedParticle::get_phi(RecoMuonsAbove10)")
            .Define("RecoMuons_e",      "ReconstructedParticle::get_e(RecoMuonsAbove10)")
            .Define("RecoMuons_charge", "ReconstructedParticle::get_charge(RecoMuonsAbove10)")
            .Define("RecoMuons_theta",  "ReconstructedParticle::get_theta(RecoMuonsAbove10)")

            # Visible energy
            .Define("RecoVisibleEnergy",
                "float visE = 0.;"
                "for (size_t i = 0; i < n_RecoJets; i++) { visE += RecoJet_e[i]; }"
                "for (size_t i = 0; i < n_RecoElectrons; i++) { visE += RecoElectrons_e[i]; }"
                "for (size_t i = 0; i < n_RecoMuons; i++) { visE += RecoMuons_e[i]; }"
                "return visE;"
            )
            .Define("RecoMissingEnergy3D", "240. - RecoVisibleEnergy")
        
            # .Define("RecoTauTracks", "VertexingUtils::get_tracksInJets(Jet, _EFlowTrack_trackStates, Jet_to_Track_indices, 0)")
            # .Define("RecoTauDecayVertexObject", "VertexFitterSimple::VertexFitter_Tk(2, RecoTauTracks)")
            # .Define("RecoTauDecayVertex", "VertexingUtils::get_VertexData(RecoTauDecayVertexObject)")
            # .Define("RecoTau_Lxy", "sqrt(RecoTauDecayVertex.position.x*RecoTauDecayVertex.position.x + RecoTauDecayVertex.position.y*RecoTauDecayVertex.position.y)")
            # .Define("RecoTau_Lxyz", "sqrt(RecoTauDecayVertex.position.x*RecoTauDecayVertex.position.x + RecoTauDecayVertex.position.y*RecoTauDecayVertex.position.y + RecoTauDecayVertex.position.z*RecoTauDecayVertex.position.z)")
        )

        return df2

# Mandatory: output function
    def output(self):
        '''
        Output branches to be saved in the final ROOT file.
        '''
        branch_list = [
            # Gen-level stau
            # "debug_xpos",
            # "debug_ypos",
            # "debug_zpos",
            "GenStau",
            "n_GenStau",
            "GenStau_status",
            "GenStau_vx",
            "GenStau_vy",
            "GenStau_vz",
            "GenStau_Lxy",
            "GenStau_e",
            "GenStau_phi",
            "GenStau_Lxyz",

            # # Stau daughters
            "GenTau",
            "GenGravitino",
            "n_GenGravitino",
            "GenGravitino_status",
            "GenGravitino_e",
            # Tau kinematics
            "GenTau_e",
            "GenTau_pt",
            "GenTau_eta",
            "GenTau_phi",
            "GenTau_px",
            "GenTau_py",
            "GenTau_pz",
            "GenTau_charge",

            # Tau vertex
            "GenTau_vx",
            "GenTau_vy",
            "GenTau_vz",
            "decayLengthTau",

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

            # Neutrinos
            "FSGenNeutrino",
            "n_FSGenNeutrino",
            "FSGenNeutrino_e",
            "FSGenNeutrino_p",
            "FSGenNeutrino_pt",
            "FSGenNeutrino_px",
            "FSGenNeutrino_py",
            "FSGenNeutrino_pz",
            "FSGenNeutrino_eta",
            "FSGenNeutrino_theta",
            "FSGenNeutrino_phi",
            "FSGenNeutrino_charge",

            # Track information
            "n_RecoedPrimaryTracks", 
            "n_AcceptedTracks",
            "PrimaryVertex_ntracks",
            "n_RecoTracks",
            "sel_tracks_pt_DV",
            "sel_tracks_D0_DV",
            "sel_tracks_Z0_DV",
            "n_nonprimary_tracks",
            "Recoed_Primary_D0sig",
            "sel_tracks_D0sig_DV",
            
            ## DV information
            "nDisplaced_Vertices",
            "nTracks_DV",
            "nDisplacedVertices_failInnerHitVeto",

            "nTracks_DV_failInnerHitVeto",
            "invMass_seltracks_DVs",
            "DV_evt_seltracks_chi2",
            "DV_evt_seltracks_normchi2",
            "Reco_seltracks_DVs_Lxy",
            "Reco_seltracks_DVs_Lxyz",

            # Hit pattern information
            "RecoParticles_firstHitLoc",
            "RecoParticles_lastHitLoc",
            "RecoParticles_nHits",
            "RecoParticles_nDriftChamberHits",

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
            "RecoElectrons_p_NEW",
            "RecoElectrons_pt",
            "RecoElectrons_px",
            "RecoElectrons_py",
            "RecoElectrons_pz",
            "RecoElectrons_eta",
            "RecoElectrons_theta",
            "RecoElectrons_phi",
            "RecoElectrons_charge",
            # "Reco_ee_energy",
            # "Reco_ee_px",
            # "Reco_ee_py",
            # "Reco_ee_pz",
            # "Reco_ee_invMass",

            # Reco Muons
            "RecoMuons",
            "n_RecoMuons",
            "RecoMuons_e",
            "RecoMuons_p",
            "RecoMuons_p_NEW",
            "RecoMuons_pt",
            "RecoMuons_px",
            "RecoMuons_py",
            "RecoMuons_pz",
            "RecoMuons_eta",
            "RecoMuons_theta",
            "RecoMuons_phi",
            "RecoMuons_charge",
            # "Reco_mumu_energy",
            # "Reco_mumu_px",
            # "Reco_mumu_py",
            # "Reco_mumu_pz",
            # "Reco_mumu_invMass",
            # "muon_electron_overlap",

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

            # Time variables
            "GenStau_time",
            "GenTau_time",
            "FSGenElectron_time",
            "FSGenMuon_time",
            "GenTau_status",
            "GenTau_cTau",
            "GenStau_theta",
            "GenTau_theta",

            "RecoVisibleEnergy",
            "RecoMissingEnergy3D",
            "RecoedPrimaryTracks_charge",
            "sel_tracks_charge",
            "RecoedPrimaryTracks_d0",
            "RecoedPrimaryTracks_phi",
            # "RecoedPrimaryTracks_theta",
            "RecoedPrimaryTracks_p",
            # "GenStau_daughters",
            # "GenTau_daughters",
            ## Kinked candidates :
         
            "KinkCandidates",
            "KinkCandidates_passInnerHitVeto",
            "nKinkCandidates_passVeto",
            "KinkVertex_invMass",
            "nKinkVertices",
            "KinkVertex_SV",
            "KinkVertex_ntracks",
            "KinkAngle",
            "KinkVertex_dxy",
            "KinkVertex_d3d",

        ]

        return branch_list
