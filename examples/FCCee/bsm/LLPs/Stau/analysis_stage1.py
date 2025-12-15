# Took inspiration from the Alps - different energies- fcc analyses - stage 1
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
            'FCC_100stau_240com_000_events'  : {'fraction': 1.0},
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
        self.input_dir = '/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/'
        self.output_dir = '/eos/user/s/svashish/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/output'

        self.tree_name = "Delphes"
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
        # '''
        # df = df.Alias("MCRecoAssociations0", "_RecoMCLink_from.index")
        # df = df.Alias("MCRecoAssociations1", "_RecoMCLink_to.index")

        df2 = (
            df
            # --------------------------
            # Generator-level particles
            # --------------------------
            .Alias("Particle1", "_Particle_daughters.index")
            .Define("GenStau", "FCCAnalyses::MCParticle::sel_pdgID(1000015,true)(Particle)")
            .Define("n_GenStau", "FCCAnalyses::MCParticle::get_n(GenStau)")
            .Define("GenStau_vx", "FCCAnalyses::MCParticle::get_vertex_x(GenStau)")
            .Define("GenStau_vy", "FCCAnalyses::MCParticle::get_vertex_y(GenStau)")
            .Define("GenStau_vz", "FCCAnalyses::MCParticle::get_vertex_z(GenStau)")
            .Define("GenStau_Lxy", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy)")
            .Define("GenStau_Lxyz", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy + GenStau_vz*GenStau_vz)")

            # Stau daughters
            .Define("StauDaughters", "FCCAnalyses::MCParticle::get_daughters(GenStau)")
            .Define("GenTau", "FCCAnalyses::MCParticle::sel_pdgID(15,true)(StauDaughters)")
            .Define("GenGravitino", "FCCAnalyses::MCParticle::sel_pdgID(1000049,false)(StauDaughters)")

            # Tau kinematics and vertex
            .Define("GenTau_px", "FCCAnalyses::MCParticle::get_px(GenTau)")
            .Define("GenTau_py", "FCCAnalyses::MCParticle::get_py(GenTau)")
            .Define("GenTau_pz", "FCCAnalyses::MCParticle::get_pz(GenTau)")
            .Define("GenTau_pt", "FCCAnalyses::MCParticle::get_pt(GenTau)")
            .Define("GenTau_eta", "FCCAnalyses::MCParticle::get_eta(GenTau)")
            .Define("GenTau_phi", "FCCAnalyses::MCParticle::get_phi(GenTau)")
            .Define("GenTau_e", "FCCAnalyses::MCParticle::get_e(GenTau)")
            .Define("GenTau_charge", "FCCAnalyses::MCParticle::get_charge(GenTau)")
            .Define("GenTau_vx", "FCCAnalyses::MCParticle::get_vertex_x(GenTau)")
            .Define("GenTau_vy", "FCCAnalyses::MCParticle::get_vertex_y(GenTau)")
            .Define("GenTau_vz", "FCCAnalyses::MCParticle::get_vertex_z(GenTau)")

            # Final-state generator particles (status 1)
            .Define("FSGenElectron", "FCCAnalyses::MCParticle::sel_genStatus(1)(FCCAnalyses::MCParticle::sel_pdgID(11,true)(ParticleMC))")
            .Define("FSGenMuon",     "FCCAnalyses::MCParticle::sel_genStatus(1)(FCCAnalyses::MCParticle::sel_pdgID(13,true)(ParticleMC))")
            .Define("FSGenTau",      "FCCAnalyses::MCParticle::sel_genStatus(1)(FCCAnalyses::MCParticle::sel_pdgID(15,true)(ParticleMC))")
            .Define("FSGenPion",     "FCCAnalyses::MCParticle::sel_genStatus(1)(FCCAnalyses::MCParticle::sel_pdgID(211,true)(ParticleMC))")

            # Kinematics helper for FSGen
            .Define("FSGenElectron_px","FCCAnalyses::MCParticle::get_px(FSGenElectron)")
            .Define("FSGenElectron_py","FCCAnalyses::MCParticle::get_py(FSGenElectron)")
            .Define("FSGenElectron_pz","FCCAnalyses::MCParticle::get_pz(FSGenElectron)")
            .Define("FSGenElectron_pt","FCCAnalyses::MCParticle::get_pt(FSGenElectron)")
            .Define("FSGenElectron_eta","FCCAnalyses::MCParticle::get_eta(FSGenElectron)")
            .Define("FSGenElectron_phi","FCCAnalyses::MCParticle::get_phi(FSGenElectron)")
            .Define("FSGenElectron_e","FCCAnalyses::MCParticle::get_e(FSGenElectron)")
            .Define("FSGenElectron_charge","FCCAnalyses::MCParticle::get_charge(FSGenElectron)")
            .Define("FSGenElectron_vx","FCCAnalyses::MCParticle::get_vertex_x(FSGenElectron)")
            .Define("FSGenElectron_vy","FCCAnalyses::MCParticle::get_vertex_y(FSGenElectron)")
            .Define("FSGenElectron_vz","FCCAnalyses::MCParticle::get_vertex_z(FSGenElectron)")

            .Define("FSGenMuon_px","FCCAnalyses::MCParticle::get_px(FSGenMuon)")
            .Define("FSGenMuon_py","FCCAnalyses::MCParticle::get_py(FSGenMuon)")
            .Define("FSGenMuon_pz","FCCAnalyses::MCParticle::get_pz(FSGenMuon)")
            .Define("FSGenMuon_pt","FCCAnalyses::MCParticle::get_pt(FSGenMuon)")
            .Define("FSGenMuon_eta","FCCAnalyses::MCParticle::get_eta(FSGenMuon)")
            .Define("FSGenMuon_phi","FCCAnalyses::MCParticle::get_phi(FSGenMuon)")
            .Define("FSGenMuon_e","FCCAnalyses::MCParticle::get_e(FSGenMuon)")
            .Define("FSGenMuon_charge","FCCAnalyses::MCParticle::get_charge(FSGenMuon)")
            .Define("FSGenMuon_vx","FCCAnalyses::MCParticle::get_vertex_x(FSGenMuon)")
            .Define("FSGenMuon_vy","FCCAnalyses::MCParticle::get_vertex_y(FSGenMuon)")
            .Define("FSGenMuon_vz","FCCAnalyses::MCParticle::get_vertex_z(FSGenMuon)")

            .Define("FSGenTau_px","FCCAnalyses::MCParticle::get_px(FSGenTau)")
            .Define("FSGenTau_py","FCCAnalyses::MCParticle::get_py(FSGenTau)")
            .Define("FSGenTau_pz","FCCAnalyses::MCParticle::get_pz(FSGenTau)")
            .Define("FSGenTau_pt","FCCAnalyses::MCParticle::get_pt(FSGenTau)")
            .Define("FSGenTau_eta","FCCAnalyses::MCParticle::get_eta(FSGenTau)")
            .Define("FSGenTau_phi","FCCAnalyses::MCParticle::get_phi(FSGenTau)")
            .Define("FSGenTau_e","FCCAnalyses::MCParticle::get_e(FSGenTau)")
            .Define("FSGenTau_charge","FCCAnalyses::MCParticle::get_charge(FSGenTau)")
            .Define("FSGenTau_vx","FCCAnalyses::MCParticle::get_vertex_x(FSGenTau)")
            .Define("FSGenTau_vy","FCCAnalyses::MCParticle::get_vertex_y(FSGenTau)")
            .Define("FSGenTau_vz","FCCAnalyses::MCParticle::get_vertex_z(FSGenTau)")

            .Define("FSGenPion_px","FCCAnalyses::MCParticle::get_px(FSGenPion)")
            .Define("FSGenPion_py","FCCAnalyses::MCParticle::get_py(FSGenPion)")
            .Define("FSGenPion_pz","FCCAnalyses::MCParticle::get_pz(FSGenPion)")
            .Define("FSGenPion_pt","FCCAnalyses::MCParticle::get_pt(FSGenPion)")
            .Define("FSGenPion_eta","FCCAnalyses::MCParticle::get_eta(FSGenPion)")
            .Define("FSGenPion_phi","FCCAnalyses::MCParticle::get_phi(FSGenPion)")
            .Define("FSGenPion_e","FCCAnalyses::MCParticle::get_e(FSGenPion)")
            .Define("FSGenPion_charge","FCCAnalyses::MCParticle::get_charge(FSGenPion)")
            .Define("FSGenPion_vx","FCCAnalyses::MCParticle::get_vertex_x(FSGenPion)")
            .Define("FSGenPion_vy","FCCAnalyses::MCParticle::get_vertex_y(FSGenPion)")
            .Define("FSGenPion_vz","FCCAnalyses::MCParticle::get_vertex_z(FSGenPion)")

            # --------------------------
            # Reconstructed particles
            # --------------------------
            # Electrons
            .Alias("Electron0", "Electron#0.index")
            .Define("RecoElectrons_px", "ReconstructedParticle::get_px(Electron)")
            .Define("RecoElectrons_py", "ReconstructedParticle::get_py(Electron)")
            .Define("RecoElectrons_pz", "ReconstructedParticle::get_pz(Electron)")
            .Define("RecoElectrons_pt", "ReconstructedParticle::get_pt(Electron)")
            .Define("RecoElectrons_eta", "ReconstructedParticle::get_eta(Electron)")
            .Define("RecoElectrons_phi", "ReconstructedParticle::get_phi(Electron)")
            .Define("RecoElectrons_e", "ReconstructedParticle::get_e(Electron)")
            .Define("RecoElectrons_charge", "ReconstructedParticle::get_charge(Electron)")
            .Define("RecoElectrons_vx", "ReconstructedParticle::get_vertex_x(Electron)")
            .Define("RecoElectrons_vy", "ReconstructedParticle::get_vertex_y(Electron)")
            .Define("RecoElectrons_vz", "ReconstructedParticle::get_vertex_z(Electron)")

            # Muons
            .Alias("Muon0", "Muon#0.index")
            .Define("RecoMuons_px", "ReconstructedParticle::get_px(Muon)")
            .Define("RecoMuons_py", "ReconstructedParticle::get_py(Muon)")
            .Define("RecoMuons_pz", "ReconstructedParticle::get_pz(Muon)")
            .Define("RecoMuons_pt", "ReconstructedParticle::get_pt(Muon)")
            .Define("RecoMuons_eta", "ReconstructedParticle::get_eta(Muon)")
            .Define("RecoMuons_phi", "ReconstructedParticle::get_phi(Muon)")
            .Define("RecoMuons_e", "ReconstructedParticle::get_e(Muon)")
            .Define("RecoMuons_charge", "ReconstructedParticle::get_charge(Muon)")
            .Define("RecoMuons_vx", "ReconstructedParticle::get_vertex_x(Muon)")
            .Define("RecoMuons_vy", "ReconstructedParticle::get_vertex_y(Muon)")
            .Define("RecoMuons_vz", "ReconstructedParticle::get_vertex_z(Muon)")

            # Taus
            .Alias("Tau0", "Tau#0.index")
            .Define("RecoTaus_px", "ReconstructedParticle::get_px(Tau)")
            .Define("RecoTaus_py", "ReconstructedParticle::get_py(Tau)")
            .Define("RecoTaus_pz", "ReconstructedParticle::get_pz(Tau)")
            .Define("RecoTaus_pt", "ReconstructedParticle::get_pt(Tau)")
            .Define("RecoTaus_eta", "ReconstructedParticle::get_eta(Tau)")
            .Define("RecoTaus_phi", "ReconstructedParticle::get_phi(Tau)")
            .Define("RecoTaus_e", "ReconstructedParticle::get_e(Tau)")
            .Define("RecoTaus_charge", "ReconstructedParticle::get_charge(Tau)")
            .Define("RecoTaus_vx", "ReconstructedParticle::get_vertex_x(Tau)")
            .Define("RecoTaus_vy", "ReconstructedParticle::get_vertex_y(Tau)")
            .Define("RecoTaus_vz", "ReconstructedParticle::get_vertex_z(Tau)")

            # Pions
            .Define("RecoPions_px", "ReconstructedParticle::get_px(ParticleFlowCandidate)")
            .Define("RecoPions_py", "ReconstructedParticle::get_py(ParticleFlowCandidate)")
            .Define("RecoPions_pz", "ReconstructedParticle::get_pz(ParticleFlowCandidate)")
            .Define("RecoPions_pt", "ReconstructedParticle::get_pt(ParticleFlowCandidate)")
            .Define("RecoPions_eta", "ReconstructedParticle::get_eta(ParticleFlowCandidate)")
            .Define("RecoPions_phi", "ReconstructedParticle::get_phi(ParticleFlowCandidate)")
            .Define("RecoPions_e", "ReconstructedParticle::get_e(ParticleFlowCandidate)")
            .Define("RecoPions_charge", "ReconstructedParticle::get_charge(ParticleFlowCandidate)")
            .Define("RecoPions_vx", "ReconstructedParticle::get_vertex_x(ParticleFlowCandidate)")
            .Define("RecoPions_vy", "ReconstructedParticle::get_vertex_y(ParticleFlowCandidate)")
            .Define("RecoPions_vz", "ReconstructedParticle::get_vertex_z(ParticleFlowCandidate)")
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
            "GenStau_observed_lifetime_xyz",
            "GenStau_vx",
            "GenStau_vy",
            "GenStau_vz",
            "GenStau_Lxy",
            "GenStau_Lxyz",

            # Stau daughters
            "StauDaughters",
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

            # Tau vertex
            "GenTau_vx",
            "GenTau_vy",
            "GenTau_vz",

            # Final-state pions
            "FSGenPion",
            "n_FSGenPion",
            "FSGenPion_e",
            "FSGenPion_p",
            "FSGenPion_pt",
            "FSGenPion_px",
            "FSGenPion_py",
            "FSGenPion_pz",
            "FSGenPion_eta",
            "FSGenPion_theta",
            "FSGenPion_phi",
            "FSGenPion_charge",
            "FSGenPion_vertex_x",
            "FSGenPion_vertex_y",
            "FSGenPion_vertex_z",

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
            "FSGenMuon_vertex_x",
            "FSGenMuon_vertex_y",
            "FSGenMuon_vertex_z",

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
            "FSGenElectron_vertex_x",
            "FSGenElectron_vertex_y",
            "FSGenElectron_vertex_z",

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

            # Reco Electrons
            "RecoElectrons",
            "n_RecoElectrons",
            "RecoElectron_e",
            "RecoElectron_p",
            "RecoElectron_pt",
            "RecoElectron_px",
            "RecoElectron_py",
            "RecoElectron_pz",
            "RecoElectron_eta",
            "RecoElectron_theta",
            "RecoElectron_phi",
            "RecoElectron_charge",

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
