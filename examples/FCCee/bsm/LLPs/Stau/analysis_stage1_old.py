'''
Analysis for ALP study: stau⁺stau⁻ → tau + gravitino → hadronic final states
'''
from argparse import ArgumentParser

import ROOT

ROOT.gInterpreter.Declare("""
template<typename T>
auto SafeAlias(T df, const std::string &alias, const std::string &target) {
    auto cols = df.GetColumnNames();
    auto pos = target.find(".");
    std::string base = (pos == std::string::npos ? target : target.substr(0, pos));

    if (cols.contains(base)) {
        return df.Alias(alias, target);
    } else {
        std::cout << "Warning: skipping alias " << alias
                  << " because " << base << " does not exist" << std::endl;
        return df;
    }
}
""");


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
        '''
        df = df.Alias("MCRecoAssociations0", "_RecoMCLink_from.index")
        df = df.Alias("MCRecoAssociations1", "_RecoMCLink_to.index")

        df2 = (
            df
            .Alias("Particle1", "_Particle_daughters.index")
            .Define("GenStau","FCCAnalyses::MCParticle::sel_pdgID(1000015,true)(Particle)")
            .Define("n_GenStau",   "FCCAnalyses::MCParticle::get_n(GenStau)") #number of staus
            # .Define("GenStau_observed_lifetime_xyz", "{ if(n_GenStau>0) return ((FCCAnalyses::MCParticle::get_p(GenStau) / FCCAnalyses::MCParticle::get_mass(GenStau)) * FCCAnalyses::MCParticle::get_lifetime_xyz(GenStau)); else return -2.0f; }")
            .Define("GenStau_observed_lifetime_xyz","(GenStau_p / GenStau_mass) * GenStau_lifetime_xyz")


            
            #gen level decay vertex of stau
            .Define("GenStau_vx", "FCCAnalyses::MCParticle::get_vertex_x(GenStau)")
            .Define("GenStau_vy", "FCCAnalyses::MCParticle::get_vertex_y(GenStau)")
            .Define("GenStau_vz", "FCCAnalyses::MCParticle::get_vertex_z(GenStau)")

            .Define("GenStau_Lxy",  "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy)")
            .Define("GenStau_Lxyz", "sqrt(GenStau_vx*GenStau_vx + GenStau_vy*GenStau_vy + GenStau_vz*GenStau_vz)")

            #select daughters of stau - tau and a gravitino
            .Define("StauDaughters", "FCCAnalyses::MCParticle::get_daughters(GenStau)")  
            .Define("GenTau",        "FCCAnalyses::MCParticle::sel_pdgID(15, true)(StauDaughters)")
            .Define("GenGravitino",  "FCCAnalyses::MCParticle::sel_pdgID(1000049, false)(StauDaughters)")

            #Tau 4-vectors
            .Define("GenTau_e",  "FCCAnalyses::MCParticle::get_e(GenTau)")
            .Define("GenTau_pt", "FCCAnalyses::MCParticle::get_pt(GenTau)")
            .Define("GenTau_eta","FCCAnalyses::MCParticle::get_eta(GenTau)")
            .Define("GenTau_phi","FCCAnalyses::MCParticle::get_phi(GenTau)")
            .Define("GenTau_px", "FCCAnalyses::MCParticle::get_px(GenTau)")
            .Define("GenTau_py", "FCCAnalyses::MCParticle::get_py(GenTau)")
            .Define("GenTau_pz", "FCCAnalyses::MCParticle::get_pz(GenTau)")

            #Tau decay vertex
            .Define("GenTau_vx", "FCCAnalyses::MCParticle::get_vertex_x(GenTau)")
            .Define("GenTau_vy", "FCCAnalyses::MCParticle::get_vertex_y(GenTau)")
            .Define("GenTau_vz", "FCCAnalyses::MCParticle::get_vertex_z(GenTau)")

            # final state particles: genStatus(1)
            # 1. HADRONIC decay
            .Define("GenPion_PID", "FCCAnalyses::MCParticle::sel_pdgID(211, true)(Particle)")
            .Define("FSGenPion",   "FCCAnalyses::MCParticle::sel_genStatus(1)(GenPion_PID)")

            # 2. MUONS and ELECTRONS
            .Define("GenMuon_PID", "FCCAnalyses::MCParticle::sel_pdgID(13, true)(Particle)")
            .Define("FSGenMuon",   "FCCAnalyses::MCParticle::sel_genStatus(1)(GenMuon_PID)")
            .Define("GenElectron_PID", "FCCAnalyses::MCParticle::sel_pdgID(11, true)(Particle)")
            .Define("FSGenElectron",   "FCCAnalyses::MCParticle::sel_genStatus(1)(GenElectron_PID)")

            # 3. NEUTRINOS
            .Define("GenNuE",  "FCCAnalyses::MCParticle::sel_pdgID(12,true)(Particle)")
            .Define("GenNuMu", "FCCAnalyses::MCParticle::sel_pdgID(14,true)(Particle)")
            .Define("GenNuTau", "FCCAnalyses::MCParticle::sel_pdgID(16,true)(Particle)")
            .Define("GenNeutrinos", "Concatenate(GenNuE, GenNuMu, GenNuTau)")

            # Kinematics of final state particles - taken directly from DisplacedHNL
            # ELECTRON:
            .Define("n_FSGenElectron", "FCCAnalyses::MCParticle::get_n(FSGenElectron)")
            .Define("FSGenElectron_e",      "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_e(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_p",      "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_p(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_pt",     "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_pt(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_px",     "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_px(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_py",     "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_py(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_pz",     "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_pz(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_eta",    "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_eta(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_theta",  "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_theta(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_phi",    "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_phi(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_charge", "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_charge(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_vertex_x", "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_vertex_x(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_vertex_y", "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_vertex_y(FSGenElectron); else return -2.0f; }")
            .Define("FSGenElectron_vertex_z", "{ if(n_FSGenElectron>0) return FCCAnalyses::MCParticle::get_vertex_z(FSGenElectron); else return -2.0f; }")

            # MUON:
            .Define("n_FSGenMuon", "FCCAnalyses::MCParticle::get_n(FSGenMuon)")
            .Define("FSGenMuon_e",      "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_e(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_p",      "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_p(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_pt",     "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_pt(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_px",     "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_px(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_py",     "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_py(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_pz",     "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_pz(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_eta",    "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_eta(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_theta",  "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_theta(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_phi",    "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_phi(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_charge", "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_charge(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_vertex_x", "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_vertex_x(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_vertex_y", "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_vertex_y(FSGenMuon); else return -2.0f; }")
            .Define("FSGenMuon_vertex_z", "{ if(n_FSGenMuon>0) return FCCAnalyses::MCParticle::get_vertex_z(FSGenMuon); else return -2.0f; }")

            # PION:
            .Define("n_FSGenPion", "FCCAnalyses::MCParticle::get_n(FSGenPion)")
            .Define("FSGenPion_e",      "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_e(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_p",      "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_p(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_pt",     "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_pt(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_px",     "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_px(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_py",     "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_py(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_pz",     "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_pz(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_eta",    "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_eta(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_theta",  "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_theta(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_phi",    "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_phi(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_charge", "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_charge(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_vertex_x", "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_vertex_x(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_vertex_y", "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_vertex_y(FSGenPion); else return -2.0f; }")
            .Define("FSGenPion_vertex_z", "{ if(n_FSGenPion>0) return FCCAnalyses::MCParticle::get_vertex_z(FSGenPion); else return -2.0f; }")

            # Reconstructed tracks of final state particles
            # JETS:
            .Define("n_RecoJets", "ReconstructedParticle::get_n(Jet)") # count how many jets are in the event in total
            .Define("RecoJet_e",      "{ if(n_RecoJets>0) return ReconstructedParticle::get_e(Jet); else return -2.0f; }")
            .Define("RecoJet_p",      "{ if(n_RecoJets>0) return ReconstructedParticle::get_p(Jet); else return -2.0f; }")
            .Define("RecoJet_pt",     "{ if(n_RecoJets>0) return ReconstructedParticle::get_pt(Jet); else return -2.0f; }")
            .Define("RecoJet_px",     "{ if(n_RecoJets>0) return ReconstructedParticle::get_px(Jet); else return -2.0f; }")
            .Define("RecoJet_py",     "{ if(n_RecoJets>0) return ReconstructedParticle::get_py(Jet); else return -2.0f; }")
            .Define("RecoJet_pz",     "{ if(n_RecoJets>0) return ReconstructedParticle::get_pz(Jet); else return -2.0f; }")
            .Define("RecoJet_eta",    "{ if(n_RecoJets>0) return ReconstructedParticle::get_eta(Jet); else return -2.0f; }")
            .Define("RecoJet_theta",  "{ if(n_RecoJets>0) return ReconstructedParticle::get_theta(Jet); else return -2.0f; }")
            .Define("RecoJet_phi",    "{ if(n_RecoJets>0) return ReconstructedParticle::get_phi(Jet); else return -2.0f; }")
            .Define("RecoJet_charge", "{ if(n_RecoJets>0) return ReconstructedParticle::get_charge(Jet); else return -2.0f; }")

            # ELECTRON:
            .Alias("Electron0", "Electron#0.index")            
            .Define("n_RecoElectrons", "ReconstructedParticle::get_n(RecoElectrons)")
            .Define("RecoElectrons", "{ if(n_RecoElectrons>0) return ReconstructedParticle::get(Electron0, ReconstructedParticles); else return -2.0f; }")
            .Define("RecoElectron_e",      "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_e(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_p",      "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_p(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_pt",     "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_pt(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_px",     "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_px(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_py",     "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_py(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_pz",     "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_pz(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_eta",    "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_eta(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_theta",  "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_theta(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_phi",    "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_phi(RecoElectrons); else return -2.0f; }")
            .Define("RecoElectron_charge", "{ if(n_RecoElectrons>0) return ReconstructedParticle::get_charge(RecoElectrons); else return -2.0f; }")

            # MUON:
            .Alias("Muon0", "Muon#0.index")
            .Define("n_RecoMuons", "ReconstructedParticle::get_n(RecoMuons)")
            .Define("RecoMuons", "{ if(n_RecoMuons>0) return ReconstructedParticle::get(Muon0, ReconstructedParticles); else return -2.0f; }")
            .Define("RecoMuons_e",      "{ if(n_RecoMuons>0) return ReconstructedParticle::get_e(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_p",      "{ if(n_RecoMuons>0) return ReconstructedParticle::get_p(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_pt",     "{ if(n_RecoMuons>0) return ReconstructedParticle::get_pt(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_px",     "{ if(n_RecoMuons>0) return ReconstructedParticle::get_px(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_py",     "{ if(n_RecoMuons>0) return ReconstructedParticle::get_py(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_pz",     "{ if(n_RecoMuons>0) return ReconstructedParticle::get_pz(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_eta",    "{ if(n_RecoMuons>0) return ReconstructedParticle::get_eta(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_theta",  "{ if(n_RecoMuons>0) return ReconstructedParticle::get_theta(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_phi",    "{ if(n_RecoMuons>0) return ReconstructedParticle::get_phi(RecoMuons); else return -2.0f; }")
            .Define("RecoMuons_charge", "{ if(n_RecoMuons>0) return ReconstructedParticle::get_charge(RecoMuons); else return -2.0f; }")

            # MET:
            .Define("RecoMissingEnergy_e", "ReconstructedParticle::get_e(MissingET)") 
            .Define("RecoMissingEnergy_p", "ReconstructedParticle::get_p(MissingET)")
            .Define("RecoMissingEnergy_pt", "ReconstructedParticle::get_pt(MissingET)")
            .Define("RecoMissingEnergy_px", "ReconstructedParticle::get_px(MissingET)") 
            .Define("RecoMissingEnergy_py", "ReconstructedParticle::get_py(MissingET)") 
            .Define("RecoMissingEnergy_pz", "ReconstructedParticle::get_pz(MissingET)")
            .Define("RecoMissingEnergy_eta", "ReconstructedParticle::get_eta(MissingET)")
            .Define("RecoMissingEnergy_theta", "ReconstructedParticle::get_theta(MissingET)")
            .Define("RecoMissingEnergy_phi", "ReconstructedParticle::get_phi(MissingET)")

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
