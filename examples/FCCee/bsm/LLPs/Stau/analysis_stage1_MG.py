# import ROOT
from analysis_stage1 import Analysis as OriginalAnalysis
class Analysis(OriginalAnalysis):
    def __init__(self, cmdline_args):
        super().__init__(cmdline_args)
        self.process_list = {
            # uncomment the samples you want to run over
            #######################################################
            #               CME: 240 GeV (ZH)                     #
            #######################################################
            # 'FCCee_100_stau_10mm_ctau_ecm_240'  : {'fraction': 1.0},
            # 'FCCee_110_stau_10mm_ctau_ecm_240'  : {'fraction': 1.0},

            # "FCCee_110_stau_20cm_ctau_ecm_240"  : {'fraction': 1.0},
            # "FCCee_110_stau_1p5m_ctau_ecm_240"  : {'fraction': 1.0},
            "FCCee_110_stau_3m_ctau_ecm_240_changed_delphes"  : {'fraction': 1.0},
            
            #######################################################
            #               SPRING 2021                           #
            #######################################################
            # 'p8_ee_WW_mumu_ecm240':{'fraction': 0.5,'chunks':100},

            #######################################################
            #               WINTER 2023                           #
            #######################################################
            # 'p8_ee_WW_ecm240': {'fraction': 0.01,'chunks':100},
            # "p8_ee_ZZ_ecm240": {'fraction': 0.1,'chunks':100},
            # "mgp8_ee_zh_ecm240_hbb": {'fraction': 0.5,'chunks':100},
            # 'wzp6_ee_nuenueH_Htautau_ecm240': {'fraction': 0.5,'chunks':100},
            # 'wzp6_ee_bbH_Htautau_ecm240': {'fraction': 0.5,'chunks':100},


            # 'p8_ee_WW_ee_ecm240': {'fraction': 0.1,'chunks':100},
            # 'p8_ee_Zbb_ecm91_EvtGen_Bd2KstarTauTau' : {'fraction': 0.5},
            } 

        # prodTag     = "FCCee/winter2023/IDEA/"
        # prodTag     = "FCCee/spring2021/IDEA/"

        # Input/output directories
        # uncomment the input directory you want to use

        # backgrounds are stored here:
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA'
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023_training/IDEA/'
        # self.input_dir = '/eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/'

        self.input_dir = '/home/goblirsc/Code/FccStudies/Delphes_dev/out'
        self.output_dir = '/home/goblirsc/Code/FccStudies/Delphes_dev/output_stage1'

        # Optional: read the input files with podio::DataSource
        # self.use_data_source = True
    def analyzers(self, df):
        df2 = super().analyzers(df)
        df2 = (df2 
            .Define("GenStau_p", "MCParticle::get_p(GenStau)")
            .Define("GenStau_m", "MCParticle::get_mass(GenStau)")
            .Define("GenStau_beta", "GenStau_p/GenStau_e")
            .Define("GenStau_gamma", "GenStau_e/GenStau_m")
            .Define("ParticleIndices", "Utils::index_range(Particle)")
            .Define("RecoIndices", "Utils::index_range(ReconstructedParticles)")
            .Define("GenStauIndices","MCParticle::sel_pdgID(1000015,true)(ParticleIndices,Particle)")
            .Define("GenStauIndicesStat22","MCParticle::sel_genStatus(22)(GenStauIndices,Particle)")
            .Define("StauPDG","Particle.PDG[GenStauIndicesStat22]")
            .Define("StauChildrenIndices1","MCParticle::get_list_of_stable_particles_from_decay(GenStauIndicesStat22[0], Particle, ParticleIndices)")
            .Define("StauChildrenIndices2","MCParticle::get_list_of_stable_particles_from_decay(GenStauIndicesStat22[1], Particle, ParticleIndices)")
            .Define("ChargedStauChildrenIndices1","MCParticle::sel_charged()(StauChildrenIndices1, Particle)")
            .Define("ChargedStauChildrenIndices2","MCParticle::sel_charged()(StauChildrenIndices2, Particle)")
            .Define("ChargedStauChildren1","Utils::sel_byIndex(ChargedStauChildrenIndices1, Particle)")
            .Define("ChargedStauChildren2","Utils::sel_byIndex(ChargedStauChildrenIndices2, Particle)")
            # .Define("ChargedStauChildren",
                    # "MCParticle::broadCast(MCParticle::sel_byIndex)(GenStauIndicesStat22, Particle)")
                    # "MCParticle::broadCast<edm4hep::MCParticleData,int,ROOT::VecOps::RVec<edm4hep::MCParticleData>>(MCParticle::sel_byIndex)(GenStauIndicesStat22, Particle)")
            .Define("ChargedStauChildren1_PDG","MCParticle::get_pdg(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_PDG","MCParticle::get_pdg(ChargedStauChildren2)")
            .Define("ChargedStauChildren1_q","MCParticle::get_charge(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_q","MCParticle::get_charge(ChargedStauChildren2)")
            .Define("ChargedStauChildren1_vx","MCParticle::get_vertex_x(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_vx","MCParticle::get_vertex_x(ChargedStauChildren2)")
            .Define("ChargedStauChildren1_vy","MCParticle::get_vertex_y(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_vy","MCParticle::get_vertex_y(ChargedStauChildren2)")
            .Define("ChargedStauChildren1_vz","MCParticle::get_vertex_z(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_vz","MCParticle::get_vertex_z(ChargedStauChildren2)")
            .Define("ChargedStauChildren1_n","Utils::getsize(ChargedStauChildren1)")
            .Define("ChargedStauChildren2_n","Utils::getsize(ChargedStauChildren2)")      
            .Define("Stau_nChildren","Utils::merge<unsigned long>({ChargedStauChildren1_n},{ChargedStauChildren2_n})")
            .Define("ChargedStauChildren1_recoMatch_ix","ReconstructedParticle2MC::selRP_indices_matched_to_list(ChargedStauChildrenIndices1,MCRecoAssociations0, MCRecoAssociations1, Particle, false )")
            .Define("ChargedStauChildren2_recoMatch_ix","ReconstructedParticle2MC::selRP_indices_matched_to_list(ChargedStauChildrenIndices2,MCRecoAssociations0, MCRecoAssociations1, Particle, false )")
            .Define("ChargedStauChildren1_recoMatch_n","Utils::count_valid_indices(ChargedStauChildren1_recoMatch_ix, ReconstructedParticles)")
            .Define("ChargedStauChildren2_recoMatch_n","Utils::count_valid_indices(ChargedStauChildren2_recoMatch_ix, ReconstructedParticles)")
            .Define("Stau_nRecoTracks","Utils::merge<int>({ChargedStauChildren1_recoMatch_n},{ChargedStauChildren2_recoMatch_n})")
            .Define("ChargedStauChildren1_vxMatch_ix","VertexingUtils::getVertex_matching_recoParticles(DV_evt_seltracks, ChargedStauChildren1_recoMatch_ix , ReconstructedParticles)")
            .Define("ChargedStauChildren2_vxMatch_ix","VertexingUtils::getVertex_matching_recoParticles(DV_evt_seltracks, ChargedStauChildren2_recoMatch_ix , ReconstructedParticles)")
            .Define("Stau_DVmatch_index","Utils::merge<int>({ChargedStauChildren1_vxMatch_ix},{ChargedStauChildren2_vxMatch_ix})")
            .Define("Stau1_recoMatch_ix","ReconstructedParticle2MC::selRP_indices_matched_to_list({GenStauIndicesStat22[0]}, MCRecoAssociations0, MCRecoAssociations1, Particle, false )")
            .Define("Stau2_recoMatch_ix","ReconstructedParticle2MC::selRP_indices_matched_to_list({GenStauIndicesStat22[1]},MCRecoAssociations0, MCRecoAssociations1, Particle, false )")

            .Define ("Stau1_reco_for_KinkVx_ix","Utils::merge(Stau1_recoMatch_ix, ChargedStauChildren1_recoMatch_ix)")
            .Define ("Stau2_reco_for_KinkVx_ix","Utils::merge(Stau2_recoMatch_ix, ChargedStauChildren2_recoMatch_ix)")
            .Define("Stau1_kinkVertexMatch_ix","VertexingUtils::getVertex_matching_recoParticles(KinkCandidates_VertexObject, Stau1_reco_for_KinkVx_ix , ReconstructedParticles)")
            .Define("Stau2_kinkVertexMatch_ix","VertexingUtils::getVertex_matching_recoParticles(KinkCandidates_VertexObject, Stau2_reco_for_KinkVx_ix , ReconstructedParticles)")
            .Define("Stau_KVmatch_index","Utils::merge<int>({Stau1_kinkVertexMatch_ix},{Stau2_kinkVertexMatch_ix})")

            .Define("ChargedStauChildrenIndices","Utils::merge(ChargedStauChildrenIndices1, ChargedStauChildrenIndices2)")
            .Define("ChargedStauChildren","Utils::merge(ChargedStauChildren1, ChargedStauChildren2)")
            .Define("ChargedStauChildren_recoMatch_ix","Utils::merge(ChargedStauChildren1_recoMatch_ix, ChargedStauChildren2_recoMatch_ix)")
            .Define("ChargedStauChildren_pt","MCParticle::get_pt(ChargedStauChildren)")
            .Define("ChargedStauChildren_p","MCParticle::get_p(ChargedStauChildren)")
            .Define("ChargedStauChildren_eta","MCParticle::get_eta(ChargedStauChildren)")
            .Define("ChargedStauChildren_phi","MCParticle::get_phi(ChargedStauChildren)")
            .Define("ChargedStauChildren_pdg","MCParticle::get_pdg(ChargedStauChildren)")
            .Define("ChargedStauChildren_charge","MCParticle::get_charge(ChargedStauChildren)")
            .Define("ChargedStauChildren_vertex_x","MCParticle::get_vertex_x(ChargedStauChildren)")
            .Define("ChargedStauChildren_vertex_y","MCParticle::get_vertex_y(ChargedStauChildren)")
            .Define("ChargedStauChildren_vertex_z","MCParticle::get_vertex_z(ChargedStauChildren)")

        )
        for prop in ["status",
                        "e",
                        "time",
                        "theta",
                        "beta",
                        "gamma",
                        "phi" ]:
            df2 = df2.Define(f"Stau_{prop}",f"GenStau_{prop}[GenStau_status==22]")


        # Stau lifetime 
        df2 = df2.Define("Stau1_Lxy","sqrt(ChargedStauChildren1_vx[0]*ChargedStauChildren1_vx[0] + ChargedStauChildren1_vy[0]*ChargedStauChildren1_vy[0])")
        df2 = df2.Define("Stau1_Lxyz","sqrt(Stau1_Lxy*Stau1_Lxy + ChargedStauChildren1_vz[0]*ChargedStauChildren1_vz[0])")
        df2 = df2.Define("Stau2_Lxy","sqrt(ChargedStauChildren2_vx[0]*ChargedStauChildren2_vx[0] + ChargedStauChildren2_vy[0]*ChargedStauChildren2_vy[0])")
        df2 = df2.Define("Stau2_Lxyz","sqrt(Stau2_Lxy*Stau2_Lxy + ChargedStauChildren2_vz[0]*ChargedStauChildren2_vz[0])")
        df2 = df2.Define("Stau_Lxy","Utils::merge<float>({Stau1_Lxy},{Stau2_Lxy})")
        df2 = df2.Define("Stau_Lxyz","Utils::merge<float>({Stau1_Lxyz},{Stau2_Lxyz})")

        return df2
    
    def output(self):
        out = super().output() 
        out += ["GenStauIndices"]
        out += ["GenStauIndicesStat22"]
        out += ["ChargedStauChildrenIndices1"]
        out += ["ChargedStauChildrenIndices2"]
        out += ["ChargedStauChildren1_PDG"]
        out += ["ChargedStauChildren2_PDG"]
        out += ["ChargedStauChildren1_q"]
        out += ["ChargedStauChildren2_q"]
        out += ["ChargedStauChildren1_n"]
        out += ["ChargedStauChildren2_n"]
        out += ["ChargedStauChildren1_recoMatch_ix"]
        out += ["ChargedStauChildren2_recoMatch_ix"]
        out += ["ChargedStauChildren1_recoMatch_n"]
        out += ["ChargedStauChildren2_recoMatch_n"]
        out += ["ChargedStauChildren1_vxMatch_ix"]
        out += ["ChargedStauChildren2_vxMatch_ix"]
        out += ["ChargedStauChildrenIndices"]
        out += ["ChargedStauChildren_recoMatch_ix"]
        out += ["ChargedStauChildren_pt"]
        out += ["ChargedStauChildren_p"]
        out += ["ChargedStauChildren_eta"]
        out += ["ChargedStauChildren_phi"]
        out += ["ChargedStauChildren_pdg"]
        out += ["ChargedStauChildren_charge"]
        out += ["ChargedStauChildren_vertex_x"]
        out += ["ChargedStauChildren_vertex_y"]
        out += ["ChargedStauChildren_vertex_z"]
        out += ["Stau1_recoMatch_ix"]
        out += ["Stau2_recoMatch_ix"]
        out += ["Stau1_reco_for_KinkVx_ix"]
        out += ["Stau2_reco_for_KinkVx_ix"]
        out += ["Stau1_kinkVertexMatch_ix"]
        out += ["Stau2_kinkVertexMatch_ix"]
        out += ["Stau_nChildren"]
        out += ["Stau_nRecoTracks"]
        out += ["Stau_DVmatch_index"]
        out += ["Stau_KVmatch_index"]
        out += [f"Stau_{prop}" for prop in ["status",
                        "e",
                        "time",
                        "theta",
                        "beta",
                        "gamma",
                        "phi",
                         "Lxy","Lxyz" ]]
        # out += ["StauPDG"]
        return out 
