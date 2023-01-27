//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 24 09:23:45 2022 by ROOT version 6.12/06
// from TTree stv_tree/STV analysis tree
// found on file: /uboone/data/users/chilgenb/post_process_feb2022/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2_postprocess.root
//////////////////////////////////////////////////////////

#ifndef kinrecoana_h
#define kinrecoana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TVector3.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class kinrecoana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          is_mc;
   Bool_t          mc_neutrino_is_numu;
   Bool_t          mc_vertex_in_FV;
   Bool_t          mc_muon_in_mom_range;
   Bool_t          mc_lead_p_in_mom_range;
   Bool_t          mc_no_fs_pi0;
   Bool_t          mc_no_charged_pi_above_threshold;
   Bool_t          mc_no_fs_mesons;
   Bool_t          mc_is_signal;
   Int_t           category;
   Float_t         spline_weight;
   Float_t         tuned_cv_weight;
   vector<double>  *weight_All_UBGenie;
   vector<double>  *weight_AxFFCCQEshape_UBGenie;
   vector<double>  *weight_DecayAngMEC_UBGenie;
   vector<double>  *weight_NormCCCOH_UBGenie;
   vector<double>  *weight_NormNCCOH_UBGenie;
   vector<double>  *weight_RPA_CCQE_UBGenie;
   vector<double>  *weight_RootinoFix_UBGenie;
   vector<double>  *weight_ThetaDelta2NRad_UBGenie;
   vector<double>  *weight_Theta_Delta2Npi_UBGenie;
   vector<double>  *weight_TunedCentralValue_UBGenie;
   vector<double>  *weight_VecFFCCQEshape_UBGenie;
   vector<double>  *weight_XSecShape_CCMEC_UBGenie;
   vector<double>  *weight_flux_all;
   vector<double>  *weight_reint_all;
   vector<double>  *weight_splines_general_Spline;
   vector<double>  *weight_xsr_scc_Fa3_SCC;
   vector<double>  *weight_xsr_scc_Fv3_SCC;
   Float_t         nu_completeness_from_pfp;
   Float_t         nu_purity_from_pfp;
   Int_t           nslice;
   Bool_t          sel_nu_mu_cc;
   Bool_t          sel_reco_vertex_in_FV;
   Bool_t          sel_topo_cut_passed;
   Bool_t          sel_cosmic_ip_cut_passed;
   Bool_t          sel_pfp_starts_in_PCV;
   Bool_t          sel_no_reco_showers;
   Bool_t          sel_has_muon_candidate;
   Bool_t          sel_muon_passed_mom_cuts;
   Bool_t          sel_has_p_candidate;
   Bool_t          sel_passed_proton_pid_cut;
   Bool_t          sel_protons_contained;
   Bool_t          sel_lead_p_passed_mom_cuts;
   Bool_t          sel_CCNp0pi;
   Int_t           muon_candidate_idx;
   Int_t           lead_p_candidate_idx;
   TVector3        *p3_mu;
   TVector3        *p3_lead_p;
   vector<TVector3> *p3_p_vec;
   TVector3        *mc_p3_mu;
   TVector3        *mc_p3_lead_p;
   vector<TVector3> *mc_p3_p_vec;
   Float_t         delta_pT;
   Float_t         delta_phiT;
   Float_t         delta_alphaT;
   Float_t         delta_pL;
   Float_t         pn;
   Float_t         delta_pTx;
   Float_t         delta_pTy;
   Float_t         mc_delta_pT;
   Float_t         mc_delta_phiT;
   Float_t         mc_delta_alphaT;
   Float_t         mc_delta_pL;
   Float_t         mc_pn;
   Float_t         mc_delta_pTx;
   Float_t         mc_delta_pTy;
   Float_t         topological_score;
   Float_t         CosmicIP;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           mc_nu_pdg;
   Float_t         mc_nu_vtx_x;
   Float_t         mc_nu_vtx_y;
   Float_t         mc_nu_vtx_z;
   Float_t         mc_nu_energy;
   Int_t           mc_ccnc;
   Int_t           mc_interaction;
   vector<unsigned int> *pfp_generation_v;
   vector<unsigned int> *pfp_trk_daughters_v;
   vector<unsigned int> *pfp_shr_daughters_v;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *shr_start_x_v;
   vector<float>   *shr_start_y_v;
   vector<float>   *shr_start_z_v;
   vector<float>   *shr_dist_v;
   vector<float>   *trk_len_v;
   vector<float>   *trk_sce_start_x_v;
   vector<float>   *trk_sce_start_y_v;
   vector<float>   *trk_sce_start_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_sce_end_x_v;
   vector<float>   *trk_sce_end_y_v;
   vector<float>   *trk_sce_end_z_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_range_muon_mom_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_llr_pid_v;
   vector<float>   *trk_llr_pid_u_v;
   vector<float>   *trk_llr_pid_v_v;
   vector<float>   *trk_llr_pid_y_v;
   vector<float>   *trk_llr_pid_score_v;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   Float_t         kin_reco_enu;
   Float_t         kin_reco_enu_diff;
   Float_t         kin_reco_enu_frac;
   Float_t         kin_reco_pmu;
   Float_t         kin_reco_costhetamu;
   Float_t         kin_reco_pmu_diff;
   Float_t         kin_reco_pmu_frac;
   Float_t         kin_reco_costhetamu_diff;
   Int_t           np;
   Int_t           nn;
   Float_t         mc_kin_reco_enu;
   Float_t         mc_kin_reco_enu_diff;
   Float_t         mc_kin_reco_enu_frac;
   Float_t         mc_kin_reco_pmu;
   Float_t         mc_kin_reco_costhetamu;
   Float_t         mc_kin_reco_pmu_diff;
   Float_t         mc_kin_reco_pmu_frac;
   Float_t         mc_kin_reco_costhetamu_diff;
   Int_t           mc_np;
   Int_t           mc_nn;
   TLorentzVector  *mc_p4_had;
   TLorentzVector  *p4_had;

   Bool_t sel_muon_contained;
   Bool_t mc_sublead_p_in_mom_range;
   Bool_t sel_sublead_p_passed_mom_cuts;

   vector<string>  mc_var_names  = {"mc_kin_reco_enu", "mc_kin_reco_pmu", "mc_kin_reco_costhetamu", "mc_nu_energy", "mc_kin_reco_pmu_diff", 
                                    "mc_kin_reco_pmu_frac", "mc_kin_reco_enu_diff", "mc_kin_reco_enu_frac" };
   vector<Float_t> mc_vars_low   = {-3., -3., -1., 0. , -3., -3., -3., -3. };
   vector<Float_t> mc_vars_high  = {3. , 3. , 1. , 3. , 3., 3., 3., 3. };
   vector<Int_t>   mc_vars_nbin  = {101, 101, 101, 101, 101, 101, 101, 101 };
   vector<string>  mc_var_titles = {"E_{#nu}^{kin-truth} [GeV]","p_{#mu}^{kin-truth} [GeV/c]","cos#theta_{#mu}^{kin-truth}", "E_{#nu}^{true} [GeV]", 
                                    "p_{#mu}^{kin-truth} - p_{#mu}^{true} [GeV/c]", "(p_{#mu}^{kin-truth} - p_{#mu}^{true})/p_{#mu}^{true}", 
                                    "E_{#nu}^{kin-truth} - E_{#nu}^{true} [GeV]", "(E_{#nu}^{kin-truth} - E_{#nu}^{true})/E_{#nu}^{true}"};

   vector<string>  reco_var_names  = {"kin_reco_enu", "kin_reco_pmu", "kin_reco_costhetamu","kin_reco_pmu_diff", 
                                      "kin_reco_pmu_frac", "kin_reco_enu_diff", "kin_reco_enu_frac"  };
   vector<Float_t> reco_vars_low   = {-3., -3., -1., -3, -3, -3, -3};
   vector<Float_t> reco_vars_high  = {3., 3., 1., 3., 3., 3., 3.};
   vector<Int_t>   reco_vars_nbin  = {101, 101, 101, 101, 101, 101, 101};
   vector<string>  reco_var_titles = {"E_{#nu}^{kin-reco} [GeV]","p_{#mu}^{kin-reco} [GeV/c]","cos#theta_{#mu}^{kin-reco}",
                                      "p_{#mu}^{kin-reco} - p_{#mu}^{pndr} [GeV/c]", "(p_{#mu}^{kin-reco} - p_{#mu}^{pndr})/p_{#mu}^{pndr}",
                                      "E_{#nu}^{kin-reco} - E_{#nu}^{true} [GeV]", "(E_{#nu}^{kin-reco} - E_{#nu}^{true})/E_{#nu}^{true}"};

   // List of branches
   TBranch        *b_is_mc;   //!
   TBranch        *b_mc_neutrino_is_numu;   //!
   TBranch        *b_mc_vertex_in_FV;   //!
   TBranch        *b_mc_muon_in_mom_range;   //!
   TBranch        *b_mc_lead_p_in_mom_range;   //!
   TBranch        *b_mc_no_fs_pi0;   //!
   TBranch        *b_mc_no_charged_pi_above_threshold;   //!
   TBranch        *b_mc_no_fs_mesons;   //!
   TBranch        *b_mc_is_signal;   //!
   TBranch        *b_category;   //!
   TBranch        *b_spline_weight;   //!
   TBranch        *b_tuned_cv_weight;   //!
   TBranch        *b_weight_All_UBGenie;   //!
   TBranch        *b_weight_AxFFCCQEshape_UBGenie;   //!
   TBranch        *b_weight_DecayAngMEC_UBGenie;   //!
   TBranch        *b_weight_NormCCCOH_UBGenie;   //!
   TBranch        *b_weight_NormNCCOH_UBGenie;   //!
   TBranch        *b_weight_RPA_CCQE_UBGenie;   //!
   TBranch        *b_weight_RootinoFix_UBGenie;   //!
   TBranch        *b_weight_ThetaDelta2NRad_UBGenie;   //!
   TBranch        *b_weight_Theta_Delta2Npi_UBGenie;   //!
   TBranch        *b_weight_TunedCentralValue_UBGenie;   //!
   TBranch        *b_weight_VecFFCCQEshape_UBGenie;   //!
   TBranch        *b_weight_XSecShape_CCMEC_UBGenie;   //!
   TBranch        *b_weight_flux_all;   //!
   TBranch        *b_weight_reint_all;   //!
   TBranch        *b_weight_splines_general_Spline;   //!
   TBranch        *b_weight_xsr_scc_Fa3_SCC;   //!
   TBranch        *b_weight_xsr_scc_Fv3_SCC;   //!
   TBranch        *b_nu_completeness_from_pfp;   //!
   TBranch        *b_nu_purity_from_pfp;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_sel_nu_mu_cc;   //!
   TBranch        *b_sel_reco_vertex_in_FV;   //!
   TBranch        *b_sel_topo_cut_passed;   //!
   TBranch        *b_sel_cosmic_ip_cut_passed;   //!
   TBranch        *b_sel_pfp_starts_in_PCV;   //!
   TBranch        *b_sel_no_reco_showers;   //!
   TBranch        *b_sel_has_muon_candidate;   //!
   TBranch        *b_sel_muon_passed_mom_cuts;   //!
   TBranch        *b_sel_has_p_candidate;   //!
   TBranch        *b_sel_passed_proton_pid_cut;   //!
   TBranch        *b_sel_protons_contained;   //!
   TBranch        *b_sel_lead_p_passed_mom_cuts;   //!
   TBranch        *b_sel_CCNp0pi;   //!
   TBranch        *b_muon_candidate_idx;   //!
   TBranch        *b_lead_p_candidate_idx;   //!
   TBranch        *b_p3_mu;   //!
   TBranch        *b_p3_lead_p;   //!
   TBranch        *b_p3_p_vec;   //!
   TBranch        *b_mc_p3_mu;   //!
   TBranch        *b_mc_p3_lead_p;   //!
   TBranch        *b_mc_p3_p_vec;   //!
   TBranch        *b_delta_pT;   //!
   TBranch        *b_delta_phiT;   //!
   TBranch        *b_delta_alphaT;   //!
   TBranch        *b_delta_pL;   //!
   TBranch        *b_pn;   //!
   TBranch        *b_delta_pTx;   //!
   TBranch        *b_delta_pTy;   //!
   TBranch        *b_mc_delta_pT;   //!
   TBranch        *b_mc_delta_phiT;   //!
   TBranch        *b_mc_delta_alphaT;   //!
   TBranch        *b_mc_delta_pL;   //!
   TBranch        *b_mc_pn;   //!
   TBranch        *b_mc_delta_pTx;   //!
   TBranch        *b_mc_delta_pTy;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_CosmicIP;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_mc_nu_pdg;   //!
   TBranch        *b_mc_nu_vtx_x;   //!
   TBranch        *b_mc_nu_vtx_y;   //!
   TBranch        *b_mc_nu_vtx_z;   //!
   TBranch        *b_mc_nu_energy;   //!
   TBranch        *b_mc_ccnc;   //!
   TBranch        *b_mc_interaction;   //!
   TBranch        *b_pfp_generation_v;   //!
   TBranch        *b_pfp_trk_daughters_v;   //!
   TBranch        *b_pfp_shr_daughters_v;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_shr_start_x_v;   //!
   TBranch        *b_shr_start_y_v;   //!
   TBranch        *b_shr_start_z_v;   //!
   TBranch        *b_shr_dist_v;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_sce_start_x_v;   //!
   TBranch        *b_trk_sce_start_y_v;   //!
   TBranch        *b_trk_sce_start_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_sce_end_x_v;   //!
   TBranch        *b_trk_sce_end_y_v;   //!
   TBranch        *b_trk_sce_end_z_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_range_muon_mom_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_llr_pid_v;   //!
   TBranch        *b_trk_llr_pid_u_v;   //!
   TBranch        *b_trk_llr_pid_v_v;   //!
   TBranch        *b_trk_llr_pid_y_v;   //!
   TBranch        *b_trk_llr_pid_score_v;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_kin_reco_enu;   //!
   TBranch        *b_kin_reco_enu_diff;   //!
   TBranch        *b_kin_reco_enu_frac;   //!
   TBranch        *b_kin_reco_pmu;   //!
   TBranch        *b_kin_reco_costhetamu;   //!
   TBranch        *b_kin_reco_pmu_diff;   //!
   TBranch        *b_kin_reco_pmu_frac;   //!
   TBranch        *b_kin_reco_costhetamu_diff;   //!
   TBranch        *b_np;
   TBranch        *b_nn;
   TBranch        *b_mc_kin_reco_enu;   //!
   TBranch        *b_mc_kin_reco_enu_diff;   //!
   TBranch        *b_mc_kin_reco_enu_frac;   //!
   TBranch        *b_mc_kin_reco_pmu;   //!
   TBranch        *b_mc_kin_reco_costhetamu;   //!
   TBranch        *b_mc_kin_reco_pmu_diff;   //!
   TBranch        *b_mc_kin_reco_pmu_frac;   //!
   TBranch        *b_mc_kin_reco_costhetamu_diff;   //!
   TBranch        *b_mc_np;
   TBranch        *b_mc_nn;
   TBranch        *b_mc_p4_had;
   TBranch        *b_p4_had;

   TBranch *b_sel_muon_contained;
   TBranch *b_mc_sublead_p_in_mom_range;
   TBranch *b_sel_sublead_p_passed_mom_cuts;

   kinrecoana(TTree *tree=0);
   kinrecoana(string infile="");
   virtual ~kinrecoana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void SetList(bool useMC=false, string opts="");
   Int_t GetHadMult(const Int_t& pdg) const;
   void EnuEmu();
   void cc1p();
   void ccnp();
   void DeactivateBranches();
   void BwrdPs();
   void Binning();
   void Resolution();
   void ProtonMult();
   void Plot2p();
   void SignalDef();
};

#endif

#ifdef kinrecoana_cxx
kinrecoana::kinrecoana(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uboone/data/users/chilgenb/post_process_feb2022/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_postprocess_add_phad.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uboone/data/users/chilgenb/post_process_feb2022/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_postprocess_add_phad.root");
      }
      f->GetObject("stv_tree",tree);

   }
   Init(tree);
}

kinrecoana::kinrecoana(string infile) : fChain(0)
{
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(infile.c_str());
      }
      TTree *tree;
      f->GetObject("stv_tree",tree);

   Init(tree);
}



kinrecoana::~kinrecoana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kinrecoana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);

   /*mc_vars.clear();
   mc_vars.push_back(mc_kin_reco_enu);
   mc_vars.push_back(mc_kin_reco_pmu);
   mc_vars.push_back(mc_kin_reco_costhetamu);
   mc_vars.push_back(mc_nu_energy);

   reco_vars.clear();
   reco_vars.push_back(kin_reco_enu);
   reco_vars.push_back(kin_reco_pmu);
   reco_vars.push_back(kin_reco_costhetamu);*/
}
Long64_t kinrecoana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void kinrecoana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   weight_All_UBGenie = 0;
   weight_AxFFCCQEshape_UBGenie = 0;
   weight_DecayAngMEC_UBGenie = 0;
   weight_NormCCCOH_UBGenie = 0;
   weight_NormNCCOH_UBGenie = 0;
   weight_RPA_CCQE_UBGenie = 0;
   weight_RootinoFix_UBGenie = 0;
   weight_ThetaDelta2NRad_UBGenie = 0;
   weight_Theta_Delta2Npi_UBGenie = 0;
   weight_TunedCentralValue_UBGenie = 0;
   weight_VecFFCCQEshape_UBGenie = 0;
   weight_XSecShape_CCMEC_UBGenie = 0;
   weight_flux_all = 0;
   weight_reint_all = 0;
   weight_splines_general_Spline = 0;
   weight_xsr_scc_Fa3_SCC = 0;
   weight_xsr_scc_Fv3_SCC = 0;
   p3_mu = 0;
   p3_lead_p = 0;
   p3_p_vec = 0;
   mc_p3_mu = 0;
   mc_p3_lead_p = 0;
   mc_p3_p_vec = 0;
   pfp_generation_v = 0;
   pfp_trk_daughters_v = 0;
   pfp_shr_daughters_v = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   shr_start_x_v = 0;
   shr_start_y_v = 0;
   shr_start_z_v = 0;
   shr_dist_v = 0;
   trk_len_v = 0;
   trk_sce_start_x_v = 0;
   trk_sce_start_y_v = 0;
   trk_sce_start_z_v = 0;
   trk_distance_v = 0;
   trk_sce_end_x_v = 0;
   trk_sce_end_y_v = 0;
   trk_sce_end_z_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_energy_proton_v = 0;
   trk_range_muon_mom_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_pid_chipr_v = 0;
   trk_llr_pid_v = 0;
   trk_llr_pid_u_v = 0;
   trk_llr_pid_v_v = 0;
   trk_llr_pid_y_v = 0;
   trk_llr_pid_score_v = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;

   mc_p4_had=0;
   p4_had=0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("is_mc", &is_mc, &b_is_mc);
   fChain->SetBranchAddress("mc_neutrino_is_numu", &mc_neutrino_is_numu, &b_mc_neutrino_is_numu);
   fChain->SetBranchAddress("mc_vertex_in_FV", &mc_vertex_in_FV, &b_mc_vertex_in_FV);
   fChain->SetBranchAddress("mc_muon_in_mom_range", &mc_muon_in_mom_range, &b_mc_muon_in_mom_range);
   fChain->SetBranchAddress("mc_lead_p_in_mom_range", &mc_lead_p_in_mom_range, &b_mc_lead_p_in_mom_range);
   fChain->SetBranchAddress("mc_no_fs_pi0", &mc_no_fs_pi0, &b_mc_no_fs_pi0);
   fChain->SetBranchAddress("mc_no_charged_pi_above_threshold", &mc_no_charged_pi_above_threshold, &b_mc_no_charged_pi_above_threshold);
   fChain->SetBranchAddress("mc_no_fs_mesons", &mc_no_fs_mesons, &b_mc_no_fs_mesons);
   fChain->SetBranchAddress("mc_is_signal", &mc_is_signal, &b_mc_is_signal);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("spline_weight", &spline_weight, &b_spline_weight);
   fChain->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight, &b_tuned_cv_weight);
   fChain->SetBranchAddress("weight_All_UBGenie", &weight_All_UBGenie, &b_weight_All_UBGenie);
   fChain->SetBranchAddress("weight_AxFFCCQEshape_UBGenie", &weight_AxFFCCQEshape_UBGenie, &b_weight_AxFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("weight_DecayAngMEC_UBGenie", &weight_DecayAngMEC_UBGenie, &b_weight_DecayAngMEC_UBGenie);
   fChain->SetBranchAddress("weight_NormCCCOH_UBGenie", &weight_NormCCCOH_UBGenie, &b_weight_NormCCCOH_UBGenie);
   fChain->SetBranchAddress("weight_NormNCCOH_UBGenie", &weight_NormNCCOH_UBGenie, &b_weight_NormNCCOH_UBGenie);
   fChain->SetBranchAddress("weight_RPA_CCQE_UBGenie", &weight_RPA_CCQE_UBGenie, &b_weight_RPA_CCQE_UBGenie);
   fChain->SetBranchAddress("weight_RootinoFix_UBGenie", &weight_RootinoFix_UBGenie, &b_weight_RootinoFix_UBGenie);
   fChain->SetBranchAddress("weight_ThetaDelta2NRad_UBGenie", &weight_ThetaDelta2NRad_UBGenie, &b_weight_ThetaDelta2NRad_UBGenie);
   fChain->SetBranchAddress("weight_Theta_Delta2Npi_UBGenie", &weight_Theta_Delta2Npi_UBGenie, &b_weight_Theta_Delta2Npi_UBGenie);
   fChain->SetBranchAddress("weight_TunedCentralValue_UBGenie", &weight_TunedCentralValue_UBGenie, &b_weight_TunedCentralValue_UBGenie);
   fChain->SetBranchAddress("weight_VecFFCCQEshape_UBGenie", &weight_VecFFCCQEshape_UBGenie, &b_weight_VecFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("weight_XSecShape_CCMEC_UBGenie", &weight_XSecShape_CCMEC_UBGenie, &b_weight_XSecShape_CCMEC_UBGenie);
   fChain->SetBranchAddress("weight_flux_all", &weight_flux_all, &b_weight_flux_all);
   fChain->SetBranchAddress("weight_reint_all", &weight_reint_all, &b_weight_reint_all);
   fChain->SetBranchAddress("weight_splines_general_Spline", &weight_splines_general_Spline, &b_weight_splines_general_Spline);
   fChain->SetBranchAddress("weight_xsr_scc_Fa3_SCC", &weight_xsr_scc_Fa3_SCC, &b_weight_xsr_scc_Fa3_SCC);
   fChain->SetBranchAddress("weight_xsr_scc_Fv3_SCC", &weight_xsr_scc_Fv3_SCC, &b_weight_xsr_scc_Fv3_SCC);
   fChain->SetBranchAddress("nu_completeness_from_pfp", &nu_completeness_from_pfp, &b_nu_completeness_from_pfp);
   fChain->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp, &b_nu_purity_from_pfp);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("sel_nu_mu_cc", &sel_nu_mu_cc, &b_sel_nu_mu_cc);
   fChain->SetBranchAddress("sel_reco_vertex_in_FV", &sel_reco_vertex_in_FV, &b_sel_reco_vertex_in_FV);
   fChain->SetBranchAddress("sel_topo_cut_passed", &sel_topo_cut_passed, &b_sel_topo_cut_passed);
   fChain->SetBranchAddress("sel_cosmic_ip_cut_passed", &sel_cosmic_ip_cut_passed, &b_sel_cosmic_ip_cut_passed);
   fChain->SetBranchAddress("sel_pfp_starts_in_PCV", &sel_pfp_starts_in_PCV, &b_sel_pfp_starts_in_PCV);
   fChain->SetBranchAddress("sel_no_reco_showers", &sel_no_reco_showers, &b_sel_no_reco_showers);
   fChain->SetBranchAddress("sel_has_muon_candidate", &sel_has_muon_candidate, &b_sel_has_muon_candidate);
   fChain->SetBranchAddress("sel_muon_passed_mom_cuts", &sel_muon_passed_mom_cuts, &b_sel_muon_passed_mom_cuts);
   fChain->SetBranchAddress("sel_has_p_candidate", &sel_has_p_candidate, &b_sel_has_p_candidate);
   fChain->SetBranchAddress("sel_passed_proton_pid_cut", &sel_passed_proton_pid_cut, &b_sel_passed_proton_pid_cut);
   fChain->SetBranchAddress("sel_protons_contained", &sel_protons_contained, &b_sel_protons_contained);
   fChain->SetBranchAddress("sel_lead_p_passed_mom_cuts", &sel_lead_p_passed_mom_cuts, &b_sel_lead_p_passed_mom_cuts);
   fChain->SetBranchAddress("sel_CCNp0pi", &sel_CCNp0pi, &b_sel_CCNp0pi);
   fChain->SetBranchAddress("muon_candidate_idx", &muon_candidate_idx, &b_muon_candidate_idx);
   fChain->SetBranchAddress("lead_p_candidate_idx", &lead_p_candidate_idx, &b_lead_p_candidate_idx);
   fChain->SetBranchAddress("p3_mu", &p3_mu, &b_p3_mu);
   fChain->SetBranchAddress("p3_lead_p", &p3_lead_p, &b_p3_lead_p);
   fChain->SetBranchAddress("p3_p_vec", &p3_p_vec, &b_p3_p_vec);
   fChain->SetBranchAddress("mc_p3_mu", &mc_p3_mu, &b_mc_p3_mu);
   fChain->SetBranchAddress("mc_p3_lead_p", &mc_p3_lead_p, &b_mc_p3_lead_p);
   fChain->SetBranchAddress("mc_p3_p_vec", &mc_p3_p_vec, &b_mc_p3_p_vec);
   fChain->SetBranchAddress("delta_pT", &delta_pT, &b_delta_pT);
   fChain->SetBranchAddress("delta_phiT", &delta_phiT, &b_delta_phiT);
   fChain->SetBranchAddress("delta_alphaT", &delta_alphaT, &b_delta_alphaT);
   fChain->SetBranchAddress("delta_pL", &delta_pL, &b_delta_pL);
   fChain->SetBranchAddress("pn", &pn, &b_pn);
   fChain->SetBranchAddress("delta_pTx", &delta_pTx, &b_delta_pTx);
   fChain->SetBranchAddress("delta_pTy", &delta_pTy, &b_delta_pTy);
   fChain->SetBranchAddress("mc_delta_pT", &mc_delta_pT, &b_mc_delta_pT);
   fChain->SetBranchAddress("mc_delta_phiT", &mc_delta_phiT, &b_mc_delta_phiT);
   fChain->SetBranchAddress("mc_delta_alphaT", &mc_delta_alphaT, &b_mc_delta_alphaT);
   fChain->SetBranchAddress("mc_delta_pL", &mc_delta_pL, &b_mc_delta_pL);
   fChain->SetBranchAddress("mc_pn", &mc_pn, &b_mc_pn);
   fChain->SetBranchAddress("mc_delta_pTx", &mc_delta_pTx, &b_mc_delta_pTx);
   fChain->SetBranchAddress("mc_delta_pTy", &mc_delta_pTy, &b_mc_delta_pTy);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("CosmicIP", &CosmicIP, &b_CosmicIP);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("mc_nu_pdg", &mc_nu_pdg, &b_mc_nu_pdg);
   fChain->SetBranchAddress("mc_nu_vtx_x", &mc_nu_vtx_x, &b_mc_nu_vtx_x);
   fChain->SetBranchAddress("mc_nu_vtx_y", &mc_nu_vtx_y, &b_mc_nu_vtx_y);
   fChain->SetBranchAddress("mc_nu_vtx_z", &mc_nu_vtx_z, &b_mc_nu_vtx_z);
   fChain->SetBranchAddress("mc_nu_energy", &mc_nu_energy, &b_mc_nu_energy);
   fChain->SetBranchAddress("mc_ccnc", &mc_ccnc, &b_mc_ccnc);
   fChain->SetBranchAddress("mc_interaction", &mc_interaction, &b_mc_interaction);
   fChain->SetBranchAddress("pfp_generation_v", &pfp_generation_v, &b_pfp_generation_v);
   fChain->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v, &b_pfp_trk_daughters_v);
   fChain->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v, &b_pfp_shr_daughters_v);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("shr_start_x_v", &shr_start_x_v, &b_shr_start_x_v);
   fChain->SetBranchAddress("shr_start_y_v", &shr_start_y_v, &b_shr_start_y_v);
   fChain->SetBranchAddress("shr_start_z_v", &shr_start_z_v, &b_shr_start_z_v);
   fChain->SetBranchAddress("shr_dist_v", &shr_dist_v, &b_shr_dist_v);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v, &b_trk_sce_start_x_v);
   fChain->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v, &b_trk_sce_start_y_v);
   fChain->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v, &b_trk_sce_start_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v, &b_trk_sce_end_x_v);
   fChain->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v, &b_trk_sce_end_y_v);
   fChain->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v, &b_trk_sce_end_z_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v, &b_trk_range_muon_mom_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_llr_pid_v", &trk_llr_pid_v, &b_trk_llr_pid_v);
   fChain->SetBranchAddress("trk_llr_pid_u_v", &trk_llr_pid_u_v, &b_trk_llr_pid_u_v);
   fChain->SetBranchAddress("trk_llr_pid_v_v", &trk_llr_pid_v_v, &b_trk_llr_pid_v_v);
   fChain->SetBranchAddress("trk_llr_pid_y_v", &trk_llr_pid_y_v, &b_trk_llr_pid_y_v);
   fChain->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v, &b_trk_llr_pid_score_v);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("kin_reco_enu", &kin_reco_enu, &b_kin_reco_enu);
   fChain->SetBranchAddress("kin_reco_enu_diff", &kin_reco_enu_diff, &b_kin_reco_enu_diff);
   fChain->SetBranchAddress("kin_reco_enu_frac", &kin_reco_enu_frac, &b_kin_reco_enu_frac);
   fChain->SetBranchAddress("kin_reco_pmu", &kin_reco_pmu, &b_kin_reco_pmu);
   fChain->SetBranchAddress("kin_reco_costhetamu", &kin_reco_costhetamu, &b_kin_reco_costhetamu);
   fChain->SetBranchAddress("kin_reco_pmu_diff", &kin_reco_pmu_diff, &b_kin_reco_pmu_diff);
   fChain->SetBranchAddress("kin_reco_pmu_frac", &kin_reco_pmu_frac, &b_kin_reco_pmu_frac);
   fChain->SetBranchAddress("kin_reco_costhetamu_diff", &kin_reco_costhetamu_diff, &b_kin_reco_costhetamu_diff);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("nn", &nn, &b_nn);
   fChain->SetBranchAddress("mc_kin_reco_enu", &mc_kin_reco_enu, &b_mc_kin_reco_enu);
   fChain->SetBranchAddress("mc_kin_reco_enu_diff", &mc_kin_reco_enu_diff, &b_mc_kin_reco_enu_diff);
   fChain->SetBranchAddress("mc_kin_reco_enu_frac", &mc_kin_reco_enu_frac, &b_mc_kin_reco_enu_frac);
   fChain->SetBranchAddress("mc_kin_reco_pmu", &mc_kin_reco_pmu, &b_mc_kin_reco_pmu);
   fChain->SetBranchAddress("mc_kin_reco_costhetamu", &mc_kin_reco_costhetamu, &b_mc_kin_reco_costhetamu);
   fChain->SetBranchAddress("mc_kin_reco_pmu_diff", &mc_kin_reco_pmu_diff, &b_mc_kin_reco_pmu_diff);
   fChain->SetBranchAddress("mc_kin_reco_pmu_frac", &mc_kin_reco_pmu_frac, &b_mc_kin_reco_pmu_frac);
   fChain->SetBranchAddress("mc_kin_reco_costhetamu_diff", &mc_kin_reco_costhetamu_diff, &b_mc_kin_reco_costhetamu_diff);
   fChain->SetBranchAddress("mc_np", &mc_np, &b_mc_np);
   fChain->SetBranchAddress("mc_nn", &mc_nn, &b_mc_nn);
   fChain->SetBranchAddress("mc_p4_had", &mc_p4_had, &b_mc_p4_had);
   fChain->SetBranchAddress("p4_had", &p4_had, &b_p4_had);

   fChain->SetBranchAddress("sel_muon_contained", &sel_muon_contained, &b_sel_muon_contained);
   fChain->SetBranchAddress("mc_sublead_p_in_mom_range", &mc_sublead_p_in_mom_range, &b_mc_sublead_p_in_mom_range);
   fChain->SetBranchAddress("sel_sublead_p_passed_mom_cuts", &sel_sublead_p_passed_mom_cuts, &b_sel_sublead_p_passed_mom_cuts);

   Notify();
}

Bool_t kinrecoana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kinrecoana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kinrecoana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef kinrecoana_cxx
