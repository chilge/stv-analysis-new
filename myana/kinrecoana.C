#define kinrecoana_cxx
#include "kinrecoana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "stv-analysis-new/KinematicReco.hh"

void Rotate(const Float_t& x, const Float_t&y, Float_t& xprime, Float_t&yprime, const Float_t& slope){
  float theta = atan(1./slope);
  xprime = x*cos(theta) - y*sin(theta);
  yprime = y*sin(theta) + x*cos(theta);
} 

void kinrecoana::SetList(bool useMC, string opts) {

  this->DeactivateBranches();
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  if (elist) delete elist;

  std::string sel = "sel_CCNp0pi";
  if(useMC) sel = "mc_is_signal"; // FIXME: should probably use both if doing reco to avoid introducing backgrounds (for now)
  if(opts!="") sel+="&&"+opts;
  std::cout << "using selection: " << sel << std::endl;
  fChain->Draw(">>elist",sel.c_str(),"entrylist");
  elist = (TEntryList*)gDirectory->Get("elist");

  cout << "total tree entries: " << fChain->GetEntries() << endl;
  cout << "list entries: " << elist->GetN() << endl;
  if(elist->GetN()==0) {
    std::cout << "ERROR in SetList() - entry list is empty!" << std::endl;
    return;
  }
  fChain->SetEntryList(elist);
}

// deactivate all branches and enable only the ones we need (massive speed boost)
void kinrecoana::DeactivateBranches(){

  fChain->SetBranchStatus("*","0"); //deactivate all branches

  // pandora
  fChain->SetBranchStatus("trk_range_muon_mom_v",1);
  fChain->SetBranchStatus("trk_mcs_muon_mom_v",1);
  fChain->SetBranchStatus("muon_candidate_idx",1);
  fChain->SetBranchStatus("sel_muon_contained",1);

  // truth
  fChain->SetBranchStatus("mc_is_signal",1);
  fChain->SetBranchStatus("mc_np",1);
  fChain->SetBranchStatus("mc_nn",1);
  fChain->SetBranchStatus("mc_p3_mu",1);
  fChain->SetBranchStatus("mc_p4_had",1);
  fChain->SetBranchStatus("mc_p3_p_vec",1);
  fChain->SetBranchStatus("mc_kin_reco_enu",1);
  fChain->SetBranchStatus("mc_kin_reco_enu_diff",1);
  fChain->SetBranchStatus("mc_kin_reco_enu_frac",1);
  fChain->SetBranchStatus("mc_kin_reco_emu",1);
  fChain->SetBranchStatus("mc_kin_reco_pmu",1);
  fChain->SetBranchStatus("mc_kin_reco_pmu_diff",1);
  fChain->SetBranchStatus("mc_kin_reco_pmu_frac",1);
  fChain->SetBranchStatus("mc_kin_reco_costhetamu",1);
  fChain->SetBranchStatus("mc_kin_reco_costhetamu_diff",1);

  // reco
  fChain->SetBranchStatus("sel_CCNp0pi",1);
  fChain->SetBranchStatus("np",1);
  fChain->SetBranchStatus("nn",1);
  fChain->SetBranchStatus("p3_mu",1);
  fChain->SetBranchStatus("p4_had",1);
  fChain->SetBranchStatus("kin_reco_enu",1);
  fChain->SetBranchStatus("kin_reco_enu_diff",1);
  fChain->SetBranchStatus("kin_reco_enu_frac",1);
  fChain->SetBranchStatus("kin_reco_emu",1);
  fChain->SetBranchStatus("kin_reco_pmu",1);
  fChain->SetBranchStatus("kin_reco_pmu_diff",1);
  fChain->SetBranchStatus("kin_reco_pmu_frac",1);
  fChain->SetBranchStatus("kin_reco_costhetamu",1);
  fChain->SetBranchStatus("kin_reco_costhetamu_diff",1);
}

void kinrecoana::Test() {
  this->SetList();
  TH1F* h = new TH1F("h","",100,-1,1);

  for(size_t ientry=0; ientry<1000; ientry++) {
  //for(size_t ientry=0; ientry<fChain->GetEntries(); ientry++) {

    Int_t entryNumber = fChain->GetEntryNumber(ientry);
    if(entryNumber < 0) break;
    fChain->GetEntry(entryNumber);

    TLorentzVector hadvec = GetHadVec(*p3_p_vec);
    //double enup = KinRecoEnu(hadvec);
    double enu = KinRecoEnu(kin_reco_emu,hadvec.E());
     
    h->Fill(enu-kin_reco_enu);
    //h->Fill(enu-enup);
  }

  h->Draw();
}

Int_t kinrecoana::GetHadMult(const int& pdg) const {

  Int_t n=0;
  for(int const& ipdg : *mc_pdg){
    if(pdg==ipdg)
      n++;
  }

  return n;
}

void kinrecoana::EnuEmu() {

  bool usemc = false;
  bool viewneg = true;
  const size_t nmults=5;
  const float maxe = 3;
  this->SetList(usemc);

  string htitle = "MC ";
  if(usemc) htitle += "truth - CCNp0#pi;N_{p}";
  else htitle += "reco - CCNp0#pi;N_{p}";
  TH1F* hmult = new TH1F("hmult",htitle.c_str(),24,-0.5,23.5);
  if(usemc) fChain->Draw("mc_np>>hmult");
  else fChain->Draw("np>>hmult");

  vector<TGraph*> gs;
  for(int i=1; i<nmults+1; i++){

    string cut;
    string var1 = "kin_reco_enu", var2 = "kin_reco_emu";
    if(usemc){
      var1 = "mc_" + var1;
      var2 = "mc_" + var2;
    }
    string draw = var1 + ":" + var2;

    if(i<nmults) cut="abs("+var1+")<" + to_string(maxe) + " && abs(" + var2 + ")<" + to_string(maxe) + " && mc_np=="+to_string(i);
    else cut="abs("+var1+")<" + to_string(maxe) + " && abs("+var2+")<" + to_string(maxe) +" && mc_np>="+to_string(nmults);

    fChain->Draw(draw.c_str(),cut.c_str());
    gs.push_back((TGraph*)gPad->FindObject("Graph")->Clone());
  }

  //NC1p0pi only
  gs[0]->SetMarkerColor(kBlue);
  gs[0]->SetLineWidth(3);
  gs[0]->SetTitle("Kinematic reco with CC1p0#pi");
  gs[0]->GetXaxis()->SetTitle("E_{#mu} [GeV]");
  gs[0]->GetYaxis()->SetTitle("E_{#nu} [GeV]");

  gs[0]->GetHistogram()->SetMaximum(maxe);
  //conflicting ranges when !viewneg && usemc?
  if(!viewneg){
    gs[0]->GetXaxis()->SetLimits(0,maxe);
    gs[0]->GetHistogram()->SetMinimum(0);
  }

  // use with proton mass correction
  if(usemc){
    gs[0]->GetXaxis()->SetLimits(-maxe,maxe);
    gs[0]->GetHistogram()->SetMinimum(-maxe);
  }

  TCanvas *c1p = new TCanvas("c1p","");
  gs[0]->Draw("ap");
  if(usemc)
    c1p->SaveAs("cc1p0pi_mc_kin_enu_emu.gif"); 
  else
    c1p->SaveAs("cc1p0pi_kin_enu_emu.gif"); 

  cout << "loop over multiplicities..." << endl;
  for(size_t i=1; i<nmults; i++){
    
    //string title = "N_{n} ";
    string title = "N_{p} ";
    //if(i<2) title+=" = "+to_string(i); //neutron
    if(i<nmults-1) title+=" = "+to_string(i+1); //proton
    else title += " #geq"+to_string(nmults);
    gs[i]->SetTitle(title.c_str());
    gs[i]->GetXaxis()->SetTitle("E_{#mu} [GeV]");
    gs[i]->GetYaxis()->SetTitle("E_{#nu} [GeV]");
    gs[i]->SetLineWidth(3);

    //conflicting ranges when !viewneg && usemc?
    gs[i]->GetHistogram()->SetMaximum(maxe);
    if(!viewneg){
      gs[i]->GetXaxis()->SetLimits(0,maxe);
      gs[i]->GetHistogram()->SetMinimum(0);
    }
    else{
      gs[i]->GetXaxis()->SetLimits(-maxe,maxe);
      gs[i]->GetHistogram()->SetMinimum(-maxe);
    }


    // use with nucleon mass correction
    /*if(usemc){
      gs[i]->GetXaxis()->SetLimits(0.1,1.5);
      gs[i]->GetHistogram()->SetMaximum(1);
      gs[i]->GetHistogram()->SetMinimum(-0.7);
    }*/
  }

  // styling/drawing
  gs[0]->SetTitle("Kinematic reco with CCNp0#pi");
  if(nmults>=2) gs[1]->SetMarkerColor(kCyan);
  if(nmults>=3) gs[2]->SetMarkerColor(kGreen-2);
  if(nmults>=4) gs[3]->SetMarkerColor(kOrange);
  if(nmults>=5) gs[4]->SetMarkerColor(kMagenta);
  if(nmults>=6) gs[5]->SetMarkerColor(kRed);
  if(nmults>=7) gs[6]->SetMarkerColor(46);
  for(size_t i=0; i<gs.size(); i++)
    gs[i]->SetLineColor(gs[i]->GetMarkerColor());

  // drawing //
  string outname;
  TCanvas *cmult = new TCanvas("cmult","");
  hmult->SetLineWidth(3);
  hmult->Draw("e0hist");
  cmult->SetLogy();
  outname = "ccnp0pi_proton_mult.gif";
  if(usemc)
    outname = "mc_" + outname;
  cmult->SaveAs(outname.c_str());

  TCanvas *call = new TCanvas("call","");
  gs[0]->Draw("ap");
  //for(int i=nmults-1; i>-1; i--)
  for(size_t i=0; i<nmults-1; i++)
    gs[i]->Draw("samep");

  TLegend *lall = new TLegend(0.7,0.2,0.85,0.5);
  lall->SetBorderSize(0);
  for(size_t i=0; i<nmults-1; i++) {
    string txt = to_string(i+1) + "p";
    lall->AddEntry(gs[i],txt.c_str(), "l");
  }
  string geqlabel="#geq"+to_string(nmults)+"p";
  lall->AddEntry(gs[nmults-1],geqlabel.c_str(), "l");
  lall->Draw();

  outname = "ccnp0pi_kin_enu_emu_by_prot_mult_overlaid.gif";
  if(usemc)
    outname = "mc_" + outname;
  call->SaveAs(outname.c_str());

  TCanvas *cbymult = new TCanvas("cbymult","by proton multiplicity",1400,800);
  if(nmults==7) cbymult->Divide(3,2);
  if(nmults==5) cbymult->Divide(2,2);
  cbymult->Modified();
  cbymult->Update();
  for(size_t i=1; i<nmults; i++){
    cbymult->cd(i);
    gs.at(i)->Draw("ap");
  }
  cbymult->Modified();

  outname = "ccnp0pi_kin_enu_emu_by_prot_mult.gif";
  if(usemc)
    outname = "mc_" + outname;
  cbymult->SaveAs(outname.c_str());

}// end EnuEmu()

void kinrecoana::ccnp() {

  bool usemc=true;
  this->SetList(true); //pass true for mc selection

  // 1D kin reco vars
  TH1F *hkin_emu_neg = new TH1F("hkin_emu_neg","Kin. Reco.; E_{#mu} [GeV]",  100,-3,0);
  TH1F *hkin_emu     = new TH1F("hkin_emu",    "Kin. Reco.; E_{#mu} [GeV]",  100,0,3);
  TH1F *hkin_enu_neg = new TH1F("hkin_enu_neg","Kin. Reco.; E_{#nu} [GeV]",  100,-3,0);
  TH1F *hkin_enu     = new TH1F("hkin_enu",    "Kin. Reco.; E_{#nu} [GeV]",  100,0,3);
  TH1F *hkin_cosmu   = new TH1F("hkin_cosmu",  "Kin. Reco.; cos#theta_{#mu}",101,-1,1);
  TH1F *hkin_cosmuhad = new TH1F("hkin_cosmuhad","Kin. Reco.; cos#theta_{#mu-h}",101,-1,1);

  // 2D kin reco vars
  TH2F* h2_emu_enu   = new TH2F("h2_emu_enu",  "Kin. Reco.;E_{#mu} [GeV]; E_{#nu} [GeV]",101,-3,3,101,-3,3); 
  TH2F* h2_emu_cosmu = new TH2F("h2_emu_cosmu","Kin. Reco.;E_{#mu} [GeV]; cos#theta_{#mu}",51,-3,3,51,-1,1); 
  TH2F* h2_cosmu_enu = new TH2F("h2_cosmu_enu","Kin. Reco.;cos#theta_{#mu}; E_{#nu} [GeV]",51,-1,1,51,-3,3); 
  TH2F* h2_cosmu_nu  = new TH2F("h2_cosmu_nu", "Kin. Reco.;cos#theta_{#mu}; #nu [GeV]",51,-1,1,51,-0.1,3); 

  // 1D resolution
  TH1F *hpmu_diff    = new TH1F("hpmu_diff",  ";p_{#mu}^{kin} - p_{#mu}^{pndr} [GeV]",           101,-2,2);
  TH1F *hpmu_frac    = new TH1F("hpmu_frac",  ";(p_{#mu}^{kin} - p_{#mu}^{pndr})/p_{#mu}^{pndr}",101,-1,2);
  TH1F *hcosmu_diff  = new TH1F("hcosmu_diff",";cos#theta_{#mu}^{kin} - cos#theta_{#mu}^{pndr}", 101,-2,2);
  TH1F *henu_diff    = new TH1F("henu_diff",  ";E_{#nu}^{kin} - E_{#nu}^{pndr} [GeV]",           101,-2,2);
  TH1F *henu_frac    = new TH1F("henu_frac",  ";(E_{#nu}^{kin} - E_{#nu}^{pndr})/E_{#nu}^{pndr}",101,-1,1);

  // 2D resolution
  TH2F *h2_kin_pndr_pmu   = new TH2F("h2_kin_pndr_pmu",  ";p_{#mu}^{pndr} [GeV];p_{#mu}^{kin} [GeV]",    5,0,3,51,-2,2); 
  TH2F *h2_kin_pndr_cosmu = new TH2F("h2_kin_pndr_cosmu",";cos#theta_{#mu}^{pndr};cos#theta_{#mu}^{kin}",51,-1,1,51,-1,1);
  TH2F *h2_kin_pndr_enu   = new TH2F("h2_kin_pndr_enu",  ";E_{#nu}^{pndr} [GeV];E_{#nu}^{kin} [GeV]",    50,0,3,51,-2,2);

  // fill histos
  TCanvas* c = new TCanvas("c","");
  string outname = "ccNp0pi_plots.pdf";
  if(usemc) outname = "mctruth_" +outname;
  string outname_tmp = outname + "[";
  c->Print(outname_tmp.c_str());

  gStyle->SetOptStat(0);

  // 1D kin reco vars
  string var = "kin_reco_emu";
  if(usemc) var = "mc_" + var;
  string draw = var + ">>hkin_emu_neg";
  string cut = var + "<0";
  fChain->Draw(draw.c_str(),cut.c_str());
  c->Print(outname.c_str());

  draw = var + ">>hkin_emu";
  cut = var + ">0";
  fChain->Draw(draw.c_str(),cut.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_enu";
  if(usemc) var = "mc_" + var;
  draw = var + ">>hkin_enu_neg";
  cut = var + "<0";
  fChain->Draw(draw.c_str(),cut.c_str());
  c->Print(outname.c_str());

  draw = var + ">>hkin_enu";
  cut = var + ">0";
  fChain->Draw(draw.c_str(),cut.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_costhetamu";
  if(usemc) var = "mc_" + var;
  draw = var + ">>hkin_cosmu";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  // 2D kin reco vars
  var = "kin_reco_emu";
  string var2 = "kin_reco_enu";
  if(usemc) {
    var = "mc_" + var;
    var2 = "mc_" + var2;
  }
  draw = var2 + ":" + var + ">>h2_emu_enu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  var = "kin_reco_emu";
  var2 = "kin_reco_costhetamu";
  if(usemc) {
    var = "mc_" + var;
    var2 = "mc_" + var2;
  }
  draw = var2 + ":" + var + ">>h2_emu_cosmu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  var2 = "kin_reco_enu";
  var = "kin_reco_costhetamu";
  if(usemc) {
    var = "mc_" + var;
    var2 = "mc_" + var2;
  }
  draw = var2 + ":" + var + ">>h2_cosmu_enu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  string var3 = "kin_reco_emu";
  var2 = "kin_reco_enu";
  var = "kin_reco_costhetamu";
  if(usemc) {
    var = "mc_" + var;
    var2 = "mc_" + var2;
    var3 = "mc_" + var3;
  }
  draw = var2 + "-" + var3 + ":" + var + ">>h2_cosmu_nu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  // 1D resolution plots
  var = "kin_reco_pmu_diff";
  if(usemc) var = "mc_" + var;
  draw = var + ">>hpmu_diff";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_pmu_frac";
  if(usemc) 
    var = "mc_" + var;
  draw = var + ">>hpmu_frac";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_enu_diff";
  if(usemc) var = "mc_" + var;
  draw = var + ">>henu_diff";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_enu_frac";
  if(usemc) var = "mc_" + var;
  draw = var + ">>henu_frac";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  var = "kin_reco_costhetamu_diff";
  if(usemc) var = "mc_" + var;
  draw = var + ">>hcosmu_diff";
  fChain->Draw(draw.c_str());
  c->Print(outname.c_str());

  // 2D resolution
  var = "kin_reco_pmu";
  var2 = "p3_mu";
  if(usemc) {
    var = "mc_" + var;
    var2 = "mc_" + var2;
  }
  //var2 = "sqrt(" + var2 + "->Mag2()+0.105658^2)";
  draw = var2 + ":" + var + ">>h2_kin_pndr_pmu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  var2 = "kin_reco_costhetamu";
  var = "p3_mu->CosTheta()";
  if(usemc) {
    var2 = "mc_" + var2;
  }
  draw = var2 + ":" + var + ">>h2_kin_pndr_cosmu";
  fChain->Draw(draw.c_str(),"","colz");
  c->Print(outname.c_str());

  var = "kin_reco_enu";
  var2 = "mc_nu_energy";
  if(usemc) 
    var = "mc_" + var;
  draw = var + ":" + var2 + ">>h2_kin_pndr_enu";
  fChain->Draw(draw.c_str(),"","colz");

  outname_tmp = outname + "]";
  c->Print(outname.c_str());
  c->Clear();
  c->Print(outname_tmp.c_str());
}

int CountBwrdPs(vector<TVector3>* ps, float pthresh=0.){

  int n=0;
  for(auto const& vec : *ps){
    if(vec.Mag()<=pthresh) continue;
    if(vec.CosTheta()<0) n++;
  }
  return n;
}

TVector3 GetPh(vector<TVector3>* ps, float pthresh=0.){
  TVector3 ph;
  for(auto const& vec : *ps){
    if(vec.Mag()<=pthresh) continue;
    ph+=vec;
  }
  return ph;
}

void SetStyle(TH1F* h, float scale=1.0){
	h->Sumw2();
	h->SetLineWidth(2);
	h->Scale(1.0/scale);
}

void kinrecoana::BwrdPs() {

	this->SetList(true); //pass true for mc selection  

	TH1F* hmult_np     = new TH1F("hbwrdpmult_np",    "backward proton multiplicity (mc_is_signal);N_{p}^{true}",8,-0.5,7.5);
	TH1F* hmult_1ptrue = new TH1F("hbwrdpmult_1ptrue","backward proton multiplicity (mc_is_signal);N_{p}^{true}",8,-0.5,7.5);
	TH1F* hmult_1preco = new TH1F("hbwrdpmult_1preco","backward proton multiplicity (mc_is_signal);N_{p}^{true}",8,-0.5,7.5);

	TH1F* hkinenutrue_np     = new TH1F("hbkwdp_kinetrue_np",    "reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-true} [GeV]",151,-3,3);
	TH1F* hkinenutrue_1ptrue = new TH1F("hbkwdp_kinetrue_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-true} [GeV]",151,-3,3);
	TH1F* hkinenutrue_1preco = new TH1F("hbkwdp_kinetrue_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-true} [GeV]",151,-3,3);

	TH1F* hkinenureco_np     = new TH1F("hbkwdp_kinereco_np",    "reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);
	TH1F* hkinenureco_1ptrue = new TH1F("hbkwdp_kinereco_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);
	TH1F* hkinenureco_1preco = new TH1F("hbkwdp_kinereco_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{true}=0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);

	TH1F* hkinenutrue_bkwd_np     = new TH1F("hbkwdp_kinetrue_bkwd_np",    "reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-true} [GeV]",151,-3,3);
	TH1F* hkinenutrue_bkwd_1ptrue = new TH1F("hbkwdp_kinetrue_bkwd_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-true} [GeV]",151,-3,3);
	TH1F* hkinenutrue_bkwd_1preco = new TH1F("hbkwdp_kinetrue_bkwd_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-true} [GeV]",151,-3,3);

	TH1F* hkinenureco_bkwd_np     = new TH1F("hbkwdp_kinereco_bkwd_np",    "reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);
	TH1F* hkinenureco_bkwd_1ptrue = new TH1F("hbkwdp_kinereco_bkwd_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);
	TH1F* hkinenureco_bkwd_1preco = new TH1F("hbkwdp_kinereco_bkwd_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{true}>0);E_{#nu}^{kin-reco} [GeV]",151,-3,3);

	TH2F* h2kinenutrue_np     = new TH2F("h2bkwdp_kinenutrue_np",    "reco. impact from backward protons (mc_is_signal);N_{p}^{true};E_{#nu}^{kin-true} [GeV]",               6,-0.5,5.5,151,-3,3);
	TH2F* h2kinenutrue_1ptrue = new TH2F("h2bkwdp_kinenutrue_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=1);N_{p}^{true};E_{#nu}^{kin-true} [GeV]",6,-0.5,5.5,151,-3,3);
	TH2F* h2kinenutrue_1preco = new TH2F("h2bkwdp_kinenutrue_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{reco}=1);N_{p}^{true};E_{#nu}^{kin-true} [GeV]",6,-0.5,5.5,151,-3,3);

	TH2F* h2kinenureco_np     = new TH2F("h2bkwdp_kinenureco_np",    "reco. impact from backward protons (mc_is_signal);N_{p}^{true};E_{#nu}^{kin-reco} [GeV]",               6,-0.5,7.5,151,-3,3);
	TH2F* h2kinenureco_1ptrue = new TH2F("h2bkwdp_kinenureco_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=1);N_{p}^{true};E_{#nu}^{kin-reco} [GeV]",6,-0.5,7.5,151,-3,3);
	TH2F* h2kinenureco_1preco = new TH2F("h2bkwdp_kinenureco_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{reco}=1);N_{p}^{true};E_{#nu}^{kin-reco} [GeV]",6,-0.5,7.5,151,-3,3);

	TH2F* h2kinenutrue_cosh_np     = new TH2F("h2bkwdp_kinenutrue_cosh_np",    "reco. impact from backward protons (mc_is_signal);cos#theta_{h}^{true};E_{#nu}^{kin-true} [GeV]",               151,-1,1,151,-3,3);
	TH2F* h2kinenutrue_cosh_1ptrue = new TH2F("h2bkwdp_kinenutrue_cosh_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=1);cos#theta_{h}^{true};E_{#nu}^{kin-true} [GeV]",151,-1,1,151,-3,3);
	TH2F* h2kinenutrue_cosh_1preco = new TH2F("h2bkwdp_kinenutrue_cosh_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{reco}=1);cos#theta_{h}^{true};E_{#nu}^{kin-true} [GeV]",151,-1,1,151,-3,3);

	TH2F* h2kinenureco_cosh_np     = new TH2F("h2bkwdp_kinenureco_cosh_np",    "reco. impact from backward protons (mc_is_signal);cos#theta_{h}^{true};E_{#nu}^{kin-reco} [GeV]",               151,-1,1,151,-3,3);
	TH2F* h2kinenureco_cosh_1ptrue = new TH2F("h2bkwdp_kinenureco_cosh_1ptrue","reco. impact from backward protons (mc_is_signal N_{p}^{true}=1);cos#theta_{h}^{true};E_{#nu}^{kin-reco} [GeV]",151,-1,1,151,-3,3);
	TH2F* h2kinenureco_cosh_1preco = new TH2F("h2bkwdp_kinenureco_cosh_1preco","reco. impact from backward protons (mc_is_signal N_{p}^{reco}=1);cos#theta_{h}^{true};E_{#nu}^{kin-reco} [GeV]",151,-1,1,151,-3,3);

	size_t ientry=0;
	for(; ientry<fChain->GetEntries(); ientry++){

		Int_t entryNumber = fChain->GetEntryNumber(ientry);
		if(entryNumber < 0) break;
		fChain->GetEntry(entryNumber);

		int mult=CountBwrdPs(mc_p3_p_vec,0.25); //0.25 GeV momentum threshold
		float cosh = GetPh(mc_p3_p_vec,0.25).CosTheta();

		hmult_np->Fill(mult);
		if(mult==0){
			hkinenutrue_np->Fill(mc_kin_reco_enu);
			hkinenureco_np->Fill(kin_reco_enu);
		}
		else{
			hkinenutrue_bkwd_np->Fill(mc_kin_reco_enu);
			hkinenureco_bkwd_np->Fill(kin_reco_enu);
		}
		h2kinenutrue_np->Fill(mult,mc_kin_reco_enu);
		h2kinenureco_np->Fill(mult,kin_reco_enu);
		h2kinenutrue_cosh_np->Fill(cosh,mc_kin_reco_enu);
		h2kinenureco_cosh_np->Fill(cosh,kin_reco_enu);

		if(mc_np==1){ // true number of protons
			hmult_1ptrue->Fill(mult);
			if(mult==0){
				hkinenutrue_1ptrue->Fill(mc_kin_reco_enu);
				hkinenureco_1ptrue->Fill(kin_reco_enu);
			}
			else{
				hkinenutrue_bkwd_1ptrue->Fill(mc_kin_reco_enu);
				hkinenureco_bkwd_1ptrue->Fill(kin_reco_enu);
			}
			h2kinenutrue_1ptrue->Fill(mult,mc_kin_reco_enu);
			h2kinenureco_1ptrue->Fill(mult,kin_reco_enu);
			h2kinenutrue_cosh_1ptrue->Fill(cosh,mc_kin_reco_enu);
			h2kinenureco_cosh_1ptrue->Fill(cosh,kin_reco_enu);
		}

		if(np==1){ // reco number of protons
			hmult_1preco->Fill(mult);
			if(mult==0){
				hkinenutrue_1preco->Fill(mc_kin_reco_enu);
				hkinenureco_1preco->Fill(kin_reco_enu);
			}
			else {
				hkinenutrue_bkwd_1preco->Fill(mc_kin_reco_enu);
				hkinenureco_bkwd_1preco->Fill(kin_reco_enu);
			}
			h2kinenutrue_1preco->Fill(mult,mc_kin_reco_enu);
			h2kinenureco_1preco->Fill(mult,kin_reco_enu);
			h2kinenutrue_cosh_1preco->Fill(cosh,mc_kin_reco_enu);
			h2kinenureco_cosh_1preco->Fill(cosh,kin_reco_enu);
		}
	}

	SetStyle(hmult_np,ientry);
	SetStyle(hmult_1ptrue,ientry);
	SetStyle(hmult_1preco,ientry);
	SetStyle(hkinenutrue_np,ientry);
	SetStyle(hkinenutrue_1ptrue,ientry);
	SetStyle(hkinenutrue_1preco,ientry);
	SetStyle(hkinenureco_np,ientry);
	SetStyle(hkinenureco_1ptrue,ientry);
	SetStyle(hkinenureco_1preco,ientry);
	SetStyle(hkinenutrue_bkwd_np,ientry);
	SetStyle(hkinenutrue_bkwd_1ptrue,ientry);
	SetStyle(hkinenutrue_bkwd_1preco,ientry);
	SetStyle(hkinenureco_bkwd_np,ientry);
	SetStyle(hkinenureco_bkwd_1ptrue,ientry);
	SetStyle(hkinenureco_bkwd_1preco,ientry);

	hmult_np->SetLineColor(kGreen-2);
	hkinenutrue_np->SetLineColor(kGreen-2);
	hkinenureco_np->SetLineColor(kGreen-2);
	hkinenutrue_bkwd_np->SetLineColor(kGreen-2);
	hkinenureco_bkwd_np->SetLineColor(kGreen-2);

	hmult_1ptrue->SetLineStyle(kDashed);
	hkinenutrue_1ptrue->SetLineStyle(kDashed);
	hkinenureco_1ptrue->SetLineStyle(kDashed);
	hkinenutrue_bkwd_1ptrue->SetLineStyle(kDashed);
	hkinenureco_bkwd_1ptrue->SetLineStyle(kDashed);

	h2kinenutrue_np->Scale(1.0/ientry);
	h2kinenutrue_1ptrue->Scale(1.0/ientry);
	h2kinenutrue_1preco->Scale(1.0/ientry);

	h2kinenureco_np->Scale(1.0/ientry);
	h2kinenureco_1ptrue->Scale(1.0/ientry);
	h2kinenureco_1preco->Scale(1.0/ientry);

	h2kinenutrue_cosh_np->Scale(1.0/ientry);
	h2kinenutrue_cosh_1ptrue->Scale(1.0/ientry);
	h2kinenutrue_cosh_1preco->Scale(1.0/ientry);

	h2kinenureco_cosh_np->Scale(1.0/ientry);
	h2kinenureco_cosh_1ptrue->Scale(1.0/ientry);
	h2kinenureco_cosh_1preco->Scale(1.0/ientry);

	gStyle->SetOptStat(0);

	TCanvas *cmult = new TCanvas();
	cmult->SetLogy();
	hmult_np->Draw("ehist");
	hmult_1ptrue->Draw("ehistsame");
	hmult_1preco->Draw("ehistsame");

	TLegend *leg = new TLegend(0.6,0.6,0.85,0.8);
	leg->SetBorderSize(0);
	leg->AddEntry("hbwrdpmult_np",    "Np",       "l");
	leg->AddEntry("hbwrdpmult_1ptrue","1p - true","l");
	leg->AddEntry("hbwrdpmult_1preco","1p - reco","l");
	leg->Draw();

	cmult->SaveAs("backward_proton_multiplicity_truepmom_mcissignal.png");

	TCanvas *cenukintrue = new TCanvas();
	hkinenutrue_np->Draw("ehist");
	hkinenutrue_1ptrue->Draw("ehistsame");
	hkinenutrue_1preco->Draw("ehistsame");

	leg->SetX1NDC(0.2);
	leg->SetX2NDC(0.45);
	leg->Draw();

	cenukintrue->SaveAs("backward_proton_kinenutrue_truepmom_mcissignal.png");

	TCanvas *cenukinreco_bkwd = new TCanvas();
	hkinenureco_bkwd_np->Draw("ehist");
	hkinenureco_bkwd_1ptrue->Draw("ehistsame");
	hkinenureco_bkwd_1preco->Draw("ehistsame");
	leg->Draw();
	cenukinreco_bkwd->SaveAs("backward_proton_kinenureco_bkwd_truepmom_mcissignal.png");

	TCanvas *cenukintrue_bkwd = new TCanvas();
	hkinenutrue_bkwd_np->Draw("ehist");
	hkinenutrue_bkwd_1ptrue->Draw("ehistsame");
	hkinenutrue_bkwd_1preco->Draw("ehistsame");
	leg->Draw();
	cenukintrue_bkwd->SaveAs("backward_proton_kinenutrue_bkwd_truepmom_mcissignal.png");

	TCanvas* c2enukintrue_np = new TCanvas();
	c2enukintrue_np->SetLogz();
	h2kinenutrue_np->Draw("colz");
	c2enukintrue_np->SaveAs("backward_proton_kinenutrue_np_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_np = new TCanvas();
	c2enukinreco_np->SetLogz();
	h2kinenureco_np->Draw("colz");
	c2enukinreco_np->SaveAs("backward_proton_kinenureco_np_truepmom_mcissignal.png");

	TCanvas* c2enukintrue_1ptrue = new TCanvas();
	c2enukintrue_1ptrue->SetLogz();
	h2kinenutrue_1ptrue->Draw("colz");
	c2enukintrue_1ptrue->SaveAs("backward_proton_kinenutrue_1ptrue_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_1ptrue = new TCanvas();
	c2enukinreco_1ptrue->SetLogz();
	h2kinenureco_1ptrue->Draw("colz");
	c2enukinreco_1ptrue->SaveAs("backward_proton_kinenureco_1ptrue_truepmom_mcissignal.png");

	TCanvas* c2enukintrue_1preco = new TCanvas();
	c2enukintrue_1preco->SetLogz();
	h2kinenutrue_1ptrue->Draw("colz");
	c2enukintrue_1ptrue->SaveAs("backward_proton_kinenutrue_1preco_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_1preco = new TCanvas();
	c2enukinreco_1preco->SetLogz();
	h2kinenureco_1preco->Draw("colz");
	c2enukinreco_1preco->SaveAs("backward_proton_kinenureco_1preco_truepmom_mcissignal.png");


	TCanvas* c2enukintrue_cosh_np = new TCanvas();
	c2enukintrue_cosh_np->SetLogz();
	h2kinenutrue_cosh_np->Draw("colz");
	c2enukintrue_cosh_np->SaveAs("backward_proton_kinenutrue_cosh_np_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_cosh_np = new TCanvas();
	c2enukinreco_cosh_np->SetLogz();
	h2kinenureco_cosh_np->Draw("colz");
	c2enukinreco_cosh_np->SaveAs("backward_proton_kinenureco_cosh_np_truepmom_mcissignal.png");

	TCanvas* c2enukintrue_cosh_1ptrue = new TCanvas();
	c2enukintrue_cosh_1ptrue->SetLogz();
	h2kinenutrue_cosh_1ptrue->Draw("colz");
	c2enukintrue_cosh_1ptrue->SaveAs("backward_proton_kinenutrue_cosh_1ptrue_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_cosh_1ptrue = new TCanvas();
	c2enukinreco_cosh_1ptrue->SetLogz();
	h2kinenureco_cosh_1ptrue->Draw("colz");
	c2enukinreco_cosh_1ptrue->SaveAs("backward_proton_kinenureco_cosh_1ptrue_truepmom_mcissignal.png");

	TCanvas* c2enukintrue_cosh_1preco = new TCanvas();
	c2enukintrue_cosh_1preco->SetLogz();
	h2kinenutrue_cosh_1ptrue->Draw("colz");
	c2enukintrue_cosh_1ptrue->SaveAs("backward_proton_kinenutrue_cosh_1preco_truepmom_mcissignal.png");

	TCanvas* c2enukinreco_cosh_1preco = new TCanvas();
	c2enukinreco_cosh_1preco->SetLogz();
	h2kinenureco_cosh_1preco->Draw("colz");
	c2enukinreco_cosh_1preco->SaveAs("backward_proton_kinenureco_cosh_1preco_truepmom_mcissignal.png");

}

void kinrecoana::cc1p() {

	// probably should use same selection for both truth and reco
	// either mc_is_signal && select_CCNp0pi  or
	// mc_is_signal or select_CCNP0pi
	this->SetList(false,"mc_np==1&&mc_is_signal"); //pass true for mc selection

	vector<TH1F*> hists1d;
	vector<TH2F*> hists2d;

	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	string outname = "cc1p0pi_mcnp_plots_0bkgrnd.pdf";
	string outname_tmp = outname + "[";
	c->Print(outname_tmp.c_str()); //open file for writing

	// mc truth only (no backgrounds)
	for(size_t i=0; i<mc_var_names.size(); i++){

		string name = "h_" + mc_var_names[i];
		string title = ";" + mc_var_titles[i];
		hists1d.push_back( new TH1F( name.c_str(), title.c_str(), mc_vars_nbin[i], mc_vars_low[i], mc_vars_high[i]));
		hists1d.back()->SetLineWidth(2);

		string draw = mc_var_names[i] + ">>" + name;
		fChain->Draw(draw.c_str(),"","ehist");
		c->Print(outname.c_str());

		for(size_t j=i+1; j<mc_var_names.size(); j++){

			name = "h_"+ mc_var_names[i] + "_" + mc_var_names[j];
			title = ";" + mc_var_titles[i] + ";" + mc_var_titles[j];
			hists2d.push_back( new TH2F(name.c_str(), title.c_str(), 
						mc_vars_nbin[i], mc_vars_low[i], mc_vars_high[i],
						mc_vars_nbin[j], mc_vars_low[j], mc_vars_high[j]) );
			draw = mc_var_names[j] + ":"+ mc_var_names[i] + ">>" + name;
			fChain->Draw(draw.c_str(),"","colz");
			c->Print(outname.c_str());
		}
	}

	//this->SetList(false,"mc_np==1");

	// reco only ( with backgrounds)
	for(size_t i=0; i<reco_var_names.size(); i++){

		string name = "h_reco_" + reco_var_names[i];
		string title = ";" + reco_var_titles[i];
		hists1d.push_back( new TH1F( name.c_str(), title.c_str(), reco_vars_nbin[i], reco_vars_low[i], reco_vars_high[i]));
		hists1d.back()->SetLineWidth(2);

		string draw = reco_var_names[i] + ">>" + name;
		fChain->Draw(draw.c_str(),"","ehist");
		c->Print(outname.c_str());

		for(size_t j=i+1; j<reco_var_names.size(); j++){

			name = "h_reco_"+ reco_var_names[i] + "_" + reco_var_names[j];
			title = ";" + reco_var_titles[i] + ";" + reco_var_titles[j];
			hists2d.push_back( new TH2F(name.c_str(), title.c_str(),
						reco_vars_nbin[i], reco_vars_low[i], reco_vars_high[i],
						reco_vars_nbin[j], reco_vars_low[j], reco_vars_high[j]) );
			draw = reco_var_names[j] + ":"+ reco_var_names[i] + ">>" + name;
			fChain->Draw(draw.c_str(),"","colz");
			c->Print(outname.c_str());
		}
	}

	// reco vs. mc truth
	for(size_t i=0; i<mc_var_names.size(); i++){

		for(size_t j=0; j<reco_var_names.size(); j++){

			string name = "h_recotruth_"+ mc_var_names[i] + "_" + reco_var_names[j];
			string title = ";" + mc_var_titles[i] + ";" + reco_var_titles[j];
			hists2d.push_back( new TH2F(name.c_str(), title.c_str(),
						mc_vars_nbin[i], mc_vars_low[i], mc_vars_high[i],
						reco_vars_nbin[j], reco_vars_low[j], reco_vars_high[j]) );
			string draw = reco_var_names[j] + ":"+ mc_var_names[i] + ">>" + name;
			fChain->Draw(draw.c_str(),"","colz");
			c->Print(outname.c_str());
		}
	}

	outname_tmp = outname + "]";
	c->Print(outname_tmp.c_str()); //close file

	// write histo object to root file
	TFile* fout = new TFile("cc1p_histos.root","RECREATE");
	for(auto const& h : hists1d){
		h->Write();
	}
	for(auto const& h : hists2d){
		h->Write();
	}
	cout << hists1d.size()+hists2d.size() << " histograms written to file, cc1p_histos.root" << endl;
	fout->Close();

}

// normalize each row of bins independently 
// to the integral counts of that row
// bin occupancy is M_ij = mu_ij/sum_k mu_kj 
// where mu_ij is bin counts in bin with truth index j, reco index i
// error in M_ij = sqrt(mu_ij*sum_k,k!=i mu_kj / (sum_k mu_kj)^3 )
void RowNormalize(TH2D* h){

        cout << "row normalizing..." << endl;
	TH1D* hreco = h->ProjectionY("hreco");

	for(int ireco=1; ireco<h->GetNbinsY()+1; ireco++) {

		int nrow = hreco->GetBinContent(ireco); 
		if(nrow==0) continue;

                float rowfrac=0.;
		for(int itrue=1; itrue<h->GetNbinsX()+1; itrue++){
                        float a = h->GetBinContent(itrue,ireco); 
                        float b = nrow - a;
			h->SetBinContent(itrue,ireco,a/nrow);
                        h->SetBinError(itrue,ireco,sqrt(a*b/pow(a+b,3)));
                        rowfrac+=h->GetBinContent(itrue,ireco);
		}//for true bins
                cout << "total fraction in row " << ireco << ": " << rowfrac << endl;
	}//for reco bins
}

// want >~ 2/3 entries in each reco assigned to correct truth bin
// want sytematic uncertainties ~ statistical uncertainties
void kinrecoana::Binning() {

	this->SetList(false,"kin_reco_enu>0&&mc_kin_reco_enu>0&&kin_reco_pmu_diff!=9999&&mc_kin_reco_pmu_diff!=9999&&mc_kin_reco_pmu<2.5&&kin_reco_pmu<2.5&&np==1&&mc_np==1&&mc_is_signal");

	// delta_p in momentum space
	const size_t nedge = 6;
	const size_t nbins = nedge-1;
	//float bins_dp[nedge] = {-1.8,-0.8,-0.22,0.18,0.78,1.8};
	float bins_dp[nedge] = {-1.5,-0.54,-0.1,0.1,0.6,1.8};
	TH2D* hdp_mom = new TH2D("hdp_mom",";#Deltap_{true} [GeV/c];#Deltap_{reco} [GeV/c]",nbins,bins_dp,nbins,bins_dp);

	// delta_p in bin space
	TH2D* hdp_bin = new TH2D("hdp_bin",";#Deltap_{true} [bin number];#Deltap_{reco} [bin number]",nbins,0.5,nbins+0.5,nbins,0.5,nbins+0.5);
        //TH2D* hdp_bin = new TH2D("hdp_bin",";#Deltap_{true} [bin number];#Deltap_{reco} [bin number]",nbins+2,-0.5,nbins+1.5,nbins+2,-0.5,nbins+1.5);

	fChain->Draw("kin_reco_pmu_diff:mc_kin_reco_pmu_diff>>hdp_mom");
	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("8.3f");

	new TCanvas();
	TH1D* hdp_mom_reco = hdp_mom->ProjectionY("hdp_mom_reco"); //ProjectionY includes under/overflow bins by default
	hdp_mom_reco->Sumw2();
	for(int i=1; i<hdp_mom_reco->GetNbinsX()+1; i++)
		hdp_mom_reco->SetBinContent(i,1.0*hdp_mom_reco->GetBinContent(i)/hdp_mom_reco->GetBinWidth(i));
	hdp_mom_reco->Draw("ehist");

	RowNormalize(hdp_mom);
	for(int i=1; i<nbins+1; i++){ 
        //for(int i=0; i<nbins+2; i++){
		for(int j=1; j<nbins+1; j++){
                //for(int j=0; j<nbins+2; j++){
			hdp_bin->SetBinContent(i,j,hdp_mom->GetBinContent(i,j));
		}
	}  

	hdp_bin->GetXaxis()->SetNdivisions(nbins);
	hdp_bin->GetYaxis()->SetNdivisions(nbins);
	TCanvas *c = new TCanvas();
	hdp_bin->Draw("colztexte");

	// delta_costheta

}

void SetResolutionStyle(TH1D* h, bool norm=true){

	h->Sumw2();
	h->SetLineWidth(2);
	if(norm) h->Scale(1.0/h->Integral()); //0,-1));

}

void kinrecoana::Resolution(){

	// set entry list
	string mucuts="kin_reco_enu>0&&mc_kin_reco_enu>0&&kin_reco_emu>0.106&&mc_kin_reco_emu>0.106"; //apply this only to muon stuff, trying to select only physical values
	this->SetList(false,"kin_reco_pmu_diff!=9999&&mc_kin_reco_pmu_diff!=9999&&mc_is_signal"); //sel_CCNp0pi
	fChain->SetMakeClass(0);

	///// LOTS of histograms, avoiding duplication where possible using TH1::Add and TH2::ProjectionX/Y
	///// Create 1 set of histos for each subsample: CC0pi1p and CC0piNp, N>1
	///// N.B. histogram binning must match between subsamples
	// muon kinematics
	TH2D* hdp_diff   = new TH2D("hdp_diff",  "CC0#piNp (N#geq2);p_{#mu}^{kin-true}-p_{#mu}^{true} [GeV/c];p_{#mu}^{kin-reco}-p_{#mu}^{pndr} [GeV/c]",                       251,-2,3, 251,-2,3);
	TH2D* hdp_frac   = new TH2D("hdp_frac",  "CC0#piNp (N#geq2);#Deltap_{#mu}^{true}/p_{#mu}^{true};#Deltap_{#mu}^{reco}/p_{#mu}^{pndr}",                                   151,-1.1,5, 151,-1.1,5);
	TH2D* hdcos_diff = new TH2D("hdcos_diff","CC0#piNp (N#geq2);cos#theta_{#mu}^{kin-true}-cos#theta_{#mu}^{true};cos#theta_{#mu}^{kin-reco}-cos#theta_{#mu}^{pndr}",       251,-2.1,2.1, 251,-2.1,2.1);

	TH2D* hdp_1p_diff   = new TH2D("hdp_1p_diff",  "CC0#pi1p;p_{#mu}^{kin-true}-p_{#mu}^{true} [GeV/c];p_{#mu}^{kin-reco}-p_{#mu}^{pndr} [GeV/c]",                          251,-2,3, 251,-2,3);
	TH2D* hdp_1p_frac   = new TH2D("hdp_1p_frac",  "CC0#pi1p;#Deltap_{#mu}^{true}/p_{#mu}^{true};#Deltap_{#mu}^{reco}/p_{#mu}^{pndr}",                                      151,-1.1,5, 151,-1.1,5);
	TH2D* hdcos_1p_diff = new TH2D("hdcos_1p_diff","CC0#pi1p;cos#theta_{#mu}^{kin-true}-cos#theta_{#mu}^{true};cos#theta_{#mu}^{kin-reco}-cos#theta_{#mu}^{pndr}",          251,-2.1,2.1, 251,-2.1,2.1);

	// delta delta is a combination of kin reco and conventional resolutions kin reco resolution 
	TH1D* hddp_diff   = new TH1D("hddp_diff",  ";#Deltap_{#mu}^{reco}-#Deltap_{#mu}^{true} [GeV/c]", 301,-3,3);
	TH1D* hddp_frac   = new TH1D("hddp_frac",  ";(#Deltap_{#mu}^{reco}-#Deltap_{#mu}^{true})/#Deltap_{#mu}^{true}", 201,-5,5);
	TH1D* hddcos_diff = new TH1D("hddcos_diff",";#Deltacos#theta_{#mu}^{reco}-#Deltacos#theta_{#mu}^{true}", 401,-2.1,2.1);

	TH1D* hddp_1p_diff   = new TH1D("hddp_1p_diff", ";#Deltap_{#mu}^{reco}-#Deltap_{#mu}^{true} [GeV/c]", 301,-3,3);
	TH1D* hddp_1p_frac   = new TH1D("hddp_1p_frac", ";(#Deltap_{#mu}^{reco}-#Deltap_{#mu}^{true})/#Deltap_{#mu}^{true}", 201,-5,5);
	TH1D* hddcos_1p_diff = new TH1D("hddcos_1p_diff",";#Deltacos#theta_{#mu}^{reco}-#Deltacos#theta_{#mu}^{true}", 401,-2.1,2.1);

	// neutrino energy (separate positive and negative components)
	TH1D* henu_diff_true = new TH1D("henu_diff_true",";E_{#nu}^{kin-true}-E_{#nu}^{true} [GeV]",               301,-4,4); 
	TH1D* henu_frac_true = new TH1D("henu_frac_true",";(E_{#nu}^{kin-true}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);
	TH1D* henu_diff_reco = new TH1D("henu_diff_reco",";E_{#nu}^{kin-reco}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_frac_reco = new TH1D("henu_frac_reco",";(E_{#nu}^{kin-reco}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);

	TH1D* henu_1p_diff_true = new TH1D("henu_1p_diff_true",";E_{#nu}^{kin-true}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_1p_frac_true = new TH1D("henu_1p_frac_true",";(E_{#nu}^{kin-true}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);
	TH1D* henu_1p_diff_reco = new TH1D("henu_1p_diff_reco",";E_{#nu}^{kin-reco}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_1p_frac_reco = new TH1D("henu_1p_frac_reco",";(E_{#nu}^{kin-reco}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);

	TH1D* henu_diff_true_neg = new TH1D("henu_diff_true_neg",";E_{#nu}^{kin-true}-E_{#nu}^{true} [GeV]",               301,-4,4); 
	TH1D* henu_frac_true_neg = new TH1D("henu_frac_true_neg",";(E_{#nu}^{kin-true}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);
	TH1D* henu_diff_reco_neg = new TH1D("henu_diff_reco_neg",";E_{#nu}^{kin-reco}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_frac_reco_neg = new TH1D("henu_frac_reco_neg",";(E_{#nu}^{kin-reco}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);

	TH1D* henu_1p_diff_true_neg = new TH1D("henu_1p_diff_true_neg",";E_{#nu}^{kin-true}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_1p_frac_true_neg = new TH1D("henu_1p_frac_true_neg",";(E_{#nu}^{kin-true}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);
	TH1D* henu_1p_diff_reco_neg = new TH1D("henu_1p_diff_reco_neg",";E_{#nu}^{kin-reco}-E_{#nu}^{true} [GeV]",               301,-4,4);
	TH1D* henu_1p_frac_reco_neg = new TH1D("henu_1p_frac_reco_neg",";(E_{#nu}^{kin-reco}-E_{#nu}^{true})/E_{#nu}^{true}",    201,-1.1,4);

	// hadronic system kinematics
	TH2D* hphad         = new TH2D("hphad",        "CC0#piNp (N#geq2);p_{had}^{true} [GeV/c]; p_{had}^{reco} [GeV/c]",301,0,2.5,301,0,2.5);
	TH2D* hphad_1p      = new TH2D("hphad_1p",     "CC0#pi1p;p_{had}^{true} [GeV/c]; p_{had}^{reco} [GeV/c]",        301,0,1.5,301,0,1.5);
	TH1D* hphad_diff    = new TH1D("hphad_diff",   ";p_{had}^{reco} - p_{had}^{true} [GeV/c]",       401,-1,1.5);
	TH1D* hphad_1p_diff = new TH1D("hphad_1p_diff",";p_{had}^{reco} - p_{had}^{true} [GeV/c]",       401,-1,1.5);

	TH2D* hcoshad         = new TH2D("hcoshad",        "CC0#piNp (N#geq2);cos#theta_{had}^{true}; cos#theta_{had}^{reco}", 301,-1.1,1.1,301,-1.1,1.1);
	TH2D* hcoshad_1p      = new TH2D("hcoshad_1p",     "CC0#pi1p;cos#theta_{had}^{true}; cos#theta_{had}^{reco}", 301,-1.1,1.1,301,-1.1,1.1);
	TH1D* hcoshad_diff    = new TH1D("hcoshad_diff",   ";cos#theta_{had}^{reco} - cos#theta_{had}^{true}",301,-1.1,1.1);
	TH1D* hcoshad_1p_diff = new TH1D("hcoshad_1p_diff",";cos#theta_{had}^{reco} - cos#theta_{had}^{true}",301,-1.1,1.1);

	TH2D* hwhad         = new TH2D("hwhad",         "CC0#piNp (N#geq2);W_{true} [GeV];W_{reco} [GeV]", 201,-2,2,201,-3,1.5);
	TH2D* hwhad_1p      = new TH2D("hwhad_1p",      "CC0#pi1p;W_{true} [GeV];W_{reco} [GeV]",         201,0,1.5,201,0,1.5);
	TH1D* hwhad_diff    = new TH1D("hwhad_diff",    ";W_{reco} - W_{true} [GeV]",     401,-4,3);
	TH1D* hwhad_1p_diff = new TH1D("hwhad_1p_diff", ";W_{reco} - W_{true} [GeV]",     401,-4,3);

	TH2D* hwphad_true      = new TH2D("hwphad_true",   "CC0#piNp (N#geq2);p_{had}^{true} [GeV/c];W^{true} [GeV]",301,0,3,301,-3,2);
	TH2D* hwphad_reco      = new TH2D("hwphad_reco",   "CC0#piNp (N#geq2);p_{had}^{reco} [GeV/c];W^{reco} [GeV]",301,0,3,301,-3,2);
	TH2D* hwphad_1p_true   = new TH2D("hwphad_1p_true","CC0#pi1p;p_{had}^{true} [GeV/c];W^{true} [GeV]",         301,0,3,301,0,2);
	TH2D* hwphad_1p_reco   = new TH2D("hwphad_1p_reco","CC0#pi1p;p_{had}^{reco} [GeV/c];W^{reco} [GeV]",         301,0,3,301,0,2);

	TH2D* hwcoshad_true    = new TH2D("hwcoshad_true",   "CC0#piNp (N#geq2);cos#theta_{had}^{true} [GeV/c];W^{true} [GeV]",301,-1.1,1.1,301,-3,2);
	TH2D* hwcoshad_reco    = new TH2D("hwcoshad_reco",   "CC0#piNp (N#geq2);cos#theta_{had}^{reco} [GeV/c];W^{reco} [GeV]",301,-1.1,1.1,301,-3,2);
	TH2D* hwcoshad_1p_true = new TH2D("hwcoshad_1p_true","CC0#pi1p;cos#theta_{had}^{true} [GeV/c];W^{true} [GeV]",        301,-1.1,1.1,301,0,2);
	TH2D* hwcoshad_1p_reco = new TH2D("hwcoshad_1p_reco","CC0#pi1p;cos#theta_{had}^{reco} [GeV/c];W^{reco} [GeV]",        301,-1.1,1.1,301,0,2);

	// fill histograms
	// Np (N>1) first
	string mucuts_np = mucuts + "&&mc_np>1&&np>1";
	fChain->Draw("kin_reco_pmu_diff:mc_kin_reco_pmu_diff>>hdp_diff",mucuts_np.c_str());
	fChain->Draw("kin_reco_pmu_frac:mc_kin_reco_pmu_frac>>hdp_frac",mucuts_np.c_str());
	fChain->Draw("kin_reco_costhetamu_diff:mc_kin_reco_costhetamu_diff>>hdcos_diff",mucuts_np.c_str());

	fChain->Draw("kin_reco_pmu_diff-mc_kin_reco_pmu_diff>>hddp_diff",mucuts_np.c_str());
	fChain->Draw("(kin_reco_pmu_diff-mc_kin_reco_pmu_diff)/mc_kin_reco_pmu_diff>>hddp_frac",mucuts_np.c_str());
	fChain->Draw("kin_reco_costhetamu_diff-mc_kin_reco_costhetamu_diff>>hddcos_diff",mucuts_np.c_str());

	fChain->Draw("mc_kin_reco_enu_diff>>henu_diff_true",mucuts_np.c_str());
	fChain->Draw("mc_kin_reco_enu_frac>>henu_frac_true",mucuts_np.c_str());
	fChain->Draw("kin_reco_enu_diff>>henu_diff_reco",   mucuts_np.c_str());
	fChain->Draw("kin_reco_enu_frac>>henu_frac_reco",   mucuts_np.c_str());

	fChain->Draw("-1*mc_kin_reco_enu-mc_nu_energy>>henu_diff_true_neg",               "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np>1&&np>1");
	fChain->Draw("(-1*mc_kin_reco_enu-mc_nu_energy)/mc_nu_energy>>henu_frac_true_neg","kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np>1&&np>1");
	fChain->Draw("-1*kin_reco_enu-mc_nu_energy>>henu_diff_reco_neg",                  "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np>1&&np>1");
	fChain->Draw("(-1*kin_reco_enu-mc_nu_energy)/mc_nu_energy>>henu_frac_reco_neg",   "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np>1&&np>1");

	fChain->Draw("p4_had->P():mc_p4_had->P()>>hphad",                     "mc_np>1&&np>1");
	fChain->Draw("p4_had->P()-mc_p4_had->P()>>hphad_diff",                "mc_np>1&&np>1");
	fChain->Draw("p4_had->CosTheta():mc_p4_had->CosTheta()>>hcoshad",     "mc_np>1&&np>1");
	fChain->Draw("p4_had->CosTheta()-mc_p4_had->CosTheta()>>hcoshad_diff","mc_np>1&&np>1");
	fChain->Draw("p4_had->M():mc_p4_had->M()>>hwhad",                     "mc_np>1&&np>1");
	fChain->Draw("p4_had->M()-mc_p4_had->M()>>hwhad_diff",                "mc_np>1&&np>1");

	fChain->Draw("mc_p4_had->M():mc_p4_had->P()>>hwphad_true",            "mc_np>1&&np>1");
	fChain->Draw("p4_had->M():p4_had->P()>>hwphad_reco",                  "mc_np>1&&np>1");
	fChain->Draw("mc_p4_had->M():mc_p4_had->CosTheta()>>hwcoshad_true",   "mc_np>1&&np>1");
	fChain->Draw("p4_had->M():p4_had->CosTheta()>>hwcoshad_reco",         "mc_np>1&&np>1");

	// now repeat for 1p only
	string mucuts_1p = mucuts + "&&mc_np==1&&np==1";
	fChain->Draw("kin_reco_pmu_diff:mc_kin_reco_pmu_diff>>hdp_1p_diff",                mucuts_1p.c_str());
	fChain->Draw("kin_reco_pmu_frac:mc_kin_reco_pmu_frac>>hdp_1p_frac",                mucuts_1p.c_str());
	fChain->Draw("kin_reco_costhetamu_diff:mc_kin_reco_costhetamu_diff>>hdcos_1p_diff",mucuts_1p.c_str());

	fChain->Draw("kin_reco_pmu_diff-mc_kin_reco_pmu_diff>>hddp_1p_diff",                       mucuts_1p.c_str());
	fChain->Draw("(kin_reco_pmu_diff-mc_kin_reco_pmu_diff)/mc_kin_reco_pmu_diff>>hddp_1p_frac",mucuts_1p.c_str());
	fChain->Draw("kin_reco_costhetamu_diff-mc_kin_reco_costhetamu_diff>>hddcos_1p_diff",       mucuts_1p.c_str());

	fChain->Draw("mc_kin_reco_enu_diff>>henu_1p_diff_true",mucuts_1p.c_str());
	fChain->Draw("mc_kin_reco_enu_frac>>henu_1p_frac_true",mucuts_1p.c_str());
	fChain->Draw("kin_reco_enu_diff>>henu_1p_diff_reco",   mucuts_1p.c_str());
	fChain->Draw("kin_reco_enu_frac>>henu_1p_frac_reco",   mucuts_1p.c_str());

	fChain->Draw("-1*mc_kin_reco_enu-mc_nu_energy>>henu_1p_diff_true_neg",               "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np==1&&np==1");
	fChain->Draw("(-1*mc_kin_reco_enu-mc_nu_energy)/mc_nu_energy>>henu_1p_frac_true_neg","kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np==1&&np==1");
	fChain->Draw("-1*kin_reco_enu-mc_nu_energy>>henu_1p_diff_reco_neg",                  "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np==1&&np==1");
	fChain->Draw("(-1*kin_reco_enu-mc_nu_energy)/mc_nu_energy>>henu_1p_frac_reco_neg",   "kin_reco_enu<0&&mc_kin_reco_enu<0&&mc_np==1&&np==1");

	fChain->Draw("p4_had->P():mc_p4_had->P()>>hphad_1p",                      "mc_np==1&&np==1" );
	fChain->Draw("p4_had->P()-mc_p4_had->P()>>hphad_1p_diff",                 "mc_np==1&&np==1" );
	fChain->Draw("p4_had->CosTheta():mc_p4_had->CosTheta()>>hcoshad_1p",      "mc_np==1&&np==1" );
	fChain->Draw("p4_had->CosTheta()-mc_p4_had->CosTheta()>>hcoshad_1p_diff", "mc_np==1&&np==1" );
	fChain->Draw("p4_had->M():mc_p4_had->M()>>hwhad_1p",                      "mc_np==1&&np==1" );
	fChain->Draw("p4_had->M()-mc_p4_had->M()>>hwhad_1p_diff",                 "mc_np==1&&np==1" );

	fChain->Draw("mc_p4_had->M():mc_p4_had->P()>>hwphad_1p_true",             "mc_np==1&&np==1" );
	fChain->Draw("p4_had->M():p4_had->P()>>hwphad_1p_reco",                   "mc_np==1&&np==1" );
	fChain->Draw("mc_p4_had->M():mc_p4_had->CosTheta()>>hwcoshad_1p_true",    "mc_np==1&&np==1" );
	fChain->Draw("p4_had->M():p4_had->CosTheta()>>hwcoshad_1p_reco",          "mc_np==1&&np==1" );


	// extract 1D histos from 2D ones and make x-axis label truth-reco agnostic
	TH1D* hdp_diff_true   = (TH1D*)hdp_diff->ProjectionX("hdp_diff_true",  0, -1, "e"); 
	TH1D* hdp_frac_true   = (TH1D*)hdp_frac->ProjectionX("hdp_frac_true",  0, -1, "e"); 
	TH1D* hdcos_diff_true = (TH1D*)hdcos_diff->ProjectionX("hdcos_diff_true",0, -1, "e"); 

	hdp_diff_true->GetXaxis()->SetTitle("p_{#mu}^{kin}-p_{#mu}^{ref} [GeV/c]");  
	hdp_frac_true->GetXaxis()->SetTitle("(p_{#mu}^{kin}-p_{#mu}^{ref})/p_{#mu}^{ref}"); 
	hdcos_diff_true->GetXaxis()->SetTitle("cos#theta_{#mu}^{kin}-cos#theta_{#mu}^{ref}");

	TH1D* hdp_diff_reco   = (TH1D*)hdp_diff->ProjectionY("hdp_diff_reco",  0, -1, "e"); 
	TH1D* hdp_frac_reco   = (TH1D*)hdp_frac->ProjectionY("hdp_frac_reco",  0, -1, "e"); 
	TH1D* hdcos_diff_reco = (TH1D*)hdcos_diff->ProjectionY("hdcos_diff_reco",0, -1, "e");

	TH1D* hdp_1p_diff_true   = (TH1D*)hdp_1p_diff->ProjectionX("hdp_1p_diff_true",  0, -1, "e"); 
	TH1D* hdp_1p_frac_true   = (TH1D*)hdp_1p_frac->ProjectionX("hdp_1p_frac_true",  0, -1, "e"); 
	TH1D* hdcos_1p_diff_true = (TH1D*)hdcos_1p_diff->ProjectionX("hdcos_1p_diff_true",0, -1, "e");

	TH1D* hdp_1p_diff_reco   = (TH1D*)hdp_1p_diff->ProjectionY("hdp_1p_diff_reco",  0, -1, "e");
	TH1D* hdp_1p_frac_reco   = (TH1D*)hdp_1p_frac->ProjectionY("hdp_1p_frac_reco",  0, -1, "e");
	TH1D* hdcos_1p_diff_reco = (TH1D*)hdcos_1p_diff->ProjectionY("hdcos_1p_diff_reco",0, -1, "e");

	TH1D* hphad_true   = (TH1D*)hphad->ProjectionX("hphad_true",0,-1,"e");  
	TH1D* hphad_reco   = (TH1D*)hphad->ProjectionY("hphad_reco",0,-1,"e"); 
	TH1D* hcoshad_true = (TH1D*)hcoshad->ProjectionX("hcoshad_true",0,-1,"e"); 
	TH1D* hcoshad_reco = (TH1D*)hcoshad->ProjectionY("hcoshad_reco",0,-1,"e"); 
	TH1D* hwhad_true   = (TH1D*)hwhad->ProjectionX("hwhad_true",0,-1,"e"); 
	TH1D* hwhad_reco   = (TH1D*)hwhad->ProjectionY("hwhad_reco",0,-1,"e"); 

	TH1D* hphad_1p_true   = (TH1D*)hphad_1p->ProjectionX("hphad_1p_true",0,-1,"e");
	TH1D* hphad_1p_reco   = (TH1D*)hphad_1p->ProjectionY("hphad_1p_reco",0,-1,"e");
	TH1D* hcoshad_1p_true = (TH1D*)hcoshad_1p->ProjectionX("hcoshad_1p_true",0,-1,"e");
	TH1D* hcoshad_1p_reco = (TH1D*)hcoshad_1p->ProjectionY("hcoshad_1p_reco",0,-1,"e");
	TH1D* hwhad_1p_true   = (TH1D*)hwhad_1p->ProjectionX("hwhad_1p_true",0,-1,"e");
	TH1D* hwhad_1p_reco   = (TH1D*)hwhad_1p->ProjectionY("hwhad_1p_reco",0,-1,"e");

	hphad_1p_true  ->GetXaxis()->SetTitle("p_{had} [GeV/c]");
	hcoshad_1p_true->GetXaxis()->SetTitle("cos#theta_{had}");
	hwhad_1p_true  ->GetXaxis()->SetTitle("W [GeV]");

	// style
	int color_np_true = kGreen-4;
	int color_np_reco = kBlue;
	int color_1p_true = kGreen-2;
	int color_1p_reco = kBlue+2;

	hdp_diff->GetYaxis()->SetTitleOffset(1.4); 
	hdp_frac->GetYaxis()->SetTitleOffset(1.4);
	hdcos_diff->GetYaxis()->SetTitleOffset(1.4);

	hdp_1p_diff->GetYaxis()->SetTitleOffset(1.4);
	hdp_1p_frac->GetYaxis()->SetTitleOffset(1.4);
	hdcos_1p_diff->GetYaxis()->SetTitleOffset(1.4);

	hphad->GetYaxis()->SetTitleOffset(1.4);
	hphad_1p->GetYaxis()->SetTitleOffset(1.4);
	hcoshad->GetYaxis()->SetTitleOffset(1.4);
	hcoshad_1p->GetYaxis()->SetTitleOffset(1.4);
	hwhad->GetYaxis()->SetTitleOffset(1.4);
	hwhad_1p->GetYaxis()->SetTitleOffset(1.4);

	hwphad_true->GetYaxis()->SetTitleOffset(1.4);     
	hwphad_reco     ->GetYaxis()->SetTitleOffset(1.4);
	hwphad_1p_true  ->GetYaxis()->SetTitleOffset(1.4);
	hwphad_1p_reco  ->GetYaxis()->SetTitleOffset(1.4);

	hwcoshad_true   ->GetYaxis()->SetTitleOffset(1.4);
	hwcoshad_reco   ->GetYaxis()->SetTitleOffset(1.4);
	hwcoshad_1p_true->GetYaxis()->SetTitleOffset(1.4);
	hwcoshad_1p_reco->GetYaxis()->SetTitleOffset(1.4);

	SetResolutionStyle(hdp_diff_true);   hdp_diff_true->SetLineColor(color_np_true); 
	SetResolutionStyle(hdp_frac_true);   hdp_frac_true->SetLineColor(color_np_true);
	SetResolutionStyle(hdcos_diff_true); hdcos_diff_true->SetLineColor(color_np_true);

	SetResolutionStyle(hdp_1p_diff_true);   hdp_1p_diff_true->SetLineColor(color_1p_true);
	SetResolutionStyle(hdp_1p_frac_true);   hdp_1p_frac_true->SetLineColor(color_1p_true);
	SetResolutionStyle(hdcos_1p_diff_true); hdcos_1p_diff_true->SetLineColor(color_1p_true);

	SetResolutionStyle(hdp_diff_reco);   hdp_diff_reco->SetLineColor(color_np_reco); 
	SetResolutionStyle(hdp_frac_reco);   hdp_frac_reco->SetLineColor(color_np_reco);
	SetResolutionStyle(hdcos_diff_reco); hdcos_diff_reco->SetLineColor(color_np_reco);

	SetResolutionStyle(hdp_1p_diff_reco);   hdp_1p_diff_reco->SetLineColor(color_1p_reco);
	SetResolutionStyle(hdp_1p_frac_reco);   hdp_1p_frac_reco->SetLineColor(color_1p_reco);
	SetResolutionStyle(hdcos_1p_diff_reco); hdcos_1p_diff_reco->SetLineColor(color_1p_reco);

	SetResolutionStyle(hddp_diff);   hddp_diff->SetLineColor(color_np_reco);
	SetResolutionStyle(hddp_frac);   hddp_frac->SetLineColor(color_np_reco);
	SetResolutionStyle(hddcos_diff); hddcos_diff->SetLineColor(color_np_reco);

	SetResolutionStyle(hddp_1p_diff);   hddp_1p_diff->SetLineColor(color_1p_reco);
	SetResolutionStyle(hddp_1p_frac);   hddp_1p_frac->SetLineColor(color_1p_reco);
	SetResolutionStyle(hddcos_1p_diff); hddcos_1p_diff->SetLineColor(color_1p_reco);

	SetResolutionStyle(henu_diff_true_neg,false); henu_diff_true_neg->SetLineColor(color_np_true); henu_diff_true_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_frac_true_neg,false); henu_frac_true_neg->SetLineColor(color_np_true); henu_frac_true_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_diff_reco_neg,false); henu_diff_reco_neg->SetLineColor(color_np_reco); henu_diff_reco_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_frac_reco_neg,false); henu_frac_reco_neg->SetLineColor(color_np_reco); henu_frac_reco_neg->SetLineStyle(kDashed);

	SetResolutionStyle(henu_1p_diff_true_neg,false); henu_1p_diff_true_neg->SetLineColor(color_1p_true);  henu_1p_diff_true_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_1p_frac_true_neg,false); henu_1p_frac_true_neg->SetLineColor(color_1p_true);  henu_1p_frac_true_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_1p_diff_reco_neg,false); henu_1p_diff_reco_neg->SetLineColor(color_1p_reco); henu_1p_diff_reco_neg->SetLineStyle(kDashed);
	SetResolutionStyle(henu_1p_frac_reco_neg,false); henu_1p_frac_reco_neg->SetLineColor(color_1p_reco); henu_1p_frac_reco_neg->SetLineStyle(kDashed);

	henu_diff_true_neg->Scale(1.0/henu_diff_true->Integral()); //0,-1));
	henu_frac_true_neg->Scale(1.0/henu_frac_true->Integral()); //0,-1));
	henu_diff_reco_neg->Scale(1.0/henu_diff_reco->Integral()); //0,-1));
	henu_frac_reco_neg->Scale(1.0/henu_frac_reco->Integral()); //0,-1));

	henu_1p_diff_true_neg->Scale(1.0/henu_1p_diff_true->Integral()); //0,-1));
	henu_1p_frac_true_neg->Scale(1.0/henu_1p_frac_true->Integral()); //0,-1));
	henu_1p_diff_reco_neg->Scale(1.0/henu_1p_diff_reco->Integral()); //0,-1));
	henu_1p_frac_reco_neg->Scale(1.0/henu_1p_frac_reco->Integral()); //0,-1));

	SetResolutionStyle(henu_diff_true); henu_diff_true->SetLineColor(color_np_true);
	SetResolutionStyle(henu_frac_true); henu_frac_true->SetLineColor(color_np_true);
	SetResolutionStyle(henu_diff_reco); henu_diff_reco->SetLineColor(color_np_reco);
	SetResolutionStyle(henu_frac_reco); henu_frac_reco->SetLineColor(color_np_reco);

  SetResolutionStyle(henu_1p_diff_true); henu_1p_diff_true->SetLineColor(color_1p_true);
  SetResolutionStyle(henu_1p_frac_true); henu_1p_frac_true->SetLineColor(color_1p_true);
  SetResolutionStyle(henu_1p_diff_reco); henu_1p_diff_reco->SetLineColor(color_1p_reco);
  SetResolutionStyle(henu_1p_frac_reco); henu_1p_frac_reco->SetLineColor(color_1p_reco);

  henu_1p_diff_true->GetXaxis()->SetTitle("E_{#nu}^{kin}-E_{#nu}^{true} [GeV]");
  henu_1p_frac_true->GetXaxis()->SetTitle("(E_{#nu}^{kin}-E_{#nu}^{true})/E_{#nu}^{true}");

  SetResolutionStyle(hphad_diff  );  hphad_diff  ->SetLineColor(color_np_reco);
  SetResolutionStyle(hcoshad_diff);  hcoshad_diff->SetLineColor(color_np_reco);
  SetResolutionStyle(hwhad_diff  );  hwhad_diff  ->SetLineColor(color_np_reco);

  SetResolutionStyle(hphad_1p_diff  ); hphad_1p_diff  ->SetLineColor(color_1p_reco);
  SetResolutionStyle(hcoshad_1p_diff); hcoshad_1p_diff->SetLineColor(color_1p_reco);
  SetResolutionStyle(hwhad_1p_diff  ); hwhad_1p_diff  ->SetLineColor(color_1p_reco);

  SetResolutionStyle(hphad_true  );  hphad_true  ->SetLineColor(color_np_true);   
  SetResolutionStyle(hphad_reco  );  hphad_reco  ->SetLineColor(color_np_reco);   
  SetResolutionStyle(hcoshad_true);  hcoshad_true->SetLineColor(color_np_true);     
  SetResolutionStyle(hcoshad_reco);  hcoshad_reco->SetLineColor(color_np_reco);    
  SetResolutionStyle(hwhad_true  );  hwhad_true  ->SetLineColor(color_np_true);    
  SetResolutionStyle(hwhad_reco  );  hwhad_reco  ->SetLineColor(color_np_reco);   

  SetResolutionStyle(hphad_1p_true  );  hphad_1p_true  ->SetLineColor(color_1p_true);
  SetResolutionStyle(hphad_1p_reco  );  hphad_1p_reco  ->SetLineColor(color_1p_reco);
  SetResolutionStyle(hcoshad_1p_true);  hcoshad_1p_true->SetLineColor(color_1p_true);
  SetResolutionStyle(hcoshad_1p_reco);  hcoshad_1p_reco->SetLineColor(color_1p_reco);
  SetResolutionStyle(hwhad_1p_true  );  hwhad_1p_true  ->SetLineColor(color_1p_true);
  SetResolutionStyle(hwhad_1p_reco  );  hwhad_1p_reco  ->SetLineColor(color_1p_reco);

  gStyle->SetOptStat(0);

  // draw
  auto cdmu2d = new TCanvas("cdmu2d","",1800,1000);
  cdmu2d->Divide(3,2);
  cdmu2d->cd(1);
  hdp_diff->Draw("colz");
  cdmu2d->cd(2);
  hdp_frac->Draw("colz");
  cdmu2d->cd(3);
  hdcos_diff->Draw("colz");

  cdmu2d->cd(4);
  hdp_1p_diff->Draw("colz");
  cdmu2d->cd(5);
  hdp_1p_frac->Draw("colz");
  cdmu2d->cd(6);
  hdcos_1p_diff->Draw("colz");
  cdmu2d->SaveAs("deltas2d.png");

  auto cdmu = new TCanvas("cdmu","",1800,600);
  cdmu->Divide(3,1);
  cdmu->cd(1);
  hdp_1p_diff_true->Draw("ehist");
  hdp_1p_diff_reco->Draw("ehistsame");
  hdp_diff_true->Draw("ehistsame");
  hdp_diff_reco->Draw("ehistsame");

  TLegend* legmu = new TLegend(0.6,0.4,0.9,0.9);
  legmu->SetBorderSize(0);
  legmu->AddEntry("hdp_diff_true","Np w/truth input","l");
  legmu->AddEntry("hdp_diff_reco","Np w/reco input","l");
  legmu->AddEntry("hdp_1p_diff_true","1p w/truth input","l");
  legmu->AddEntry("hdp_1p_diff_reco","1p w/reco input","l");

  cdmu->cd(2);
  hdp_1p_frac_true->Draw("ehist");
  hdp_1p_frac_reco->Draw("ehistsame");
  hdp_frac_true->Draw("ehistsame");
  hdp_frac_reco->Draw("ehistsame");

  cdmu->cd(3);
  hdcos_1p_diff_true->Draw("ehist");
  hdcos_1p_diff_reco->Draw("ehistsame");
  hdcos_diff_true->Draw("ehistsame");
  hdcos_diff_reco->Draw("ehistsame");
  legmu->Draw();

  cdmu->SaveAs("deltas.png");

  auto cddmu = new TCanvas("cddmu","",1800,600);
  cddmu->Divide(3,1);
  cddmu->cd(1);
  hddp_1p_diff->Draw("e0hist");
  hddp_1p_diff->GetXaxis()->SetRangeUser(-2,2);
  hddp_diff->Draw("e0histsame");
  cddmu->cd(2);
  hddp_1p_frac->Draw("e0hist");
  hddp_frac->Draw("e0histsame");
  cddmu->cd(3);
  hddcos_1p_diff->Draw("e0hist");
  hddcos_diff->Draw("e0histsame");
  cddmu->SaveAs("deltadeltas.png");

  auto cnu = new TCanvas("cnu","",1200,600);
  cnu->Divide(2,1);
  cnu->cd(1);
  henu_1p_diff_true->Draw("e0hist");
  henu_1p_diff_reco->Draw("e0histsame");
  henu_1p_diff_true_neg->Draw("e0histsame");
  henu_1p_diff_reco_neg->Draw("e0histsame");

  henu_diff_true->Draw("e0histsame");
  henu_diff_reco->Draw("e0histsame");
  henu_diff_true_neg->Draw("e0histsame");
  henu_diff_reco_neg->Draw("e0histsame");

  TLegend* legnu = new TLegend(0.35,0.5,0.88,0.88);
  legnu->SetBorderSize(0);
  legnu->AddEntry("henu_diff_true",     "Np w/truth inputs", "l");
  legnu->AddEntry("henu_diff_reco",     "Np w/reco inputs", "l");
  legnu->AddEntry("henu_diff_true_neg", "Np w/truth inputs, E_{#nu}^{kin}<0 and E_{#nu}^{kin}#rightarrow-E_{#nu}^{kin}", "l");
  legnu->AddEntry("henu_diff_reco_neg", "Np w/reco inputs, E_{#nu}^{kin}<0 and E_{#nu}^{kin}#rightarrow-E_{#nu}^{kin}", "l");
  legnu->AddEntry("henu_1p_diff_true",     "1p w/truth inputs", "l");
  legnu->AddEntry("henu_1p_diff_reco",     "1p w/reco inputs", "l");
  legnu->AddEntry("henu_1p_diff_true_neg", "1p w/truth inputs, E_{#nu}^{kin}<0 and E_{#nu}^{kin}#rightarrow-E_{#nu}^{kin}", "l");
  legnu->AddEntry("henu_1p_diff_reco_neg", "1p w/reco inputs, E_{#nu}^{kin}<0 and E_{#nu}^{kin}#rightarrow-E_{#nu}^{kin}", "l");

  cnu->cd(2);
  henu_1p_frac_true->Draw("e0hist");
  henu_1p_frac_reco->Draw("e0histsame");
  henu_1p_frac_true_neg->Draw("e0histsame");
  henu_1p_frac_reco_neg->Draw("e0histsame");
  henu_frac_true->Draw("e0histsame");
  henu_frac_reco->Draw("e0histsame");
  henu_frac_true_neg->Draw("e0histsame");
  henu_frac_reco_neg->Draw("e0histsame");
  legnu->Draw();

  cnu->SaveAs("enu_eres.png");

  TCanvas* cphad = new TCanvas("cphad","",1800,600);
  cphad->Divide(4,1);
  cphad->cd(1);
  hphad->Draw("colz");
  cphad->cd(2);
  hphad_1p->Draw("colz");
  cphad->cd(3);
  hphad_1p_true->Draw("ehist");
  hphad_1p_reco->Draw("ehistsame");
  hphad_true->Draw("ehistsame");
  hphad_reco->Draw("ehistsame");
  legmu->Draw();
  cphad->cd(4);
  hphad_1p_diff->Draw("ehist");
  hphad_diff->Draw("ehistsame");
  gPad->SetLogy();
  cphad->SaveAs("had_mom.png");

  TCanvas* ccoshad = new TCanvas("ccoshad","",1800,600);
  ccoshad->Divide(4,1);
  ccoshad->cd(1);
  hcoshad->Draw("colz");
  ccoshad->cd(2);
  hcoshad_1p->Draw("colz");
  ccoshad->cd(3);
  hcoshad_1p_true->Draw("ehist");
  hcoshad_1p_reco->Draw("ehistsame");
  hcoshad_true->Draw("ehistsame");
  hcoshad_reco->Draw("ehistsame");

  TLegend* legmuclone = (TLegend*)legmu->Clone();
  legmuclone->SetX1NDC(0.15);
  legmuclone->SetX2NDC(0.45);
  legmuclone->Draw();

  ccoshad->cd(4);
  hcoshad_1p_diff->Draw("ehist");
  hcoshad_diff->Draw("ehistsame");
  ccoshad->SaveAs("had_cos.png");

  TCanvas* cwhad = new TCanvas("cwhad","",1800,600);
  cwhad->Divide(4,1);
  cwhad->cd(1);
  hwhad->Draw("colz");
  cwhad->cd(2);
  hwhad_1p->Draw("colz");
  cwhad->cd(3);
  hwhad_1p_true->Draw("ehist");
  hwhad_1p_reco->Draw("ehistsame");
  hwhad_true->Draw("ehistsame");
  hwhad_reco->Draw("ehistsame");
  legmu->Draw();
  cwhad->cd(4);
  hwhad_1p_diff->Draw("ehist");
  hwhad_diff->Draw("ehistsame");
  cwhad->SaveAs("had_w.png");

  TCanvas* cwvshad = new TCanvas("cwvshad","",1400,1200);
  cwvshad->Divide(2,2);
  cwvshad->cd(1);
  hwphad_true->Draw("colz");
  cwvshad->cd(2);
  hwphad_reco->Draw("colz");
  cwvshad->cd(3);
  hwcoshad_true->Draw("colz");
  cwvshad->cd(4);
  hwcoshad_reco->Draw("colz");
  cwvshad->SaveAs("had_w_vs.png");

  TCanvas* cwvshad_1p = new TCanvas("cwvshad_1p","",1400,1200);
  cwvshad_1p->Divide(2,2);
  cwvshad_1p->cd(1);
  hwphad_1p_true->Draw("colz");
  cwvshad_1p->cd(2);
  hwphad_1p_reco->Draw("colz");
  cwvshad_1p->cd(3);
  hwcoshad_1p_true->Draw("colz");
  cwvshad_1p->cd(4);
  hwcoshad_1p_reco->Draw("colz");
  cwvshad_1p->SaveAs("had_w_vs_1p.png");

}

void kinrecoana::ProtonMult() {

  this->SetList(false,"mc_is_signal"); //sel_CCNp0pi
  fChain->SetMakeClass(0);


  TH2F* h = new TH2F("h","proton multiplicity;N_{true};N_{reco}",21,-0.5,20.5,21,-0.5,20.5);
  fChain->Draw("np:mc_np>>h");

  h->GetXaxis()->SetNdivisions(211);
  h->GetYaxis()->SetNdivisions(211);
  auto htrue = h->ProjectionX("htrue");
  auto hreco = h->ProjectionY("hreco");

  hreco->GetYaxis()->SetTitle("");
  hreco->GetXaxis()->SetTitle("N_{p}");
  hreco->SetLineColor(kBlue);
  hreco->SetLineWidth(2);
  htrue->SetLineColor(kGreen-2);
  htrue->SetLineWidth(2);

  gStyle->SetOptStat(0);

  auto c = new TCanvas();
  h->Draw("colz");
  h->GetYaxis()->SetNdivisions(208);
  h->GetYaxis()->SetRangeUser(-0.5,8.5);
  c->SetLogz();
  c->SaveAs("proton_mult_2d_ntrue_nreco.png");

  auto c1 = new TCanvas();
  hreco->Draw("ehist");
  htrue->Draw("sameehist");
  c1->SetLogy();

  TLegend* leg = new TLegend(0.7,0.7,0.88,0.88);
  leg->SetBorderSize(0);
  leg->AddEntry("htrue","true","el");
  leg->AddEntry("hreco","reco","el");
  leg->Draw();

  c1->SaveAs("proton_mult_1d_ntrue_nreco.png");
}
