#include "core/test.hpp"

#include <algorithm>
#include <bitset>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include <cmath>
#include "Math/Vector4D.h"
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/xingchen.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  //   Process::Type data =  Process::Type::data;

  //ROOT file path
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/");
  string sfolder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/");
  string mc_path("/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/mc/");
  string sig_path("/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/signal/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  //NamedFunc trigs(el_trigs || mu_trigs);


  auto proc_smzg  = Process::MakeShared<Baby_nano>("SM Z#gamma",back,TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, "1");
  //auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back,
  //                                               TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_nano>("DY",               back,
                                                   TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  "1");
//auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back,
  //                                                 TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  auto proc_ttg   = Process::MakeShared<Baby_nano>("ttbar",            back,
                                                 TColor::GetColor("#ED702D"),{mc_path+"*TTJets*"}, "1");
  //auto proc_WJtoLNu  = Process::MakeShared<Baby_nano>("WJets to LNu",back, kBlue,{mc_path+"WJetsToLNu*"}, "1");
  auto proc_ewk  = Process::MakeShared<Baby_nano>("EWK(WW,WZ,ZZ)",back, kRed,{mc_path+"WW*", mc_path+"WZ*",mc_path+"ZZ*"}, "1");
  auto proc_hzg   = Process::MakeShared<Baby_nano>("HToZ#gamma", sig,
                                                   kRed     ,{sig_path+"*.root"},   "1");
  auto proc_hzg_vbf = Process::MakeShared<Baby_nano>("HToZ#gamma VBF", sig,
                          kMagenta                   ,{sig_path+"VBF*.root"},   "1");
  auto proc_hzg_gg  = Process::MakeShared<Baby_nano>("HToZ#gamma ggF", sig,
                                                     kBlue ,{sig_path+"GluGlu*.root"},   "1");
  auto proc_hzg_tth  = Process::MakeShared<Baby_nano>("HToZ#gamma ttH", sig,
                                                     kGreen ,{sig_path+"ttH*.root"},   "1");

  auto proc_hzg_vh =  Process::MakeShared<Baby_nano>("HToZ#gamma VH", sig,
                                                     kBlack                  ,{sig_path+"*WplusH*.root",sig_path+"*WminusH*.root",sig_path+"*ZH*.root"},   "1");
  proc_smzg->SetLineWidth (1); proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); proc_ewk->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  //proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(3);
  proc_hzg_vbf->SetLineWidth(3);
  proc_hzg_tth->SetLineWidth(3);
  //vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg_vbf, proc_hzg_gg };
  //vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_gg};
  vector<shared_ptr<Process>> sigprocs = {proc_hzg_vbf, proc_hzg_gg, proc_hzg_tth, proc_hzg_vh };
  vector<shared_ptr<Process>> bkgprocs = { proc_ttg,  proc_ewk};

  /*  const NamedFunc mu_sig("mu_sig",[](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> mu_sig_;
      for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
	bool is_sig = true;
	if (b.Muon_pt()->at(imu) <= 5.0) is_sig = false;
	if (abs(b.Muon_eta()->at(imu)) > 2.4) is_sig = false;
	if (abs(b.Muon_dz()->at(imu))>1.0)  is_sig = false;
	if (abs(b.Muon_dxy()->at(imu))>0.5) is_sig = false; 
	if (!((b.Muon_looseId()->at(imu) || (b.Muon_pt()->at(imu) > 200 && b.Muon_highPtId()->at(imu))) && 
	      b.Muon_pfRelIso03_all()->at(imu) < 0.35 &&
	      b.Muon_sip3d()->at(imu) < 4)) is_sig = false;
	if (is_sig)
	  mu_sig_.push_back(1.0);
	else
	  mu_sig_.push_back(0.0);
      }
      return mu_sig_;
      });*/
  const NamedFunc el_sig("el_sig",[](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> el_sig_;
      for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
	bool is_sig = true;
	if (b.Electron_pt()->at(iel) <= 7.0) is_sig = false;
	float etasc = b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel);
	if (abs(etasc) > 2.5) is_sig = false;
	if (abs(b.Electron_dz()->at(iel))>1.0) is_sig = false;
	if (abs(b.Electron_dxy()->at(iel))>0.5) is_sig = false; 
	/*double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
			   { 0.604775,  0.628743,   0.896462}};
	int ipt(1), ieta(2);
	if(b.Electron_pt()->at(iel)>10.0) ipt = 0;
	if(fabs(etasc) < 0.8) ieta = 0;
	else if(fabs(etasc) < 1.479) ieta = 1;*/
	//double mva = b.Electron_mvaFall17V2Iso()->at(iel);
	bool mva = b.Electron_mvaFall17V2Iso_WPL()->at(iel);
	if (!mva)//wp[ipt][ieta])
	  is_sig = false;
	if (is_sig)
	  el_sig_.push_back(1.0);
	else
	  el_sig_.push_back(0.0);
      }
      return el_sig_;
    });
  const NamedFunc z_decay_pdgid("z_decay_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
      for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
	int abs_pdgid = abs(b.GenPart_pdgId()->at(imc));
	if (abs_pdgid==11 || abs_pdgid == 13 || abs_pdgid == 15) {
	  int mom_idx = b.GenPart_genPartIdxMother()->at(imc);
	  if (mom_idx != -1) {
	    if ((b.GenPart_pdgId()->at(mom_idx))==23) {
	      return abs_pdgid;
	    }
	  }
	}
      }
      return 0;
    });
  const NamedFunc offline_nmu("offline_nmu",[/*mu_sig*/](const Baby &b) -> NamedFunc::ScalarType{
      int nmu_sig = 0;
      std::vector<double> mu_sig_ = mu_sig.GetVector(b);
      for (double imu_sig_ : mu_sig_) {
	if (imu_sig_ > 0.5)
	  nmu_sig++;
      }
      return nmu_sig;
    });
  const NamedFunc offline_nel("offline_nel",[el_sig](const Baby &b) -> NamedFunc::ScalarType{
      int nel_sig = 0;
      std::vector<double> el_sig_ = el_sig.GetVector(b);
      for (double iel_sig_ : el_sig_) {
	if (iel_sig_ > 0.5)
	  nel_sig++;
      }
      return nel_sig;
    });
 
  /*  const NamedFunc ph_sig("ph_sig",[](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> ph_sig_;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
	if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15) ph_sig_.push_back(1.0);
	else ph_sig_.push_back(0.0);
      }
      return ph_sig_;
      });*/

  const NamedFunc ph_sig("ph_sig",[](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> ph_sig_;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if (b.Photon_pt()->at(iph) <= 15. || abs(b.Photon_eta()->at(iph)) > 2.5 || !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) ph_sig_.push_back(0.0);
        //else if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0) ph_sig_.push_back(1.0);
	else if (((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58)) && b.Photon_electronVeto()->at(iph)) ph_sig_.push_back(1.0);
	else ph_sig_.push_back(0.0);
      }
      return ph_sig_;
    });

 
 NamedFunc ele("nElectron >= 2 && nPhoton >= 1");
 NamedFunc muon("nMuon >= 2 && nPhoton >= 1");
 NamedFunc muon_mll("muon_mll", [](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector m1, m2, mm;
     m1.SetPtEtaPhiM(b.Muon_pt()->at(0), b.Muon_eta()->at(0), b.Muon_phi()->at(0), b.Muon_mass()->at(0));
     m2.SetPtEtaPhiM(b.Muon_pt()->at(1), b.Muon_eta()->at(1), b.Muon_phi()->at(1), b.Muon_mass()->at(1));
     mm = m1 + m2;
     mass = mm.M();
     return mass>50;
   });

 NamedFunc ele_mll("ele_mll", [](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector e1, e2, ee;
     e1.SetPtEtaPhiM(b.Electron_pt()->at(0), b.Electron_eta()->at(0), b.Electron_phi()->at(0), b.Electron_mass()->at(0));
     e2.SetPtEtaPhiM(b.Electron_pt()->at(1), b.Electron_eta()->at(1), b.Electron_phi()->at(1), b.Electron_mass()->at(1));
     ee = e1 + e2;
     mass = ee.M();
     return mass>50;
   });

 NamedFunc off_muon_mll("off_muon_mll", [/*mu_sig*/](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector m1, m2, mm;
     int nmu = 0;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
	 m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	 nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
	 m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	 break;
       }
     }
     mm = m1 + m2;
     mass = mm.M();
     return mass>50;
   });
 NamedFunc off_muon_mll_cut("off_muon_mll_cut", [/*mu_sig*/](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector m1, m2, mm;
     int nmu = 0;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     mm = m1 + m2;
     mass = mm.M();
     return mass>81.2 && mass<101.2;
   });

 NamedFunc off_ele_mll("off_ele_mll", [el_sig](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector e1, e2, ee;
     int nel = 0;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
	 e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	 nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
	 e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	 break;
       }
     }
     ee = e1 + e2;
     mass = ee.M();
     return mass>50;
   });

 NamedFunc off_ele_mll_cut("off_ele_mll_cut", [el_sig](const Baby &b) -> NamedFunc::ScalarType{
     float mass;
     TLorentzVector e1, e2, ee;
     int nel = 0;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     ee = e1 + e2;
     mass = ee.M();
     return mass>81.2 && mass<101.2;
   });

 NamedFunc off_ele_dR_list("off_ele_dR_list", [el_sig, ph_sig](const Baby &b) -> NamedFunc::VectorType{
     TLorentzVector e, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     std::vector<double> dRlist;
     double mindr = 9999.0;
     for (unsigned j(0); j<ph_sig_.size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else{
         dRlist.push_back(0.0);
         continue;
       }
       for (unsigned i(0); i<el_sig_.size();i++) {
         if (el_sig_[i]>0.5)  e.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         else continue;
         double tempdr = e.DeltaR(p);
         if (tempdr<mindr) mindr = tempdr;
       }
       if (mindr>0.4 && mindr<9000.) dRlist.push_back(1.0);
       else dRlist.push_back(0.0);
     }
     return dRlist;
   });

 NamedFunc off_muon_dR_list("off_muon_dR_list", [/*mu_sig, */ph_sig](const Baby &b) -> NamedFunc::VectorType{
     TLorentzVector m, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     std::vector<double> dRlist;
     double mindr = 9999.0;
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else {
         dRlist.push_back(0.0);
         continue;
       }
       for (unsigned i(0); i<mu_sig_.size();i++) {
         if (mu_sig_[i]>0.5)  m.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         else continue;
         double tempdr = m.DeltaR(p);
         if (tempdr<mindr) mindr = tempdr;
       }
       if (mindr>0.4 && mindr<9000.) dRlist.push_back(1.0);
       else dRlist.push_back(0.0);
     }
     return dRlist;
   });

 const NamedFunc offline_nph("offline_nph",[off_muon_dR_list, off_ele_dR_list](const Baby &b) -> NamedFunc::ScalarType{
     int nph_sig = 0;
     std::vector<double> mudR = off_muon_dR_list.GetVector(b);
     std::vector<double> eledR = off_ele_dR_list.GetVector(b);
     for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
       if (b.Photon_pt()->at(iph) <= 15.0) continue;
       if (abs(b.Photon_eta()->at(iph)) > 2.5) continue;
       if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) continue;
       // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
       //if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0)
       if (((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58)) && b.Photon_electronVeto()->at(iph) && (mudR[iph]>0.5 || eledR[iph]>0.5))
	 nph_sig++;
     }
     return nph_sig;
   });

 /*NamedFunc off_ele_dR("off_ele_dR", [el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector e, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     bool cut = false;
     double mindr = 9999.0; 
     for (unsigned j(0); j<ph_sig_.size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       for (unsigned i(0); i<el_sig_.size();i++) {
	 if (el_sig_[i]>0.5)  e.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	 else continue;
	 double tempdr = e.DeltaR(p);
	 if (tempdr<mindr) mindr = tempdr;
       }
       if (mindr>0.4 && mindr<9000.){
	 cut = true;
	 break;
       }
     }
     return cut;
     });*/

 NamedFunc off_ele_sumcut("off_ele_sumcut", [el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector e1, e2, p, mll, mllg;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nel = 0;
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4){
	 mll = e1 + e2;
	 mllg = e1 + e2 + p;
	 return (mll.M() + mllg.M())>185; 
       }
     }
     return false;
   });

 /* NamedFunc off_muon_dR("off_muon_dR", [mu_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector m, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     bool cut = false;
     double mindr = 9999.0;
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       for (unsigned i(0); i<mu_sig_.size();i++) {
         if (mu_sig_[i]>0.5)  m.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         else continue;
         double tempdr = m.DeltaR(p);
         if (tempdr<mindr) mindr = tempdr;
       }
       if (mindr>0.4 && mindr<9000.){
	 cut = true;
	 break;
       }     
     }
     return cut;
     });*/

 NamedFunc off_muon_sumcut("off_muon_sumcut", [/*mu_sig, */ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector m1, m2, p, mll, mllg;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nmu = 0;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4) {
	 mll = m1 + m2;
	 mllg = m1 + m2 +p;
	 return (mllg.M() + mll.M())>185;
       }
     }
     return false;
   });

 NamedFunc off_muon_mll_mass("off_muon_mll_mass", [/*mu_sig*/](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector m1, m2, mll;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     double closemass = 0.0;
     double closediff = 999.;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       for (unsigned j(i+1); j<mu_sig_.size();j++) {
	 if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
	   m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	   m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
	   mll = m1 + m2;	 
	   double massdiff = abs(mll.M()-91.1876);
	   if (massdiff < closediff) {
	     closemass = mll.M();
	     closediff = massdiff;
	   }
	 }
	 else continue;
       }
     }
     return closemass;
   });

 NamedFunc off_muon_mllg_mass("off_muon_mllg_mass", [/*mu_sig, */ph_sig,off_muon_dR_list](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector m, m1, m2, mll, m1c, m2c, p, mllg;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     std::vector<double> mudR = off_muon_dR_list.GetVector(b);
     //std::vector<double> eledR = off_ele_dR_list.GetVector(b);
     double closediff = 999.;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       for (unsigned j(i+1); j<mu_sig_.size();j++) {
         if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
           m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
           m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
           mll = m1 + m2;
           double massdiff = abs(mll.M()-91.1876);
           if (massdiff < closediff) {
             closediff = massdiff;
	     m1c = m1;
	     m2c = m2;
           }
         }
         else continue;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5 && mudR[j]>0.5){
	 p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
	 mllg = m1c + m2c + p;
	 break;
       }
     }
     return mllg.M();
   });

 NamedFunc off_ele_mll_mass("off_ele_mll_mass", [el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector e1, e2, mll;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     double closemass = 0.0;
     double closediff = 999.;
     for (unsigned i(0); i<el_sig_.size();i++) {
       for (unsigned j(i+1); j<el_sig_.size();j++) {      
	 if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
	   e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	   e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
	   mll = e1 + e2;
	   double massdiff = abs(mll.M()-91.1876);
	   if (massdiff < closediff){
	     closemass = mll.M();
	     closediff = massdiff;
	   }
	 }
	 else continue;
       }
     }
     return closemass;
   });

 NamedFunc off_ele_mllg_mass("off_ele_mllg_mass", [el_sig, ph_sig, off_ele_dR_list](const Baby &b) -> NamedFunc::ScalarType{
     TLorentzVector e, e1, e2,mll,e1c,e2c, p, mllg;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     std::vector<double> eledR = off_ele_dR_list.GetVector(b);
     double closediff = 999.;
     for (unsigned i(0); i<el_sig_.size();i++) {
       for (unsigned j(i+1); j<el_sig_.size();j++) {
         if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
           e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
           e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
           mll = e1 + e2;
           double massdiff = abs(mll.M()-91.1876);
           if (massdiff < closediff){
             closediff = massdiff;
	     e1c = e1;
	     e2c = e2;
           }
         }
         else continue;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5 && eledR[j]>0.5){
         p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
         mllg = e1c + e2c + p;
         break;
       }
     }
     return mllg.M();
   });

 NamedFunc ele_pdgid("ele_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector e1, e2, p;
     e1.SetPtEtaPhiM(b.Electron_pt()->at(0), b.Electron_eta()->at(0), b.Electron_phi()->at(0), b.Electron_mass()->at(0));
     e2.SetPtEtaPhiM(b.Electron_pt()->at(1), b.Electron_eta()->at(1), b.Electron_phi()->at(1), b.Electron_mass()->at(1));
     for (unsigned i(0); i<b.Photon_pt()->size(); i++){
       p.SetPtEtaPhiM(b.Photon_pt()->at(i), b.Photon_eta()->at(i), b.Photon_phi()->at(i), b.Photon_mass()->at(i));
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4){
	 index = i;
	 break;
       }
     }
     //cout <<"Photon Index: " << index <<endl;
     //cout << "GenSize: " << b.GenPart_pdgId()->size() << " GenIndex: " << b.Photon_genPartIdx()->at(index) << endl; 
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_pdgId()->at(b.Photon_genPartIdx()->at(index));
   });

 NamedFunc muon_pdgid("muon_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector m1, m2, p;
     m1.SetPtEtaPhiM(b.Muon_pt()->at(0), b.Muon_eta()->at(0), b.Muon_phi()->at(0), b.Muon_mass()->at(0));
     m2.SetPtEtaPhiM(b.Muon_pt()->at(1), b.Muon_eta()->at(1), b.Muon_phi()->at(1), b.Muon_mass()->at(1));
     for (unsigned i(0); i<b.Photon_pt()->size(); i++){
       p.SetPtEtaPhiM(b.Photon_pt()->at(i), b.Photon_eta()->at(i), b.Photon_phi()->at(i), b.Photon_mass()->at(i));
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4){
         index = i;
         break;
       }
     }
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_pdgId()->at(b.Photon_genPartIdx()->at(index));
   });

 NamedFunc off_ele_pdgid("off_ele_pdgid",[el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector e1, e2, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nel = 0;
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4) {
	 index = j;
	 break;
       }
     }
     //cout <<"Photon Index: " << index <<endl;
     //cout << "GenSize: " << b.GenPart_pdgId()->size() << " GenIndex: " << b.Photon_genPartIdx()->at(index) << endl;
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_pdgId()->at(b.Photon_genPartIdx()->at(index));
   });

 NamedFunc off_muon_pdgid("off_muon_pdgid",[/*mu_sig, */ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector m1, m2, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nmu = 0;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4){
	 index = j;
	 break;
       }
     }
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_pdgId()->at(b.Photon_genPartIdx()->at(index));
   });

 NamedFunc off_ele_motherpdgid("off_ele_motherpdgid",[el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector e1, e2, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nel = 0;
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4) {
         index = j;
         break;
       }
     }
     //cout <<"Photon Index: " << index <<endl;
     //cout << "GenSize: " << b.GenPart_pdgId()->size() << " GenIndex: " << b.Photon_genPartIdx()->at(index) << endl;
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else if (b.GenPart_genPartIdxMother()->at(b.Photon_genPartIdx()->at(index))<0) return -2;
     else return abs(b.GenPart_pdgId()->at(b.GenPart_genPartIdxMother()->at(b.Photon_genPartIdx()->at(index))));
   });

 NamedFunc off_muon_motherpdgid("off_muon_motherpdgid",[/*mu_sig, */ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector m1, m2, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nmu = 0;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4){
         index = j;
         break;
       }
     }
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else if (b.GenPart_genPartIdxMother()->at(b.Photon_genPartIdx()->at(index))<0) return -2;
     else return abs(b.GenPart_pdgId()->at(b.GenPart_genPartIdxMother()->at(b.Photon_genPartIdx()->at(index))));
   });

 NamedFunc off_ele_genstatus("off_ele_genstatus",[el_sig, ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector e1, e2, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nel = 0;
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4) {
         index = j;
         break;
       }
     }
     //cout <<"Photon Index: " << index <<endl;
     //cout << "GenSize: " << b.GenPart_pdgId()->size() << " GenIndex: " << b.Photon_genPartIdx()->at(index) << endl;
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_status()->at(b.Photon_genPartIdx()->at(index));
   });

 NamedFunc off_muon_genstatus("off_muon_genstatus",[/*mu_sig, */ph_sig](const Baby &b) -> NamedFunc::ScalarType{
     int index = 0;
     TLorentzVector m1, m2, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nmu = 0;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4){
         index = j;
         break;
       }
     }
     if ( b.Photon_genPartIdx()->at(index)<0) return -1;
     else return b.GenPart_status()->at(b.Photon_genPartIdx()->at(index));
   });

 const NamedFunc off_ele_flags("off_ele_flags",[el_sig, ph_sig](const Baby &b) -> NamedFunc::VectorType{
     std::vector<double> flag;
     std::bitset<15> bit;
     int index = 0;
     TLorentzVector e1, e2, p;
     std::vector<double> el_sig_ = el_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nel = 0;
     for (unsigned i(0); i<el_sig_.size();i++) {
       if (el_sig_[i] > 0.5 && nel==0){
         e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         nel += 1;
       }
       else if (el_sig_[i] > 0.5 && nel==1) {
         e2.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (e1.DeltaR(p)>0.4 && e2.DeltaR(p)>0.4) {
         index = j;
         break;
       }
     }
     //cout <<"Photon Index: " << index <<endl;
     //cout << "GenSize: " << b.GenPart_pdgId()->size() << " GenIndex: " << b.Photon_genPartIdx()->at(index) << endl;
     if ( b.Photon_genPartIdx()->at(index)<0) return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     else {
       bit = std::bitset<15>(b.GenPart_statusFlags()->at(b.Photon_genPartIdx()->at(index)));
       for(int i(0); i<15; i++){
	 if (bit[i]==1) flag.push_back(1.0);
	 else flag.push_back(0.0);
       }
       return flag;
     }
   });

 const NamedFunc off_muon_flags("off_muon_flags",[/*mu_sig,*/ ph_sig](const Baby &b) -> NamedFunc::VectorType{
     std::vector<double> flag;
     std::bitset<15> bit;
     int index = 0;
     TLorentzVector m1, m2, p;
     std::vector<double> mu_sig_ = mu_sig.GetVector(b);
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     int nmu = 0;
     for (unsigned i(0); i<mu_sig_.size();i++) {
       if (mu_sig_[i] > 0.5 && nmu==0){
         m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         nmu += 1;
       }
       else if (mu_sig_[i] > 0.5 && nmu==1) {
         m2.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
         break;
       }
     }
     for (unsigned j(0); j<b.Photon_pt()->size(); j++){
       if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
       else continue;
       if (m1.DeltaR(p)>0.4 && m2.DeltaR(p)>0.4){
         index = j;
         break;
       }
     }
     int gind;
     gind =  b.Photon_genPartIdx()->at(index);
     if (gind<0) return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     else {
       bit = std::bitset<15>(b.GenPart_statusFlags()->at(gind));
       for(int i(0); i<15; i++){
         if (bit[i] == 1) flag.push_back(1.0);
	 else flag.push_back(0.0);
       }
       return flag;
     }
   });

 const NamedFunc off_muon_flags0("off_muon_flags0",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[0] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags1("off_muon_flags1",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[1] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags2("off_muon_flags2",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[2] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags3("off_muon_flags3",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[3] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags4("off_muon_flags4",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[4] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags5("off_muon_flags5",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[5] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags6("off_muon_flags6",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[6] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags7("off_muon_flags7",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[7] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags8("off_muon_flags8",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[8] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags9("off_muon_flags9",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[9] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags10("off_muon_flags10",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[10] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags11("off_muon_flags11",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[11] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags12("off_muon_flags12",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[12] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags13("off_muon_flags13",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[13] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_muon_flags14("off_muon_flags14",[off_muon_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_muon_flags.GetVector(b);
     if (flag_[14] > 0.5) return true;
     else return false;
   });


 const NamedFunc off_ele_flags0("off_ele_flags0",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[0] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags1("off_ele_flags1",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[1] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags2("off_ele_flags2",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[2] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags3("off_ele_flags3",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[3] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags4("off_ele_flags4",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[4] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags5("off_ele_flags5",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[5] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags6("off_ele_flags6",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[6] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags7("off_ele_flags7",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[7] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags8("off_ele_flags8",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[8] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags9("off_ele_flags9",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[9] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags10("off_ele_flags10",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[10] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags11("off_ele_flas11",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[11] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags12("off_ele_flags12",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[12] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags13("off_ele_flags13",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[13] > 0.5) return true;
     else return false;
   });
 const NamedFunc off_ele_flags14("off_ele_flags14",[off_ele_flags](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> flag_ = off_ele_flags.GetVector(b);
     if (flag_[14] > 0.5) return true;
     else return false;
   });


 NamedFunc stitch_cut ("stitch_cut", [ph_sig, off_muon_dR_list, off_ele_dR_list](const Baby &b) -> NamedFunc::ScalarType{
     std::vector<double> ph_sig_ = ph_sig.GetVector(b);
     std::vector<double> eledR = off_ele_dR_list.GetVector(b);
     std::vector<double> muondR = off_muon_dR_list.GetVector(b);
     for (int imc(0); imc<b.nGenPart(); ++imc) {
       bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
       int mc_id = b.GenPart_pdgId()->at(imc);
       if(mc_id == 22 && b.GenPart_status()->at(imc) == 1)
	 if(mc_statusFlags[0] || mc_statusFlags[8])
	   for (int igamma(0); igamma<b.nPhoton(); ++igamma)
	     if(ph_sig_[igamma] > 0.5 && (eledR[igamma]>0.5 || muondR[igamma]>0.5)){
	       TVector3 genPhoton, recoPhoton;
	     genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc), 
			     b.GenPart_eta()->at(imc), 
			     b.GenPart_phi()->at(imc));
	     recoPhoton.SetPtEtaPhi(b.Photon_pt()->at(igamma), 
			      b.Photon_eta()->at(igamma), 
			      b.Photon_phi()->at(igamma));
	     if(genPhoton.DeltaR(recoPhoton) < 0.1)
	       return  false;
	     }	 
     }
     return true;
   });
 NamedFunc stitch_cut_more ("stitch_cut_more", [](const Baby &b) -> NamedFunc::ScalarType{
     for (int imc(0); imc<b.nGenPart(); ++imc) {
       bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
       int mc_id = b.GenPart_pdgId()->at(imc);
       if(mc_id == 22 && b.GenPart_status()->at(imc) == 1)
         if(mc_statusFlags[0] || mc_statusFlags[8]) return  false;
     }
     return true;
   });
 auto proc_dy_nos= Process::MakeShared<Baby_nano>("DY stitch",back,TColor::GetColor("#ffb400"),{mc_path+"*DYJetsToLL_M-50*049988F3*"}, "1");
 bkgprocs = {proc_smzg, proc_dy_nos};
 vector<shared_ptr<Process>> procs = {proc_smzg, proc_dy_nos,proc_hzg_vbf, proc_hzg_gg, proc_hzg_tth, proc_hzg_vh }; 
 proc_dy_nos->SetLineWidth(1);
 vector<shared_ptr<Process>> stitch = {proc_dy_nos};

 PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
 log_lumi.Title(TitleType::info)
      .YAxis(YAxisType::log)
      .Stack(StackType::shapes)
      //           .Overflow(OverflowType::none)
      .YTitleOffset(1.)
      .AutoYAxis(false)
      .UseCMYK(false)
      .LeftMargin(0.17)
      .LegendColumns(1)
      //           .CanvasWidth(1077)
      .FileExtensions({"pdf"});

    //Plot options (linear/log stack or not)
 PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
 PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
 PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);

 vector<PlotOpt> ops1 = {log_stack};
 vector<PlotOpt> ops2 = {lin_stack};
vector<PlotOpt> ops3 = {lin_lumi, lin_stack};
    //vector<NamedFunc> loose_lep = {"ll_lepid[0] == 11","ll_lepid[0] == 13"};
    //NamedFunc baseline("nll >= 1 && nphoton > 0 && ll_m[0] > 50");// && ll_m[0] > 50");

    //Weight
 const NamedFunc weight("weight",[](const Baby &b) -> NamedFunc::ScalarType{
     //extra factor for weighting is xs/nevts_effective where effective events are calculated including negative weights
     //0.007519850838359999 GluGluH->ZG->llG xs (pb)
     //967054
     float w_lumi = 1.0;
     float w_year = 1.0;
     if (b.Generator_weight()<0) w_lumi = -1.0;
     if (b.SampleType()==2016) w_year = 36.32264;
     else if (b.SampleType()==2017) w_year = 41.52756;
     else if (b.SampleType()==2018) w_year = 59.67377;
     float corr_factor = 112.0/104.0; //temp, just to match scale with AN
     return w_lumi*w_year*7.7760403e-09*1000.0*corr_factor;
    });

 const NamedFunc weight_1("weight_1",[](const Baby &b) -> NamedFunc::ScalarType{
     float w = 1.0;
     if (b.Generator_weight()<0) return w;
     else return w;
   });

 //Set up plot maker
 string tag;
 PlotOpt style("txt/plot_styles.txt","Scatter");
 vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
			     .LogMinimum(1)
			     .CanvasWidth(600)
			     .LabelSize(0.04)
			     .YAxis(YAxisType::log)
			     .Title(TitleType::info)};
 PlotMaker pm;
 //vector<PlotOpt> ops;
 tag = "Nano_";
   // ele && ele_mll && ele_dR
   //pm.Push<Hist1D>(Axis(41,-1,40, off_ele_pdgid,    "Photon PDG ID"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_");
   //pm.Push<Hist1D>(Axis(41,-1,40, off_muon_pdgid,    "Photon PDG ID"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_");
   
 //pm.Push<Hist1D>(Axis(42,-2,40, off_ele_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR, bkgprocs, ops2).Weight(weight).Tag(tag);
 //pm.Push<Hist1D>(Axis(42,-2,40, off_muon_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR, bkgprocs, ops2).Weight(weight).Tag(tag);
 //pm.Push<Hist1D>(Axis(42,-2,40, off_ele_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR &&  off_ele_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag+"_sumcut_");
 //pm.Push<Hist1D>(Axis(42,-2,40, off_muon_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR &&  off_muon_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag+"_sumcut_");
 
// pm.Push<Hist1D>(Axis(42,-2,40, off_ele_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_");
 //pm.Push<Hist1D>(Axis(42,-2,40, off_muon_motherpdgid,    "Mother PDG ID"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_");
   //pm.Push<Hist1D>(Axis(100,0,350, "MET_pt",    "MET_pt"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_mll_cut, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_mcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "MET_pt",    "MET_pt"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_mll_cut, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_mcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "MET_pt",    "MET_pt"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_sumcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "MET_pt",    "MET_pt"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_sumcut");
   //pm.Push<Hist1D>(Axis(80,0,80, "Pileup_nPU",    "Pileup"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_mll_cut, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_mcut");
   //pm.Push<Hist1D>(Axis(80,0,80, "Pileup_nPU",    "Pileup"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_mll_cut, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_mcut");
   //pm.Push<Hist1D>(Axis(80,0,80, "Pileup_nPU",    "Pileup"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_sumcut");
   //pm.Push<Hist1D>(Axis(80,0,80, "Pileup_nPU",    "Pileup"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_sumcut");
   
   //pm.Push<Hist1D>(Axis(100,0,350, "GenMET_pt",    "GenMET_pt"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_mll_cut, procs, ops2).Weight(weight).Tag(tag + "ele_mcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "GenMET_pt",    "GenMET_pt"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_mll_cut, procs, ops2).Weight(weight).Tag(tag + "muon_mcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "GenMET_pt",    "GenMET_pt"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut, procs, ops2).Weight(weight).Tag(tag + "ele_sumcut");
   //pm.Push<Hist1D>(Axis(100,0,350, "GenMET_pt",    "GenMET_pt"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut, procs, ops2).Weight(weight).Tag(tag + "muon_sumcut");
   //pm.Push<Hist1D>(Axis(101,-1,100, off_ele_genstatus,    "Generator particle status"      ,{}), (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "ele_bkg_");
   //pm.Push<Hist1D>(Axis(101,-1,100, off_muon_genstatus,    "Generator particle status"      ,{}), (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR, bkgprocs, ops2).Weight(weight).Tag(tag + "muon_bkg_");
   //pm.Push<Hist1D>(Axis(10,0,10, nE,    "nE"      ,{}),  ele && ele_mll && ele_dR, sigprocs, ops2).Weight(weight).Tag(tag + "test_");
 /*
 pm.Push<Hist2D>(Axis(80,50,130,  off_muon_mll_mass,"m_{#mu#mu} [GeV]",{}), Axis(200,50,250,  off_muon_mllg_mass, "m_{#mu#mu#gamma} [GeV]",{100,80}),
 		 (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == 13), stitch, bkg_hist).Tag(tag + "mu_dy_l_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_ele_mll_mass,"m_{ee} [GeV]",{}), Axis(200,50,250,  off_ele_mllg_mass, "m_{ee#gamma} [GeV]",{100,80}),
 		 (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == 11), stitch, bkg_hist).Tag(tag + "el_dy_l_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_muon_mll_mass,"m_{#mu#mu} [GeV]",{}), Axis(200,50,250,  off_muon_mllg_mass, "m_{#mu#mu#gamma} [GeV]",{100,80}),
                 (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == -2), stitch, bkg_hist).Tag(tag + "mu_dy_nomom_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_ele_mll_mass,"m_{ee} [GeV]",{}), Axis(200,50,250,  off_ele_mllg_mass, "m_{ee#gamma} [GeV]",{100,80}),
                 (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == -2), stitch, bkg_hist).Tag(tag + "el_dy_nomom_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_muon_mll_mass,"m_{#mu#mu} [GeV]",{}), Axis(200,50,250,  off_muon_mllg_mass, "m_{#mu#mu#gamma} [GeV]",{100,80}),
                 (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == -1), stitch, bkg_hist).Tag(tag + "mu_dy_nop_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_ele_mll_mass,"m_{ee} [GeV]",{}), Axis(200,50,250,  off_ele_mllg_mass, "m_{ee#gamma} [GeV]",{100,80}),
                 (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == -1), stitch, bkg_hist).Tag(tag + "el_dy_nop_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_muon_mll_mass,"m_{#mu#mu} [GeV]",{}), Axis(200,50,250,  off_muon_mllg_mass, "m_{#mu#mu#gamma} [GeV]",{100,80}),
                 (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid >0.5 && off_muon_motherpdgid <9), stitch, bkg_hist).Tag(tag + "mu_dy_quark_").Weight(weight);
 pm.Push<Hist2D>(Axis(80,50,130,  off_ele_mll_mass,"m_{ee} [GeV]",{}), Axis(200,50,250,  off_ele_mllg_mass, "m_{ee#gamma} [GeV]",{100,80}),
                 (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_muon_motherpdgid >0.5 && off_muon_motherpdgid <9), stitch, bkg_hist).Tag(tag + "el_dy_quark_").Weight(weight);
pm.Push<Hist2D>(Axis(80,50,130,  off_muon_mll_mass,"m_{#mu#mu} [GeV]",{}), Axis(200,50,250,  off_muon_mllg_mass, "m_{#mu#mu#gamma} [GeV]",{100,80}),
		(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid >=40), stitch, bkg_hist).Tag(tag + "mu_dy_other_").Weight(weight);
pm.Push<Hist2D>(Axis(80,50,130,  off_ele_mll_mass,"m_{ee} [GeV]",{}), Axis(200,50,250,  off_ele_mllg_mass, "m_{ee#gamma} [GeV]",{100,80}),
(offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_muon_motherpdgid >=40), stitch, bkg_hist).Tag(tag + "el_dy_other_").Weight(weight);*/


//pm.Push<Hist1D>(Axis(100,15,200, lpospt,    "Positive l pt"      ,{}), cut, procs, ops).Weight(wgt).Tag	(tag);
/*pm.Push<Hist1D>(Axis(150,50,200, off_muon_mllg_mass,    "m_{ll#gamma} [GeV] (Nano)"      ,{}), stitch_cut && (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll_mass > 50 && off_muon_dR && off_ele_dR && (off_muon_mllg_mass <= 200), stitch, ops3).Weight(weight_1).Tag(tag + "Check");
  pm.Push<Hist1D>(Axis(150,50,200, off_ele_mllg_mass,    "m_{ll#gamma} [GeV] (Nano)"      ,{}), stitch_cut && (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll_mass > 50 && off_ele_dR && off_muon_dR && (off_ele_mllg_mass <= 200), stitch, ops3).Weight(weight_1).Tag(tag + "Check");
  pm.Push<Hist1D>(Axis(150,50,200, off_muon_mllg_mass,    "m_{ll#gamma} [GeV] (Nano)"      ,{}), stitch_cut_more && (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll_mass > 50 && off_muon_dR && off_ele_dR && (off_muon_mllg_mass <= 200), stitch, ops3).Weight(weight_1).Tag(tag + "Check_more");
  pm.Push<Hist1D>(Axis(150,50,200, off_ele_mllg_mass,    "m_{ll#gamma} [GeV] (Nano)"      ,{}), stitch_cut_more && (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll_mass > 50 && off_ele_dR && off_muon_dR && (off_ele_mllg_mass <= 200), stitch, ops3).Weight(weight_1).Tag(tag + "Check_more");
*/
  /*pm.Push<Table>("Muon_bkg", vector<TableRow>{
     TableRow("Prompt", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags0,0,0,weight),
       TableRow("DecayedLeptonHadron",(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags1,0,0,weight),
       TableRow("TauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags2,0,0,weight),
		TableRow("PromptTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags3,0,0,weight),
			 TableRow("DirectTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags4,0,0,weight),
				  TableRow("DirectPromptTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags5,0,0,weight),
					   TableRow("DirectHadronDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags6,0,0,weight),
						    TableRow("HardProcess", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags7,0,0,weight),
							     TableRow("FromHardProcess",(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags8,0,0,weight),
							     TableRow("HardProcessTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags9,0,0,weight),
								      TableRow("DirectHardProcessTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags10,0,0,weight),
									       TableRow("FromHardProcessBeforeFSR", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags11,0,0,weight),
											TableRow("FirstCopy", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags12,0,0,weight),
												 TableRow("LastCopy", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags13,0,0,weight),
													  TableRow("LastCopyBeforeFSR", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR   && off_muon_flags14,0,0,weight)},bkgprocs,false);
 pm.Push<Table>("Electron_bkg", vector<TableRow>{
     TableRow("Prompt", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags0,0,0,weight),
       TableRow("DecayedLeptonHadron", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags1,0,0,weight),
       TableRow("TauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags2,0,0,weight),
       TableRow("PromptTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags3,0,0,weight),
       TableRow("DirectTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags4,0,0,weight),
       TableRow("DirectPromptTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags5,0,0,weight),
       TableRow("DirectHadronDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags6,0,0,weight),
       TableRow("HardProcess", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags7,0,0,weight),
       TableRow("FromHardProcess", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags8,0,0,weight),
       TableRow("HardProcessTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags9,0,0,weight),
       TableRow("DirectHardProcessTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags10,0,0,weight),
       TableRow("FromHardProcessBeforeFSR", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags11,0,0,weight),
       TableRow("FirstCopy", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags12,0,0,weight),
       TableRow("LastCopy", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags13,0,0,weight),
       TableRow("LastCopyBeforeFSR", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_flags14,0,0,weight)}, bkgprocs, false);
 
 pm.Push<Table>("Muon_bkg_sumcut", vector<TableRow>{
     TableRow("Prompt", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags0,0,0,weight),
       TableRow("DecayedLeptonHadron",(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags1,0,0,weight),
       TableRow("TauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags2,0,0,weight),
       TableRow("PromptTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags3,0,0,weight),
       TableRow("DirectTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags4,0,0,weight),
       TableRow("DirectPromptTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags5,0,0,weight),
       TableRow("DirectHadronDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags6,0,0,weight),
       TableRow("HardProcess", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags7,0,0,weight),
       TableRow("FromHardProcess",(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags8,0,0,weight),
       TableRow("HardProcessTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags9,0,0,weight),
       TableRow("DirectHardProcessTauDecayProduct", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags10,0,0,weight),
       TableRow("FromHardProcessBeforeFSR", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags11,0,0,weight),
       TableRow("FirstCopy", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags12,0,0,weight),
       TableRow("LastCopy", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags13,0,0,weight),
       TableRow("LastCopyBeforeFSR", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && off_muon_sumcut && off_muon_flags14,0,0,weight)}\
   ,bkgprocs,false);
 pm.Push<Table>("Electron_bkg_sumcut", vector<TableRow>{
     TableRow("Prompt", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags0,0,0,weight),
       TableRow("DecayedLeptonHadron", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags1,0,0,weight),
       TableRow("TauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags2,0,0,weight),
       TableRow("PromptTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags3,0,0,weight),
       TableRow("DirectTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags4,0,0,weight),
       TableRow("DirectPromptTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags5,0,0,weight),
       TableRow("DirectHadronDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags6,0,0,weight),
       TableRow("HardProcess", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags7,0,0,weight),
       TableRow("FromHardProcess", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags8,0,0,weight),
       TableRow("HardProcessTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags9,0,0,weight),
       TableRow("DirectHardProcessTauDecayProduct", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags10,0,0,weight),
       TableRow("FromHardProcessBeforeFSR", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags11,0,0,weight),
       TableRow("FirstCopy", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags12,0,0,weight),
       TableRow("LastCopy", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags13,0,0,weight),
       TableRow("LastCopyBeforeFSR", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && off_ele_sumcut && off_ele_flags14,0,0,weight)}, bkgprocs, false);
  */
    
/*pm.Push<Table>("Electron mother", vector<TableRow>{
     TableRow("Photon w/o PDG ID", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == -1),0,0,weight),
       TableRow("Photon w/o Mother", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == -2),0,0,weight),
       TableRow("Electron Mom", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == 11),0,0,weight),
       TableRow("Photon Mom", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == 22),0,0,weight),
       TableRow("Z Mom", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == 23),0,0,weight),
       TableRow("Pion Mom", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid == 111 || off_ele_motherpdgid == 211),0,0,weight),
       TableRow("Quark Mom", (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll && off_ele_dR && (off_ele_motherpdgid >0.5 && off_ele_motherpdgid < 9),0,0,weight)}, bkgprocs, false);
 pm.Push<Table>("Muon mother", vector<TableRow>{
     TableRow("Photon w/o PDG ID", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == -1),0,0,weight),
       TableRow("Photon w/o Mother", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == -2),0,0,weight),
       TableRow("Muon Mom", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == 13),0,0,weight),
       TableRow("Photon Mom", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == 22),0,0,weight),
       TableRow("Z Mom", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == 23),0,0,weight),
       TableRow("Pion Mom", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid == 111 || off_muon_motherpdgid == 211),0,0,weight),
       TableRow("Quark Mom", (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll && off_muon_dR && (off_muon_motherpdgid >0.5 && off_muon_motherpdgid < 9),0,0,weight)}, bkgprocs, false);*/
  pm.Push<Table>("Cutflow electron DY one file", vector<TableRow>{
    TableRow("Total","1" ,0,0,weight_1),
      TableRow("Stitch",stitch_cut,0,0,weight_1),
      TableRow("Stricter stitch",stitch_cut_more ,0,0,weight_1),
      TableRow("offline n lepton >= 2", offline_nel>=2,0,0,weight_1),
      TableRow("offline n photon >= 1", offline_nph>=1,0,0,weight_1),
      TableRow("Basecut",(offline_nel>=2)&&(offline_nph>=1) && off_ele_mll_mass > 50,0,0,weight_1),
      TableRow("Basecut and stitch", stitch_cut && (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll_mass > 50,0,0,weight_1),
      TableRow("Basecut and stricter stitch", stitch_cut_more && (offline_nel>=2)&&(offline_nph>=1) && off_ele_mll_mass > 50,0,0,weight_1)}, stitch, false);
 
  pm.Push<Table>("Cutflow muon DY one file", vector<TableRow>{
      TableRow("Total","1" ,0,0,weight_1),
	TableRow("offline n lepton >= 2", offline_nmu>=2,0,0,weight_1),
	//TableRow("offline dR >0.4", off_ele_dR && off_muon_dR,0,0,weight_1),
	TableRow("n lepton >= 2", "nMuon >= 2",0,0,weight_1),
	TableRow("Basecut",(offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll_mass > 50,0,0,weight_1),
	TableRow("Basecut and stitch", stitch_cut && (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll_mass > 50,0,0,weight_1),
	TableRow("Basecut and stricter stitch", stitch_cut_more && (offline_nmu>=2)&&(offline_nph>=1) && off_muon_mll_mass > 50,0,0,weight_1)}, stitch, false);

  pm.min_print_ = true;
  pm.MakePlots(1);
}
