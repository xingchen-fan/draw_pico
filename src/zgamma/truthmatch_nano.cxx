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
  string mc_path("/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/");
  string sig_path("/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/signal/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");//HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  NamedFunc trigs(el_trigs || mu_trigs);

  /*
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

  auto proc_hzg_vh_wplus =  Process::MakeShared<Baby_nano>("HToZ#gamma VH (W+)", sig,
                                                     kBlack                  ,{sig_path+"*WplusH*.root"},   "1");
  auto proc_hzg_vh_wminus =  Process::MakeShared<Baby_nano>("HToZ#gamma VH (W-)", sig,
                                                     kBlack                  ,{sig_path+"*WminusH*.root"},   "1");
  auto proc_hzg_vh_ZH =  Process::MakeShared<Baby_nano>("HToZ#gamma VH (ZH)", sig,
							    kBlack                  ,{sig_path+"*ZH*.root"},   "1");
  proc_smzg->SetLineWidth (1); proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); proc_ewk->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  //proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(3);
  proc_hzg_vbf->SetLineWidth(3);
  proc_hzg_tth->SetLineWidth(3);
  //vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg_vbf, proc_hzg_gg };
  //vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_gg};
  vector<shared_ptr<Process>> sigprocs = {proc_hzg, proc_hzg_vbf, proc_hzg_gg, proc_hzg_tth, proc_hzg_vh_wplus,  proc_hzg_vh_wminus, proc_hzg_vh_ZH };
  vector<shared_ptr<Process>> bkgprocs = { proc_ttg,  proc_ewk};
  */
  const NamedFunc nmu("nmu",fnmu);
  const NamedFunc nel("nel", fnel);
  const NamedFunc nph("nph",fnph);
  const NamedFunc jk1nph("jk1nph",fjk1nph);
  const NamedFunc jk2nph("jk2nph",fjk2nph);
  const NamedFunc ele_dR_list("ele_dR_list", fele_dR_list);
  const NamedFunc muon_dR_list("muon_dR_list", fmuon_dR_list);
  const NamedFunc pdgid("pdgid", fpdgid);
  const NamedFunc fc_pdgid("fc_pdgid", ffc_pdgid);
  const NamedFunc mother_pdgid("mother_pdgid", fmother_pdgid);
  const NamedFunc fc_mother_pdgid("fc_mother_pdgid", ffc_mother_pdgid);
  const NamedFunc muon_mll("muon_mll", fmuon_mll);
  const NamedFunc muon_mll_more("muon_mll_more", fmuon_mll_more);  
  const NamedFunc ele_mll("ele_mll", fele_mll);
  const NamedFunc ele_mll_more("ele_mll_more", fele_mll_more);
  const NamedFunc muon_mllg("muon_mllg", fmuon_mllg);
  const NamedFunc muon_mllg_more("muon_mllg_more", fmuon_mllg_more);
  const NamedFunc ele_mllg("ele_mllg", fele_mllg);
  const NamedFunc ele_mllg_more("ele_mllg_more", fele_mllg_more);
  const NamedFunc photon_truth_dR ("photon_truth_dR",fphoton_truth_dR);
  const NamedFunc photon_flags0("photon_flags0", fphoton_flags0);
  const NamedFunc photon_flags1("photon_flags1", fphoton_flags1);
  const NamedFunc photon_flags2("photon_flags2", fphoton_flags2);
  const NamedFunc photon_flags3("photon_flags3", fphoton_flags3);
  const NamedFunc photon_flags4("photon_flags4", fphoton_flags4);
  const NamedFunc photon_flags5("photon_flags5", fphoton_flags5);
  const NamedFunc photon_flags6("photon_flags6", fphoton_flags6);
  const NamedFunc photon_flags7("photon_flags7", fphoton_flags7);
  const NamedFunc photon_flags8("photon_flags8", fphoton_flags8);
  const NamedFunc photon_flags9("photon_flags9", fphoton_flags9);
  const NamedFunc photon_flags10("photon_flags10", fphoton_flags10);
  const NamedFunc photon_flags11("photon_flags11", fphoton_flags11);
  const NamedFunc photon_flags12("photon_flags12", fphoton_flags12);
  const NamedFunc photon_flags13("photon_flags13", fphoton_flags13);
  const NamedFunc photon_flags14("photon_flags14", fphoton_flags14);
  const NamedFunc photon_cone("photon_cone", fphoton_cone);
  const NamedFunc mindR_v("mindR_v", fmindR); //No dR cut, min dR varaible
//const NmaedFunc
  NamedFunc el_pt_cut ("el_pt_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> el_sig_ = fel_sig(b);
      int num1 = 0;
      int num2 = 0;
      for (unsigned i = 0; i < el_sig_.size(); i++){
	if (el_sig_[i] > 2.5) num1 +=1;
	if (el_sig_[i] > 1.5) num2 +=1;
      }
      return (num1 >=1 && num2 >=2);
    });

  NamedFunc mu_pt_cut ("mu_pt_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> mu_sig_ = fmu_sig(b);
      int num1 = 0;
      int num2 = 0;
      for (unsigned i = 0; i < mu_sig_.size(); i++){
        if (mu_sig_[i] > 2.5) num1 +=1;
        if (mu_sig_[i] > 1.5) num2 +=1;
      }
      return (num1 >=1 && num2 >=2);
    });

  NamedFunc photon_pt ("photon_pt", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double pt = 0.;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
	if(ph_sig_[iph]>0.5 && edR[iph] >0.5 && mdR[iph] > 0.5){
	  pt = b.Photon_pt()->at(iph);
	  break;
	}
      }
      return pt;
    });

  NamedFunc photon_mva ("photon_mva", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mva = 0.;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if(ph_sig_[iph]>0.5 && edR[iph] >0.5 && mdR[iph] > 0.5){
          mva = b.Photon_mvaID()->at(iph);
          break;
        }
      }
      return mva;
    });

  NamedFunc truth_part_pt("truth_part_pt", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
	if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
	  p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
	  for (int i(0); i<b.nGenPart();i++){
	    pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	    double tempdr = abs(pg.DeltaR(p));
	    if (tempdr<mindr) {
	      mindr = tempdr;
	      index = i;
	    }
	  }
	  if (mindr > 0.1) return -1;
	  break;
	}
      }
      return abs(b.GenPart_pt()->at(index));
    });

  NamedFunc fc_truth_part_pt("fc_truth_part_pt", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            bitset<15> mc_statusFlags0(b.GenPart_statusFlags()->at(i));
	    if (mc_statusFlags0[12]){
	      pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	      double tempdr = abs(pg.DeltaR(p));
	      if (tempdr<mindr) {
		mindr = tempdr;
		index = i;
	      }
	    }
	  }
          if (mindr > 0.1) return -1;
          break;
        }
      }
      return abs(b.GenPart_pt()->at(index));
    });

  NamedFunc truth_part_eta("truth_part_eta", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
            double tempdr = abs(pg.DeltaR(p));
            if (tempdr<mindr) {
              mindr = tempdr;
              index = i;
            }
          }
          if (mindr > 0.1) return -1;
          break;
        }
      }
      return abs(abs(b.GenPart_eta()->at(index)));
    });

  NamedFunc fc_truth_part_eta("fc_truth_part_eta", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            bitset<15> mc_statusFlags0(b.GenPart_statusFlags()->at(i));
	    if (mc_statusFlags0[12]){
	      pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	      double tempdr = abs(pg.DeltaR(p));
	      if (tempdr<mindr) {
		mindr = tempdr;
		index = i;
	      }
	    }
	  }
          if (mindr > 0.1) return -1;
          break;
        }
      }
      return abs(abs(b.GenPart_eta()->at(index)));
    });

  NamedFunc truth_part_mindR("truth_part_mindR", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
	    pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	    double tempdr = abs(pg.DeltaR(p));
	    if (tempdr<mindr) {
	      mindr = tempdr;
	      index = i;
	      
	    }
          }
          if (mindr > 0.1) return -0.1;
          break;
        }
      }
      TLorentzVector tpg;
      double tmindr = 999.;
      pg.SetPtEtaPhiM(b.GenPart_pt()->at(index), b.GenPart_eta()->at(index),b.GenPart_phi()->at(index),b.GenPart_mass()->at(index));
      for (int i(0); i<b.nGenPart();i++){
	bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(i));
	if(mc_statusFlags[8] && b.GenPart_pdgId()->at(i) != 22 && b.GenPart_pt()->at(i) > 5.0){
	  tpg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	  double tempdr = abs(tpg.DeltaR(pg));
	  if (tempdr<tmindr) {
	    tmindr = tempdr;
	  }
	}
      }
      return tmindr;
    });

  NamedFunc fc_truth_part_mindR("fc_truth_part_mindR", [](const Baby &b) -> NamedFunc::ScalarType{
      int index = -1;
      TLorentzVector pg, p;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double mindr = 999.;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            bitset<15> mc_statusFlags0(b.GenPart_statusFlags()->at(i));
	    if (mc_statusFlags0[12]){
	      pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	      double tempdr = abs(pg.DeltaR(p));
	      if (tempdr<mindr) {
		mindr = tempdr;
		index = i;
	      }
	    }
	  }
          if (mindr > 0.1) return -0.1;
          break;
        }
      }
      TLorentzVector tpg;
      double tmindr = 999.;
      pg.SetPtEtaPhiM(b.GenPart_pt()->at(index), b.GenPart_eta()->at(index),b.GenPart_phi()->at(index),b.GenPart_mass()->at(index));
      for (int i(0); i<b.nGenPart();i++){
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(i));
        if(mc_statusFlags[12] && i != index){
          tpg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
          double tempdr = abs(tpg.DeltaR(pg));
          if (tempdr<tmindr) {
            tmindr = tempdr;
          }
        }
      }
      return tmindr;
    });

  NamedFunc photon_iso("photon_iso", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      int index = -1;
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
	if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
	  index = j;
	  break;
	}
      }
      return abs(b.Photon_pfRelIso03_all()->at(index));
    });

  NamedFunc quark_check("quark_check", [](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector pg, p;
      bool pi = false;
      bool ph = false;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
	  p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
            double tempdr = abs(pg.DeltaR(p));
	    if (tempdr < 0.1){
	      if (abs(b.GenPart_pdgId()->at(i)) == 22) ph = true;
	      if (abs(b.GenPart_pdgId()->at(i)) == 111) pi = true;
	    }          
	  }
	  break;
	}
      }
      if (ph && pi) return 1;
      else if (ph && !pi) return 2;
      else if (!ph && pi) return 3;
      else return 4;
    });

  NamedFunc pion_check("pion_check", [](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector pg, p;
      bool ph = false;
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      for (unsigned j(0); j<b.Photon_pt()->size(); j++){
        if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
          p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
          for (int i(0); i<b.nGenPart();i++){
            pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
            double tempdr = abs(pg.DeltaR(p));
            if (tempdr < 0.1){
              if (abs(b.GenPart_pdgId()->at(i)) == 22) ph = true;
	    }
          }
          break;
        }
      }
      if (ph) return 1;
      else return 2;
    });

  NamedFunc npileup("npileup", [](const Baby &b) -> NamedFunc::ScalarType{
      return b.Pileup_nPU();
    });

  NamedFunc stitch_cut ("stitch_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
	bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
	int mc_id = b.GenPart_pdgId()->at(imc);
	//int mom_index = b.GenPart_genPartIdxMother()->at(imc);
	//int mom_id = -1;
	//if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index); 
	if(mc_id == 22){
	  bool find_part = false;
	  if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
	    genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc),
				  b.GenPart_eta()->at(imc),
				  b.GenPart_phi()->at(imc));
	    for (int imc2(0); imc2<b.nGenPart(); ++imc2){
	      bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
	      if(imc2 != imc){
		genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
				      b.GenPart_eta()->at(imc2),
				      b.GenPart_phi()->at(imc2));
		//int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
		if(genPhoton.DeltaR(genPart) < 0.05 && genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22) find_part = true;
	      }
	    }
	    if (!find_part) return false;
	  }
	}
      }
      return true;
    });

  NamedFunc dRcheck_cut ("dRcheck_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      bool find_one = false;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
	  if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
	    find_one = true;
	    break;
	  }
        }
      }
      return find_one;
    });

  NamedFunc dRcheck ("dRcheck", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      double mindR = 999.;
      double tempdR = 100.;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
          if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
            genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc),
                                  b.GenPart_eta()->at(imc),
                                  b.GenPart_phi()->at(imc));
            for (int imc2(0); imc2<b.nGenPart(); ++imc2){
              bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
              if(imc2 != imc){
                genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
				    b.GenPart_eta()->at(imc2),
				    b.GenPart_phi()->at(imc2));
                //int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
                if(genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22){
		  tempdR = genPhoton.DeltaR(genPart);
		  if (tempdR < mindR) mindR =  tempdR;
		}
              }
            }
           
          }
        }
      }
      return mindR;
    });


  NamedFunc dRcheck_IDs ("dRcheck_IDs", [](const Baby &b) -> NamedFunc::VectorType{
      TVector3 genPhoton, genPart;
      double mindR = 999.;
      double tempdR = 100.;
      vector<double> IDs = {0.,0.};
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
          if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
            genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc),
                                  b.GenPart_eta()->at(imc),
                                  b.GenPart_phi()->at(imc));
            for (int imc2(0); imc2<b.nGenPart(); ++imc2){
              bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
              if(imc2 != imc){
                genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
                                    b.GenPart_eta()->at(imc2),
                                    b.GenPart_phi()->at(imc2));
                //int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
                if(genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22){
                  tempdR = genPhoton.DeltaR(genPart);
                  if (tempdR < mindR) {
		    mindR =  tempdR;
		    IDs[0] = imc;
		    IDs[1] = imc2;
		  }
                }
              }
            }

          }
        }
      }
      return IDs;
    });

  NamedFunc dRcheck_photonpT ("dRcheck_photonpT", [dRcheck_IDs](const Baby &b) -> NamedFunc::ScalarType{
      int pimc = dRcheck_IDs.GetVector(b)[0];
      return b.GenPart_pt()->at(pimc);
    });

  NamedFunc dRcheck_pdg ("dRcheck_pdg", [dRcheck_IDs](const Baby &b) -> NamedFunc::ScalarType{
      int pimc = dRcheck_IDs.GetVector(b)[1];
      return  abs(b.GenPart_pdgId()->at(pimc));
    });

  NamedFunc dRcheck_lead ("dRcheck_lead", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      double mindR = 999.;
      double tempdR = 100.;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
          if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
            genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc),
                                  b.GenPart_eta()->at(imc),
                                  b.GenPart_phi()->at(imc));
            for (int imc2(0); imc2<b.nGenPart(); ++imc2){
              bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
              if(imc2 != imc){
                genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
                                    b.GenPart_eta()->at(imc2),
                                    b.GenPart_phi()->at(imc2));
                //int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
                if(genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22){
                  tempdR = genPhoton.DeltaR(genPart);
                  if (tempdR < mindR) mindR =  tempdR;
                }
              }
            }
	    break;
          }
        }
      }
      return mindR;
    });

  NamedFunc dRcheck_lead2 ("dRcheck_lead", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      double mindR = 999.;
      double tempdR = 100.;
      int lead_mc = 0;
      bool found = false;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
          if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
            found = true;
	    if (b.GenPart_pt()->at(imc) > b.GenPart_pt()->at(lead_mc)) lead_mc = imc;
	  }
	}
      }
      if (found){
	genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(lead_mc),
			      b.GenPart_eta()->at(lead_mc),
			      b.GenPart_phi()->at(lead_mc));
	for (int imc2(0); imc2<b.nGenPart(); ++imc2){
	  bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
	  if(imc2 != lead_mc){
	    genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
				b.GenPart_eta()->at(imc2),
				b.GenPart_phi()->at(imc2));
	    //int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
	    if(genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22){
	      tempdR = genPhoton.DeltaR(genPart);
	      if (tempdR < mindR) mindR =  tempdR;
	    }
	  }
	}
	return mindR;
      }
      else return -1;
    });

  NamedFunc dRchecklead_IDs ("dRchecklead_IDs", [](const Baby &b) -> NamedFunc::VectorType{
      TVector3 genPhoton, genPart;
      double mindR = 999.;
      double tempdR = 100.;
      vector<double> IDs = {0.,0.};
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        //int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
        if(mc_id == 22){
          if((mc_statusFlags[0] || mc_statusFlags[8]) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
            genPhoton.SetPtEtaPhi(b.GenPart_pt()->at(imc),
                                  b.GenPart_eta()->at(imc),
                                  b.GenPart_phi()->at(imc));
            for (int imc2(0); imc2<b.nGenPart(); ++imc2){
              bitset<15> mc_statusFlags2(b.GenPart_statusFlags()->at(imc2));
              if(imc2 != imc){
                genPart.SetPtEtaPhi(b.GenPart_pt()->at(imc2),
                                    b.GenPart_eta()->at(imc2),
                                    b.GenPart_phi()->at(imc2));
                //int mom_index2 = b.GenPart_genPartIdxMother()->at(imc2);
                if(genPart.Pt() > 5.0 && mc_statusFlags2[8] && b.GenPart_pdgId()->at(imc2) != 22){
                  tempdR = genPhoton.DeltaR(genPart);
                  if (tempdR < mindR) {
                    mindR =  tempdR;
                    IDs[0] = imc;
                    IDs[1] = imc2;
                  }
                }
              }
            }
	    break;
          }
        }
      }
      return IDs;
    });

  NamedFunc dRchecklead_pdg ("dRchecklead_pdg", [dRchecklead_IDs](const Baby &b) -> NamedFunc::ScalarType{
      int pimc = dRchecklead_IDs.GetVector(b)[1];
      return  abs(b.GenPart_pdgId()->at(pimc));
    });

  NamedFunc npnh ("npnh", [](const Baby &b) -> NamedFunc::ScalarType{
      bool npnh_ = true;
      int num = 0;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
	int mc_id = b.GenPart_pdgId()->at(imc);
	//int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        //int mom_id = -1;
        //if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
	if(mc_id == 22 && mc_statusFlags[12]){
	  num += 1;
	  if(mc_statusFlags[0] || mc_statusFlags[8]) npnh_ = false;
	}
      }
      if (num == 0) npnh_ = false;
      return npnh_;
    });

  NamedFunc old_stitch_cut ("old_stitch_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> eledR = fele_dR_list(b);
      std::vector<double> muondR = fmuon_dR_list(b);
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
        int mc_id = b.GenPart_pdgId()->at(imc);
        if(mc_id == 22 && b.GenPart_status()->at(imc) == 1)
          if(mc_statusFlags[0] || mc_statusFlags[8])
            for (int igamma(0); igamma<b.nPhoton(); ++igamma)
              if(ph_sig_[igamma] > 0.5 && eledR[igamma]>0.5 && muondR[igamma]>0.5){
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


  NamedFunc selection_el = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >100 && ele_mllg_more <180 && (ele_mllg_more + ele_mll_more) > 185 && photon_pt/ele_mllg_more > (15./110) && stitch_cut;
  NamedFunc selection_mu = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >100 && muon_mllg_more < 180 && (muon_mllg_more + muon_mll_more) > 185 && photon_pt/muon_mllg_more > (15./110) && stitch_cut;
  NamedFunc selection_el_dR = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >100 && ele_mllg_more <180 && (ele_mllg_more + ele_mll_more) > 185 && photon_pt/ele_mllg_more > (15./110) && dRcheck_cut;
  NamedFunc selection_mu_dR = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >100 && muon_mllg_more < 180 && (muon_mllg_more + muon_mll_more) > 185 && photon_pt/muon_mllg_more > (15./110) && dRcheck_cut;
  NamedFunc old_selection_el = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >100 && ele_mllg_more <180 && (ele_mllg_more + ele_mll_more) > 185 && photon_pt/ele_mllg_more > (15./110) && old_stitch_cut;
  NamedFunc old_selection_mu = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >100 && muon_mllg_more < 180 && (muon_mllg_more + muon_mll_more) > 185 && photon_pt/muon_mllg_more > (15./110) && old_stitch_cut;
  NamedFunc selection_el_nos = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >100 && ele_mllg_more <180 && (ele_mllg_more + ele_mll_more) > 185 && photon_pt/ele_mllg_more > (15./110);
  NamedFunc selection_mu_nos = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >100 && muon_mllg_more < 180 && (muon_mllg_more + muon_mll_more) > 185 && photon_pt/muon_mllg_more > (15./110);
  vector<NamedFunc> selection = {selection_el, selection_mu};
  vector<NamedFunc> selection_nos = {selection_el_nos, selection_mu_nos};
  vector<NamedFunc> old_selection = {old_selection_el, old_selection_mu};
  vector<NamedFunc> selection_dR = {selection_el_dR, selection_mu_dR};
  vector<NamedFunc> mllg = {ele_mllg_more, muon_mllg_more};


  NamedFunc print_photon_pt ("print_photon_pt", [el_pt_cut, photon_pt, stitch_cut](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> ph_sig_ = fph_sig(b);
      std::vector<double> edR = fele_dR_list(b);
      std::vector<double> mdR = fmuon_dR_list(b);
      double pt = 0.;
      for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
        if(ph_sig_[iph]>0.5 && edR[iph] >0.5 && mdR[iph] > 0.5){
          pt = b.Photon_pt()->at(iph);
          break;
        }
      }
      if( fnel(b)>=2 && fnph(b)>=1 && el_pt_cut.GetScalar(b) && fele_mll_more(b)>50 && fele_mllg_more(b) >100 && fele_mllg_more(b) <180 && (fele_mllg_more(b) + fele_mll_more(b)) > 185 && photon_pt.GetScalar(b)/fele_mllg_more(b) > (15/110) && stitch_cut.GetScalar(b) && ffc_pdgid(b) == -1 && fpdgid(b) != 22 && fpdgid(b) != -1) 
	std::cout <<"Run: " << b.run() <<  "Luminosity block: " << b.luminosityBlock() << " Event: " << b.event() << std::endl;
      return pt;
    });


  auto proc_dy_new= Process::MakeShared<Baby_nano>("DY + Jets",back,TColor::GetColor("#ffb400"),{mc_path+"*DYJetsToLL_M-50*"}, trigs && stitch_cut);//changed
  auto proc_sm= Process::MakeShared<Baby_nano>("SM Z#gamma",back,TColor::GetColor("#16bac5"),{mc_path+"ZGToLLG_*"}, trigs);//changed
  auto proc_dy_old= Process::MakeShared<Baby_nano>("DY + Jets Old",sig,kRed,{mc_path+"*DYJetsToLL_M-50*"}, trigs && old_stitch_cut);
  auto proc_ggf_sig= Process::MakeShared<Baby_nano>("GGF",back,TColor::GetColor("#ffb400"),{sig_path+"*GluGlu*"}, trigs); 
  //vector<shared_ptr<Process>> procs = {proc_smzg, proc_dy_nos,proc_hzg_vbf, proc_hzg_gg, proc_hzg_tth, proc_hzg_vh };
  proc_dy_new->SetLineWidth(1);
  proc_dy_old->SetLineWidth(1);
  proc_sm->SetLineWidth(1); 
  proc_ggf_sig->SetLineWidth(3);

  vector<shared_ptr<Process>> dy = {proc_dy_new};
  vector<shared_ptr<Process>> sm = {proc_sm};
  vector<shared_ptr<Process>> stitchcom = {proc_dy_new, proc_dy_old};
  vector<shared_ptr<Process>> ggf = {proc_ggf_sig};
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
    .YAxis(YAxisType::log)
    .Stack(StackType::shapes)
    .Overflow(OverflowType::none)
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
  vector<PlotOpt> op_per = {lin_lumi().Stack(StackType::shapes).Bottom(BottomType::diff)}; 
  const NamedFunc weight_p("weight_p",[](const Baby &b) -> NamedFunc::ScalarType{
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

  const NamedFunc weight_tot("weight_tot",[](const Baby &b) -> NamedFunc::ScalarType{
      float w = 1.0;
      if (b.Generator_weight()<0) return -w;
      else return w;
    });



  const NamedFunc weight("weight",[](const Baby &b) -> NamedFunc::ScalarType{
      float w = 0.0;
      float xsec = -999.;
      float HToZG(0.001533), ZToLL(0.100974);
      int tot = 1;
      float sign =  b.Generator_weight() > 0? 1.: -1.;
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("ZGToLLG")) {  xsec = 1000*55.48; tot = 18727540;}
        if(file.Contains("DYJetsToLL_M-50")) {xsec = 1000* 6077.22; tot = 102486448;}

        if(file.Contains("GluGluHToZG"))          {xsec = 1000 * 1000 * HToZG * ZToLL * 48.58 ; tot = 400000;}
        if(file.Contains("VBFHToZG"))             {xsec = 1000 * 1000 * HToZG * ZToLL * 3.782 ; tot = 200000;}
        if(file.Contains("WplusH_HToZG"))         {xsec = 1000 * 1000 * HToZG * 0.831 ; tot = 299978;}
        if(file.Contains("WminusH_HToZG"))        {xsec = 1000 * 1000 * HToZG * 0.527 ; tot = 299276;}
        if(file.Contains("ZH_HToZG"))             {xsec = 1000 * 1000 * HToZG * 0.8839; tot = 297389;}
        if(file.Contains("ttHToZG"))              {xsec = 1000 * 1000 * HToZG * 0.5071; tot = 200000;}
        w = xsec * 41.5 * sign/tot;


      }
      return w;
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

  vector<PlotOpt> bkg_hist_lin = {style().Stack(StackType::data_norm)
                              .LogMinimum(1)
                              .CanvasWidth(600)
                              .LabelSize(0.04)
                              .YAxis(YAxisType::linear)
                              .Title(TitleType::info)};
  PlotMaker pm;
  tag = "Nano_";
  vector<string> lepton ={"ele_", "muon_"};
  for (int i(0); i<1; i++){
    //    pm.Push<Hist1D>(Axis(80,100,180, mllg[i],    "m_{ll#gamma}"      ,{}), selection_nos.at(i), stitch, op_per).Weight(weight).Tag(tag + "stitch_comp_" + lepton[i]);  
    //pm.Push<Hist1D>(Axis(80,-0.6,1, photon_mva,    "MVA"      ,{}), selection_nos.at(i) && photon_mva > -0.4, stitchcom, ops2).Weight(weight_1).Tag(tag + "stitch_compmva_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(80,100,180, mllg[i],    "m_{ll#gamma}"      ,{}), nel>=2 && nph>=1 && el_pt_cut, stitch, ops2).Weight(weight).Tag(tag + "weight_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(80,100,180, mllg[i],    "m_{ll#gamma}"      ,{}), nel>=2 && nph>=1 && el_pt_cut, stitch, ops2).Weight(weight_1).Tag(tag + "weight1_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(90,10,100, photon_cone,    "no mom photon pt"      ,{}), nel>=2 && nph>=1 && ele_mll>50 && ele_mllg >100 && ele_mllg<180 && stitch_cut && mother_pdgid==-2, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_nomom_ele_");
  //pm.Push<Hist1D>(Axis(90,10,100, photon_cone,    "no mom photon pt"      ,{}), nmu>=2 && nph>=1 && muon_mll>50 && muon_mllg >100 && muon_mllg<180 && stitch_cut && mother_pdgid==-2,stitch, ops2).Weight(weight_1).Tag(tag + "stitch_nomom_muon_");
  /*pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), nel>=2 && nph>=1 && ele_mll>50 && stitch_cut_more, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_pdg_ele_");
  pm.Push<Hist1D>(Axis(41,-2,40, pdgid,    "Photon PDG ID"      ,{}), nmu>=2 && nph>=1 && muon_mll>50 && stitch_cut_more,stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_pdg_muon_");
  pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), nel>=2 && nph>=1 && ele_mll>50 && stitch_cut_more, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_mpdg_ele_");
  pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), nmu>=2 && nph>=1 && muon_mll>50 && stitch_cut_more,stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_mpdg_muon_");
  pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), nel>=2 && nph>=1 && ele_mll>50 && stitch_cut_more, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_true_ele_");
  pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), nmu>=2 && nph>=1 && muon_mll>50 && stitch_cut_more,stitch, ops2).Weight(weight_1).Tag(tag + "stitch_more_true_muon_");*/
  /*  
  pm.Push<Table>("Electron_stitch_check_", vector<TableRow>{
      TableRow("No cut","1",0,0,weight_1),
	TableRow("Lepton",nel>=2,0,0,weight_1),
	TableRow("Photon",nel>=2 && nph>=1,0,0,weight_1),
	TableRow("Photon (only edR)",nel>=2 && jk1nph>=1,0,0,weight_1),
	TableRow("Photon (only mdR)",nel>=2 && jk2nph>=1,0,0,weight_1),
	TableRow("basic",nel>=2 && nph>=1 && ele_mll>50,0,0,weight_1),
	TableRow("All",nel>=2 && nph>=1 && ele_mll>50 && stitch_cut,0,0,weight_1)}, stitch, false);
  pm.Push<Table>("Muon_stitch_check_", vector<TableRow>{
      TableRow("No cut","1",0,0,weight_1),
        TableRow("Lepton",nmu>=2,0,0,weight_1),
        TableRow("Photon",nmu>=2 && nph>=1,0,0,weight_1),
	TableRow("Photon (only edR)",nmu>=2 && jk1nph>=1,0,0,weight_1),
        TableRow("Photon (only mdR)",nmu>=2 && jk2nph>=1,0,0,weight_1),
	TableRow("basic",nmu>=2 && nph>=1 && muon_mll>50,0,0,weight_1),
        TableRow("All",nmu>=2 && nph>=1 && muon_mll>50 && stitch_cut,0,0,weight_1)}, stitch, false);
  */    
  /*
    pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_pdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_mpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_true_" + lepton[i]);

    pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_pdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_mpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_true_" + lepton[i]);

    pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_pdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_mpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_true_" + lepton[i]);

    pm.Push<Table>("Stitch_" + lepton[i], vector<TableRow>{
      TableRow("Prompt", selection.at(i) && photon_flags0,0,0,weight_1),
	TableRow("DecayedLeptonHadron", selection.at(i) && photon_flags1,0,0,weight_1),
	TableRow("TauDecayProduct", selection.at(i) && photon_flags2,0,0,weight_1),
	TableRow("PromptTauDecayProduct", selection.at(i) && photon_flags3,0,0,weight_1),
	TableRow("DirectTauDecayProduct", selection.at(i) && photon_flags4,0,0,weight_1),
	TableRow("DirectPromptTauDecayProduct", selection.at(i) && photon_flags5,0,0,weight_1),
	TableRow("DirectHadronDecayProduct", selection.at(i) && photon_flags6,0,0,weight_1),
	TableRow("HardProcess", selection.at(i) && photon_flags7,0,0,weight_1),
	TableRow("FromHardProcess", selection.at(i) && photon_flags8,0,0,weight_1),
	TableRow("HardProcessTauDecayProduct", selection.at(i) && photon_flags9,0,0,weight_1),
	TableRow("DirectHardProcessTauDecayProduct", selection.at(i) && photon_flags10,0,0,weight_1),
	TableRow("FromHardProcessBeforeFSR", selection.at(i) && photon_flags11,0,0,weight_1),
	TableRow("FirstCopy", selection.at(i) && photon_flags12,0,0,weight_1),
	TableRow("LastCopy", selection.at(i) && photon_flags13,0,0,weight_1),
	TableRow("LastCopyBeforeFSR", selection.at(i) && photon_flags14,0,0,weight_1)}, stitch, false);
    
    pm.Push<Table>("Old_Stitch_" + lepton[i], vector<TableRow>{
	TableRow("Prompt", old_selection.at(i) && photon_flags0,0,0,weight_1),
	  TableRow("DecayedLeptonHadron", old_selection.at(i) && photon_flags1,0,0,weight_1),
	  TableRow("TauDecayProduct", old_selection.at(i) && photon_flags2,0,0,weight_1),
	  TableRow("PromptTauDecayProduct", old_selection.at(i) && photon_flags3,0,0,weight_1),
	  TableRow("DirectTauDecayProduct", old_selection.at(i) && photon_flags4,0,0,weight_1),
	  TableRow("DirectPromptTauDecayProduct", old_selection.at(i) && photon_flags5,0,0,weight_1),
	  TableRow("DirectHadronDecayProduct", old_selection.at(i) && photon_flags6,0,0,weight_1),
	  TableRow("HardProcess", old_selection.at(i) && photon_flags7,0,0,weight_1),
	  TableRow("FromHardProcess", old_selection.at(i) && photon_flags8,0,0,weight_1),
	  TableRow("HardProcessTauDecayProduct", old_selection.at(i) && photon_flags9,0,0,weight_1),
	  TableRow("DirectHardProcessTauDecayProduct", old_selection.at(i) && photon_flags10,0,0,weight_1),
	  TableRow("FromHardProcessBeforeFSR", old_selection.at(i) && photon_flags11,0,0,weight_1),
	  TableRow("FirstCopy", old_selection.at(i) && photon_flags12,0,0,weight_1),
	  TableRow("LastCopy", old_selection.at(i) && photon_flags13,0,0,weight_1),
	  TableRow("LastCopyBeforeFSR", old_selection.at(i) && photon_flags14,0,0,weight_1)}, stitch, false);

    pm.Push<Table>("Stitch_npnh" + lepton[i], vector<TableRow>{
	TableRow("Prompt", selection.at(i) && photon_flags0 && npnh,0,0,weight_1),
	  TableRow("DecayedLeptonHadron", selection.at(i) && photon_flags1 && npnh,0,0,weight_1),
	  TableRow("TauDecayProduct", selection.at(i) && photon_flags2 && npnh,0,0,weight_1),
	  TableRow("PromptTauDecayProduct", selection.at(i) && photon_flags3 && npnh,0,0,weight_1),
	  TableRow("DirectTauDecayProduct", selection.at(i) && photon_flags4 && npnh,0,0,weight_1),
	  TableRow("DirectPromptTauDecayProduct", selection.at(i) && photon_flags5 && npnh,0,0,weight_1),
	  TableRow("DirectHadronDecayProduct", selection.at(i) && photon_flags6 && npnh,0,0,weight_1),
	  TableRow("HardProcess", selection.at(i) && photon_flags7 && npnh,0,0,weight_1),
	  TableRow("FromHardProcess", selection.at(i) && photon_flags8 && npnh,0,0,weight_1),
	  TableRow("HardProcessTauDecayProduct", selection.at(i) && photon_flags9 && npnh,0,0,weight_1),
	  TableRow("DirectHardProcessTauDecayProduct", selection.at(i) && photon_flags10 && npnh,0,0,weight_1),
	  TableRow("FromHardProcessBeforeFSR", selection.at(i) && photon_flags11 && npnh,0,0,weight_1),
	  TableRow("FirstCopy", selection.at(i) && photon_flags12 && npnh,0,0,weight_1),
	  TableRow("LastCopy", selection.at(i) && photon_flags13 && npnh,0,0,weight_1),
	  TableRow("LastCopyBeforeFSR", selection.at(i) && photon_flags14 && npnh,0,0,weight_1)}, stitch, false);
  */
    //pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Motherless Photon PDG ID"      ,{}), selection.at(i) && mother_pdgid==-2, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_mlpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(41,-1,40, fc_pdgid,    "First Copy Photon PDG ID"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fcpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(42,-2,40, fc_mother_pdgid,    "First Copy Photon Mother PDG ID"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fcmpdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(51,-1,50, truth_part_pt,    "Truth photon pt"      ,{}), selection.at(i) && pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_truthpt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(51,-1,50, truth_part_pt,    "Truth photon pt"      ,{}), old_selection.at(i) && pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_truthpt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(52,0,2.6, truth_part_eta,    "Truth photon eta"      ,{}), selection.at(i) && pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_trutheta_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(52,0,2.6, truth_part_eta,    "Truth photon eta"      ,{}), old_selection.at(i) && pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_trutheta_" + lepton[i]);

    //pm.Push<Hist1D>(Axis(200,0,0.2, truth_part_mindR,    "Truth photon mindR"      ,{}), selection.at(i) && pdgid==22 && truth_part_mindR<900, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_truthdr_900_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(200,0,2, photon_iso,    "Isolation"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_photoniso_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(200,0,0.2, truth_part_mindR,    "Truth photon mindR"      ,{}), old_selection.at(i) && pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_truthdr_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(51,-1,50, fc_truth_part_pt,    "First Copy Truth photon pt"      ,{}), selection.at(i) && fc_pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fctruthpt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(52,0,2.6, fc_truth_part_eta,    "First Copy Truth photon eta"      ,{}), selection.at(i) && fc_pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fctrutheta_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(200,0,0.2, fc_truth_part_mindR,    "First Copy Truth photon mindR"      ,{}), selection.at(i) && fc_pdgid==22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fctruthdr_" + lepton[i]);
/*
    pm.Push<Table>("Stitch_motherless_" + lepton[i], vector<TableRow>{
        TableRow("Prompt", selection.at(i) && photon_flags0 && mother_pdgid==-2,0,0,weight_1),
          TableRow("DecayedLeptonHadron", selection.at(i) && photon_flags1 && mother_pdgid==-2,0,0,weight_1),
          TableRow("TauDecayProduct", selection.at(i) && photon_flags2 && mother_pdgid==-2,0,0,weight_1),
          TableRow("PromptTauDecayProduct", selection.at(i) && photon_flags3 && mother_pdgid==-2,0,0,weight_1),
          TableRow("DirectTauDecayProduct", selection.at(i) && photon_flags4 && mother_pdgid==-2,0,0,weight_1),
          TableRow("DirectPromptTauDecayProduct", selection.at(i) && photon_flags5 && mother_pdgid==-2,0,0,weight_1),
          TableRow("DirectHadronDecayProduct", selection.at(i) && photon_flags6 && mother_pdgid==-2,0,0,weight_1),
          TableRow("HardProcess", selection.at(i) && photon_flags7 && mother_pdgid==-2,0,0,weight_1),
          TableRow("FromHardProcess", selection.at(i) && photon_flags8 && mother_pdgid==-2,0,0,weight_1),
          TableRow("HardProcessTauDecayProduct", selection.at(i) && photon_flags9 && mother_pdgid==-2,0,0,weight_1),
          TableRow("DirectHardProcessTauDecayProduct", selection.at(i) && photon_flags10 && mother_pdgid==-2,0,0,weight_1),
          TableRow("FromHardProcessBeforeFSR", selection.at(i) && photon_flags11 && mother_pdgid==-2,0,0,weight_1),
          TableRow("FirstCopy", selection.at(i) && photon_flags12 && mother_pdgid==-2,0,0,weight_1),
          TableRow("LastCopy", selection.at(i) && photon_flags13 && mother_pdgid==-2,0,0,weight_1),
          TableRow("LastCopyBeforeFSR", selection.at(i) && photon_flags14 && mother_pdgid==-2,0,0,weight_1)}, stitch, false);
*/
/*    pm.Push<Table>("Stitch_pdgid_" + lepton[i], vector<TableRow>{
	TableRow("Total", selection.at(i),0,0,weight_1),
	  TableRow("No truth information", selection.at(i) && pdgid == -1,0,0,weight_1),
	  TableRow("Quarks", selection.at(i) && pdgid >= 1 && pdgid <= 8,0,0,weight_1),
	  TableRow("Electron", selection.at(i) && pdgid == 11,0,0,weight_1),
	  TableRow("Muon", selection.at(i) && pdgid == 13,0,0,weight_1),
	  TableRow("Photon", selection.at(i) && pdgid == 22,0,0,weight_1),
	  TableRow("Pions 0", selection.at(i) && pdgid == 111,0,0,weight_1),
	  TableRow("Pions +/-", selection.at(i) && pdgid == 211,0,0,weight_1)},stitch, false);
    
        pm.Push<Table>("Stitch_mpdgid_quarks_" + lepton[i], vector<TableRow>{
	TableRow("Total", selection.at(i) && pdgid >= 1 && pdgid <= 8,0,0,weight_1),
	TableRow("No mother information", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid == 111,0,0,weight_1),
	  TableRow("Others(including pions)", selection.at(i) && pdgid >= 1 && pdgid <= 8 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);
   
    pm.Push<Table>("Stitch_mpdgid_electron_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) &&  pdgid == 11,0,0,weight_1),
	  TableRow("No mother information", selection.at(i) && pdgid == 11 && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid == 11 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid == 11 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && pdgid == 11 && mother_pdgid == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid == 11 && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid == 11 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid == 11 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && pdgid == 11 && mother_pdgid == 111,0,0,weight_1),
	  TableRow("Others(including pions)", selection.at(i) && pdgid == 11 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);
    
    pm.Push<Table>("Stitch_mpdgid_photon_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) &&  pdgid == 22,0,0,weight_1),
          TableRow("No mother information", selection.at(i) && pdgid == 22 && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid == 22 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid == 22 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && pdgid == 22 && mother_pdgid == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid == 22 && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid == 22 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid == 22 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && pdgid == 22 && mother_pdgid == 111,0,0,weight_1),
	  TableRow("Others(including pions)", selection.at(i) && pdgid == 22 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);
   
    pm.Push<Table>("Stitch_mpdgid_p0_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) &&  pdgid == 111,0,0,weight_1),
          TableRow("No mother information", selection.at(i) && pdgid == 111 && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid == 111 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid == 111 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && pdgid == 111 && mother_pdgid == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid == 111 && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid == 111 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid == 111 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && pdgid == 111 && mother_pdgid == 111,0,0,weight_1),
	  TableRow("Others(including pions)", selection.at(i) && pdgid == 111 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);
*/ 
    /* pm.Push<Table>("Stitch_mpdgid_ppm_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) &&  pdgid == 211,0,0,weight_1),
          TableRow("No mother information", selection.at(i) && pdgid == 211 && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid == 211 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid == 211 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && pdgid == 211 && mother_pdgid == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid == 211 && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid == 211 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid == 211 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && pdgid == 211 && mother_pdgid == 111,0,0,weight_1),
          TableRow("Pions +/-", selection.at(i) && pdgid == 211 && mother_pdgid == 211,0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && pdgid == 211 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);*/
    /*pm.Push<Table>("Old_stitch_pdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", old_selection.at(i),0,0,weight_1),
          TableRow("No truth particle", old_selection.at(i) && pdgid == -1,0,0,weight_1),
          TableRow("Quarks", old_selection.at(i) && pdgid >= 1 && pdgid <= 8,0,0,weight_1),
          TableRow("Electron", old_selection.at(i) && pdgid == 11,0,0,weight_1),
          TableRow("Photon", old_selection.at(i) && pdgid == 22,0,0,weight_1),
          TableRow("Pions", old_selection.at(i) && (pdgid == 111 || pdgid == 211),0,0,weight_1),
	  TableRow("Others(including pions)", old_selection.at(i) && pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);
    */
    /*    pm.Push<Table>("Stitch_mpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i),0,0,weight_1),
          TableRow("No mother information", selection.at(i) && mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && mother_pdgid == 11,0,0,weight_1),
	  TableRow("Muon", selection.at(i) && mother_pdgid == 13,0,0,weight_1),
	  TableRow("Photon", selection.at(i) && mother_pdgid == 22,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && mother_pdgid == 21,0,0,weight_1),
	  TableRow("Z", selection.at(i) && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", selection.at(i) && (mother_pdgid == 111 || mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && mother_pdgid >= 40,0,0,weight_1)},stitch, false);*/
    /*
    pm.Push<Table>("Stitch_fcpdgid_" + lepton[i], vector<TableRow>{
	TableRow("Total", selection.at(i),0,0,weight_1),
	  TableRow("No truth particle", selection.at(i) && fc_pdgid == -1,0,0,weight_1),
	  TableRow("Quarks", selection.at(i) && fc_pdgid >= 1 && fc_pdgid <= 8,0,0,weight_1),
	  TableRow("Electron", selection.at(i) && fc_pdgid == 11,0,0,weight_1),
	  TableRow("Photon", selection.at(i) && fc_pdgid == 22,0,0,weight_1),
	  TableRow("Pions", selection.at(i) && (fc_pdgid == 111 || fc_pdgid == 211),0,0,weight_1),
	  TableRow("Others(including pions)", selection.at(i) && fc_pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);
    pm.Push<Table>("Stitch_fcmpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i),0,0,weight_1),
	  TableRow("No mother particle", selection.at(i) && fc_mother_pdgid == -2,0,0,weight_1),
          TableRow("No truth particle", selection.at(i) && fc_mother_pdgid == -1,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && fc_mother_pdgid >= 1 && fc_mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && fc_mother_pdgid == 11,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && fc_mother_pdgid == 21,0,0,weight_1),
	  TableRow("Z", selection.at(i) && fc_mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", selection.at(i) && (fc_mother_pdgid == 111 || fc_mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && fc_mother_pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);
    pm.Push<Table>("Stitch_truth_photon_mpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) && pdgid == 22,0,0,weight_1),
          TableRow("No mother particle", selection.at(i) && pdgid == 22 &&  mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && pdgid == 22 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && pdgid == 22 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && pdgid == 22 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && pdgid == 22 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", selection.at(i) && pdgid == 22 && (mother_pdgid == 111 || mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && pdgid == 22 && mother_pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);
    pm.Push<Table>("Stitch_fc_truth_photon_mpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) && fc_pdgid == 22,0,0,weight_1),
          TableRow("No mother particle", selection.at(i) && fc_pdgid == 22 &&  fc_mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && fc_pdgid == 22 && fc_mother_pdgid >= 1 && fc_mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && fc_pdgid == 22 && fc_mother_pdgid == 11,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && fc_pdgid == 22 && fc_mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && fc_pdgid == 22 && fc_mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", selection.at(i) && fc_pdgid == 22 && (fc_mother_pdgid == 111 || fc_mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && fc_pdgid == 22 && fc_mother_pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);

    pm.Push<Table>("Old_stitch_truth_photon_mpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", old_selection.at(i) && pdgid == 22,0,0,weight_1),
          TableRow("No mother particle", old_selection.at(i) && pdgid == 22 &&  mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", old_selection.at(i) && pdgid == 22 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", old_selection.at(i) && pdgid == 22 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Gluon", old_selection.at(i) && pdgid == 22 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", old_selection.at(i) && pdgid == 22 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", old_selection.at(i) && pdgid == 22 && (mother_pdgid == 111 || mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", old_selection.at(i) && pdgid == 22 && mother_pdgid >= 40,0,0,weight_1)},stitch, false, true, true, true, false, false);
    */

    /* pm.Push<Table>("Stitch_fc_notruth_mpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid != -2,0,0,weight_1),
          //TableRow("No mother particle", selection.at(i) && fc_pdgid == -1 && pdgid = 22 &&  fc_mother_pdgid == -2,0,0,weight_1),
          TableRow("Quarks", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid >= 1 && mother_pdgid <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid == 11,0,0,weight_1),
          TableRow("Gluon", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid == 21,0,0,weight_1),
          TableRow("Z", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid == 23,0,0,weight_1),
          TableRow("Pions", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && (mother_pdgid == 111 || mother_pdgid == 211),0,0,weight_1),
          TableRow("Others(including pions)", selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid >= 40,0,0,weight_1)},stitch, false);
    */

    //pm.Push<Hist1D>(Axis(50,0,50, photon_pt,    "Signal Photon Pt (GeV)"      ,{}), selection.at(i) && pdgid == 22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_truthphotonpt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(100,0,100, photon_pt,    "Signal Photon Pt (GeV)"      ,{}), selection.at(i) && pdgid == 11, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_truthelept_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(200,0,0.2, photon_truth_dR,    "Photon and Truth Part dR"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_dR_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(100,0,0.1, photon_truth_dR,    "Photon and Truth Part dR"      ,{}), selection.at(i) && pdgid == 22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_truthphotondR_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(100,0,0.1, photon_truth_dR,    "Photon and Truth Part dR"      ,{}), selection.at(i) && pdgid == 11, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_trutheledR_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(50,0,50, print_photon_pt,    "Signal Photon Pt (GeV)"      ,{}), selection.at(i) && fc_pdgid == -1 && pdgid == 22 && mother_pdgid != -2, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_103pt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(50,0,50, photon_pt,    "Signal Photon Pt (GeV)"      ,{}), selection.at(i) && fc_pdgid == -1 && pdgid != 22 && pdgid != -1, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_68pt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(50,0,50, print_photon_pt,    "Signal Photon Pt (GeV)"      ,{}), selection.at(i) && fc_pdgid == 22, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_fctruthpt_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(100,0,100, "Pileup_nPU",    "Pile up number"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npu_" + lepton[i]);

    //pm.Push<Hist1D>(Axis(40,0,40, dRchecklead_pdg,    "PDG ID"      ,{}), selection_dR.at(i), sm, ops2).Weight(weight_1).Tag(tag + "sm_pdg_" + lepton[i]);
    // pm.Push<Hist1D>(Axis(40,0,40, dRchecklead_pdg,    "PDG ID"      ,{}), selection_dR.at(i) && dRcheck_lead <= 0.05, sm, ops2).Weight(weight_1).Tag(tag + "sm_005pdg_" + lepton[i]);
    //pm.Push<Hist1D>(Axis(50,0,0.5, dRcheck_lead,    "mindR"      ,{}), selection_dR.at(i) && dRcheck_cut && dRcheck_lead <= 0.5, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_mindR_" );
    //pm.Push<Hist1D>(Axis(40,0,40, dRchecklead_pdg,    "PDG ID"      ,{}), selection_dR.at(i) && dRcheck_cut && dRcheck_lead <= 0.5, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_pdg_" );
    //pm.Push<Hist1D>(Axis(50,50,150, mllg[i],    "m_{ll#gamma}"      ,{}), nel>=2 && nph>=1 && el_pt_cut && , stitch, ops2).Weight(weight_1).Tag(tag + "mllg_" + lepton[i] );
    //pm.Push<Hist1D>(Axis(100,0,100, dRcheck_photonpT,    "Truth Photon pT"      ,{}), selection_dR.at(i) && dRcheck > 0.01, sm, ops2).Weight(weight_1).Tag(tag + "sm_lppT_" + lepton[i]);
    /*pm.Push<Table>("SM_dRpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection_dR.at(i),0,0,weight_1),
	  TableRow("Quarks", selection.at(i) && dRchecklead_pdg >= 1 && dRchecklead_pdg <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && dRchecklead_pdg == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && dRchecklead_pdg == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && dRchecklead_pdg == 22,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && dRchecklead_pdg == 111,0,0,weight_1),
          TableRow("Pions +/-", selection.at(i) && dRchecklead_pdg == 211,0,0,weight_1)},sm, false);
    pm.Push<Table>("SM_005dRpdgid_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection_dR.at(i),0,0,weight_1),
          TableRow("Quarks", selection.at(i) && dRcheck_lead <= 0.05 && dRchecklead_pdg >= 1 && dRchecklead_pdg <= 8,0,0,weight_1),
          TableRow("Electron", selection.at(i) && dRcheck_lead <= 0.05 && dRchecklead_pdg == 11,0,0,weight_1),
          TableRow("Muon", selection.at(i) && dRcheck_lead <= 0.05 && dRchecklead_pdg == 13,0,0,weight_1),
          TableRow("Photon", selection.at(i) && dRchecklead_pdg == 22 && dRcheck_lead <= 0.05,0,0,weight_1),
          TableRow("Pions 0", selection.at(i) && dRcheck_lead <= 0.05 && dRchecklead_pdg == 111,0,0,weight_1),
          TableRow("Pions +/-", selection.at(i) && dRcheck_lead <= 0.05 && dRchecklead_pdg == 211,0,0,weight_1)},sm, false);
    */
    /*pm.Push<Table>("Stitch_quark_pion_" + lepton[i], vector<TableRow>{
        TableRow("Total Quarks", selection.at(i) && pdgid >= 1 && pdgid <= 8 ,0,0,weight_1),
          TableRow("Photon & Pi0", selection.at(i) && pdgid >= 1 && pdgid <= 8 && quark_check == 1,0,0,weight_1),
          TableRow("Only Photon", selection.at(i) && pdgid >= 1 && pdgid <= 8 && quark_check == 2,0,0,weight_1),
          TableRow("Only Pi0", selection.at(i) && pdgid >= 1 && pdgid <= 8 && quark_check == 3,0,0,weight_1),
          TableRow("No Photon or Pi0", selection.at(i) && pdgid >= 1 && pdgid <= 8 && quark_check == 4,0,0,weight_1),
          TableRow("Total Pion", selection.at(i) && pdgid == 111,0,0,weight_1),
          TableRow("Photon", selection.at(i) && pdgid == 111 && pion_check == 1,0,0,weight_1),
          TableRow("No Photon", selection.at(i) && pdgid == 111 && pion_check == 2,0,0,weight_1)},stitch, false);
    
    pm.Push<Table>("Stitch_pileup_" + lepton[i], vector<TableRow>{
        TableRow("Total", selection.at(i),0,0,weight_1),
	  TableRow("PU <= 20", selection.at(i) && npileup <= 20,0,0,weight_1),
          TableRow("PU <= 20 No truth info", selection.at(i) && pdgid == -1 && npileup <= 20,0,0,weight_1),
          TableRow("PU < 20 <= 30", selection.at(i) && npileup > 20 && npileup <= 30,0,0,weight_1),
	  TableRow("PU < 20 <= 30 No truth info", selection.at(i) && pdgid == -1 && npileup > 20 && npileup <= 30,0,0,weight_1),
	  TableRow("PU < 30 <= 40", selection.at(i) && npileup > 30 && npileup <= 40,0,0,weight_1),
	  TableRow("PU < 30 <= 40 No truth info", selection.at(i) && pdgid == -1 && npileup > 30 && npileup <= 40,0,0,weight_1),
	  TableRow("PU < 40 <= 50", selection.at(i) && npileup > 40 && npileup <= 50,0,0,weight_1),
	  TableRow("PU < 40 <= 50 No truth info", selection.at(i) && pdgid == -1 && npileup > 40 && npileup <= 50,0,0,weight_1),
	  TableRow("PU > 50", selection.at(i) && npileup > 50,0,0,weight_1),
	  TableRow("PU > 50 No truth info", selection.at(i) && pdgid == -1 && npileup > 50,0,0,weight_1)},stitch, false);
    */
   
    pm.Push<Hist2D>(
                    Axis(100,15,80, photon_pt,      "Photon pT [GeV]",{}),
                    Axis(100,0.1,1.1,     mindR_v, "min dR SMZ#gamma",{}), 
                    selection_nos[i], sm, bkg_hist).Tag(tag + "mindR_sm_" + lepton[i]).Weight(weight_1);
    //pm.Push<Hist2D>(
    //              Axis(100,15,80, photon_pt,      "Photon pT [GeV]",{}),
    //              Axis(100,0.1,1.1,     mindR_v, "min dR DY",{}),
    //              selection[i], dy, bkg_hist).Tag(tag + "mindR_dy_" + lepton[i]).Weight(weight_1);
    pm.Push<Hist2D>(
                    Axis(100,15,80, photon_pt,      "Photon pT [GeV]",{}),
                    Axis(100,0.1,1.1,     mindR_v, "min dR GGF Sig",{}),
                    selection_nos[i], ggf, bkg_hist_lin).Tag(tag + "mindR_ggf_" + lepton[i]).Weight(weight_1);
  }
  //pm.Push<Hist1D>(Axis(50,0,0.5, dRcheck_lead,    "mindR"      ,{}), dRcheck_cut && dRcheck_lead <= 0.5, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_bare_" );
  //pm.Push<Hist1D>(Axis(40,0,40, dRchecklead_pdg,    "PDG ID"      ,{}), dRcheck_cut && dRcheck_lead <= 0.5, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_bare_pdg_");
  //pm.Push<Hist1D>(Axis(100,50,150, mllg[0],    "m_{ll#gamma}"      ,{}), nel>=2 && nph>=1 && ele_mll_more>50 && el_pt_cut && mllg[0] < 150 && mllg[0] > 50, sm, ops2).Weight(weight_1).Tag(tag + "mllg_" + lepton[0] );
  pm.min_print_ = true;
  pm.MakePlots(1);
}
