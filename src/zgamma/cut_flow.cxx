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
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");//HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  NamedFunc trigs(el_trigs || mu_trigs);

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

  const NamedFunc nmu("nmu",fnmu);
  const NamedFunc nel("nel", fnel);
  const NamedFunc nph("nph",fnph);
  const NamedFunc jk1nph("jk1nph",fjk1nph);
  const NamedFunc jk2nph("jk2nph",fjk2nph);
  const NamedFunc ele_dR_list("ele_dR_list", fele_dR_list);
  const NamedFunc muon_dR_list("muon_dR_list", fmuon_dR_list);
  const NamedFunc pdgid("pdgid", fpdgid);
  const NamedFunc mother_pdgid("mother_pdgid", fmother_pdgid);
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
  // Cut flow
  const NamedFunc ph_sig_pt("ph_sig_pt",  fph_sig_pt);
  const NamedFunc ph_sig_eta("ph_sig_eta",  fph_sig_eta);
  const NamedFunc ph_sig_barrel("ph_sig_barrel",  fph_sig_barrel);
  const NamedFunc ph_sig_endcap("ph_sig_endcap",  fph_sig_endcap);
  const NamedFunc ph_sig_eleveto("ph_sig_eleveto",  fph_sig_eleveto);
  const NamedFunc el_lead("el_lead", fel_lead);
  const NamedFunc el_sublead("el_sublead", fel_sublead);
  const NamedFunc mu_lead("mu_lead", fmu_lead);
  const NamedFunc mu_sublead("mu_sublead", fmu_sublead);
  
//const NmaedFunc
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

  NamedFunc stitch_cut ("stitch_cut", [](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 genPhoton, genPart;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
	bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
	int mc_id = b.GenPart_pdgId()->at(imc);
	//int mom_index = b.GenPart_genPartIdxMother()->at(imc);
	//int mom_id = -1;
	//if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index); 
	if(mc_id == 22 &&  mc_statusFlags[12]){
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
		if(genPhoton.DeltaR(genPart) < 0.05 && genPart.Pt() > 5.0 && mc_statusFlags2[8] && mc_statusFlags2[12] &&  b.GenPart_pdgId()->at(imc2) != 22) find_part = true;
	      }
	    }
	    if (!find_part) return false;
	  }
	}
      }
      return true;
    });

  NamedFunc npnh ("npnh", [](const Baby &b) -> NamedFunc::ScalarType{
      bool npnh_ = true;
      for (int imc(0); imc<b.nGenPart(); ++imc) {
        bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(imc));
	int mc_id = b.GenPart_pdgId()->at(imc);
	int mom_index = b.GenPart_genPartIdxMother()->at(imc);
        int mom_id = -1;
        if (mom_index >= 0) mom_id = b.GenPart_pdgId()->at(mom_index);
	if(mc_id == 22 && (mc_statusFlags[12] || (mc_statusFlags[13] && (mom_id == 22 || mom_id == -1))) && b.GenPart_pt()->at(imc)>15.0 && abs(b.GenPart_eta()->at(imc))<2.6){
	  if(mc_statusFlags[0] || mc_statusFlags[8]) npnh_ = false;
	}
      }
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


  /*NamedFunc selection_el = nel>=2 && nph>=1 && ele_mll>50 && ele_mllg >100 && ele_mllg <180 && (ele_mllg + ele_mll) > 185 && photon_pt/ele_mllg > 0.14 && stitch_cut;
  NamedFunc selection_mu = nmu>=2 && nph>=1 && muon_mll>50 && muon_mllg >100 && muon_mllg < 180 && (muon_mllg + muon_mll) > 185 && photon_pt/muon_mllg > 0.14 && stitch_cut;
  NamedFunc old_selection_el = nel>=2 && nph>=1 && ele_mll>50 && ele_mllg >100 && ele_mllg <180 && (ele_mllg + ele_mll) > 185 && photon_pt/ele_mllg > 0.14 && old_stitch_cut;
  NamedFunc old_selection_mu = nmu>=2 && nph>=1 && muon_mll>50 && muon_mllg >100 && muon_mllg < 180 && (muon_mllg + muon_mll) > 185 && photon_pt/muon_mllg > 0.14 && old_stitch_cut;
  vector<NamedFunc> selection = {selection_el, selection_mu};
  vector<NamedFunc> old_selection = {old_selection_el, old_selection_mu};
  */

  vector<NamedFunc> selection_el = {"1",
				    trigs,
				    trigs && nel>=2, 
				    trigs && nel>=2 && ph_sig_pt, 
				    trigs && nel>=2 && ph_sig_eta, 
				    trigs && nel>=2 && ph_sig_barrel,  
				    trigs && nel>=2 && ph_sig_endcap, 
				    trigs && nel>=2 && ph_sig_eleveto,
				    trigs && nel>=2 && ph_sig_eleveto && el_lead, 
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead, 
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >=1 , 
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50, 
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50 && ele_mllg_more >100 && ele_mllg_more <180,
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50 && ele_mllg_more >100 && ele_mllg_more <180 && photon_pt/ele_mllg_more > (15/110),
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50 && ele_mllg_more >100 && ele_mllg_more <180 && photon_pt/ele_mllg_more > (15/110) && (ele_mllg_more + ele_mll_more)>185,
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50 && ele_mllg_more >100 && ele_mllg_more <180 && photon_pt/ele_mllg_more > (15/110) && (ele_mllg_more + ele_mll_more)>185&& stitch_cut,
				    trigs && nel>=2 && ph_sig_eleveto && el_sublead && nph >= 1 && ele_mll_more > 50 && ele_mllg_more >100 && ele_mllg_more <180 && photon_pt/ele_mllg_more > (15/110) && (ele_mllg_more + ele_mll_more)>185&& old_stitch_cut}; 
  vector<NamedFunc> selection_mu = {"1",
				    trigs, 
				    trigs && nmu>=2, 
				    trigs && nmu>=2 && ph_sig_pt, 
				    trigs && nmu>=2 && ph_sig_eta, 
				    trigs && nmu>=2 && ph_sig_barrel,  
				    trigs && nmu>=2 && ph_sig_endcap, 
				    trigs && nmu>=2 && ph_sig_eleveto,
                                    trigs && nmu>=2 && ph_sig_eleveto && mu_lead, 
				    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead, 
				    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1,
                                    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50,
                                    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50 && muon_mllg_more >100 && muon_mllg_more <180,
                                    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50 && muon_mllg_more >100 && muon_mllg_more <180 && photon_pt/muon_mllg_more > (15/110),
				    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50 && muon_mllg_more >100 && muon_mllg_more <180 && photon_pt/muon_mllg_more > (15/110) && (muon_mllg_more + muon_mll_more)>185,
				    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50 && muon_mllg_more >100 && muon_mllg_more <180 && photon_pt/muon_mllg_more > (15/110) && (muon_mllg_more + muon_mll_more)>185 && stitch_cut,
				    trigs && nmu>=2 && ph_sig_eleveto && mu_sublead && nph >= 1 && muon_mll_more > 50 && muon_mllg_more >100 && muon_mllg_more <180 && photon_pt/muon_mllg_more > (15/110) && (muon_mllg_more + muon_mll_more)>185 && old_stitch_cut};
  vector<vector<NamedFunc>> selection {selection_el, selection_mu};


  auto proc_dy_nos= Process::MakeShared<Baby_nano>("DY stitch",back,TColor::GetColor("#ffb400"),{mc_path+"*DYJetsToLL_M-50*"}, "1");
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
  tag = "Nano_";
  vector<string> lepton ={"ele_", "muon_"};
  for (int i(0); i<2; i++){
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
    pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "stitch_true_" + lepton[i]);

    pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_pdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_mpdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), old_selection.at(i), stitch, ops2).Weight(weight_1).Tag(tag + "oldstitch_true_" + lepton[i]);

    pm.Push<Hist1D>(Axis(41,-1,40, pdgid,    "Photon PDG ID"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_pdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(42,-2,40, mother_pdgid,    "Photon Mother PDG ID"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_mpdg_" + lepton[i]);
    pm.Push<Hist1D>(Axis(10,-0.1,0.4, photon_truth_dR,    "Photon and true photon dR"      ,{}), selection.at(i) && npnh, stitch, ops2).Weight(weight_1).Tag(tag + "stitch_npnh_true_" + lepton[i]);
  */
    pm.Push<Table>("Stitch_cut_flow0-4" + lepton[i], vector<TableRow>{
      TableRow("Total", selection[i][0],0,0,weight_1),
	TableRow("Triggers", selection[i][1],0,0,weight_1),
	TableRow("n lepton>=2", selection[i][2],0,0,weight_1),
	TableRow("photon Pt > 15", selection[i][3],0,0,weight_1),
        TableRow("Eta < 2.5", selection[i][4],0,0,weight_1),
        TableRow("Barrel", selection[i][5],0,0,weight_1),
	TableRow("Endcap", selection[i][6],0,0,weight_1),
        TableRow("Electron veto", selection[i][7],0,0,weight_1),
        TableRow("Leading > 25(15)", selection[i][8],0,0,weight_1),
        TableRow("Subleading > 15(10)", selection[i][9],0,0,weight_1),
        TableRow("dR > 0.4", selection[i][10],0,0,weight_1),
        TableRow("m_ll > 50", selection[i][11],0,0,weight_1),
	TableRow("100 < m_llg < 180", selection[i][12],0,0,weight_1),
        TableRow("gamma pt/m_llg > 15/110", selection[i][13],0,0,weight_1),
	TableRow("m_ll + m_llg >185", selection[i][14],0,0,weight_1),        
	TableRow("New stitch", selection[i][15],0,0,weight_1),
       	TableRow("Old Stitch", selection[i][16],0,0,weight_1)}, stitch, false);

  }
  pm.min_print_ = true;
  pm.MakePlots(1);
}
