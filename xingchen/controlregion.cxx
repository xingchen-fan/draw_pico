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

NamedFunc NminusOne(vector<NamedFunc> cuts, int i) {
  NamedFunc total("1");
  for(int j(0); j < static_cast<int>(cuts.size()); j++) 
    if(j != i)
      total = total && cuts.at(j);
  return total;
}

bool checkBit(int i, int n) {
  return((i%static_cast<int>(pow(2,n+1)))/static_cast<int>(pow(2,n)));
}

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/");
  //string mc_path( bfolder+"zgamma_mc_ul/2017/mc/merged_zgmc_llg/");
  //string mc_path("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v0/2017/mc/skim_llg/");
  //string sig_path(bfolder+"zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string mc_path("/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/");
  string sig_path("/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/signal/");
  string mc2016_path ("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_llg/");
  //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  string rawpico("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/unskimmed/");
  NamedFunc trigs(el_trigs || mu_trigs);

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .YTitleOffset(1.)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Stack(StackType::signal_overlay);;
  vector<PlotOpt> ops = {lin_lumi};

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


  NamedFunc selection_el_1 = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >50 && ele_mllg_more <100 /*&& (ele_mllg_more + ele_mll_more) > 185*/ && photon_pt/ele_mllg_more >(15./110);
  NamedFunc selection_mu_1 = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >50 && muon_mllg_more < 100 /*&& (muon_mllg_more + muon_mll_more) > 185*/ && photon_pt/muon_mllg_more > (15./110);
  NamedFunc selection_el_2 = nel>=2 && nph>=1 && el_pt_cut && ele_mll_more>50 && ele_mllg_more >100 && ele_mllg_more <180 && /*(ele_mllg_more + ele_mll_more) > 185 &&*/ photon_pt/ele_mllg_more < (15./110);
  
  NamedFunc selection_mu_2 = nmu>=2 && nph>=1 && mu_pt_cut && muon_mll_more>50 && muon_mllg_more >100 && muon_mllg_more < 180 && /*(muon_mllg_more + muon_mll_more) > 185 &&*/ photon_pt/muon_mllg_more < (15./110);
  vector<NamedFunc> selection1 = {selection_el_1, selection_mu_1};
  vector<NamedFunc> selection2 = {selection_el_2, selection_mu_2};
  vector<NamedFunc> mll = {ele_mll_more, muon_mll_more};
  vector<NamedFunc> mllg = {ele_mllg_more, muon_mllg_more};

  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 1000*weight;
    else if(b.type() == 6200)
      return weight/1.664;
    return        weight;
  });

  const NamedFunc weight_1("weight_1",[](const Baby &b) -> NamedFunc::ScalarType{
      float w = 1.0;
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
	TString name = *it;
	if(name.Contains("GluGluHToZG")||name.Contains("ttHToZG")||name.Contains("VBFHToZG")||name.Contains("WminusH")||name.Contains("WplusH")||name.Contains("ZH_HToZG_ZToAll")){
	  w = 1.;
	  break;
	}
	else w = 1.0;
      }
      return w;
    });
  
  const NamedFunc weight("weight",[](const Baby &b) -> NamedFunc::ScalarType{
      float w = 1.0;
      float xsec = -999.;
      float HToZG(0.001533), ZToLL(0.100974);
      int tot = 1;
      float sign =  b.Generator_weight() > 0? 1.: -1.;
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
	if(file.Contains("ZGToLLG")) {  xsec = 1000 * 55.48; tot = 18727540;}
	if(file.Contains("DYJetsToLL_M-50_TuneCP5")) {xsec = 1000 * 6077.22; tot = 102486448;}

	if(file.Contains("GluGluHToZG"))          {xsec = 1000 * 1000 * HToZG * ZToLL * 48.58 ; tot = 395990;}
	if(file.Contains("VBFHToZG"))             {xsec = 1000 * 1000 * HToZG * ZToLL * 3.782 ; tot = 199732;}
	if(file.Contains("WplusH_HToZG"))         {xsec = 1000 * 1000 * HToZG * 0.831 ; tot = 282330;}
	if(file.Contains("WminusH_HToZG"))        {xsec = 1000 * 1000 * HToZG * 0.527 ; tot = 283432;}
	if(file.Contains("ZH_HToZG"))             {xsec = 1000 * 1000 * HToZG * 0.8839; tot = 279837;}
	if(file.Contains("ttHToZG"))              {xsec = 1000 * 1000 * HToZG * 0.5071; tot = 195782;}
	w = xsec * 41.5 * sign/tot;
	

      }
      return w;
    });

  const NamedFunc total_finder("total_finder",[](const Baby &b) -> NamedFunc::ScalarType{
      double ind = 0.;
     
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("ZGToLLG")) {ind = 1.;}
        if(file.Contains("DYJetsToLL_M-50_TuneCP5")) {ind = 2.;}

        if(file.Contains("GluGluHToZG"))         ind = 3.;
        if(file.Contains("VBFHToZG"))           ind = 4.;
        if(file.Contains("WplusH_HToZG"))        ind = 5.;
        if(file.Contains("WminusH_HToZG"))       ind = 6.;
        if(file.Contains("ZH_HToZG"))            ind = 7.;
        if(file.Contains("ZH_ZToAll_HToZG"))     ind = 8.;
        if(file.Contains("ttHToZG"))             ind = 9.;

      }
      return ind;
    });

  auto proc_smzg  = Process::MakeShared<Baby_nano>("SM Z#gamma",       back,
						   TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG_*"}, trigs && !stitch_cut);
  //auto proc_ewkzg = Process::MakeShared<Baby_nano>("EWK Z#gamma",      back,						   TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_nano>("DY + Jets",               back,
						   TColor::GetColor("#ffb400"),{mc_path+"*DYJetsToLL_M-50*"}, trigs && stitch_cut);
  //auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back,						   TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  //auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back,						   TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_nano>("HToZ#gamma (X1000)", sig,
						   kRed                      ,{sig_path+"*.root"},    trigs );
  //auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig,						     kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  //auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig,						     kBlue ,{sig_path+"*GluGlu*.root"},   trigs);
  //auto raw_proc_dy    = Process::MakeShared<Baby_pico>("DY",               back,						       TColor::GetColor("#ffb400"),{rawpico + "*DYJetsToLL*"}, "1");
  ///net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/raw_pico
  //auto new_proc_dy = Process::MakeShared<Baby_pico>("DY", back, TColor::GetColor("#ffb400"),{"/homes/abarzdukas/nano2pico_test/nano2pico/out/unskimmed/*DYJetsToLL*"}, "1");
  proc_smzg->SetLineWidth (1); //proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); //proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(4);
  //proc_dy_nos->SetLineWidth(1);
  //proc_hzg_gg->SetLineWidth(4);
  //proc_hzg_vbf->SetLineWidth(4);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> procssig = {proc_hzg};
  //vector<shared_ptr<Process>> procsbkg = {proc_dy, proc_smzg,proc_ewkzg,proc_ttg};

  //vector<shared_ptr<Process>> theory_procs = {proc_smzg, proc_dy, proc_ewkzg, proc_hzg};
  //vector<shared_ptr<Process>> stitch = {new_proc_dy};

  auto proc_smzg_total  = Process::MakeShared<Baby_nano>("SM Z#gamma",       back,
                                                   TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG_*"}, "1");
  auto proc_dy_total    = Process::MakeShared<Baby_nano>("DY + Jets",               back,
                                                   TColor::GetColor("#ffb400"),{mc_path+"*DYJetsToLL_M-50*"}, "1");
  auto proc_hzg_total   = Process::MakeShared<Baby_nano>("HToZ#gamma (X100)", sig,
                                                   kRed                      ,{sig_path+"*.root"},   "1");
  
  vector<shared_ptr<Process>> procs_total = {proc_dy_total, proc_smzg_total, proc_hzg_total};
  vector<shared_ptr<Process>> proc_sm_tot = {proc_smzg_total};
  
  PlotOpt style("txt/plot_styles.txt","Scatter");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
			      .LogMinimum(1)
			      .CanvasWidth(600)
			      .LabelSize(0.04)
			      .YAxis(YAxisType::linear)
			      .Title(TitleType::info)};
  PlotMaker pm;
  vector<string> lepton = {"electron", "muon"};
  for(int i(0); i < 1; i++) {
    //pm.Push<Hist1D>(Axis(100,50,150, mll[i], "m_{ll}",{}), selection1[i] && mll[i] > 50 && mll[i] < 150, procs, ops).Weight(weight).Tag("CR_1_mll_"+lepton[i]);
    pm.Push<Hist1D>(Axis(50,50,100, mllg[i], "m_{ll#gamma}",{}), selection1[i] && mllg[i] > 50 && mllg[i] < 100, procs, ops).Weight(weight).Tag("CR_1_mllg_"+lepton[i]);
    //pm.Push<Hist1D>(Axis(100,50,150, mll[i], "m_{ll}",{}), selection2[i] && mll[i] > 50 && mll[i] < 150, procs, ops).Weight(weight).Tag("CR_2_mll_"+lepton[i]);
    pm.Push<Hist1D>(Axis(80,100,180, mllg[i], "m_{ll#gamma}",{}), selection2[i] && mllg[i] > 100 && mllg[i] < 180, procs, ops).Weight(weight).Tag("CR_2_mllg_"+lepton[i]);
    //pm.Push<Hist1D>(Axis(150,50,200, "llphoton_m[0]", "m_{ll#gamma} [GeV]",{}), selection && zoom && !(ratio_cut), procs, ops).Weight(wgt).Tag("Ratio");
    //pm.Push<Hist2D>(Axis(115,0,115, "photon_pt[0]",       "#gamma P_{t} [GeV]",{}),
    //		    selection && zoom, procssig).Weight(wgt).Tag("sig");
    //pm.Push<Hist2D>(Axis(115,0,115, "photon_pt[0]",       "#gamma P_{t} [GeV]",{}),
    //              Axis(150,50,200, "llphoton_m[0]", "m_{#mu#mu#gamma} [GeV]",{100,80}),
    //              selection && zoom, procsbkg).Weight(wgt).Tag("back");
    //pm.Push<Hist1D>(Axis(60,100,400, "llphoton_m[0]+ll_m[0]", "m_{ll}+m_{ll#gamma} [GeV]",{185}),
    //		    selection && llcus_cut, procs, ops).Weight(wgt).Tag("Masscut");OA
    //pm.Push<Hist1D>(Axis(150,50,200, "llphoton_m[0]", "m_{ll#gamma} [GeV]",{}),
    //		    "stitch_dy" && selection && zoom, stitch, ops).Weight(weight_1).Tag("Check");
    //pm.Push<Hist1D>(Axis(100,0,1, "photon_pt[0]/llphoton_m[0]", "p_{T,#gamma}/m_{ll#gamma}",{0.136}),
    //		    selection && sum_cut, procs, ops).Weight(wgt).Tag("Sumcut");
    //pm.Push<Hist1D>(Axis(20,0,2, "photon_drmin[0]", "#DeltaR_{min}(l,#gamma)",{0.4}),
    //		    baseline && lep.at(i), procs, ops).Weight(wgt).Tag("dR");
    //pm.Push<Hist1D>(Axis(50,50,150, "met", "MissingEt [GeV]",{}),
    //                selection, procs, ops).Weight(wgt).Tag("Basecut");
  //pm.Push<Hist1D>(Axis(50,50,150, "met_tru", "MissingEtTru [GeV]",{}),
    //                    selection, procs, ops).Weight(wgt).Tag("Basecut");
    //pm.Push<Hist1D>(Axis(100,-1,1, "photon_idmva[0]", "mva",{}), "nphoton > 0", stitch, ops).Weight(weight_1).Tag("MVA");
		    
    /*  pm.Push<Table>("Pico Cutflow one file" + lepton[i]+" DY", vector<TableRow>{
	TableRow("Total","1" ,0,0,weight_1),
	  TableRow("Stitch","stitch_dy",0,0,weight_1),
	  TableRow("n lepton >= 2",lep.at(i),0,0,weight_1),
	  TableRow("nphoton >= 1","nphoton > 0",0,0,weight_1),
	  TableRow("nphoton >= 1 (mva>0.2)","nphoton > 0 && photon_idmva[0]>0.2",0,0,weight_1),
	  TableRow("nphoton >= 1 with lep","nphoton > 0 && (nel > 0 ||  nmu > 0)",0,0,weight_1),
	  TableRow("nphoton >= 1 with lep (mva>0.2)","nphoton > 0 && (nel > 0 ||  nmu > 0) && photon_idmva[0]>0.2",0,0,weight_1),
	  // TableRow("dR > 0.4","nel > 0 || nmu >0" && "nphoton >0" &&  photondr,0,0,weight_1),
	  TableRow("Basecut",selection,0,0,weight_1),
	  TableRow("Basecut and stitch", "stitch_dy" && selection,0,0,weight_1)}, stitch, false);
    */
    /*pm.Push<Table>("Pico Cutflow test" + lepton[i]+" DY", vector<TableRow>{
	    TableRow("Total","1" ,0,0,weight_1),
	      TableRow("n lepton >= 2",lep.at(i),0,0,weight_1),
	      //TableRow("nphoton >= 1","nphoton > 0",0,0,weight_1),
	      TableRow("nphoton >= 1 with lep","nphoton > 0 && photon_idmva[0]>0.2" && lep.at(i),0,0,weight_1),
	      // TableRow("dR > 0.4","nel > 0 || nmu >0" && "nphoton >0" &&  photondr,0,0,weight_1),
	      TableRow("Basecut",selection,0,0,weight_1),
	      TableRow("Basecut and stitch", "stitch_dy" && selection,0,0,weight_1)}, stitch, false);
    */
  }
  //pm.Push<Hist1D>(Axis(400,-200,200, "Generator_weight", "Generator weight",{}), "1", proc_sm_tot, ops).Weight(weight_1).Tag("gw_");
  /*pm.Push<Table>("Total event number", vector<TableRow>{
      TableRow("SM Z#gamma",total_finder < 1.5 && total_finder > 0.5 ,0,0,weight_1),
	TableRow("DY + jets",total_finder < 2.5 && total_finder > 1.5,0,0,weight_1),
	TableRow("GluGlu",total_finder < 3.5 && total_finder > 2.5,0,0,weight_1),
	TableRow("VBF",total_finder < 4.5 && total_finder > 3.5,0,0,weight_1),
        TableRow("W+H",total_finder < 5.5 && total_finder > 4.5,0,0,weight_1),
	TableRow("W-H",total_finder < 6.5 && total_finder > 5.5,0,0,weight_1),
        TableRow("ZH_HToZG",total_finder < 7.5 && total_finder > 6.5,0,0,weight_1),
	TableRow("ZH_HToAll",total_finder < 8.5 && total_finder > 7.5,0,0,weight_1),
	TableRow("ttHToZG",total_finder < 9.5 && total_finder > 8.5 ,0,0,weight_1)}, procs_total, false);
  */

  pm.min_print_ = true;
  pm.MakePlots(1);
}

