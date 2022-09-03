
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
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
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

TLorentzVector AssignLead(const Baby &b, bool gen) {
  TLorentzVector lep;
  if(gen) {
    for(size_t i = 0; i < b.mc_id()->size(); i++)
      if(b.mc_id()->at(i) == -11 ||
	 b.mc_id()->at(i) == -13 ||
	 b.mc_id()->at(i) == -15 ||
	 b.mc_id()->at(i) == 11 ||
	 b.mc_id()->at(i) == 13 ||
	 b.mc_id()->at(i) == 15) {
	lep.SetPtEtaPhiM(b.mc_pt()->at(i),
			b.mc_eta()->at(i),
			b.mc_phi()->at(i),
			b.mc_mass()->at(i));
	if(b.mc_mom()->at(i) == 23) break;
      }
  }
  else {
    int il = -1;
    if(b.ll_lepid()->at(0)==11){
      if((b.el_pt()->at(b.ll_i1()->at(0))) > (b.el_pt()->at(b.ll_i2()->at(0)))) il = b.ll_i1()->at(0);
      else il = b.ll_i2()->at(0);
      lep.SetPtEtaPhiM(b.el_pt() ->at(il),
		       b.el_eta()->at(il),
		       b.el_phi()->at(il), 0.00511);
    }
    else if(b.ll_lepid()->at(0) == 13){
      if((b.mu_pt()->at(b.ll_i1()->at(0))) > (b.mu_pt()->at(b.ll_i2()->at(0)))) il = b.ll_i1()->at(0);
      else il = b.ll_i2()->at(0);
      lep.SetPtEtaPhiM(b.mu_pt() ->at(il),
                       b.mu_eta()->at(il),
                       b.mu_phi()->at(il), 0.105);
    }
  }
    return lep;
}


int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  //   Process::Type data =  Process::Type::data;

  //ROOT file path
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/");
  string sfolder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/");
  string mc_path( bfolder+"merged_zgmc_llg/");
  string sig_path(sfolder+"merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  NamedFunc trigs(el_trigs || mu_trigs);


  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",back,TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, trigs);
  //auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back,
  //                                               TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back,
                                                   TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back,
                                                   TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  //auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back,
  //                                               TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", sig,
                                                   kRed     ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x10000)", sig,
                          kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig,
                                                     kBlue ,{sig_path+"*GluGlu*.root"},   trigs);



  proc_smzg->SetLineWidth (1); //proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); //proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg_vbf, proc_hzg_gg };
  //vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_gg};

  NamedFunc lnegpt("lnegpt",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector l1;
      l1 = AssignL1(b);
      return abs(l1.Pt());
    });

  NamedFunc lpospt("lpospt",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector l2;
      l2 = AssignL2(b);
      return abs(l2.Pt());
    });
 NamedFunc leadlep("leadlep",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector lep;
      lep = AssignLead(b,false);
      return abs(lep.Pt());
    });



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
    PlotOpt lin_stack = lin_lumi().Stack(StackType::shapes);
    PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
    vector<PlotOpt> ops1 = {log_stack};
    vector<PlotOpt> ops2 = {lin_stack};

    vector<NamedFunc> loose_lep = {"ll_lepid[0] == 11","ll_lepid[0] == 13"};
    NamedFunc baseline("nll >= 1 && nphoton > 0 && ll_m[0] > 50 && llphoton_m[0] > 100 && llphoton_m[0] < 180");// && ll_m[0] > 50");

    //Weight
    NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
	double weight = b.w_lumi();
	if(b.type() == 200100)
	  return 10000*weight;
	else if(b.type() >= 200000 && b.type() <= 205000)
	  return 1000*weight;

	//     else if(b.type() == 6200)
	//       return weight/1.664;
	return        weight;
      });

    //Set up plot maker
    PlotMaker pm;
    string tag;

    for(int i(0); i < 1; i++) {

      NamedFunc cut = baseline && loose_lep.at(i);
      vector<PlotOpt> ops;
      tag = "Pt";
      if(i == 0) {tag += "_el";
	ops = ops2;}
      else       {tag += "_mu";
	ops = ops2;}
      tag += "_baseline";

      pm.Push<Hist1D>(Axis(100,15,200, leadlep,    "Leading lepton Pt"      ,{}), "ll_lepid[0] > 0", procs, ops).Weight(wgt).Tag\
	(tag);
      //pm.Push<Hist1D>(Axis(100,15,200, lpospt,    "Positive l pt"      ,{}), cut, procs, ops).Weight(wgt).Tag	(tag);
    }
    pm.min_print_ = true;
    pm.MakePlots(1);
}
