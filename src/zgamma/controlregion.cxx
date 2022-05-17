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
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist2d.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
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
  string mc_path("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v0/2017/mc/skim_llg/");
  string sig_path(bfolder+"zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string mc2016_path ("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  string rawpico("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/unskimmed/");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back, 
                       TColor::GetColor("#16bac5"),{mc_path+"*ZGToLLG*"}, trigs);
  auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back, 
                       TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"}, trigs && "stitch_dy");
  auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back, 
                       TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", sig, 
                       kRed                      ,{sig_path+"*.root"},   trigs);
  auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig, 
                       kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);
  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, 
                       kBlue ,{sig_path+"*GluGlu*.root"},   trigs);
  auto raw_proc_dy    = Process::MakeShared<Baby_pico>("DY",               back,
  						   TColor::GetColor("#ffb400"),{rawpico + "*DYJetsToLL*"}, "1");
///net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/raw_pico
  auto new_proc_dy = Process::MakeShared<Baby_pico>("DY", back, TColor::GetColor("#ffb400"),{"/homes/abarzdukas/nano2pico_test/nano2pico/out/unskimmed/*DYJetsToLL*"}, "1");
  proc_smzg->SetLineWidth (1); //proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); //proc_ewkzg->SetLineWidth(1);
  proc_hzg->SetLineWidth(4);
  proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(4);
  proc_hzg_vbf->SetLineWidth(4);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg,proc_ewkzg,proc_ttg, proc_hzg, proc_hzg_vbf, proc_hzg_gg};
  vector<shared_ptr<Process>> procssig = {proc_hzg};
  vector<shared_ptr<Process>> procsbkg = {proc_dy, proc_smzg,proc_ewkzg,proc_ttg};

  //vector<shared_ptr<Process>> theory_procs = {proc_smzg, proc_dy, proc_ewkzg, proc_hzg};
  vector<shared_ptr<Process>> stitch = {new_proc_dy};
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
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  NamedFunc photondr("photon_drmin[0] > 0.4");
  NamedFunc ratio_cut ("photon_pt[0]/llphoton_m[0] >= 15./110");
  NamedFunc sum_cut("llphoton_m[0]+ll_m[0]>185");
  NamedFunc llcus_cut(" ll_m[0] > 81.2 && ll_m[0] < 101.2");
  vector<NamedFunc> sig_sel = {"llphoton_m[0]+ll_m[0]>185",
                               "llphoton_m[0] > 100 && llphoton_m[0] < 180",
                               "photon_pt[0]/llphoton_m[0] >= 15./110",
                               "photon_drmin[0] > 0.4"};
  NamedFunc baseline("nll >=1 && nphoton > 0 && ll_m[0] > 50 && llphoton_m[0] > 100 && llphoton_m[0]<180");
  NamedFunc zoom( "llphoton_m[0] <= 200");
  vector<NamedFunc> lep = {"nel >=2",
                           "nmu >=2"};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 1000*weight;
    else if(b.type() == 6200)
      return weight/1.664;
    return        weight;
  });
  NamedFunc weight_1("weight_1",[](const Baby &b) -> NamedFunc::ScalarType{
      double weight = 1.0;
      if(b.type() >= 200000 && b.type() <= 205000)
	return 1000*weight;
      else return weight;
    });
  PlotOpt style("txt/plot_styles.txt","Scatter");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
			      .LogMinimum(1)
			      .CanvasWidth(600)
			      .LabelSize(0.04)
			      .YAxis(YAxisType::linear)
			      .Title(TitleType::info)};
  PlotMaker pm;
  vector<string> lepton = {"electron", "muon"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection =  baseline && lep.at(i)&& photondr;
    pm.Push<Hist1D>(Axis(180,-0.8,1, "photon_idmva[0]", "MVA",{}), selection, procsbkg, ops).Weight(wgt).Tag("mva_"+lepton[i]);
    pm.Push<Hist1D>(Axis(100,50,150, "ll_m[0]", "m_{ll}",{}), selection && " photon_idmva[0]<=0.2 && photon_idmva[0]>=-0.4", procsbkg, ops).Weight(wgt).Tag("cr_"+lepton[i]);
    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && " photon_idmva[0]<=0.2 && photon_idmva[0]>=-0.4", procsbkg, ops).Weight(wgt).Tag("cr_"+lepton[i]);
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
  pm.min_print_ = true;
  pm.MakePlots(1);
}

