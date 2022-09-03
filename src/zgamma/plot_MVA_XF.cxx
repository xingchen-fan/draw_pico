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

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  //Define the data samples, i.e. signal and background
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  //Process::Type data =  Process::Type::data;

  // /net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v1/2017/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_153.root 

  //Define the location of the (MC) data smaples
  //string bfolder("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/");
  //string mc_path( bfolder+"mc/merged_zgmc_llg/");
  //string sig_path(bfolder+"HToZG/merged_zgmc_llg/");

  //string bfolder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v1/2017/");
  //string mc_path( bfolder+"mc/merged_zgmc_llg/");
  //string sig_path(bfolder+"HToZG/merged_zgmc_llg/");

  string mcpath("/net/cms37/data1/xf82/BDT/OutputV5/");

  // DY_output.rootGluGlu_output.root  SMZg_output.root

  //Define the triggers
  //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  //NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  //NamedFunc trigs(el_trigs || mu_trigs);

  NamedFunc trigs("BDT_score<10.0");

  //Define the data samples
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back, 
                       TColor::GetColor("#16bac5"),{mcpath+"SMZg_output.root"}, trigs);
  //auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma",      back, 
  //                     TColor::GetColor("#39a9df"),{mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back, 
                       TColor::GetColor("#ffb400"),{mcpath+"DY_output.root"},  trigs && "stitch_dy");
  //auto proc_dy_nos= Process::MakeShared<Baby_pico>("DY",               back, 
  //                     TColor::GetColor("#ffb400"),{mc_path+"*DYJets*"},  trigs);
  //auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",            back, 
  //                     TColor::GetColor("#ED702D"),{mc_path+"*TT_Tune*"}, trigs);
  //auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZ#gamma(x1000)", sig, 
  //                     kRed     ,{sig_path+"*.root"},   trigs);
  //auto proc_hzg_vbf = Process::MakeShared<Baby_pico>("HToZ#gamma VBF(x1000)", sig, 
  //                        kMagenta                   ,{sig_path+"*VBF*.root"},   trigs);

  auto proc_hzg_gg  = Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, 
                       TColor::GetColor("#ff0000"),{mcpath+"GluGlu_output.root"},   trigs);

  //Define line widhts for plotting?
  proc_smzg->SetLineWidth (1); //proc_ttg->SetLineWidth  (1);
  proc_dy->SetLineWidth   (1); //proc_ewkzg->SetLineWidth(1);
  //proc_hzg->SetLineWidth(3);
  //proc_dy_nos->SetLineWidth(1);
  proc_hzg_gg->SetLineWidth(3);

  //vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> theory_procs = {proc_dy, proc_smzg, proc_hzg_gg};

  //Set plot options
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
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);

  //vector<PlotOpt> ops = {lin_stack, lin_lumi};
  vector<PlotOpt> ops = {lin_stack};

  //Functions to plot
  NamedFunc MVAcut("BDT_score>-0.05");

  NamedFunc MVAcut1("BDT_score>0.05");

  NamedFunc MVAcut2("BDT_score>-0.05 && BDT_score<= 0.05");

  NamedFunc MVAcut3("BDT_score<=-0.05");


  //Some cuts
  vector<NamedFunc> loose_lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 15 && el_pt[ll_i2[0]] > 15",
                                 "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 15 && mu_pt[ll_i2[0]] > 15"};
  //vector<NamedFunc> loose_lep = {"ll_lepid[0] == 11 && el_pt[0] > 15 && el_pt[1] > 15",
  //                               "ll_lepid[0] == 13 && mu_pt[0] > 15 && mu_pt[1] > 15"};
  //NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_drmin[0] > 0.4");
  NamedFunc baseline("nphoton > 0 && ll_m[0] > 50");

  NamedFunc signalregion("llphoton_m[0]>122.3 && llphoton_m[0]<128.3");


  //NamedFunc MVAcut("MVA > 0.075");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
  //vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[0] > 25 && el_pt[1] > 15",
  //                         "ll_lepid[0] == 13 && mu_pt[0] > 20 && mu_pt[1] > 10"};


  //Not following what this is doing..
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 1000*weight;
    return        weight;
  });

  NamedFunc weight("weight", [](const Baby &b)->NamedFunc::ScalarType{
      float w = b.weight();
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("HToZG")) w *= 1000;
        if(file.Contains("Glu")) w *= 1000;
      }
      return w;
    });

  PlotMaker pm;
  for(int ilep(0); ilep < 2; ilep++) {
     NamedFunc dummy = "photon_pt[0] > -1.0";
     NamedFunc cut = lep.at(ilep);
     //NamedFunc cut = lep.at(ilep) && baseline && llphoton_cuts;
     NamedFunc cutMVA1 = lep.at(ilep) && baseline && llphoton_cuts && MVAcut1;
     NamedFunc cutMVA2 = lep.at(ilep) && baseline && llphoton_cuts && MVAcut2;
     NamedFunc cutMVA3 = lep.at(ilep) && baseline && llphoton_cuts && MVAcut3;
     NamedFunc Pt1 = lep.at(ilep) && baseline && llphoton_cuts && "photon_pt[0]>=15.0 && photon_pt[0]<20.0";
     NamedFunc Pt2 = lep.at(ilep) && baseline && llphoton_cuts && "photon_pt[0]>=20.0 && photon_pt[0]<30.0";
     NamedFunc Pt3 = lep.at(ilep) && baseline && llphoton_cuts && "photon_pt[0]>=30.0 && photon_pt[0]<45.0";
     NamedFunc Pt4 = lep.at(ilep) && baseline && llphoton_cuts && "photon_pt[0]>=45.0 && photon_pt[0]<200.0";
     NamedFunc Pt1MVA = lep.at(ilep) && baseline && llphoton_cuts && MVAcut && "photon_pt[0]>=15.0 && photon_pt[0]<20.0";
     NamedFunc Pt2MVA = lep.at(ilep) && baseline && llphoton_cuts && MVAcut && "photon_pt[0]>=20.0 && photon_pt[0]<30.0";
     NamedFunc Pt3MVA = lep.at(ilep) && baseline && llphoton_cuts && MVAcut && "photon_pt[0]>=30.0 && photon_pt[0]<45.0";
     NamedFunc Pt4MVA = lep.at(ilep) && baseline && llphoton_cuts && MVAcut && "photon_pt[0]>=45.0 && photon_pt[0]<200.0";
     NamedFunc Pt1Sig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && "photon_pt[0]>=15.0 && photon_pt[0]<20.0";
     NamedFunc Pt2Sig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && "photon_pt[0]>=20.0 && photon_pt[0]<30.0";
     NamedFunc Pt3Sig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && "photon_pt[0]>=30.0 && photon_pt[0]<45.0";
     NamedFunc Pt4Sig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && "photon_pt[0]>=45.0 && photon_pt[0]<200.0";
     NamedFunc Pt1MVASig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && MVAcut && "photon_pt[0]>=15.0 && photon_pt[0]<20.0";
     NamedFunc Pt2MVASig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && MVAcut && "photon_pt[0]>=20.0 && photon_pt[0]<30.0";
     NamedFunc Pt3MVASig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && MVAcut && "photon_pt[0]>=30.0 && photon_pt[0]<45.0";
     NamedFunc Pt4MVASig = lep.at(ilep) && baseline && llphoton_cuts && signalregion && MVAcut && "photon_pt[0]>=45.0 && photon_pt[0]<200.0";
     string tag = "";
     if(ilep == 0) tag += "_el";
     else       tag += "_mu";

     pm.Push<Hist1D>(Axis(20, -0.5, 0.5, "BDT_score",         "BDT_score"             ,{}), cut, theory_procs, ops).Weight(weight).Tag("_BDT_score");
/*
     pm.Push<Hist1D>(Axis(45, 15, 60, "photon_pt[0]",         "photon_pt"             ,{}), cut, theory_procs, ops).Weight(wgt).Tag("_photon_pt");
     pm.Push<Hist1D>(Axis(45, 15, 60, "photon_pt[0]",         "photon_pt MVA signal"             ,{}), cutMVA1, theory_procs, ops).Weight(wgt).Tag("_photon_pt_HiMVA");
     pm.Push<Hist1D>(Axis(45, 15, 60, "photon_pt[0]",         "photon_pt MVA medium"             ,{}), cutMVA2, theory_procs, ops).Weight(wgt).Tag("_photon_pt_MedMVA");
     pm.Push<Hist1D>(Axis(45, 15, 60, "photon_pt[0]",         "photon_pt MVA background"         ,{}), cutMVA3, theory_procs, ops).Weight(wgt).Tag("_photon_pt_LowMVA");
     pm.Push<Hist1D>(Axis(45, 15, 60, "photon_pt[0]",         "photon_p"             ,{}), baseline, theory_procs, ops).Weight(wgt).Tag("_photon_pt_baseline");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt1"             ,{}), Pt1, theory_procs, ops).Weight(wgt).Tag("mllg_Pt1");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt2"             ,{}), Pt2, theory_procs, ops).Weight(wgt).Tag("mllg_Pt2");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt3"             ,{}), Pt3, theory_procs, ops).Weight(wgt).Tag("mllg_Pt3");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt4"             ,{}), Pt4, theory_procs, ops).Weight(wgt).Tag("mllg_Pt4");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt1Sig"             ,{}), Pt1Sig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt1Sig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt2Sig"             ,{}), Pt2Sig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt2Sig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt3Sig"             ,{}), Pt3Sig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt3Sig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt4Sig"             ,{}), Pt4Sig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt4Sig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt1MVA"             ,{}), Pt1MVA, theory_procs, ops).Weight(wgt).Tag("mllg_Pt1MVA");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt2MVA"             ,{}), Pt2MVA, theory_procs, ops).Weight(wgt).Tag("mllg_Pt2MVA");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt3MVA"             ,{}), Pt3MVA, theory_procs, ops).Weight(wgt).Tag("mllg_Pt3MVA");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt4MVA"             ,{}), Pt4MVA, theory_procs, ops).Weight(wgt).Tag("mllg_Pt4MVA");  
*/
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt1MVASig"             ,{}), dummy, theory_procs, ops).Weight(weight).Tag("FixName:mllg_Pt1MVASig"+tag); 
/* 
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt2MVASig"             ,{}), Pt2MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt2MVASig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt3MVASig"             ,{}), Pt3MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt3MVASig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt4MVASig"             ,{}), Pt4MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_Pt4MVASig");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt1MVASig"             ,{}), Pt1MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_P1MS");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt2MVASig"             ,{}), Pt2MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_P2MS");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt3MVASig"             ,{}), Pt3MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_P3MS");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llphoton_m Pt4MVASig"             ,{}), Pt4MVASig, theory_procs, ops).Weight(wgt).Tag("mllg_P4MS");  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llg mass MVA signal"             ,{}), cutMVA1, theory_procs, ops).Weight(wgt).Tag(tag);  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llg mass MVA medium"             ,{}), cutMVA2, theory_procs, ops).Weight(wgt).Tag(tag);  
     pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]",         "llg mass MVA background"             ,{}), cutMVA3, theory_procs, ops).Weight(wgt).Tag(tag);  
     //pm.Push<Hist1D>(Axis(40, 0, 200, "llphoton_pt[0]",         "llg pT"             ,{}), cut, theory_procs, ops).Weight(wgt).Tag(tag);  
     */
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}

