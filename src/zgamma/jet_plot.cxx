#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <cmath>
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

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
    .YAxis(YAxisType::log)
    .Stack(StackType::signal_overlay)
    .YTitleOffset(1.)
    .AutoYAxis(false)
    .UseCMYK(false)
    .LeftMargin(0.17)
    .LegendColumns(1)
    .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Overflow(OverflowType::none);
  vector<PlotOpt> ops = {lin_lumi};

  NamedFunc CR_photon_size("CR_photon_size", [](const Baby &b)->NamedFunc::ScalarType{
      return b.photon_pt()->size();
    });
  //vector<NamedFunc> lep = {"nel > 1 && ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25. && el_pt[ll_i2[0]] > 15. && (el_eta[ll_i1[0]] < 2.4 && el_eta[ll_i1[0]] > -2.4) && (el_eta[ll_i2[0]] < 2.4 && el_eta[ll_i2[0]] > -2.4)",
  //                   "nmu > 1 && ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20. && mu_pt[ll_i2[0]] > 10. && (mu_eta[ll_i1[0]] < 2.5 && mu_eta[ll_i1[0]] > -2.5) && (mu_eta[ll_i2[0]] < 2.5 && mu_eta[ll_i2[0]] > -2.5)"};
  
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                         "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
  
  //NamedFunc qualityCut = "(photon_idmva[0] > -0.4 && (photon_eta[0] < 1.4442 && photon_eta[0] >-1.4442) ) || (photon_idmva[0] > -0.58 && (photon_eta[0] < -1.566 || photon_eta[0] > 1.566)  && (photon_eta[0] < 2.5&& photon_eta[0] > -2.5 ))";
  //NamedFunc el_baseline = lep[0] && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.";
  //NamedFunc mu_baseline = lep[1] && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.";
  NamedFunc el_baseline = lep[0] && "nphoton > 0" && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.";
  NamedFunc mu_baseline = lep[1] && "nphoton > 0" && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.";
  vector<NamedFunc> baseline = {el_baseline, mu_baseline};


  NamedFunc w_l("w_l", [](const Baby &b)->NamedFunc::ScalarType{
      float w = b.w_lumi();
      float w_year = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
      } else if (b.SampleType()==2017) {
        w_year = 41.52756;
      }
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("GluGluHToZG")) w *= 1000;
        if(file.Contains("ttHToZG")) w *= 10000;
	if(file.Contains("VBFHToZG")) w *= 10000;
        if(file.Contains("WminusH")) w *= 10000;
	if(file.Contains("WplusH")) w *= 10000;
        if(file.Contains("ZH_HToZG")) w *= 10000;

      }
      return w * w_year;
    });

  NamedFunc mjj("mjj", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector j1, j2, jj;
      j1.SetPtEtaPhiM(b.jet_pt()->at(0), b.jet_eta()->at(0), b.jet_phi()->at(0), b.jet_m()->at(0));
      j2.SetPtEtaPhiM(b.jet_pt()->at(1), b.jet_eta()->at(1), b.jet_phi()->at(1), b.jet_m()->at(1));
      jj = j1 + j2;
      return jj.M();
    });

  NamedFunc ptjj("ptjj", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector j1, j2, jj;
      j1.SetPtEtaPhiM(b.jet_pt()->at(0), b.jet_eta()->at(0), b.jet_phi()->at(0), b.jet_m()->at(0));
      j2.SetPtEtaPhiM(b.jet_pt()->at(1), b.jet_eta()->at(1), b.jet_phi()->at(1), b.jet_m()->at(1));
      jj = j1 + j2;
      return jj.Pt();
    });
  
  NamedFunc el_trigs17("el_trigs17", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();});
  NamedFunc mu_trigs17("mu_trigs17", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8() || b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();});
  NamedFunc trigs(el_trigs17 || mu_trigs17);
  
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  
  string folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/2017/mc/skim_llg/");
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY + Fake",               back, TColor::GetColor("#ffb400"),{folder+"*DYJetsToLL_M-50*"}, "stitch_dy" &&  trigs);
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back,  TColor::GetColor("#16bac5"),{folder+"*ZGToLLG*"}, "!stitch_dy" && trigs);
  auto proc_ewk  = Process::MakeShared<Baby_pico>("EWK",       back,  6,{folder+"*EWKZ2Jets*"}, trigs);

  auto proc_ggf    = Process::MakeShared<Baby_pico>("GGF X 1000",               sig, 600,{folder+"*GluGluHToZG_ZToLL_M-125*"},  trigs);
  auto proc_vbf  = Process::MakeShared<Baby_pico>("VBF X 10000",       sig,  632,{folder+"*VBFHToZG_ZToLL_M-125*"}, trigs);
  auto proc_vh  = Process::MakeShared<Baby_pico>("VH X 10000",       sig,  1,{folder+"*WminusH_HToZG_WToAll_M-125*", folder+"*WplusH_HToZG_WToAll_M-125*",folder+"*ZH_HToZG_ZToAll_M-125*"}, trigs);
  auto proc_tth  = Process::MakeShared<Baby_pico>("ttH X 10000",       sig,  400,{folder+"*ttHToZG_M125*"}, trigs);

  proc_dy -> SetLineWidth (1);
  proc_smzg -> SetLineWidth (1);
  proc_ewk-> SetLineWidth (1);
  proc_ggf -> SetLineWidth (3);
  proc_vbf -> SetLineWidth (3);
  proc_vh -> SetLineWidth (3);
  proc_tth -> SetLineWidth (3);

  vector<shared_ptr<Process>> proc = {proc_dy, proc_smzg, proc_ewk, proc_ggf, proc_vbf, proc_vh, proc_tth};
  vector<shared_ptr<Process>> proc_sig = {proc_ggf, proc_vbf, proc_vh, proc_tth};
  vector<string> leptag = {"el_", "mu_"};
  

  PlotMaker pm;
  for (int i(0); i < 2; i++) {
    pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]", "m_{llg}",{}), baseline[i], proc, ops).Weight(w_l).Tag(leptag[i]);
    /*pm.Push<Hist1D>(Axis(5, 0, 5, "njet", "njet",{}), baseline[i], proc, ops).Weight(w_l).Tag("njet_" + leptag[i]);
    pm.Push<Hist1D>(Axis(180, 20, 200, mjj, "dijet mass",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("mjj_" + leptag[i]);
    pm.Push<Hist1D>(Axis(150, 0, 150, ptjj, "dijet pT",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("ptjj_" + leptag[i]);
    pm.Push<Hist1D>(Axis(65, 15, 80, "photon_pt[0]", "Photon pT",{}), baseline[i], proc, ops).Weight(w_l).Tag("photonpT_" + leptag[i]);
    pm.Push<Hist1D>(Axis(130, 20, 150, "jet_pt[2]", "Jet3 pT",{}), baseline[i] && "njet > 2", proc, ops).Weight(w_l).Tag("j3pt_" + leptag[i]);
    pm.Push<Hist1D>(Axis(80, 0.4, 3, "photon_drmin[0]", "mindR",{}), baseline[i], proc, ops).Weight(w_l).Tag("mindR_" + leptag[i]);

    pm.Push<Hist1D>(Axis(130, 20, 150, "jet_pt[0]", "Jet1 pT",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1pt_" + leptag[i]);
    pm.Push<Hist1D>(Axis(130, 20, 150, "jet_pt[1]", "Jet2 pT",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2pt_" + leptag[i]);
    pm.Push<Hist1D>(Axis(52, -5, 5, "jet_eta[0]", "Jet1 eta",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1eta_" + leptag[i]);
    pm.Push<Hist1D>(Axis(52, -5, 5, "jet_eta[1]", "Jet2 eta",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2eta_" + leptag[i]);
    pm.Push<Hist1D>(Axis(52, -5, 5, "jet_eta[2]", "Jet3 eta",{}), baseline[i] && "njet > 2", proc, ops).Weight(w_l).Tag("j3eta_" + leptag[i]);
    
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_isphoton[0]", "Jet1 isphoton",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1isph_" + leptag[i]);
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_isphoton[1]", "Jet2 isphoton",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2isph_" + leptag[i]);
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_isgood[0]", "Jet1 isgood",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1isgood_" + leptag[i]);
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_isgood[1]", "Jet2 isgood",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2isgood_" + leptag[i]);
    */
    //    pm.Push<Hist1D>(Axis(130, 20, 150, "jet_pt[0]", "Jet1 pT (njet > 0)",{}), baseline[i] && "njet > 0", proc, ops).Weight(w_l).Tag("j1pt_loose_" + leptag[i]);
    //    pm.Push<Hist1D>(Axis(52, -5, 5, "jet_eta[0]", "Jet1 eta (njet > 0)",{}), baseline[i] && "njet > 0", proc, ops).Weight(w_l).Tag("j1eta_loose_" + leptag[i]);

    /*
    // pm.Push<Hist1D>(Axis(100, -50, 50, "jet_csv[0]", "Jet1 csv",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1csv_" + leptag[i]);
    //pm.Push<Hist1D>(Axis(100, -50, 50, "jet_csv[1]", "Jet2 csv",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2csv_" + leptag[i]);
    pm.Push<Hist1D>(Axis(4, -2, 2, "jet_deepcsv[0]", "Jet1 deepcsv",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1deepcsv_" + leptag[i]);
    pm.Push<Hist1D>(Axis(4, -2, 2, "jet_deepcsv[1]", "Jet2 deepcsv",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2deepcsv_" + leptag[i]);
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_deepflav[0]", "Jet1 deepflav",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j1deepflav_" + leptag[i]);
    pm.Push<Hist1D>(Axis(2, 0, 2, "jet_deepflav[1]", "Jet2 deepflav",{}), baseline[i] && "njet > 1", proc, ops).Weight(w_l).Tag("j2deepflav_" + leptag[i]);
    */


  }

  pm.min_print_ = true;
  pm.MakePlots(1);
  
}
