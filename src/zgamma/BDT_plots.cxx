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
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  string folder("/net/cms37/data1/xf82/BDT/Output/");

  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma",       back,  TColor::GetColor("#16bac5"),{folder+"*SMZg*"}, "1");
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY",               back, TColor::GetColor("#ffb400"),{folder+"*DY*"}, "stitch_dy");
  auto proc_hzg_zh   = Process::MakeShared<Baby_pico>("HToZ#gamma(x100)", sig, kRed,{folder+"*ZH*"},   "1");
  auto proc_hzg_ggf   = Process::MakeShared<Baby_pico>("HToZ#gamma(ggf)(x100)", sig, kRed,{folder+"*GGF*"},   "1");

  proc_smzg->SetLineWidth (1);
  proc_dy->SetLineWidth   (1);
  proc_hzg_ggf->SetLineWidth   (3);


  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg};
  //  vector<shared_ptr<Process>> sigprocs = {proc_hzg};
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

  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_drmin[0] > 0.4");
  NamedFunc baseline("nphoton > 0 && ll_m[0] > 50");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[0] > 25 && el_pt[1] > 15",
                           "ll_lepid[0] == 13 && mu_pt[0] > 20 && mu_pt[1] > 10"};


  NamedFunc wgt("wgt", [](const Baby &b)->NamedFunc::ScalarType{
      float w = b.weight();
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
	if(file.Contains("HToZG")) w *= 100;
	if(file.Contains("GGF")) w *= 100;
      }
      return w;
    });

  NamedFunc wgt1("wgt1", [](const Baby &b)->NamedFunc::ScalarType{
      float w = 1.;
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("HToZG")) w *= 100;
      }
      return w;
    });

  PlotMaker pm;
  for(int i(0); i < 1; i++) {
    NamedFunc selection = lep.at(i) && baseline && llphoton_cuts;
    //    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && "BDT_score > -0.05" && "photon_pt[0] > 15." && "photon_pt[0] < 30.", procs, ops).Weight(wgt).Tag("pt15_30_");
    //    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && "BDT_score > -0.05" && "photon_pt[0] > 30." && "photon_pt[0] < 50.", procs, ops).Weight(wgt).Tag("pt30_50_");
    //    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && "BDT_score > -0.05" && "photon_pt[0] > 50." && "photon_pt[0] < 80.", procs, ops).Weight(wgt).Tag("pt_50_80_");
    pm.Push<Hist1D>(Axis(75,0.25,1, "BDT_score", "BDT",{}), selection, procs, ops).Weight(wgt);
  }
  pm.min_print_ = true;
  pm.MakePlots(1);
}


///net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/
