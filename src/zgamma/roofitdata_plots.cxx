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
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

int main() {
    gErrorIgnoreLevel = 6000;
    Palette colors("txt/colors.txt","default");
    Process::Type back =  Process::Type::background;
    //    Process::Type sig =  Process::Type::signal;
    string mc_path( "/net/cms37/data1/xf82/BDT/OutputNewMVA/");
    //string mc_path( "/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v1/2017/mc/raw_pico/");
    string sig_path("/net/cms37/data1/xf82/BDT/OutputNewMVA/");
    //merged_pico_llg_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_45.root
    //DYJetsToLL_M-50
    //NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
    //NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
    //NamedFunc trigs(el_trigs || mu_trigs);
    NamedFunc B_s("BDT_score");

    auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma", back, TColor::GetColor("#16bac5"),{mc_path+"*SM*"}, "!stitch_dy");
    auto proc_dy    = Process::MakeShared<Baby_pico>("DY", back,TColor::GetColor("#ffb400"),{mc_path+"*DY*"}, "stitch_dy");
    //    auto proc_ggf =  Process::MakeShared<Baby_pico>("HToZ#gamma ggF(x1000)", sig, kRed ,{sig_path+"*GluGlu*.root"},   trigs);

    proc_smzg->SetLineWidth (1);
    proc_dy->SetLineWidth   (1);
    //proc_ggf->SetLineWidth   (3);

    vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg};
    //vector<shared_ptr<Process>> procsDY = {proc_dy};

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

    //    NamedFunc B_s("B_s", [](const Baby &b)->NamedFunc::ScalarType{
    //	float BDT = b.BDT_score();
    //	return BDT <= 0.02;
    // });

    NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185. && llphoton_m[0] > 100. && llphoton_m[0] < 180.");
    NamedFunc baseline("nll > 0 && nphoton > 0 && ll_m[0] > 50. && photon_pt[0] > 15. && photon_drmin[0] > 0.4 && photon_pt[0]/llphoton_m[0] > 15./110");
    vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25. && el_pt[ll_i2[0]] > 15.","ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20. && mu_pt[ll_i2[0]] > 10."};
    //vector<NamedFunc> bins = {"photon_pt[0] > 19." && "photon_pt[0] < 35." && B_s <= 0.02, "photon_pt[0] > 19." && "photon_pt[0] < 35." && B_s >0.02, "photon_pt[0] > 35." && B_s <= -0.05, "photon_pt[0] > 35." &&  B_s > -0.05};
    //vector<NamedFunc> lep = {"ll_lepid[0] == 11","ll_lepid[0] == 13"};




    NamedFunc wgt("wgt", [](const Baby &b)->NamedFunc::ScalarType{
        float w = b.w_lumi();
        std::set<std::string> names = b.FileNames();
        std::set<std::string>::iterator it;
        for (it = names.begin(); it != names.end(); it++){
            TString file = *it;
            //if(file.Contains("HToZG")) w *= 1000;
            if(file.Contains("GluGlu")) w *= 1000;
        }
        return w;
    });
    NamedFunc wgt1("wgt1", [](const Baby &b)->NamedFunc::ScalarType{
        float w = b.weight();
        if (w>0) return 1;
	else return 1;
      });
    NamedFunc absw("absw", [](const Baby &b)->NamedFunc::ScalarType{
        float w = b.w_lumi();
        //if (w>0) return 1;
        return abs(w);
      });


    NamedFunc dPhi("dPhi",[](const Baby &b) -> NamedFunc::ScalarType{
        TLorentzVector Z, photon;
        Double_t cosPhi, sinPhi = 0.;
	Double_t angle = 0.;
        Z.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
        photon.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0.);
        if (Z.DeltaPhi(photon)<0) angle = Z.DeltaPhi(photon) + 2*TMath::Pi();
	else angle = Z.DeltaPhi(photon);
	cosPhi = TMath::Cos(angle);
	sinPhi = TMath::Sin(angle);
        return cosPhi + sinPhi;
    });

    NamedFunc dPhi_a("dPhi_a",[](const Baby &b) -> NamedFunc::ScalarType{
        TLorentzVector Z, photon;
        //Double_t cosPhi, sinPhi = 0.;
        Z.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
        photon.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0.);
        //cosPhi = TMath::Cos(Z.DeltaPhi(photon));
        //sinPhi = TMath::Sin(Z.DeltaPhi(photon));
        if (Z.DeltaPhi(photon)<0) return Z.DeltaPhi(photon) + 2*TMath::Pi();
	else return Z.DeltaPhi(photon);
      });

    NamedFunc E2("E2",[](const Baby &b) -> NamedFunc::ScalarType{
        TLorentzVector Z;
	Z.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
       
        return Z.Energy();
      });

    vector<string> leptag = {"el", "mu"};
    PlotMaker pm;
    for(int i(0); i < 1; i++) {
      NamedFunc selection = baseline && lep.at(i) && llphoton_cuts;
      //pm.Push<Hist1D>(Axis(300,0, 6.3, dPhi_a, "d#phi",{}), selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //pm.Push<Hist1D>(Axis(300,-1.42, 1.42, dPhi, "sin+cos",{}), selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //pm.Push<Hist1D>(Axis(400,0, 4, "photon_drmin[0]", "mindR",{}), selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //      pm.Push<Hist1D>(Axis(85,15,100, "photon_pt[0]", "photon_pt ",{}),  selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //      pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_llg ",{}),  selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //      pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_llg bin1",{}),B_s < 0.02 && selection, procs, ops).Weight(wgt).Tag(leptag[i] +"bin1");
      //pm.Push<Hist1D>(Axis(80,-1,1, "BDT_score", "m_llg bin1",{}),  selection, procs, ops).Weight(wgt).Tag(leptag[i] +"bin1");

      //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_llg bin2",{}),  selection && bins[1], procs, ops).Weight(wgt).Tag(leptag[i] + "bin2");
      //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_llg bin3",{}),  selection && bins[2], procs, ops).Weight(wgt).Tag(leptag[i] + "bin3");
      //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_llg bin4",{}),  selection && bins[3], procs, ops).Weight(wgt).Tag(leptag[i] + "bin4");
      pm.Push<Hist1D>(Axis(75,25,100, "el_pt[ll_i1[0]]", "leading el pT (positive weight)",{}), selection && "w_lumi > 0.", procs, ops).Weight(absw).Tag(leptag[i]);
      pm.Push<Hist1D>(Axis(75,25,100, "el_pt[ll_i1[0]]", "leading el pT (negative weight)",{}), selection && "w_lumi < 0.", procs, ops).Weight(absw).Tag(leptag[i]);
	//pm.Push<Hist1D>(Axis(160,40,200, E2, "ll Energy",{}), selection, procs, ops).Weight(wgt).Tag(leptag[i]);
      //pm.Push<Hist1D>(Axis(200,0,200, "ll_pt[0]", "ll pT",{}), selection && "ll_pt[0] < 200", procs, ops).Weight(wgt1).Tag(leptag[i]);
      //pm.Push<Hist2D>(
      //	      Axis(80,100,180, "llphoton_m[0]",      "m_{ee#gamma} [GeV]",{}),
      //	      Axis(65,15,80,     "photon_pt[0]", "p_{T,#gamma}",{}),
      //	      selection, procsDY, ops).Tag(leptag.at(i)).Weight(wgt);
    }
    pm.min_print_ = true;
    pm.MakePlots(1);


}
