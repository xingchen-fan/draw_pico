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
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Overflow(OverflowType::none);/*.Bottom(BottomType::sorb).LabelSize(0.015);*/
  vector<PlotOpt> ops = {lin_lumi};

  //  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[0]>185. && llphoton_m[0] > 100. && llphoton_m[0] < 180. && photon_drmin[0] > 0.4");
  NamedFunc llphoton_cuts("llphoton_m[0] > 100. && llphoton_m[0] < 180. && photon_drmin[0] > 0.4");
  NamedFunc baseline("nll > 0 && nphoton > 0 && ll_m[0] > 50. && photon_pt[0] > 15. && photon_drmin[0] > 0.4");
  NamedFunc CR2_baseline("nll > 0 && ll_m[0] > 50. && photon_pt[0] > 15. && photon_drmin[0] > 0.4");

  NamedFunc CR_ph_eta("abs_eta", [](const Baby &b)->NamedFunc::ScalarType{
      return abs(b.photon_eta()->at(0));
    });
  NamedFunc CR_ph_pt("CR_ph_pt", [](const Baby &b)->NamedFunc::ScalarType{
      return b.photon_pt()->at(0);
    });

  NamedFunc CR_llphoton_m("CR_llphoton_m", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector ll, ph;
      ll.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
      ph.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0);
      return (ll + ph).M();
    });
  NamedFunc CR_llphoton_m_llm("CR_llphoton_m_llm", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector ll, ph;
      ll.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
      ph.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0);
      return (ll + ph).M() + ll.M();
    });
  NamedFunc CR_photon_size("CR_photon_size", [](const Baby &b)->NamedFunc::ScalarType{
      return b.photon_pt()->size();
    });

  //Loose mindR cut//////////////////////////////////////
  NamedFunc loose_mindR_el_cut("loose_mindR_el_cut", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector e1,e2, p;
      double mindr = 9999.0;
      double max_mindr = 0.;
      bool cut = false;
     
      if (b.ll_lepid()->at(0) == 13) cut = false;
      else if(b.ll_lepid()->at(0) == 11){
	e1.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)), b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
	e2.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)), b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
	for (unsigned i(0);i < b.photon_pt()->size();i++){
	  if (((b.photon_eta()->at(i) < 1.4442 && b.photon_eta()->at(i) > -1.4442 && b.photon_idmva()->at(i) > -0.4)||( b.photon_idmva()->at(i) > -0.58 && (b.photon_eta()->at(i) < -1.566 ||b.photon_eta()->at(i) > 1.566) && b.photon_eta()->at(i) <2.5 && b.photon_eta()->at(i) >-2.5)) 
	      && b.photon_pt()->at(i) > 15. && b.photon_elveto()->at(i)){
	    mindr = 9999.0;
	    p.SetPtEtaPhiM(b.photon_pt()->at(i),b.photon_eta()->at(i),b.photon_phi()->at(i),0.);
	    if (abs(e1.DeltaR(p)) < mindr) mindr = abs(e1.DeltaR(p));
	    if (abs(e2.DeltaR(p)) < mindr) mindr = abs(e2.DeltaR(p));
	    if (mindr > max_mindr && mindr < 9000.) max_mindr = mindr;
	  }
	}
	cut =  max_mindr > 0.1;
      }
      return cut;
    });

  NamedFunc loose_mindR_el_ind("loose_mindR_el_ind", [](const Baby &b)->NamedFunc::VectorType{
      TLorentzVector e1,e2, p;
      double mindr = 9999.0;
      double max_mindr = 0.;
      int index = -1;
      vector<double> res;
      e1.SetPtEtaPhiM(b.el_pt()->at(b.ll_i1()->at(0)),b.el_eta()->at(b.ll_i1()->at(0)), b.el_phi()->at(b.ll_i1()->at(0)),0.000511);
      e2.SetPtEtaPhiM(b.el_pt()->at(b.ll_i2()->at(0)),b.el_eta()->at(b.ll_i2()->at(0)), b.el_phi()->at(b.ll_i2()->at(0)),0.000511);
      for (unsigned i(0);i < b.photon_pt()->size();i++){
        if (((b.photon_eta()->at(i) < 1.4442 && b.photon_eta()->at(i) > -1.4442 && b.photon_idmva()->at(i) > -0.4)||( b.photon_idmva()->at(i) > -0.58 && (b.photon_eta()->at(i) < -1.566 ||b.photon_eta()->at(i) > 1.566) && b.photon_eta()->at(i) <2.5 && b.photon_eta()->at(i) >-2.5)) 
	    && b.photon_pt()->at(i) > 15. && b.photon_elveto()->at(i)){
          mindr = 9999.0;
	  p.SetPtEtaPhiM(b.photon_pt()->at(i),b.photon_eta()->at(i),b.photon_phi()->at(i),0.);
          if (abs(e1.DeltaR(p)) < mindr) mindr = abs(e1.DeltaR(p));
          if (abs(e2.DeltaR(p)) < mindr) mindr = abs(e2.DeltaR(p));
	  if (mindr > max_mindr && mindr < 9000.) {max_mindr = mindr; index = i;}
	}
      }
      res.push_back(index);
      res.push_back(max_mindr);
      return res;
    });

  NamedFunc loose_mindR_mu_cut("loose_mindR_mu_cut", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector m1,m2, p;
      double mindr = 9999.0;
      double max_mindr = 0.;
      bool cut = false;
      if (b.ll_lepid()->at(0) == 11) cut = false;
      else if (b.ll_lepid()->at(0) == 13){
	m1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)), b.mu_phi()->at(b.ll_i1()->at(0)), 0.1057);
	m2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)), b.mu_phi()->at(b.ll_i2()->at(0)), 0.1057);
	for (unsigned i(0);i < b.photon_pt()->size();i++){
	  if (((b.photon_eta()->at(i) < 1.4442 && b.photon_eta()->at(i) > -1.4442 && b.photon_idmva()->at(i) > -0.4)||( b.photon_idmva()->at(i) > -0.58 && (b.photon_eta()->at(i) < -1.566 || b.photon_eta()->at(i) > 1.566) && b.photon_eta()->at(i) <2.5 && b.photon_eta()->at(i) >-2.5)) 
	      && b.photon_pt()->at(i) > 15. && b.photon_elveto()->at(i)){
	    mindr = 9999.0;
	    p.SetPtEtaPhiM(b.photon_pt()->at(i),b.photon_eta()->at(i),b.photon_phi()->at(i),0.);
	    if (abs(m1.DeltaR(p)) < mindr) mindr = abs(m1.DeltaR(p));
	    if (abs(m2.DeltaR(p)) < mindr) mindr = abs(m2.DeltaR(p));
	    if (mindr > max_mindr && mindr < 9000.) max_mindr = mindr;
	  }
	}
	cut =  max_mindr > 0.1;
      }
      return cut;
    });

  NamedFunc loose_mindR_mu_ind("loose_mindR_mu_ind", [](const Baby &b)->NamedFunc::VectorType{
      TLorentzVector m1,m2, p;
      double mindr = 9999.0;
      double max_mindr = 0.;
      int index = -1;
      vector<double> res;
      m1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)), b.mu_phi()->at(b.ll_i1()->at(0)), 0.1057);
      m2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)), b.mu_phi()->at(b.ll_i2()->at(0)), 0.1057);
      for (unsigned i(0);i < b.photon_pt()->size();i++){
        if (((b.photon_eta()->at(i) < 1.4442 && b.photon_eta()->at(i) > -1.4442 && b.photon_idmva()->at(i) > -0.4)||( b.photon_idmva()->at(i) > -0.58 && (b.photon_eta()->at(i) < -1.566 ||b.photon_eta()->at(i) > 1.566) && b.photon_eta()->at(i) <2.5 && b.photon_eta()->at(i) >-2.5)) 
	    && b.photon_pt()->at(i) > 15. && b.photon_elveto()->at(i)){
	  mindr = 9999.0;
          p.SetPtEtaPhiM(b.photon_pt()->at(i),b.photon_eta()->at(i),b.photon_phi()->at(i),0.);
	  if (abs(m1.DeltaR(p)) < mindr) mindr = abs(m1.DeltaR(p));
	  if (abs(m2.DeltaR(p)) < mindr) mindr = abs(m2.DeltaR(p));
	  if (mindr > max_mindr && mindr < 9000.) {max_mindr = mindr; index = i;}
        }
      }
      res.push_back(index);
      res.push_back(max_mindr);
      return res;
    });

  NamedFunc loose_mindr_llphoton_m("loose_mindr_llphoton_m", [loose_mindR_mu_ind, loose_mindR_el_ind](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector ll, ph;
      vector<double> vec;
      int index = -1;
      if (b.nel()>1 && b.ll_lepid()->at(0) == 11) {vec = loose_mindR_el_ind.GetVector(b); index = vec[0];}
      if (b.nmu() > 1 && b.ll_lepid()->at(0) == 13) {vec = loose_mindR_mu_ind.GetVector(b); index = vec[0];}
      ll.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
      ph.SetPtEtaPhiM(b.photon_pt()->at(index), b.photon_eta()->at(index), b.photon_phi()->at(index), 0);
      return (ll + ph).M();
    });

  NamedFunc loose_llphoton_m_llm("loose_llphoton_m_llm", [loose_mindR_mu_ind, loose_mindR_el_ind](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector ll, ph;
      vector<double> vec;
      int index = -1;
      if (b.nel()>1 && b.ll_lepid()->at(0) == 11) {vec = loose_mindR_el_ind.GetVector(b); index = vec[0];}
      if (b.nmu() > 1 && b.ll_lepid()->at(0) == 13) {vec = loose_mindR_mu_ind.GetVector(b); index = vec[0];}
      ll.SetPtEtaPhiM(b.ll_pt()->at(0), b.ll_eta()->at(0), b.ll_phi()->at(0), b.ll_m()->at(0));
      ph.SetPtEtaPhiM(b.photon_pt()->at(index), b.photon_eta()->at(index), b.photon_phi()->at(index), 0);
      return (ll + ph).M() + ll.M();
    });  

  NamedFunc loose_mindr("loose_mindr", [loose_mindR_mu_ind, loose_mindR_el_ind](const Baby &b)->NamedFunc::ScalarType{
      double mindr = 999.;
      vector<double> vec;
      if (b.nel()>1 && b.ll_lepid()->at(0) == 11) {vec = loose_mindR_el_ind.GetVector(b); mindr = vec[1];}
      if (b.nmu() > 1 && b.ll_lepid()->at(0) == 13) {vec = loose_mindR_mu_ind.GetVector(b); mindr = vec[1];}
      return mindr;
    });
  NamedFunc loose_ph_pt("loose_ph_pt", [loose_mindR_mu_ind, loose_mindR_el_ind](const Baby &b)->NamedFunc::ScalarType{
      int index = -1;
      vector<double> vec;
      if (b.nel()>1 && b.ll_lepid()->at(0) == 11) {vec = loose_mindR_el_ind.GetVector(b); index = vec[0];}
      if (b.nmu() > 1 && b.ll_lepid()->at(0) == 13) {vec = loose_mindR_mu_ind.GetVector(b); index = vec[0];}
      return b.photon_pt()->at(index);
    });

  /////////////////////////////////////////////////////////////


  NamedFunc real_mindR_mu("real_mindR_mu", [](const Baby &b)->NamedFunc::ScalarType{
      TLorentzVector m1,m2, p;
      double mindr = 9999.0;
      p.SetPtEtaPhiM(b.photon_pt()->at(0),b.photon_eta()->at(0),b.photon_phi()->at(0),0.);
      m1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(0)),b.mu_eta()->at(b.ll_i1()->at(0)), b.mu_phi()->at(b.ll_i1()->at(0)), 0.1057);
      if (abs(m1.DeltaR(p)) < mindr) mindr = abs(m1.DeltaR(p));
      m2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(0)),b.mu_eta()->at(b.ll_i2()->at(0)), b.mu_phi()->at(b.ll_i2()->at(0)), 0.1057);
      if (abs(m2.DeltaR(p)) < mindr) mindr = abs(m2.DeltaR(p));
      return mindr;
    });

  //vector<NamedFunc> lep = {"nel > 1 && ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25. && el_pt[ll_i2[0]] > 15. && (el_eta[ll_i1[0]] < 2.4 && el_eta[ll_i1[0]] > -2.4) && (el_eta[ll_i2[0]] < 2.4 && el_eta[ll_i2[0]] > -2.4)",
  //                        "nmu > 1 && ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20. && mu_pt[ll_i2[0]] > 10. && (mu_eta[ll_i1[0]] < 2.5 && mu_eta[ll_i1[0]] > -2.5) && (mu_eta[ll_i2[0]] < 2.5 && mu_eta[ll_i2[0]] > -2.5)"};


  vector<NamedFunc> lep = {"nel > 1 && ll_lepid[0] == 11 && el_pt[0] > 25. && el_pt[1] > 15. && (el_eta[0] < 2.4 && el_eta[0] > -2.4) && (el_eta[1] < 2.4 && el_eta[1] > -2.4)",
			                           "nmu > 1 && ll_lepid[0] == 13 && mu_pt[0] > 20. && mu_pt[1] > 10. && (mu_eta[0] < 2.5 && mu_eta[0] > -2.5) && (mu_eta[1] < 2.5 && mu_eta[1] > -2.5)"};

  NamedFunc qualityCut = "(photon_idmva[0] > -0.4 && (photon_eta[0] < 1.4442 && photon_eta[0] >-1.4442) ) || (photon_idmva[0] > -0.58 && (photon_eta[0] < -1.566 || photon_eta[0] > 1.566)  && (photon_eta[0] < 2.5&& photon_eta[0] > -2.5 ))";
  NamedFunc nqualityCut = "(photon_idmva[0] < -0.4 && (photon_eta[0] < 1.4442 && photon_eta[0] >-1.4442) ) || (photon_idmva[0] < -0.58 && (photon_eta[0] < -1.566 || photon_eta[0] > 1.566) && (photon_eta[0] < 2.5&& photon_eta[0] > -2.5 ) )";
  
  NamedFunc photon_all_cut =  CR_photon_size > 0.  && "photon_pt[0] > 15." &&  qualityCut && "photon_elveto[0]" &&"photon_drmin[0] > 0.4"; 
  /*  vector<NamedFunc> SR = {lep[0]||lep[1],
			  (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50.",
			  (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50." && photon_all_cut,
			  (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50." && photon_all_cut && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.",
			  (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50." && photon_all_cut && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100." && "photon_pt[0]/llphoton_m[0] > 15./110."};*/
      vector<NamedFunc> SR = {lep[0]||lep[1],
                           (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50.",
                           (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&  "nll > 0" && "ll_m[0] > 50.",
                           (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]",
                           (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]",

                           (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" && "llphoton_m[0] > 100." && "llphoton_m[0] < 180.",
			      (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100."};

      //NamedFunc all_SR_nomindR =  (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&  qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&&  CR_ph_pt/CR_llphoton_m > 15./110. && CR_llphoton_m < 180. &&  CR_llphoton_m > 100.;

  vector<NamedFunc> CR1 = {lep[0]||lep[1], 
			   (lep[0]||lep[1]) &&  "nll > 0" && "ll_m[0] > 50.",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&  "nll > 0" && "ll_m[0] > 50.",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" && "llphoton_m[0] < 100.",  
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 100."};

  vector<NamedFunc> CR2 = {"nll > 0" && "ll_m[0] > 50." && (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15.",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15."  && nqualityCut && "nll > 0"&& "photon_idmva[0] < -0.8" && "ll_m[0] > 50." && "photon_elveto[0]",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15."  && nqualityCut && "nll > 0"&& "photon_idmva[0] < -0.8" && "ll_m[0] > 50." && "photon_elveto[0]"&&"photon_drmin[0] > 0.4",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && nqualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" &&"photon_idmva[0] < -0.8" && CR_llphoton_m > 100. && CR_llphoton_m < 180.,
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && nqualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" &&"photon_idmva[0] < -0.8" && CR_ph_pt/CR_llphoton_m > 15./110. && CR_llphoton_m > 100. && CR_llphoton_m < 180.};
  
  vector<NamedFunc> CR3 = {"nll > 0" && "ll_m[0] > 50." && (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15.",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15."  && nqualityCut && "nll > 0"&& "photon_idmva[0] > -0.8" && "ll_m[0] > 50." && "photon_elveto[0]",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15."  && nqualityCut && "nll > 0"&& "photon_idmva[0] > -0.8" && "ll_m[0] > 50." && "photon_elveto[0]"&&"photon_drmin[0] > 0.4",
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && nqualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" &&"photon_idmva[0] > -0.8" && CR_llphoton_m > 100. && CR_llphoton_m < 180.,
			   (lep[0]||lep[1]) && CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && nqualityCut && "nll > 0" && "ll_m[0] > 50." && "photon_elveto[0]" &&"photon_idmva[0] > -0.8" && CR_ph_pt/CR_llphoton_m > 15./110. && CR_llphoton_m > 100. && CR_llphoton_m < 180.};

  NamedFunc stitch("stitch", [](const Baby &b)->NamedFunc::ScalarType{
      bool st = b.stitch_dy();
      bool res = true;
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("DYJetsToLL")) res = st;
        else if(file.Contains("ZGToLLG")) res =  !st;
      }
      return res;
    });


  NamedFunc wgt("wgt", [](const Baby &b)->NamedFunc::ScalarType{
      float w = b.weight();
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
	if(file.Contains("HToZG")) w *= 1000;
	if(file.Contains("GGF")) w *= 1000;
      }
      return w;
    });
  NamedFunc w_l("w_l", [](const Baby &b)->NamedFunc::ScalarType{
      float w = b.w_lumi();
      float w_year = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
      } else if (b.SampleType()==2017) {
        w_year = 41.52756;
      } 
      else if (b.SampleType()==2018){
	w_year = 59.83;
      }
      std::set<std::string> names = b.FileNames();
      std::set<std::string>::iterator it;
      for (it = names.begin(); it != names.end(); it++){
        TString file = *it;
        if(file.Contains("HToZG")) w *= 1000;
        if(file.Contains("GGF")) w *= 1000;
      }
      return w * w_year;
    });

  //NamedFunc B_s("B_s", [](const Baby &b)->NamedFunc::ScalarType{
  //   float BDT = b.BDT_score();
  //  return BDT;
  // });

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


  NamedFunc sync_w1("sync_w1", [](const Baby &b)->NamedFunc::ScalarType{
      float w_year = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
      } else if (b.SampleType()==2017) {
        w_year = 41.52756;
      }

      return w_year * b.w_lumi() * b.w_pu();
    });
  NamedFunc sync_w2("sync_w2", [](const Baby &b)->NamedFunc::ScalarType{
      float w_year = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
      } else if (b.SampleType()==2017) {
        w_year = 41.52756;
      }

      return w_year * b.w_lumi() * b.w_photon_csev() * b.w_photon_id() * b.w_pu();
    });

  NamedFunc sync_w3("sync_w3", [](const Baby &b)->NamedFunc::ScalarType{
      float w_year = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
      } else if (b.SampleType()==2017) {
        w_year = 41.52756;
      }

      return w_year * b.w_lumi() * b.w_photon_csev() * b.w_photon_id() * b.w_pu() * b.w_prefire();
    });
  vector<NamedFunc> CRs = {
    baseline && (lep[0] || lep[1]) && qualityCut && "photon_elveto[0]" && "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 100.",
    CR_photon_size > 0. && CR2_baseline && (lep[0] || lep[1]) && nqualityCut  && "photon_idmva[0] < -0.8" && "photon_elveto[0]" && CR_ph_pt/CR_llphoton_m > 15./110. && CR_llphoton_m > 100. && CR_llphoton_m < 180.,
    CR_photon_size > 0. && CR2_baseline && (lep[0] || lep[1]) && nqualityCut && "photon_idmva[0] > -0.8" && "photon_elveto[0]" && CR_ph_pt/CR_llphoton_m > 15./110. && CR_llphoton_m > 100. && CR_llphoton_m < 180.,
     baseline && (lep[0] || lep[1]) && qualityCut && "photon_elveto[0]" && "photon_pt[0]/llphoton_m[0] < 15./110." && "llphoton_m[0] > 100." && "llphoton_m[0] < 180." };

  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  //  string folder("/net/cms37/data1/xf82/BDT/OutputpT/");

  string folder16("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/unskimmed/");
  string folder17("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/unskimmed/");
  string folder18("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/unskimmed/");

  //  string folder16("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2016/mc/merged_zgmc_llg/");
  //  string folder17("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2017/mc/merged_zgmc_llg/");
  //  string folder18("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/2018/mc/merged_zgmc_llg/");
  //auto proc_smzg_fsr  = Process::MakeShared<Baby_pico>("SM Z#gamma (FSR)",       back,  TColor::GetColor("#16bac5"),{folder+"*ZGToLLG*"}, (trig1 || (trig2||trig3)) && "llphoton_m[0] + ll_m[0] < 185.");
  //  auto proc_smzg_isr  = Process::MakeShared<Baby_pico>("SM Z#gamma (ISR)",       back,  TColor::GetColor("#16bac5"),{folder+"*ZGToLLG*"}, (trig1 || (trig2||trig3)) && "llphoton_m[0] + ll_m[0] > 185.");
  NamedFunc el_trigs16("el_trigs16", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();});
  NamedFunc el_trigs17("el_trigs17", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();});
  NamedFunc el_trigs18("el_trigs18", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();});
  NamedFunc mu_trigs17("mu_trigs17", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8() || b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();});
  NamedFunc mu_trigs16("mu_trigs16", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() || b.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();});
  NamedFunc mu_trigs18("mu_trigs18", [](const Baby &b)->NamedFunc::ScalarType{return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8() || b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();});

  NamedFunc trigs16(el_trigs16 || mu_trigs16);
  NamedFunc trigs17(el_trigs17 || mu_trigs17);
  NamedFunc tot_trigs_el("tot_trigs_el", [el_trigs16, el_trigs17, el_trigs18](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (b.SampleType() == 2016) res = el_trigs16.GetScalar(b);
      else if (b.SampleType() == 2017) res = el_trigs17.GetScalar(b);
      else if (b.SampleType() == 2018) res = el_trigs18.GetScalar(b);
      return res;
    });
  NamedFunc tot_trigs_mu("tot_trigs_mu", [mu_trigs16, mu_trigs17, mu_trigs18](const Baby &b)->NamedFunc::ScalarType{
      bool res = false;
      if (b.SampleType() == 2016) res = mu_trigs16.GetScalar(b);
      else if (b.SampleType() == 2017) res = mu_trigs17.GetScalar(b);
      else if (b.SampleType() == 2018) res = mu_trigs18.GetScalar(b);
      return res;
    });


  ///////////////////////////loose mindR cut is here
  vector<NamedFunc> loose_all_cut = {tot_trigs_el && lep[0] && CR_photon_size > 0. && "nll > 0" && "ll_m[0] > 50." && loose_mindR_el_cut && loose_mindr_llphoton_m > 100. && loose_mindr_llphoton_m < 180. && (loose_ph_pt/loose_mindr_llphoton_m) > 15./110., tot_trigs_mu && lep[1] && CR_photon_size > 0. && "nll > 0" && "ll_m[0] > 50." && loose_mindR_mu_cut && loose_mindr_llphoton_m > 100. && loose_mindr_llphoton_m < 180. && (loose_ph_pt/loose_mindr_llphoton_m) > 15./110.};
  vector<NamedFunc> loose_least_cut = {tot_trigs_el && lep[0] && CR_photon_size > 0. && "nll > 0"  && loose_mindR_el_cut, tot_trigs_mu && lep[1] && CR_photon_size > 0. && "nll > 0" && loose_mindR_mu_cut};

  ////////////////////tighter SR cut
  NamedFunc tighter_SR = (tot_trigs_el && lep[0] &&  CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut  && "photon_idmva[0] > 0.4" && "nll > 0" && "ll_m[0] > 80. && ll_m[0] < 100." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.")\
    || (tot_trigs_mu && lep[1] &&  CR_photon_size > 0. && "photon_pt[0] > 15." &&"photon_drmin[0] > 0.4" && qualityCut  && "photon_idmva[0] > 0.4" && "nll > 0" && "ll_m[0] > 80. && ll_m[0] < 100." && "photon_elveto[0]"&& "photon_pt[0]/llphoton_m[0] > 15./110." && "llphoton_m[0] < 180." &&  "llphoton_m[0] > 100.");
  


  //////////////////////////////////////////////////////

  vector<NamedFunc> loose_cut_flow = {(tot_trigs_el && lep[0] && CR_photon_size > 0.)||(tot_trigs_mu && lep[1] && CR_photon_size > 0.)};

  NamedFunc sync_w("sync_w", [el_trigs16, mu_trigs16, el_trigs17, mu_trigs17](const Baby &b)->NamedFunc::ScalarType{
      float w_year = 0;
      float w_lep = 0;
      if (b.SampleType()==2016) {
        w_year = 36.32264;
	if (el_trigs16.GetScalar(b)) w_lep = b.w_el_id();
	else if (mu_trigs16.GetScalar(b)) w_lep = b.w_mu_id() * b.w_mu_iso();
      } 
      else if (b.SampleType()==2017) {
        w_year = 41.52756;
	if (el_trigs17.GetScalar(b)) w_lep = b.w_el_id();
        else if (mu_trigs17.GetScalar(b)) w_lep = b.w_mu_id() * b.w_mu_iso();
      }

      return w_year * b.w_lumi() * b.w_photon_csev() * b.w_photon_id() * b.w_pu() * b.w_prefire() * w_lep;
    });

  //auto proc_hzg_zh   = Process::MakeShared<Baby_pico>("HToZ#gamma(x100)", sig, kRed,{folder+"*ZH*"},   "1");
  //auto proc_hzg_ggf   = Process::MakeShared<Baby_pico>("HToZ#gamma(ggf)(x1000)", sig, kRed,{folder+"GGF_output.root"},   "1");
  auto proc_smzg_isr16el  = Process::MakeShared<Baby_pico>("SM Z#gamma ISR el",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, "!stitch_dy" && el_trigs16 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm > 185.);
  auto proc_smzg_fsr16el  = Process::MakeShared<Baby_pico>("SM Z#gamma FSR el",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, "!stitch_dy" && el_trigs16 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm < 185.);
  auto proc_dy16el    = Process::MakeShared<Baby_pico>("DY + Fake el",               back, TColor::GetColor("#ffb400"),{folder16+"*DYJetsToLL_M-50*"}, "stitch_dy" &&  el_trigs16);
  auto proc_smzg_isr16mu  = Process::MakeShared<Baby_pico>("SM Z#gamma ISR mu",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, "!stitch_dy" && mu_trigs16 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm > 185.);
  auto proc_smzg_fsr16mu  = Process::MakeShared<Baby_pico>("SM Z#gamma FSR mu",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, "!stitch_dy" && mu_trigs16 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm < 185.);
  auto proc_dy16mu    = Process::MakeShared<Baby_pico>("DY + Fake mu",               back, TColor::GetColor("#ffb40"),{folder16+"*DYJetsToLL_M-50*"}, "stitch_dy" &&  mu_trigs16);
  auto proc_smzg16elcutflow  = Process::MakeShared<Baby_pico>("SM Z#gamma el",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, el_trigs16);
  auto proc_smzg16mucutflow  = Process::MakeShared<Baby_pico>("SM Z#gamma mu",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, mu_trigs16);
  auto proc_dy16mucutflow    = Process::MakeShared<Baby_pico>("DY + Fake mu",               back, TColor::GetColor("#ffb400"),{folder16+"*DYJetsToLL_M-50*"},mu_trigs16);
  auto proc_dy16elcutflow    = Process::MakeShared<Baby_pico>("DY + Fake el",               back, TColor::GetColor("#ffb400"),{folder16+"*DYJetsToLL_M-50*"}, el_trigs16);

  auto proc_smzg_isr17el  = Process::MakeShared<Baby_pico>("SM Z#gamma ISR el",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, "!stitch_dy" && el_trigs17 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm > 185.);
  auto proc_smzg_fsr17el  = Process::MakeShared<Baby_pico>("SM Z#gamma FSR el",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, "!stitch_dy" && el_trigs17 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm < 185.);
  auto proc_dy17el    = Process::MakeShared<Baby_pico>("DY + Fake el",               back, TColor::GetColor("#ffb400"),{folder17+"*DYJetsToLL_M-50*"}, "stitch_dy" &&  el_trigs17);
  auto proc_smzg_isr17mu  = Process::MakeShared<Baby_pico>("SM Z#gamma ISR mu",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, "!stitch_dy" && mu_trigs17 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm > 185.);
  auto proc_smzg_fsr17mu  = Process::MakeShared<Baby_pico>("SM Z#gamma FSR mu",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, "!stitch_dy" && mu_trigs17 &&  CR_photon_size > 0. && "nll > 0." && CR_llphoton_m_llm < 185.);
  auto proc_dy17mu    = Process::MakeShared<Baby_pico>("DY + Fake mu",               back, TColor::GetColor("#ffb400"),{folder17+"*DYJetsToLL_M-50*"}, "stitch_dy" &&  mu_trigs17);
  auto proc_smzg17elcutflow  = Process::MakeShared<Baby_pico>("SM Z#gamma el",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, el_trigs17);
  auto proc_smzg17mucutflow  = Process::MakeShared<Baby_pico>("SM Z#gamm mu",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, mu_trigs17);
  auto proc_dy17mucutflow    = Process::MakeShared<Baby_pico>("DY + Fake mu",               back, TColor::GetColor("#ffb400"),{folder17+"*DYJetsToLL_M-50*"},mu_trigs17);
  auto proc_dy17elcutflow    = Process::MakeShared<Baby_pico>("DY + Fake el",               back, TColor::GetColor("#ffb400"),{folder17+"*DYJetsToLL_M-50*"}, el_trigs17);

  auto proc_smzg17tot  = Process::MakeShared<Baby_pico>("SM Z#gamma 17",       back,  TColor::GetColor("#16bac5"),{folder17+"*ZGToLLG*"}, "1");
  auto proc_smzg16tot  = Process::MakeShared<Baby_pico>("SM Z#gamma 16",       back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*"}, "1");
  auto proc_dy17tot    = Process::MakeShared<Baby_pico>("DY + Fake 17",               back, TColor::GetColor("#ffb400"),{folder17+"*DYJetsToLL_M-50*"},"1");
  auto proc_dy16tot    = Process::MakeShared<Baby_pico>("DY + Fake 16",               back, TColor::GetColor("#ffb400"),{folder16+"*DYJetsToLL_M-50*"},"1");

  auto proc_ggf17el  = Process::MakeShared<Baby_pico>("SM Z#gamma el",       back,  TColor::GetColor("#16bac5"),{folder17+"*GluGluHToZG_ZToLL_M-125*"}, el_trigs17);
  auto proc_ggf17mu  = Process::MakeShared<Baby_pico>("SM Z#gamma mu",       back,  TColor::GetColor("#16bac5"),{folder17+"*GluGluHToZG_ZToLL_M-125*"}, mu_trigs17);

  auto proc_lep_smzg = Process::MakeShared<Baby_pico>("SM Z#gamma",  back,  TColor::GetColor("#16bac5"),{folder16+"*ZGToLLG*", folder17+"*ZGToLLG*", folder18+"*ZGToLLG*"}, stitch);
  auto proc_lep_dy = Process::MakeShared<Baby_pico>("DY + Fake",   back, TColor::GetColor("#ffb400"),{folder16+"*DYJetsToLL_M-50*", folder17+"*DYJetsToLL_M-50*", folder18+"*DYJetsToLL_M-50*"}, stitch);
  auto proc_lep_ggf = Process::MakeShared<Baby_pico>("GGF (X1000)",       sig,  kRed,{folder16+"*GluGluHToZG_ZToLL_M-125*", folder17+"*GluGluHToZG_ZToLL_M-125*",folder18+"*GluGluHToZG_ZToLL_M-125*"}, "1");
  auto proc_lep_vbf = Process::MakeShared<Baby_pico>("VBF (X1000)",       sig,  kBlue,{folder16+"*VBFHToZG_ZToLL_M-125*", folder17+"*VBFHToZG_ZToLL_M-125*",folder18+"*VBFHToZG_ZToLL_M-125*"}, "1");

  //proc_smzg->SetLineWidth (1);
  //proc_dy->SetLineWidth   (1);
  //proc_hzg_ggf->SetLineWidth   (3);
  proc_lep_smzg->SetLineWidth (1);
  proc_lep_dy->SetLineWidth (1);
  vector<shared_ptr<Process>> procs16 = {proc_dy16el, proc_smzg_isr16el, proc_smzg_fsr16el, proc_dy16mu, proc_smzg_isr16mu, proc_smzg_fsr16mu};
  vector<shared_ptr<Process>> procs17 = {proc_dy17el, proc_smzg_isr17el, proc_smzg_fsr17el, proc_dy17mu, proc_smzg_isr17mu, proc_smzg_fsr17mu};
  vector<shared_ptr<Process>> procs16nosr_el = {proc_dy16elcutflow, proc_smzg16elcutflow};
  vector<shared_ptr<Process>> procs16nosr_mu = {proc_dy16mucutflow, proc_smzg16mucutflow};
  vector<shared_ptr<Process>> procs17nosr_el = {proc_dy17elcutflow, proc_smzg17elcutflow};
  vector<shared_ptr<Process>> procs17nosr_mu = {proc_dy17mucutflow, proc_smzg17mucutflow};
  vector<shared_ptr<Process>> procstot = {proc_dy16tot, proc_smzg16tot, proc_dy17tot, proc_smzg17tot};
  vector<shared_ptr<Process>> procs17sigggf_el = {proc_ggf17el};
  vector<shared_ptr<Process>> procs17sigggf_mu = {proc_ggf17mu};
  vector<shared_ptr<Process>> procs = {proc_lep_smzg, proc_lep_dy, proc_lep_ggf, proc_lep_vbf};

  vector<string> leptag = {"el_", "mu_"};
  PlotMaker pm;
  //for(int i(0); i < 1; i++) {
  //NamedFunc selection = (lep.at(0)||lep.at(1)) && baseline && llphoton_cuts;
  //  pm.Push<Hist1D>(Axis(50,-5,5, w_l, "w_lumi * lumi",{}), stitch && SR[1] && CR_photon_size > 0.&& "photon_pt[0] < 15.", procs16nosr_el, ops).Weight(wgt1).Tag("weight_prob_");
    //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= 0.004 && B_s < 0.093, procs, ops).Weight(w_l).Tag("Untagged_2");
    //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= -0.09 && B_s < 0.004, procs, ops).Weight(w_l).Tag("Untagged_3");
    //pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= -1. && B_s < 0.09, procs, ops).Weight(w_l).Tag("Untagged_4");
    // pm.Push<Hist1D>(Axis(90,-0.5,0.4, "BDT_score", "MVA",{}), selection, procs, ops).Weight(w_l);
    //	pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= -1. && B_s < -0.05, procs, ops).Weight(w_l).Tag("Unknown1");
    //	pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= -0.05 && B_s < 0.05, procs, ops).Weight(w_l).Tag("Unknown2");
    //	pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection && B_s >= 0.05 && B_s < 1, procs, ops).Weight(w_l).Tag("Unknown3");
	//    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection, procs, ops).Weight(wgt).Tag(leptag[i]);
    //    pm.Push<Hist1D>(Axis(80,100,180, "llphoton_m[0]", "m_{ll#gamma}",{}), selection, procs, ops).Weight(w_l).Tag(leptag[i] + "lumi");
    //pm.Push<Hist1D>(Axis(65,15,80, "photon_pt[0]", "#gamma pT",{}), selection, procs, ops).Weight(w_l).Tag(leptag[i] + "lumi");
    //pm.Push<Hist1D>(Axis(50,50,100, "llphoton_m[0]", "m_{ll#gamma}",{}), selection, procs, ops).Weight(wgt).Tag("CR_");
    //pm.Push<Hist1D>(Axis(80,100,180, CR_llphoton_m, "m_{ll#gamma}",{}),CR_photon_size > 0. && (lep.at(i)||lep.at(i + 1)) && CR2_baseline && CR_llphoton_m <180. && CR_llphoton_m > 100. && (( CR_ph_eta < 1.4442 && "photon_idmva[0] < -0.4 && photon_idmva[0] > -0.8") || (CR_ph_eta > 1.566 &&  CR_ph_eta <"photon_idmva[0] < -0.58 && photon_idmva[0] > -0.8")), procs, ops).Weight(w_l).Tag("CR2_");
    
  //pm.Push<Hist1D>(Axis(50, 0.1,2, loose_mindr, "mindR electron",{}), loose_least_cut[0], procs, ops).Weight(w_l).Tag("ShortName:mindR_check_least_el_");
  //pm.Push<Hist1D>(Axis(50, 0.1,2, loose_mindr, "mindR muon",{}), loose_least_cut[1], procs, ops).Weight(w_l).Tag("ShortName:mindR_check_least_mu_");
  //pm.Push<Hist1D>(Axis(50, 0.1,2, loose_mindr, "mindR electron",{}), loose_all_cut[0], procs, ops).Weight(w_l).Tag("ShortName:mindR_check_el_");
  //  pm.Push<Hist1D>(Axis(50, 0.1,2, loose_mindr, "mindR muon",{}), loose_all_cut[1], procs, ops).Weight(w_l).Tag("ShortName:mindR_check_mu_");
  pm.Push<Hist1D>(Axis(85, 15,100, "photon_pt[0]", "Signal photon p_{T} Untagged",{}), tighter_SR && "njet <= 1", procs, ops).Weight(w_l).Tag("ShortName:tighter_untagged_unskimmed");
  pm.Push<Hist1D>(Axis(85, 15,100, "photon_pt[0]", "Signal photon p_{T} Dijet",{}), tighter_SR && "njet >= 2", procs, ops).Weight(w_l).Tag("ShortName:tighter_dijet_unskimmed");


  //pm.Push<Hist1D>(Axis(80,100,180, CR_llphoton_m, "m_{ll#gamma}",{}),CR_photon_size > 0. && (lep.at(i)||lep.at(i + 1)) && CR2_baseline && CR_llphoton_m <180. && CR_llphoton_m > 100. && "photon_idmva[0] < -0.8", procs, ops).Weight(w_l).Tag("CR3_");

  /*  pm.Push<Table>("Sync_16_sync_w", vector<TableRow>{
      TableRow("negative weight", stitch && SR[1] && CR_photon_size > 0.&& "photon_pt[0] < 15." ,0,0,w_l),
	//Tablerow("positive weight", stitch && SR[1] && CR_photon_size > 0.&& "photon_pt[0] < 15." && w_l > 0.,0,0,wgt1)
	}, procs16nosr_el, false);
  */
  /*  pm.Push<Table>("Sync_16_sync_w", vector<TableRow>{
      TableRow("Z to ll: m_llg < 100 w_lumi", CRs[0],0,0,w_l),
      TableRow("Z to ll: m_llg < 100 w_pu", CRs[0],0,0,sync_w1),
	TableRow("Z to ll: m_llg < 100 w_ph_csev w_ph_id", CRs[0],0,0,sync_w2),
	TableRow("Z to ll: m_llg < 100 w_prefire", CRs[0],0,0,sync_w3),
	TableRow("Z to ll: m_llg < 100 w_lep", CRs[0],0,0,sync_w),

	TableRow("Z to ll: idmva < -0.8 w_lumi", CRs[1],0,0,w_l),
	TableRow("Z to ll: idmva < -0.8 w_pu", CRs[1],0,0,sync_w1),
        TableRow("Z to ll: idmva < -0.8 w_ph_csev w_ph_id", CRs[1],0,0,sync_w2),
        TableRow("Z to ll: idmva < -0.8 w_prefire", CRs[1],0,0,sync_w3),
        TableRow("Z to ll: idmva < -0.8 w_lep", CRs[1],0,0,sync_w),

	TableRow("Z to ll: idmva > -0.8 w_lumi", CRs[2],0,0,w_l),
        TableRow("Z to ll: idmva > -0.8 w_pu", CRs[2],0,0,sync_w1),
        TableRow("Z to ll: idmva > -0.8 w_ph_csev w_ph_id", CRs[2],0,0,sync_w2),
        TableRow("Z to ll: idmva > -0.8 w_prefire", CRs[2],0,0,sync_w3),
        TableRow("Z to ll: idmva > -0.8 w_lep", CRs[2],0,0,sync_w),
	}, procs16, false);
  pm.Push<Table>("Sync_17_sync_w", vector<TableRow>{
      TableRow("Z to ll: m_llg < 100 w_lumi", CRs[0],0,0,w_l),
	TableRow("Z to ll: m_llg < 100 w_pu", CRs[0],0,0,sync_w1),
        TableRow("Z to ll: m_llg < 100 w_ph_csev w_ph_id", CRs[0],0,0,sync_w2),
        TableRow("Z to ll: m_llg < 100 w_prefire", CRs[0],0,0,sync_w3),
 stitch &&        TableRow("Z to ll: m_llg < 100 w_lep", CRs[0],0,0,sync_w),

        TableRow("Z to ll: idmva < -0.8 w_lumi", CRs[1],0,0,w_l),
        TableRow("Z to ll: idmva < -0.8 w_pu", CRs[1],0,0,sync_w1),
        TableRow("Z to ll: idmva < -0.8 w_ph_csev w_ph_id", CRs[1],0,0,sync_w2),
        TableRow("Z to ll: idmva < -0.8 w_prefire", CRs[1],0,0,sync_w3),
        TableRow("Z to ll: idmva < -0.8 w_lep", CRs[1],0,0,sync_w),

        TableRow("Z to ll: idmva > -0.8 w_lumi", CRs[2],0,0,w_l),
        TableRow("Z to ll: idmva > -0.8 w_pu", CRs[2],0,0,sync_w1),
        TableRow("Z to ll: idmva > -0.8 w_ph_csev w_ph_id", CRs[2],0,0,sync_w2),
        TableRow("Z to ll: idmva > -0.8 w_prefire", CRs[2],0,0,sync_w3),
        TableRow("Z to ll: idmva > -0.8 w_lep", CRs[2],0,0,sync_w),
        }, procs17, false);
  */
  /*  pm.Push<Table>("Sync_16_cutflow_sig_el", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && SR[0] && lep[0],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && SR[1] && lep[0],0,0,w_l),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[0],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[0],0,0,w_l),
	//TableRow("$mindR > 0.4$", stitch && SR[4] && lep[0],0,0,w_l),
	TableRow("Photon cut", stitch && SR[4] && lep[0],0,0,w_l),
	TableRow("$m_{llphoton}$ cut", stitch && SR[5] && lep[0],0,0,w_l),
	TableRow("$photon p_T / m_{llphoton}$ cut", stitch && SR[6] && lep[0],0,0,w_l)
	//TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs16nosr_el, false).Precision(0);
  pm.Push<Table>("Sync_16_cutflow_sig_mu", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && SR[0] && lep[1],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && SR[1] && lep[1],0,0,w_l),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[1],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[1],0,0,w_l),
        TableRow("Photon cut", stitch && SR[4] && lep[1],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && SR[5] && lep[1],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && SR[6] && lep[1],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs16nosr_mu, false).Precision(0);
  //  pm.Push<Table>("Sync_tot", vector<TableRow>{TableRow("Tot", "1",0,0,w_l)}, procstot,false).Precision(0);

  pm.Push<Table>("Sync_17_cutflow_sig_el", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && SR[0] && lep[0],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && SR[1] && lep[0],0,0,w_l),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[0],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[0],0,0,w_l),
        TableRow("Photon cut", stitch && SR[4] && lep[0],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && SR[5] && lep[0],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && SR[6] && lep[0],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17nosr_el, false).Precision(0);
  pm.Push<Table>("Sync_17_cutflow_sig_mu", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && SR[0] && lep[1],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && SR[1] && lep[1],0,0,w_l),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[1],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[1],0,0,w_l),
        //TableRow("$mindR > 0.4$", stitch && SR[4] && lep[1],0,0,w_l),
        TableRow("Photon cut", stitch && SR[4] && lep[1],0,0,w_l),
	TableRow("$m_{llphoton}$ cut", stitch && SR[5] && lep[1],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && SR[6] && lep[1],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17nosr_mu, false).Precision(0);
  */
  /*pm.Push<Table>("Sync_17_cutflow_sig_ggf_el", vector<TableRow>{
      TableRow("Triggers", "1",0,0,"1"),
        //TableRow("Overlap removal", stitch,0,0,"1"),
        TableRow("Lepton cut", SR[0] && lep[0],0,0,"1"),
        TableRow("$m_{ll} > 50", SR[1] && lep[0],0,0,"1"),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[0],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[0],0,0,w_l),
        TableRow("Photon cut", SR[4] && lep[0],0,0,"1"),
        TableRow("$m_{llphoton}$ cut", SR[5] && lep[0],0,0,"1"),
        TableRow("$photon p_T / m_{llphoton}$ cut", SR[6] && lep[0],0,0,"1")
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17sigggf_el, false).Precision(0);
  pm.Push<Table>("Sync_17_cutflow_sig_ggf_mu", vector<TableRow>{
      TableRow("Triggers", "1",0,0,"1"),
        //TableRow("Overlap removal", stitch,0,0,"1"),
        TableRow("Lepton cut", SR[0] && lep[1],0,0,"1"),
        TableRow("$m_{ll} > 50", SR[1] && lep[1],0,0,"1"),
        //TableRow("$photon p_T > 15$",  stitch && SR[2] && lep[0],0,0,w_l),
        //TableRow("Photon eta and idmva cut (and el veto)", stitch && SR[3] && lep[0],0,0,w_l),
        TableRow("Photon cut", SR[4] && lep[1],0,0,"1"),
        TableRow("$m_{llphoton}$ cut", SR[5] && lep[1],0,0,"1"),
        TableRow("$photon p_T / m_{llphoton}$ cut", SR[6] && lep[1],0,0,"1")
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17sigggf_mu, false).Precision(0);*/

  
  /*  pm.Push<Table>("Sync_16_cutflow_CR1", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && CR1[0],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && CR1[1],0,0,w_l),
        TableRow("$photon p_T > 15$",  stitch && CR1[2],0,0,w_l),
        TableRow("Photon eta and idmva cut (and el veto)", stitch && CR1[3],0,0,w_l),
        TableRow("$mindR > 0.4$", stitch && CR1[4],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && CR1[5],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && CR1[6],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs16nosr, false).Precision(0);


  pm.Push<Table>("Sync_17_cutflow_CR1", vector<TableRow>{
      TableRow("Triggers", "1",0,0,w_l),
        TableRow("Overlap removal", stitch,0,0,w_l),
        TableRow("Lepton cut", stitch && CR1[0],0,0,w_l),
        TableRow("$m_{ll} > 50", stitch && CR1[1],0,0,w_l),
        TableRow("$photon p_T > 15$",  stitch && CR1[2],0,0,w_l),
        TableRow("Photon eta and idmva cut (and el veto)", stitch && CR1[3],0,0,w_l),
        TableRow("$mindR > 0.4$", stitch && CR1[4],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && CR1[5],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && CR1[6],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR1[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR1[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17nosr, false).Precision(0);
  pm.Push<Table>("Sync_16_cutflow_CR2", vector<TableRow>{
      TableRow("$photon p_T > 15$",  stitch && CR2[0],0,0,w_l),
        TableRow("Photon eta and idmva cut (and el veto)", stitch && CR2[1],0,0,w_l),
        TableRow("$mindR > 0.4$", stitch && CR2[2],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && CR2[3],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && CR2[4],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR2[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR2[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs16nosr, false).Precision(0);

  pm.Push<Table>("Sync_17_cutflow_CR2", vector<TableRow>{
      TableRow("$photon p_T > 15$",  stitch && CR2[0],0,0,w_l),
        TableRow("Photon eta and idmva cut (and el veto)", stitch && CR2[1],0,0,w_l),
        TableRow("$mindR > 0.4$", stitch && CR2[2],0,0,w_l),
        TableRow("$m_{llphoton}$ cut", stitch && CR2[3],0,0,w_l),
        TableRow("$photon p_T / m_{llphoton}$ cut", stitch && CR2[4],0,0,w_l)
        //TableRow("$m_{llg} < 100$ ISR", stitch && CR2[3] && CR_llphoton_m_llm > 185.,0,0,w_l),
        //TableRow("$m_{llg} < 100$ FSR", stitch && CR2[3] && CR_llphoton_m_llm < 185.,0,0,w_l),

        }, procs17nosr, false).Precision(0);
	*/
  pm.min_print_ = true;
  pm.MakePlots(1);
}


///net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/
