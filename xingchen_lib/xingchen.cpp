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
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
#include "core/xingchen.hpp"

NamedFunc::VectorType fmu_sig(const Baby &b){
  std::vector<double> mu_sig_;
  for (unsigned imu = 0; imu < b.Muon_pt()->size(); imu++) {
    bool is_sig = true;
    if (b.Muon_pt()->at(imu) <= 5.0) is_sig = false;
    if (abs(b.Muon_eta()->at(imu)) > 2.4) is_sig = false;
    if (abs(b.Muon_dz()->at(imu))>1.0)  is_sig = false;
    if (abs(b.Muon_dxy()->at(imu))>0.5) is_sig = false;
    if (!((b.Muon_looseId()->at(imu) || (b.Muon_pt()->at(imu) > 200 && b.Muon_highPtId()->at(imu))) &&
	  b.Muon_pfRelIso03_all()->at(imu) < 0.35 &&
	  b.Muon_sip3d()->at(imu) < 4)) is_sig = false;
    if (is_sig && b.Muon_pt()->at(imu)>20.0)
      mu_sig_.push_back(3.0);
    else if (is_sig && b.Muon_pt()->at(imu)>10.0)
      mu_sig_.push_back(2.0);
    else if (is_sig)
      mu_sig_.push_back(1.0);
    else
      mu_sig_.push_back(0.0);
  }
  return mu_sig_;
}

NamedFunc::VectorType fel_sig(const Baby &b){
  std::vector<double> el_sig_;
  for (unsigned iel = 0; iel < b.Electron_pt()->size(); iel++) {
    bool is_sig = true;
    if (b.Electron_pt()->at(iel) <= 7.0) is_sig = false;
    float etasc = b.Electron_deltaEtaSC()->at(iel) + b.Electron_eta()->at(iel);
    if (abs(etasc) > 2.5) is_sig = false;
    if (abs(b.Electron_dz()->at(iel))>1.0) is_sig = false;
    if (abs(b.Electron_dxy()->at(iel))>0.5) is_sig = false;
    /*double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
    { 0.604775,  0.628743,   0.896462}};
        int ipt(1), ieta(2);
        if(b.Electron_pt()->at(iel)>10.0) ipt = 0;
        if(fabs(etasc) < 0.8) ieta = 0;
        else if(fabs(etasc) < 1.479) ieta = 1;
	double mva = b.Electron_mvaFall17V2Iso()->at(iel);*/
    bool mva = b.Electron_mvaFall17V2Iso_WPL()->at(iel);
    if (mva)//< wp[ipt][ieta])
      is_sig = true;
    if (is_sig && b.Electron_pt()->at(iel)>25.0)
      el_sig_.push_back(3.0);
    else if(is_sig && b.Electron_pt()->at(iel)>15.0)
      el_sig_.push_back(2.0);
    else if(is_sig)
      el_sig_.push_back(1.0);
    else
      el_sig_.push_back(0.0);
  }
  return el_sig_;
}

NamedFunc::VectorType fph_sig(const Baby &b){
  std::vector<double> ph_sig_;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10. || abs(b.Photon_eta()->at(iph)) > 2.5 || !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) ph_sig_.push_back(0.0);
    //else if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0 /*&& (b.Photon_mvaID_WP90->at(iph) && !b.Photon_mvaID_WP80->at(iph))*/) ph_sig_.push_back(1.0);
    else if ((((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58))) && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0) ph_sig_.push_back(1.0);
    else ph_sig_.push_back(0.0);
  }
  return ph_sig_;
}

//Photon cut flow
NamedFunc::ScalarType fph_sig_pt(const Baby &b){
  int nph = 0;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10. || !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) nph += 0;
    else if (b.Photon_pt()->at(iph) > 15.) nph += 1;
    else nph +=0;
  }
  return (nph >=1);
}

NamedFunc::ScalarType fph_sig_eta(const Baby &b){
  int nph = 0;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10.|| !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) nph += 0;
    else if (abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_pt()->at(iph) > 15.) nph += 1;
    else nph +=0;
  }
  return (nph >=1);
}

NamedFunc::ScalarType fph_sig_barrel(const Baby &b){
  int nph = 0;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10.|| !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) nph += 0;
    else if (abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_pt()->at(iph) > 15. && abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4) nph += 1;
    else nph +=0;
  }
  return (nph >=1);
}


NamedFunc::ScalarType fph_sig_endcap(const Baby &b){
  int nph = 0;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10.|| !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) nph += 0;
    else if (abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_pt()->at(iph) > 15. && abs(b.Photon_eta()->at(iph))>1.566 && b.Photon_mvaID()->at(iph) > -0.58) nph += 1;
    else nph +=0;
  }
  return (nph >=1);
}

NamedFunc::ScalarType fph_sig_eleveto(const Baby &b){
  int nph = 0;
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 10.|| !(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) nph += 0;
    else if (abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_pt()->at(iph) > 15. && ((abs(b.Photon_eta()->at(iph))>1.566 && b.Photon_mvaID()->at(iph) > -0.58) || (abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)) &&  b.Photon_electronVeto()->at(iph)) nph += 1;
    else nph +=0;
  }
  return (nph >=1);
}

//Electron cut flow
NamedFunc::ScalarType fel_lead(const Baby &b){
  std::vector<double> el_sig_ = fel_sig(b);
  for (unsigned i(0); i<el_sig_.size();i++) {
    if(el_sig_[i]>2.5) return true;
  }
  return false;
}

NamedFunc::ScalarType fel_sublead(const Baby &b){
  std::vector<double> el_sig_ = fel_sig(b);
  int num1 = 0;
  int num2 = 0;
  for (unsigned i(0); i<el_sig_.size();i++) {
    if(el_sig_[i]>2.5) num1 += 1;
    if(el_sig_[i]>1.5) num2 += 1;
  }
  return (num1>=1 && num2>=2);
}

//Muon cut flow
NamedFunc::ScalarType fmu_lead(const Baby &b){
  std::vector<double> mu_sig_ = fmu_sig(b);
  for (unsigned i(0); i<mu_sig_.size();i++) {
    if(mu_sig_[i]>2.5) return true;
  }
  return false;
}

NamedFunc::ScalarType fmu_sublead(const Baby &b){
  std::vector<double> mu_sig_ = fmu_sig(b);
  int num1 = 0;
  int num2 = 0;
  for (unsigned i(0); i<mu_sig_.size();i++) {
    if(mu_sig_[i]>2.5) num1 += 1;
    if(mu_sig_[i]>1.5) num2 += 1;
  }
  return (num1>=1 && num2>=2);
}


NamedFunc::VectorType fele_dR_list(const Baby &b){
  TLorentzVector e, p;
  std::vector<double> el_sig_ = fel_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> dRlist;
  double mindr = 9999.0;
  for (unsigned j(0); j<ph_sig_.size(); j++){
    if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
    else{
      dRlist.push_back(0.0);
      continue;
    }
    for (unsigned i(0); i<el_sig_.size();i++) {
      if (el_sig_[i]>0.5)  e.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
      else continue;
      double tempdr = abs(e.DeltaR(p));
      if (tempdr<mindr) mindr = tempdr;
    }
    if (mindr>0.4) dRlist.push_back(1.0);
    else dRlist.push_back(0.0);
  }
  return dRlist;
}
NamedFunc::VectorType fmuon_dR_list(const Baby &b){
  TLorentzVector m, p;
  std::vector<double> mu_sig_ = fmu_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> dRlist;
  double mindr = 9999.0;
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
    else {
      dRlist.push_back(0.0);
      continue;
    }
    for (unsigned i(0); i<mu_sig_.size();i++) {
      if (mu_sig_[i]>0.5)  m.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
      else continue;
      double tempdr = abs(m.DeltaR(p));
      if (tempdr<mindr) mindr = tempdr;
    }
    if (mindr>0.4) dRlist.push_back(1.0);
    else dRlist.push_back(0.0);
  }
  return dRlist;
}

NamedFunc::VectorType fphoton_flags(const Baby &b){
  std::vector<double> flag;
  std::bitset<15> bit;
  int index = 0;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      index = j;
      break;
    }
    else continue;
  }
  int gind;
  gind =  b.Photon_genPartIdx()->at(index);
  if (gind<0) return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  else {
    bit = std::bitset<15>(b.GenPart_statusFlags()->at(gind));
    for(int i(0); i<15; i++){
      if (bit[i] == 1) flag.push_back(1.0);
      else flag.push_back(0.0);
    }
    return flag;
  }
}

NamedFunc::ScalarType fmindR(const Baby &b){
  TLorentzVector e, m, p;
  std::vector<double> el_sig_ = fel_sig(b);
  std::vector<double> mu_sig_ = fmu_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  double mindr = 9999.0;
  double tempdr = 9999.0;
  for (unsigned j(0); j<ph_sig_.size(); j++){
    if (ph_sig_[j] > 0.5) p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
    else continue;
    for (unsigned i(0); i<el_sig_.size();i++) {
      if (el_sig_[i]>0.5)  e.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
      else continue;
      tempdr = abs(e.DeltaR(p));
      if (tempdr<mindr) mindr = tempdr;
    }
    for (unsigned i(0); i<mu_sig_.size();i++) {
      if (mu_sig_[i]>0.5)  m.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
      else continue;
      tempdr = abs(m.DeltaR(p));
      if (tempdr<mindr) mindr = tempdr;
    }
  }
  return mindr;
}

NamedFunc::ScalarType fpdgid (const Baby &b) {
  int index = -1;
  TLorentzVector pg, p;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  double mindr = 999.;
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      /* index = j;
	 break;
      }
      else continue;
    }
    if ( b.Photon_genPartIdx()->at(index)<0) {*/
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
  return abs(b.GenPart_pdgId()->at(index));
}

NamedFunc::ScalarType ffc_pdgid (const Baby &b) {
  int index = -1;
  TLorentzVector pg, p;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  double mindr = 999.;
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      /* index = j;
         break;
      }
      else continue;
    }
    if ( b.Photon_genPartIdx()->at(index)<0) {*/
      for (int i(0); i<b.nGenPart();i++){
	std::bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(i));
	if (mc_statusFlags[12]){
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
  return abs(b.GenPart_pdgId()->at(index));
}

NamedFunc::ScalarType fmother_pdgid (const Baby &b){
  int index = -1;
  TLorentzVector pg, p;
  double mindr = 999.;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      /*index = j;
      break;
    }
    else continue;
  }
  if ( b.Photon_genPartIdx()->at(index)<0) return -1;*/
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      for (int i(0); i<b.nGenPart();i++){
        pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
        double tempdr = abs(pg.DeltaR(p));
        if (tempdr<mindr) {
          mindr = tempdr;
          index = i;
        }
      }
      if (mindr > 0.1) index= -1;
      break;
    }
  }
  if (index<0) return -1;
  else if (b.GenPart_genPartIdxMother()->at(index)<0) return -2;
  else return abs(b.GenPart_pdgId()->at(b.GenPart_genPartIdxMother()->at(index)));
}

NamedFunc::ScalarType ffc_mother_pdgid (const Baby &b){
  int index = -1;
  TLorentzVector pg, p;
  double mindr = 999.;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      /*index = j;
      break;
    }
    else continue;
  }
  if ( b.Photon_genPartIdx()->at(index)<0) return -1;*/
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      for (int i(0); i<b.nGenPart();i++){
	std::bitset<15> mc_statusFlags(b.GenPart_statusFlags()->at(i));
	if (mc_statusFlags[12]){
	  pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	  double tempdr = abs(pg.DeltaR(p));
	  if (tempdr<mindr) {
	    mindr = tempdr;
	    index = i;
	  }
	}
      }
      if (mindr > 0.1) index= -1;
      break;
    }
  }
  if (index<0) return -1;
  else if (b.GenPart_genPartIdxMother()->at(index)<0) return -2;
  else return abs(b.GenPart_pdgId()->at(b.GenPart_genPartIdxMother()->at(index)));
}

NamedFunc::ScalarType fphoton_truth_dR(const Baby &b){
  TLorentzVector p, pg;
  double mindr = 999.;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      /*index = j;
      break;
    }
    else continue;
  }
  int gind = b.Photon_genPartIdx()->at(index);
  if (gind < 0) return -0.1;*/
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      for (int i(0); i<b.nGenPart();i++){
	pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	double tempdr = abs(pg.DeltaR(p));
	if (tempdr<mindr) {
	  mindr = tempdr;
	}
      }
      break;
    }
  }
  return mindr;
}

NamedFunc::ScalarType fphoton_cone(const Baby &b){
  int index = -1;
  TLorentzVector pg, p, c;
  double mindr = 999.;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  
  if (fmother_pdgid(b) == -2){
    for (unsigned j(0); j<b.Photon_pt()->size(); j++){
      if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
	p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
	/* index = j;
         break;
      }
      else continue;
    }
    if ( b.Photon_genPartIdx()->at(index)<0) {*/
	for (int i(0); i<b.nGenPart();i++){
	  pg.SetPtEtaPhiM(b.GenPart_pt()->at(i), b.GenPart_eta()->at(i),b.GenPart_phi()->at(i),b.GenPart_mass()->at(i));
	  double tempdr = abs(pg.DeltaR(p));
	  if (tempdr<mindr) {
	    mindr = tempdr;
	    index = i;
	  }
	}
	break;
      }
    }
    /*
    pg.SetPtEtaPhiM(b.GenPart_pt()->at(index), b.GenPart_eta()->at(index),b.GenPart_phi()->at(index),b.GenPart_mass()->at(index));
    for (int k(0); k<b.nGenPart();k++){
      c.SetPtEtaPhiM(b.GenPart_pt()->at(k), b.GenPart_eta()->at(k),b.GenPart_phi()->at(k),b.GenPart_mass()->at(k));
      if (c.DeltaR(pg)<=0.1) std::cout << abs(b.GenPart_pdgId()->at(k)) << " ";
    }
  }
  std::cout << std::endl << std::endl;*/
  }
  return b.GenPart_pt()->at(index);
}

NamedFunc::ScalarType fphoton_flags0(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[0] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags1(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[1] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags2(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[2] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags3(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[3] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags4(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[4] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags5(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[5] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags6(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[6] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags7(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[7] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags8(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[8] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags9(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[9] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags10(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[10] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags11(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[11] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags12(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[12] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags13(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[13] > 0.5) return true;
  else return false;
}

NamedFunc::ScalarType fphoton_flags14(const Baby &b){
  std::vector<double> flag_ = fphoton_flags(b);
  if (flag_[14] > 0.5) return true;
  else return false;
}


NamedFunc::ScalarType fnph(const Baby &b){
  int nph_sig = 0;
  std::vector<double> mudR = fmuon_dR_list(b);
  std::vector<double> eledR = fele_dR_list(b);
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 15.0) continue;
    if (abs(b.Photon_eta()->at(iph)) > 2.5) continue;
    if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) continue;
    // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
    //if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0)
    if (((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58)) && b.Photon_electronVeto()->at(iph) && mudR[iph]>0.5 && eledR[iph]>0.5)
      nph_sig++;
  }
  return nph_sig;
}

NamedFunc::ScalarType fjk1nph(const Baby &b){
  int nph_sig = 0;
  std::vector<double> mudR = fmuon_dR_list(b);
  std::vector<double> eledR = fele_dR_list(b);
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 15.0) continue;
    if (abs(b.Photon_eta()->at(iph)) > 2.5) continue;
    if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) continue;
    // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
    //if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0)
    if (((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58)) && b.Photon_electronVeto()->at(iph) && eledR[iph]>0.5)
      nph_sig++;
  }
  return nph_sig;
}

NamedFunc::ScalarType fjk2nph(const Baby &b){
  int nph_sig = 0;
  std::vector<double> mudR = fmuon_dR_list(b);
  std::vector<double> eledR = fele_dR_list(b);
  for (unsigned iph = 0; iph < b.Photon_pt()->size(); iph++) {
    if (b.Photon_pt()->at(iph) <= 15.0) continue;
    if (abs(b.Photon_eta()->at(iph)) > 2.5) continue;
    if (!(b.Photon_isScEtaEB()->at(iph) || b.Photon_isScEtaEE()->at(iph))) continue;
    // if (b.Photon_pfRelIso03_all()->at(iph)==PhotonRelIsoCut) continue; // no isolation cut in 2016...?
    //if (b.Photon_mvaID()->at(iph) > 0.2 && b.Photon_electronVeto()->at(iph) && b.Photon_pt()->at(iph) > 15.0)
    if (((abs(b.Photon_eta()->at(iph))<1.4442 && b.Photon_mvaID()->at(iph) > -0.4)||(abs(b.Photon_eta()->at(iph))>1.566 && abs(b.Photon_eta()->at(iph))<2.5 && b.Photon_mvaID()->at(iph) > -0.58)) && b.Photon_electronVeto()->at(iph) && mudR[iph]>0.5)
      nph_sig++;
  }
  return nph_sig;
}

NamedFunc::ScalarType fnmu(const Baby &b){
  int nmu_sig = 0;
  std::vector<double> mu_sig_ = fmu_sig(b);
  for (double imu_sig_ : mu_sig_) {
    if (imu_sig_ > 0.5)
      nmu_sig++;
  }
  return nmu_sig;
}
NamedFunc::ScalarType fnel(const Baby &b){
  int nel_sig = 0;
  std::vector<double> el_sig_ = fel_sig(b);
  for (double iel_sig_ : el_sig_) {
    if (iel_sig_ > 0.5)
      nel_sig++;
  }
  return nel_sig;
}

NamedFunc::ScalarType fmuon_mll(const Baby &b){
  TLorentzVector m1, m2, mll;
  std::vector<double> mu_sig_ = fmu_sig(b);
  double closemass = 0.0;
  double closediff = 999.;
  for (unsigned i(0); i<mu_sig_.size();i++) {
    for (unsigned j(i+1); j<mu_sig_.size();j++) {
      if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
	m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
	mll = m1 + m2;
	double massdiff = abs(mll.M()-91.1876);
	if (massdiff < closediff) {
	  closemass = mll.M();
	  closediff = massdiff;
	}
      }
      else continue;
    }
  }
  return closemass;
}

NamedFunc::ScalarType fmuon_mll_more(const Baby &b){
  TLorentzVector m1, m2, mll;
  std::vector<double> mu_sig_ = fmu_sig(b);
  double closemass = 0.0;
  double closediff = 999.;
  for (unsigned i(0); i<mu_sig_.size();i++) {
    for (unsigned j(i+1); j<mu_sig_.size();j++) {
      if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
	if ((b.Muon_pt()->at(i)>10.0 && b.Muon_pt()->at(j)>20.0)||(b.Muon_pt()->at(i)>20.0 && b.Muon_pt()->at(j)>10.0)){
	  m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	  m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
	  mll = m1 + m2;
	  double massdiff = abs(mll.M()-91.1876);
	  if (massdiff < closediff) {
	    closemass = mll.M();
	    closediff = massdiff;
	  }
	}
	else continue;
      }
      else continue;
    }
  }
  return closemass;
}

NamedFunc::ScalarType fele_mll (const Baby &b){
    TLorentzVector e1, e2, mll;
    std::vector<double> el_sig_ = fel_sig(b);
    double closemass = 0.0;
    double closediff = 999.;
    for (unsigned i(0); i<el_sig_.size();i++) {
      for (unsigned j(i+1); j<el_sig_.size();j++) {
	if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
	  e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	  e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
	  mll = e1 + e2;
	  double massdiff = abs(mll.M()-91.1876);
	  if (massdiff < closediff){
	    closemass = mll.M();
	    closediff = massdiff;
	  }
	}
	else continue;
      }
    }
    return closemass;
}

NamedFunc::ScalarType fele_mll_more (const Baby &b){
  TLorentzVector e1, e2, mll;
  std::vector<double> el_sig_ = fel_sig(b);
  double closemass = 0.0;
  double closediff = 999.;
  for (unsigned i(0); i<el_sig_.size();i++) {
    for (unsigned j(i+1); j<el_sig_.size();j++) {
      if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
	if ((b.Electron_pt()->at(i)>15.0 && b.Electron_pt()->at(j)>25.0) || (b.Electron_pt()->at(i)>25.0 && b.Electron_pt()->at(j)>15.0)){
	  e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	  e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
	  mll = e1 + e2;
	  double massdiff = abs(mll.M()-91.1876);
	  if (massdiff < closediff){
	    closemass = mll.M();
	    closediff = massdiff;
	  }
	}
	else continue;
      }
      else continue;
    }
  }
  return closemass;
}

NamedFunc::ScalarType fmuon_mllg(const Baby &b){
  TLorentzVector m, m1, m2, mll, m1c, m2c, p, mllg;
  std::vector<double> mu_sig_ = fmu_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  std::vector<double> edR = fele_dR_list(b);
  double closediff = 999.;
  for (unsigned i(0); i<mu_sig_.size();i++) {
    for (unsigned j(i+1); j<mu_sig_.size();j++) {
      if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
	m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
	mll = m1 + m2;
	double massdiff = abs(mll.M()-91.1876);
	if (massdiff < closediff) {
	  closediff = massdiff;
	  m1c = m1;
	  m2c = m2;
	}
      }
      else continue;
    }
  }
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && mdR[j]>0.5 && edR[j]>0.5){
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      mllg = m1c + m2c + p;
      break;
    }
  }
  return mllg.M();
}

NamedFunc::ScalarType fmuon_mllg_more(const Baby &b){
  TLorentzVector m, m1, m2, mll, m1c, m2c, p, mllg;
  std::vector<double> mu_sig_ = fmu_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  std::vector<double> edR = fele_dR_list(b);
  double closediff = 999.;
  for (unsigned i(0); i<mu_sig_.size();i++) {
    for (unsigned j(i+1); j<mu_sig_.size();j++) {
      if (mu_sig_[i] > 0.5 && mu_sig_[j] > 0.5){
	if ((b.Muon_pt()->at(i)>20.0 && b.Muon_pt()->at(j)>10.0) || (b.Muon_pt()->at(i)>10.0 && b.Muon_pt()->at(j)>20.0)){
	  m1.SetPtEtaPhiM(b.Muon_pt()->at(i), b.Muon_eta()->at(i), b.Muon_phi()->at(i), b.Muon_mass()->at(i));
	  m2.SetPtEtaPhiM(b.Muon_pt()->at(j), b.Muon_eta()->at(j), b.Muon_phi()->at(j), b.Muon_mass()->at(j));
	  mll = m1 + m2;
	  double massdiff = abs(mll.M()-91.1876);
	  if (massdiff < closediff) {
	    closediff = massdiff;
	    m1c = m1;
	    m2c = m2;
	  }
	}
	else continue;
      }
      else continue;
    }
  }
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && mdR[j]>0.5 && edR[j]>0.5){
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      mllg = m1c + m2c + p;
      break;
    }
  }
  return mllg.M();
}

NamedFunc::ScalarType fele_mllg(const Baby &b){
  TLorentzVector e, e1, e2,mll,e1c,e2c, p, mllg;
  std::vector<double> el_sig_ = fel_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  double closediff = 999.;
  for (unsigned i(0); i<el_sig_.size();i++) {
    for (unsigned j(i+1); j<el_sig_.size();j++) {
      if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
	e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
	mll = e1 + e2;
	double massdiff = abs(mll.M()-91.1876);
	if (massdiff < closediff){
	  closediff = massdiff;
	  e1c = e1;
	  e2c = e2;
	}
      }
      else continue;
    }
  }
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j]>0.5 && mdR[j]>0.5){
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      mllg = e1c + e2c + p;
      break;
    }
  }
  return mllg.M();
}

NamedFunc::ScalarType fele_mllg_more(const Baby &b){
  TLorentzVector e, e1, e2,mll,e1c,e2c, p, mllg;
  std::vector<double> el_sig_ = fel_sig(b);
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  double closediff = 999.;
  for (unsigned i(0); i<el_sig_.size();i++) {
    for (unsigned j(i+1); j<el_sig_.size();j++) {
      if (el_sig_[i] > 0.5 && el_sig_[j] > 0.5 ){
	if((b.Electron_pt()->at(i)>25.0 && b.Electron_pt()->at(j)>15.0) || (b.Electron_pt()->at(i)>15.0 && b.Electron_pt()->at(j)>25.0)){ 
	    e1.SetPtEtaPhiM(b.Electron_pt()->at(i), b.Electron_eta()->at(i), b.Electron_phi()->at(i), b.Electron_mass()->at(i));
	    e2.SetPtEtaPhiM(b.Electron_pt()->at(j), b.Electron_eta()->at(j), b.Electron_phi()->at(j), b.Electron_mass()->at(j));
	    mll = e1 + e2;
	    double massdiff = abs(mll.M()-91.1876);
	    if (massdiff < closediff){
	      closediff = massdiff;
	      e1c = e1;
	      e2c = e2;
	    }
	  }
	  else continue;
      }
      else continue;
    }
  }
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j]>0.5 && mdR[j]>0.5){
      p.SetPtEtaPhiM(b.Photon_pt()->at(j), b.Photon_eta()->at(j), b.Photon_phi()->at(j), b.Photon_mass()->at(j));
      mllg = e1c + e2c + p;
      break;
    }
  }
  return mllg.M();
}
/*NamedFunc::ScalarType fphoton_truth_dR(const Baby &b){
  int index = 0;
  TLorentzVector p, pg;
  std::vector<double> ph_sig_ = fph_sig(b);
  std::vector<double> edR = fele_dR_list(b);
  std::vector<double> mdR = fmuon_dR_list(b);
  for (unsigned j(0); j<b.Photon_pt()->size(); j++){
    if (ph_sig_[j] > 0.5 && edR[j] > 0.5 && mdR[j] > 0.5) {
      index = j;
      break;
    }
    else continue;
  }
  int gind = b.Photon_genPartIdx()->at(index);
  if (gind < 0) return -0.1;
  p.SetPtEtaPhiM(b.Photon_pt()->at(index), b.Photon_eta()->at(index), b.Photon_phi()->at(index), b.Photon_mass()->at(index));
  pg.SetPtEtaPhiM(b.GenPart_pt()->at(gind), b.GenPart_eta()->at(gind), b.GenPart_phi()->at(gind),b.GenPart_mass()->at(gind));
  return abs(pg.DeltaR(p));
}
*/
