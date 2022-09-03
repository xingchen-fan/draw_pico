#ifndef XC_FUNC
#define XC_FUNC

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
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"

NamedFunc::VectorType fmu_sig(const Baby &b);
NamedFunc::VectorType fel_sig(const Baby &b);
NamedFunc::VectorType fph_sig(const Baby &b);
NamedFunc::VectorType fele_dR_list(const Baby &b);
//NamedFunc ele_dR_list("ele_dR_list", fele_dR_list);
NamedFunc::VectorType fmuon_dR_list(const Baby &b);
NamedFunc::ScalarType fmindR(const Baby &b);
//inline NamedFunc muon_dR_list("muon_dR_list", fmuon_dR_list);
NamedFunc::ScalarType fnph(const Baby &b);
NamedFunc::ScalarType fjk1nph(const Baby &b);
NamedFunc::ScalarType fjk2nph(const Baby &b);
NamedFunc::ScalarType fnmu(const Baby &b);
NamedFunc::ScalarType fnel(const Baby &b);
NamedFunc::VectorType fphoton_flags(const Baby &b);
NamedFunc::ScalarType fphoton_flags0(const Baby &b);
NamedFunc::ScalarType fphoton_flags1(const Baby &b);
NamedFunc::ScalarType fphoton_flags2(const Baby &b);
NamedFunc::ScalarType fphoton_flags3(const Baby &b);
NamedFunc::ScalarType fphoton_flags4(const Baby &b);
NamedFunc::ScalarType fphoton_flags5(const Baby &b);
NamedFunc::ScalarType fphoton_flags6(const Baby &b);
NamedFunc::ScalarType fphoton_flags7(const Baby &b);
NamedFunc::ScalarType fphoton_flags8(const Baby &b);
NamedFunc::ScalarType fphoton_flags9(const Baby &b);
NamedFunc::ScalarType fphoton_flags10(const Baby &b);
NamedFunc::ScalarType fphoton_flags11(const Baby &b);
NamedFunc::ScalarType fphoton_flags12(const Baby &b);
NamedFunc::ScalarType fphoton_flags13(const Baby &b);
NamedFunc::ScalarType fphoton_flags14(const Baby &b);

NamedFunc::ScalarType fpdgid (const Baby &b);
NamedFunc::ScalarType ffc_pdgid (const Baby &b);
NamedFunc::ScalarType fmother_pdgid (const Baby &b);
NamedFunc::ScalarType ffc_mother_pdgid (const Baby &b);
NamedFunc::ScalarType fphoton_cone(const Baby &b);
//inline NamedFunc ele_pdgid("ele_pdgid", fele_pdgid);
//inline NamedFunc muon_pdgid("muon_pdgid", fmuon_pdgid);
NamedFunc::ScalarType fmuon_mll(const Baby &b);
NamedFunc::ScalarType fele_mll(const Baby &b);
NamedFunc::ScalarType fmuon_mllg(const Baby &b);
NamedFunc::ScalarType fele_mllg(const Baby &b);

NamedFunc::ScalarType fphoton_truth_dR(const Baby &b);
NamedFunc::ScalarType fmuon_mll_more(const Baby &b);
NamedFunc::ScalarType fele_mll_more (const Baby &b);
NamedFunc::ScalarType fmuon_mllg_more(const Baby &b);
NamedFunc::ScalarType fele_mllg_more(const Baby &b);

NamedFunc::ScalarType fph_sig_pt(const Baby &b);
NamedFunc::ScalarType fph_sig_eta(const Baby &b);
NamedFunc::ScalarType fph_sig_barrel(const Baby &b);
NamedFunc::ScalarType fph_sig_endcap(const Baby &b);
NamedFunc::ScalarType fph_sig_eleveto(const Baby &b);


NamedFunc::ScalarType fel_lead(const Baby &b);
NamedFunc::ScalarType fel_sublead(const Baby &b);
NamedFunc::ScalarType fmu_lead(const Baby &b);
NamedFunc::ScalarType fmu_sublead(const Baby &b);
#endif


