//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>

//OUR OWN CLASSES TO READ THE TREE
#include "MassParticles.h"
#include "MyBaseJet.h"
#include "MyBeamSpot.h"
#include "MyCaloJet.h"
#include "MyCastorDigi.h"
#include "MyCastorJet.h"
#include "MyCastorRecHit.h"
#include "MyDiJet.h"
#include "MyElectron.h"
#include "MyEvtId.h"
#include "MyFwdGap.h"
#include "MyGenJet.h"
#include "MyGenKin.h"
#include "MyGenMet.h"
#include "MyGenPart.h"
#include "MyHLTrig.h"
#include "MyJet.h"
#include "MyL1Trig.h"
#include "MyL1TrigOld.h"
//#include "MyMITEvtSel.h"
#include "MyMet.h"
#include "MyMuon.h"
#include "MyPFCand.h"
#include "MyPFJet.h"
#include "MyPUSumInfo.h"
#include "MyPart.h"
#include "MySimVertex.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "MyFSCHit.h"
#include "MyFSCDigi.h"

// TOTEM data formats
#include "T1Event.h"
#include "T2Event.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"

#include "analysis_tools.h"
#include "rp_aperture_config.h"

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

void template_UATreeReader(vector<string> const& fileNames, string const& outputFileName = "output.root", const Int_t nevt_max = 100){
  
  bool isMC  = false;
  bool verbose = false;
  string treeName = (!isMC) ? "cms_totem" : "evt";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  if(isMC) jetCorrName = "ak5PFL2L3"; 
  double ptJetMin = 20.0;
  double etaJetMax = 4.0;
  double etaMaxThreshold = 2.0;

  bool selectBunchCrossing = false;
  bool selectVertex = true;
  bool selectJets = true;
  bool selectEtaMax = false;
  bool selectEtaMin = false;
  bool selectZeroHitsT2Plus = false;
  bool selectZeroHitsT2Minus = false;
  bool selectSingleArmRecProton = true;
  bool selectDoubleArmRecProton = false;
  bool selectElastic = false;
  bool selectNonElastic = false;

  bool selectHLTORorAND = true; // OR=true - AND=false
  vector<string> selectHLTPathNames;
  selectHLTPathNames.push_back("HLT_L1DoubleJet20part1_v1");
  selectHLTPathNames.push_back("HLT_L1DoubleJet20part2_v1");
  //selectHLTPathNames.push_back("HLT_L1Tech53_MB_1_v1");
  //selectHLTPathNames.push_back("HLT_L1Tech53_MB_2_v1");
  //selectHLTPathNames.push_back("HLT_L1Tech53_MB_3_v1");
  //selectHLTPathNames.push_back("HLT_ZeroBias_v7");
  //selectHLTPathNames.push_back("HLT_RomanPots_Tech52_v1");
 
  vector<int> bunchCrossingList;
  bunchCrossingList.push_back(27);
  bunchCrossingList.push_back(649);
  bunchCrossingList.push_back(2991);

  ThresholdsPerRegion thresholdsPFlow;
  thresholdsPFlow[Barrel] = ThresholdsPerType(); 
  thresholdsPFlow[Endcap] = ThresholdsPerType(); 
  thresholdsPFlow[Transition] = ThresholdsPerType(); 
  thresholdsPFlow[Endcap] = ThresholdsPerType(); 
  resetPFThresholds(thresholdsPFlow[Barrel]);
  resetPFThresholds(thresholdsPFlow[Endcap]);
  resetPFThresholds(thresholdsPFlow[Transition]);
  resetPFThresholds(thresholdsPFlow[Forward]);

  thresholdsPFlow[Barrel][MyPFCand::h0]            = make_pair(-1.,1.4);
  thresholdsPFlow[Barrel][MyPFCand::gamma]         = make_pair(-1.,0.9);
  thresholdsPFlow[Endcap][MyPFCand::h0]            = make_pair(-1.,2.7);
  thresholdsPFlow[Endcap][MyPFCand::gamma]         = make_pair(-1.,2.5);
  thresholdsPFlow[Transition][MyPFCand::h0]        = make_pair(-1.,3.8);
  thresholdsPFlow[Transition][MyPFCand::gamma]     = make_pair(-1.,2.5);
  thresholdsPFlow[Transition][MyPFCand::h_HF]      = make_pair(-1.,4.0);
  thresholdsPFlow[Transition][MyPFCand::egamma_HF] = make_pair(-1.,3.5);
  thresholdsPFlow[Forward][MyPFCand::h_HF]         = make_pair(-1.,4.0);
  thresholdsPFlow[Forward][MyPFCand::egamma_HF]    = make_pair(-1.,3.5);

  ThresholdsPerType::const_iterator pfThreshold = thresholdsPFlow[Barrel].begin();
  ThresholdsPerType::const_iterator pfThresholds_end = thresholdsPFlow[Barrel].end(); 
  ostringstream oss;
  oss << "Using the following PF thresholds:\n";
  for(; pfThreshold != pfThresholds_end; ++pfThreshold){
     int key = pfThreshold->first;    
     oss << "  " << key << ": "
                 << "(" << thresholdsPFlow[Barrel][key].first
                 << "," << thresholdsPFlow[Barrel][key].second << ")  "
                 << "(" << thresholdsPFlow[Endcap][key].first
                 << "," << thresholdsPFlow[Endcap][key].second << ")  "
                 << "(" << thresholdsPFlow[Transition][key].first
                 << "," << thresholdsPFlow[Transition][key].second << ")  "
                 << "(" << thresholdsPFlow[Forward][key].first
                 << "," << thresholdsPFlow[Forward][key].second << ")\n";   
  }
  cout << oss.str();
  //==============================
  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  vector<string> hltPathNames;
  hltPathNames.push_back("HLT_L1DoubleEG3_FwdVeto_v1");
  hltPathNames.push_back("HLT_L1DoubleMu0_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20_RomanPotsOR_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part1_v1");
  hltPathNames.push_back("HLT_L1DoubleJet24_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part2_v1");
  hltPathNames.push_back("HLT_L1Tech40_BPTXAND_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_1_v1");
  hltPathNames.push_back("HLT_L1Tech_HF9OR10_v1");
  hltPathNames.push_back("HLT_T1minbias_Tech55_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_2_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_3_v1");
  hltPathNames.push_back("HLT_RomanPots_Tech52_v1");
  hltPathNames.push_back("HLT_L1Tech54_ZeroBias_v1");
  hltPathNames.push_back("HLT_ZeroBias_v7");

  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
  
  vector<string> selections;
  selections.push_back("All");
  selections.push_back("BunchCrossing");
  selections.push_back("HLT");
  selections.push_back("Vertex");
  selections.push_back("Jet");
  selections.push_back("EtaMax");
  selections.push_back("EtaMin");
  selections.push_back("ZeroHitsT2Plus");
  selections.push_back("ZeroHitsT2Minus");
  selections.push_back("SingleArmRP");
  selections.push_back("DoubleArmRP");
  selections.push_back("Elastic");
  selections.push_back("NonElastic");
  int nBinsEventSelection = selections.size();
  histosTH1F["EventSelection"] = new TH1F("EventSelection","EventSelection",nBinsEventSelection,0,nBinsEventSelection);
  for(size_t k = 0; k < selections.size(); ++k)
     histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

  histosTH1F["bunchCrossingNumber"] = new TH1F("bunchCrossingNumber", "bunchCrossingNumber" , 3900 , 0 , 3900);

  histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
  histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

  int nBinsHLT = hltPathNames.size(); 
  histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);
  for(size_t k = 0; k < nBinsHLT; ++k) 
     histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1), hltPathNames[k].c_str() );

  //histosTH1F["vtx_sumpt_max"] = new TH1F("vtx_sumpt_max", "vtx max sum(pT) index" , 30 , 0 , 30);

  histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vtx_ndof"] = new TH1F("vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
  histosTH1F["vtx_chi2"] = new TH1F("vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);

  histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertices" , 30 , 0 , 30);

  histosTH1F["prim_vtx_zpos"] = new TH1F("prim_vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["prim_vtx_xpos"] = new TH1F("prim_vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["prim_vtx_ypos"] = new TH1F("prim_vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);

  histosTH1F["prim_vtx_ndof"] = new TH1F("prim_vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
  histosTH1F["prim_vtx_chi2"] = new TH1F("prim_vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);
  histosTH1F["prim_vtx_chi2n"] = new TH1F("prim_vtx_chi2n", "chi2n(vtx)" , 100 , 0. , 10.);
  histosTH1F["prim_vtx_ntracks"] = new TH1F("prim_vtx_ntracks", "n_{trk}(vtx)" , 30 , 0 , 30);
  histosTH1F["prim_vtx_sumpt"] = new TH1F("prim_vtx_sumpt", "sum(p_{T})(vtx)" , 100 , 0. , 100.);

  histosTH1F["prim_vtx_zpos_after_vtx_sel"] = new TH1F("prim_vtx_zpos_after_vtx_sel", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["prim_vtx_xpos_after_vtx_sel"] = new TH1F("prim_vtx_xpos_after_vtx_sel", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["prim_vtx_ypos_after_vtx_sel"] = new TH1F("prim_vtx_ypos_after_vtx_sel", "y(vtx)" , 150 , -1.5 , 1.5);

  histosTH1F["prim_vtx_ndof_after_vtx_sel"] = new TH1F("prim_vtx_ndof_after_vtx_sel", "ndof(vtx)" , 100 , 0. , 15.);
  histosTH1F["prim_vtx_chi2_after_vtx_sel"] = new TH1F("prim_vtx_chi2_after_vtx_sel", "chi2(vtx)" , 100 , 0. , 10.);
  histosTH1F["prim_vtx_chi2n_after_vtx_sel"] = new TH1F("prim_vtx_chi2n_after_vtx_sel", "chi2n(vtx)" , 100 , 0. , 10.);
  histosTH1F["prim_vtx_ntracks_after_vtx_sel"] = new TH1F("prim_vtx_ntracks_after_vtx_sel", "n_{trk}(vtx)" , 30 , 0 , 30);
  histosTH1F["prim_vtx_sumpt_after_vtx_sel"] = new TH1F("prim_vtx_sumpt_after_vtx_sel", "sum(p_{T})(vtx)" , 100 , 0. , 100.);

  //histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
  histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
  histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
  histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
  histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);
  
  histosTH1F["jet_pt"] = new TH1F("jet_pt", "p_{T}(jet)" , 150 , 0. , 150.);
  histosTH1F["jet_eta"] = new TH1F("jet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
  histosTH1F["jet_phi"] = new TH1F("jet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);

  histosTH1F["leadingJet_pt"] = new TH1F("leadingJet_pt", "p_{T}(jet)" , 150 , 0. , 150.);
  histosTH1F["leadingJet_eta"] = new TH1F("leadingJet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi"] = new TH1F("leadingJet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);
  histosTH1F["leadingJet_jec"] = new TH1F("leadingJet_jec", "JEC" , 100 , 0. , 5.);
  histosTH1F["leadingJet_looseJetId"] = new TH1F("leadingJet_looseJetId", "jet Id" , 2 , 0 , 2);
  histosTH1F["leadingJet_tightJetId"] = new TH1F("leadingJet_tightJetId", "jet Id" , 2 , 0 , 2);
  histosTH1F["leadingJet_pt_after_jet_sel"] = new TH1F("leadingJet_pt_after_jet_sel", "p_{T}(jet)" , 150 , 0. , 150.);
  histosTH1F["leadingJet_eta_after_jet_sel"] = new TH1F("leadingJet_eta_after_jet_sel", "#eta(jet)" , 200 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi_after_jet_sel"] = new TH1F("leadingJet_phi_after_jet_sel", "#phi(jet)" , 200 , -M_PI , M_PI);
  histosTH1F["leadingJet_jec_after_jet_sel"] = new TH1F("leadingJet_jec_after_jet_sel", "JEC" , 100 , 0. , 5.);

  histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , 150 , 0. , 150.);
  histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
  histosTH1F["secondJet_phi"] = new TH1F("secondJet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);
  histosTH1F["secondJet_jec"] = new TH1F("secondJet_jec", "JEC" , 100 , 0. , 5.);
  histosTH1F["secondJet_looseJetId"] = new TH1F("secondJet_looseJetId", "jet Id" , 2 , 0 , 2);
  histosTH1F["secondJet_tightJetId"] = new TH1F("secondJet_tightJetId", "jet Id" , 2 , 0 , 2);
  histosTH1F["secondJet_pt_after_jet_sel"] = new TH1F("secondJet_pt_after_jet_sel", "p_{T}(jet)" , 150 , 0. , 150.);
  histosTH1F["secondJet_eta_after_jet_sel"] = new TH1F("secondJet_eta_after_jet_sel", "#eta(jet)" , 200 , -5.2 , 5.2);
  histosTH1F["secondJet_phi_after_jet_sel"] = new TH1F("secondJet_phi_after_jet_sel", "#phi(jet)" , 200 , -M_PI , M_PI);
  histosTH1F["secondJet_jec_after_jet_sel"] = new TH1F("secondJet_jec_after_jet_sel", "JEC" , 100 , 0. , 5.);

  histosTH1F["pf_etaMax"] = new TH1F("pf_etaMax","#eta^{max}",82,etaBinsHCALBoundaries);
  histosTH1F["pf_etaMin"] = new TH1F("pf_etaMin","#eta^{min}",82,etaBinsHCALBoundaries);
  histosTH1F["pf_deltaEta"] = new TH1F("pf_deltaEta","#Delta#eta",100,0.,10.);
  histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",24,binningEPlusPz);
  histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",24,binningEPlusPz);
  histosTH1F["pf_xiPlus"] = new TH1F("pf_xiPlus","#xi^{+}",200,-1.,1.);
  histosTH1F["pf_xiMinus"] = new TH1F("pf_xiMinus","#xi^{-}",200,-1.,1.);
  histosTH1F["pf_logXiPlus"] = new TH1F("pf_logXiPlus","log(#xi^{+})",200,-5.,0.);
  histosTH1F["pf_logXiMinus"] = new TH1F("pf_logXiMinus","log(#xi^{-})",200,-5.,0.);

  histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
  histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);

  histosTH1F["t2_track_chi2Prob_zplus"] = new TH1F("t2_track_chi2Prob_zplus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zplus"] = new TH1F("t2_track_entryX_zplus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zplus"] = new TH1F("t2_track_entryY_zplus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zplus"] = new TH1F("t2_track_multiplicity_zplus", "n tracks" , 100 , 0 , 100);
  histosTH1F["t2_track_chi2Prob_zminus"] = new TH1F("t2_track_chi2Prob_zminus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zminus"] = new TH1F("t2_track_entryX_zminus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zminus"] = new TH1F("t2_track_entryY_zminus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zminus"] = new TH1F("t2_track_multiplicity_zminus", "n tracks" , 100 , 0 , 100);

  histosTH1F["rp_track_posx_020"] = new TH1F("rp_track_posx_020", "x(RP track)" , 200, -10., 10.);
  histosTH1F["rp_track_posy_020"] = new TH1F("rp_track_posy_020", "y(RP track)" , 500, -50., 50.);

  histosTH1F["rp_track_posx_120"] = new TH1F("rp_track_posx_120", "x(RP track)" , 200, -10., 10.);
  histosTH1F["rp_track_posy_120"] = new TH1F("rp_track_posy_120", "y(RP track)" , 500, -50., 50.);

  histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 200 , -1. , 1.);
  histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",200,-5.,0.);
  histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 200 , 0. , 5.);
  histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 100 , 0. , 100.);
  histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 200 , -1. , 1.);
  histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",200,-5.,0.);
  histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 200 , 0. , 5.);
  histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 100 , 0. , 100.);

  histosTH1F["proton_pair_right_xi"] = new TH1F("proton_pair_right_xi", "#xi" , 200 , -1. , 1.);
  histosTH1F["proton_pair_right_logXi"] = new TH1F("proton_pair_right_logXi","log(#xi)",200,-5.,0.);
  histosTH1F["proton_pair_right_t"] = new TH1F("proton_pair_right_t", "-t" , 200 , 0. , 5.);
  histosTH1F["proton_pair_left_xi"] = new TH1F("proton_pair_left_xi", "#xi" , 200 , -1. , 1.);
  histosTH1F["proton_pair_left_logXi"] = new TH1F("proton_pair_left_logXi","log(#xi)",200,-5.,0.);
  histosTH1F["proton_pair_left_t"] = new TH1F("proton_pair_left_t", "-t" , 200 , 0. , 5.);
  histosTH1F["proton_pair_chi2"] = new TH1F("proton_pair_chi2", "#chi^{2}" , 100 , 0. , 100.);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);

  histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
  histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);

  histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" , 200, 0., 1.);
  histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" , 200, 0., 5.);
  histosTH1F["thx_proton_minus"] = new TH1F("thx_proton_minus", "thx_proton_minus" , 200, -5e-4, 5e-4);
  histosTH1F["thy_proton_minus"] = new TH1F("thy_proton_minus", "thy_proton_minus" , 200, -5e-4, 5e-4);

  //FIXME
  histosTH1F["xi_proton_plus_accepted"] = new TH1F("xi_proton_plus_accepted", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted"] = new TH1F("t_proton_plus_accepted", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["xi_proton_minus_accepted"] = new TH1F("xi_proton_minus_accepted", "xi_proton_minus" , 200, 0., 1.);
  histosTH1F["t_proton_minus_accepted"] = new TH1F("t_proton_minus_accepted", "t_proton_minus" , 200, 0., 5.);

  histosTH1F["xi_proton_plus_accepted_020"] = new TH1F("xi_proton_plus_accepted_020", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted_020"] = new TH1F("t_proton_plus_accepted_020", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["posx_proton_plus_accepted_020"] = new TH1F("posx_proton_plus_accepted_020", "posx_proton_plus" , 200, -0.05, 0.05);
  histosTH1F["posy_proton_plus_accepted_020"] = new TH1F("posy_proton_plus_accepted_020", "posy_proton_plus" , 200, -0.05, 0.05);

  histosTH1F["xi_proton_plus_accepted_021"] = new TH1F("xi_proton_plus_accepted_021", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted_021"] = new TH1F("t_proton_plus_accepted_021", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["xi_proton_plus_accepted_024"] = new TH1F("xi_proton_plus_accepted_024", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted_024"] = new TH1F("t_proton_plus_accepted_024", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["xi_proton_plus_accepted_025"] = new TH1F("xi_proton_plus_accepted_025", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted_025"] = new TH1F("t_proton_plus_accepted_025", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["xi_proton_plus_accepted_120"] = new TH1F("xi_proton_plus_accepted_120", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus_accepted_120"] = new TH1F("t_proton_plus_accepted_120", "t_proton_plus" , 200, 0., 5.);

  histosTH1F["xi_proton_minus_accepted_120"] = new TH1F("xi_proton_minus_accepted_120", "xi_proton_minus" , 200, 0., 1.);
  histosTH1F["t_proton_minus_accepted_120"] = new TH1F("t_proton_minus_accepted_120", "t_proton_minus" , 200, 0., 5.);
  histosTH1F["posx_proton_minus_accepted_120"] = new TH1F("posx_proton_minus_accepted_120", "posx_proton_minus" , 200, -0.05, 0.05);
  histosTH1F["posy_proton_minus_accepted_120"] = new TH1F("posy_proton_minus_accepted_120", "posy_proton_minus" , 200, -0.05, 0.05);

  histosTH1F["xi_proton_minus_accepted_020"] = new TH1F("xi_proton_minus_accepted_020", "xi_proton_minus" , 200, 0., 1.);
  histosTH1F["t_proton_minus_accepted_020"] = new TH1F("t_proton_minus_accepted_020", "t_proton_minus" , 200, 0., 5.);

  map<string,TH2F*> histosTH2F;  
  histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
  histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"] = new TH2F("t2_track_multiplicity_vs_leadingJet_pt","t2_track_multiplicity_vs_leadingJet_pt", 150 , 0. , 150., 100 , 0 , 100);
  histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
  histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

  histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_right_xi_vs_pf_xiMinus"] = new TH2F("proton_right_xi_vs_pf_xiMinus","proton_right_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_right_xi_vs_pf_xiPlus"] = new TH2F("proton_right_xi_vs_pf_xiPlus","proton_right_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_xi_vs_pf_xiMinus"] = new TH2F("proton_left_xi_vs_pf_xiMinus","proton_left_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_xi_vs_pf_xiPlus"] = new TH2F("proton_left_xi_vs_pf_xiPlus","proton_left_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);

  histosTH2F["proton_right_t_vs_leadingJet_pt"] = new TH2F("proton_right_t_vs_leadingJet_pt","proton_right_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);
  histosTH2F["proton_left_t_vs_leadingJet_pt"] = new TH2F("proton_left_t_vs_leadingJet_pt","proton_left_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);

  histosTH2F["rp_track_pos_y_vs_x_020"] = new TH2F("rp_track_pos_y_vs_x_020", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
  histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

  histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 200, -0.05, 0.05, 200, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 200, -0.05, 0.05, 200, -0.05, 0.05);

  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;
  histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energyVsEtaAllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energyVsEtaUndefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaChargedHadron"] = new TH2F("energyVsEtaChargedHadron","energyVsEtaChargedHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energyVsEtaElectron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energyVsEtaMuon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energyVsEtaGamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energyVsEtaNeutralHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energyVsEtaHadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energyVsEtaEGammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  //===================

  //vector<TString>* vfiles = new vector<TString>(1,"merged_reduced_8372_198903_LP_Jets1_1_test_v1.root"); 
  vector<TString>* vfiles = new vector<TString>; 
  //vfiles->push_back( fileName ); 
  for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );
  
  //Declaration of tree and its branches variables
  //TTree* tree = new TTree(treeName.c_str(),"");
  TTree* tree = NULL;
  MyEvtId*           evtId         = NULL;
  MyL1TrigOld*       l1Trig        = NULL;  
  MyHLTrig*          hltTrig       = NULL;
  vector<MyTracks>*  track_coll    = NULL;
  vector<MyVertex>*  vertex_coll   = NULL;
  vector<MyPFJet>*   pfJet_coll    = NULL;
  vector<MyPFCand>*  pFlow_coll    = NULL;
  vector<MyFSCHit>*  fscHits_coll  = NULL;
  vector<MyFSCDigi>* fscDigis_coll = NULL;

  vector<MyGenPart>* genPart       = NULL;
  MyGenKin*          genKin        = NULL;
  //===================
  T2Event* t2_event = NULL; 
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
  map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
  //===================
  
  //===================
  //===================
  if(isMC) rp_aperture_config();
  //===================
  //===================

  int i_tot = 0 , nevt_tot = 0;
  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max_corr ; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    // Access TTree from current file
    tree = (TTree*) file->Get( treeName.c_str() );

    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    // Add branches to TTree ----------------------------------------------------------------------
    if(isMC){
       tree->SetBranchAddress("evtId",&evtId);
       //tree->SetBranchAddress("L1TrigOld",&l1Trig);
       //tree->SetBranchAddress("HLTrig",&hltTrig);
       tree->SetBranchAddress("generalTracks",&track_coll); 
       tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
       tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
       tree->SetBranchAddress("particleFlow",&pFlow_coll);

       tree->SetBranchAddress("genKin",&genKin);
       tree->SetBranchAddress("genPart",&genPart);
    } else{   
       tree->SetBranchAddress("cmsEvtUA",&evtId);
       tree->SetBranchAddress("cmsTrigUA",&l1Trig);
       tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
       tree->SetBranchAddress("cmsTracksUA",&track_coll); 
       tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
       //tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
       tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
       tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
       //tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
       //tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
       tree->SetBranchAddress("branchT2EV.",&t2_event);
       tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
       tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
       tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
       /*//rp_track_info[20] = NULL;
	 tree->SetBranchAddress("track_rp_20.", &rp_track_info[20]);*/
       std::vector<unsigned int> rp_list;
       rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(24); rp_list.push_back(25);
       rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(124); rp_list.push_back(125);
       char br_name[200];
       for (unsigned int a = 0; a < 2; ++a) {
	  int s = 2;
	  for (unsigned int r = 0; r < 6; r++) {
	     unsigned int id = 100 * a + 10 * s + r;
	     if( std::find(rp_list.begin(), rp_list.end(), id) == rp_list.end() ) continue;

	     sprintf(br_name, "track_rp_%u.", id);
	     std::cout << br_name << std::endl;
	     tree->SetBranchAddress(br_name, &rp_track_info[id]);
	  }
       }
       /*readRPBranches(tree,
	 rec_proton_left,
	 rec_proton_right,
	 rec_proton_pair,
	 rp_track_info,
	 rp_digi_info,
	 rp_par_patterns_info,
	 rp_nonpar_patterns_info,
	 rp_multi_track_info);*/
    }
    
    /*//Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/
  
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){
    
      //printing the % of events done every 10k evts
      if( ((i_tot+1) % 5000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
      //if( ((i_tot+1) % 100) == 0 ) cout << (i_tot+1) << " done" << endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      //continue;

      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = 1.;

      histosTH1F["EventSelection"]->Fill( "All", event_weight );

      if(!isMC){
	 //int orbitNumber = evtId->Orbit;
	 int bunchCrossingNumber = evtId->Bunch;
	 histosTH1F["bunchCrossingNumber"]->Fill( bunchCrossingNumber, event_weight );

	 if( selectBunchCrossing &&
	       find(bunchCrossingList.begin(),bunchCrossingList.end(),bunchCrossingNumber) == bunchCrossingList.end() )
	    continue;

	 histosTH1F["EventSelection"]->Fill( "BunchCrossing", event_weight );

	 for (int itrig = 0 ; itrig < 128 ; ++itrig){
	    if( l1Trig->PhysTrigWord[itrig] == 1) 
	       histosTH1F["decisionPhysTrig"]->Fill( itrig, event_weight );
	 }

	 for (int itrig = 0 ; itrig < 64 ; ++itrig){
	    if( l1Trig->TechTrigWord[itrig] == 1 )
	       histosTH1F["decisionTechTrig"]->Fill( itrig, event_weight );
	 }

	 map<string,bool>::iterator it_hlt = (*hltTrig).HLTmap.begin();
	 map<string,bool>::iterator it_hlt_end = (*hltTrig).HLTmap.end();
	 for(; it_hlt != it_hlt_end; ++it_hlt){
	    string const& hltName = it_hlt->first;
	    vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
	    if(it_pos != hltPathNames.end()){
	       /*size_t idx = it_pos - hltPathNames.begin();
		 if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( idx, event_weight );*/
	       if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );
	    }
	 }

	 bool select_HLT = false;
	 bool hlt_result = false;
	 for(vector<string>::const_iterator it_pathname_sel = selectHLTPathNames.begin();
	       it_pathname_sel != selectHLTPathNames.end(); ++it_pathname_sel){
	    /*bool accept = (*hltTrig).HLTmap[ *it_pathname_sel ];
	      if(accept) { select_HLT = true; break; }*/

	    if( it_pathname_sel == selectHLTPathNames.begin() ) hlt_result = selectHLTORorAND ? false : true;
	    bool accept = (*hltTrig).HLTmap[ *it_pathname_sel ];
	    if(selectHLTORorAND) hlt_result = hlt_result || accept;
	    else                 hlt_result = hlt_result && accept;
	 }
	 select_HLT = hlt_result;
	 if(!select_HLT) continue;

	 histosTH1F["EventSelection"]->Fill( "HLT", event_weight );
      }

      //-------------------------------------------------------------------------------------------------
      //filling pt distribution for the generated particles
      //ie those from pythia generator, without reconstruction
      /*if(isMC){
        for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
          pt_gen->Fill(p->Pt());
      }*/
      
      if(isMC){
	 // Gen. Particles
	 double genEPlusPz = 0;
	 double genEMinusPz = 0;
	 double proton_pi = 4000.;
	 double proton_pz_plus = -999.;
	 double proton_px_plus = -999.;
	 double proton_py_plus = -999.;
	 double proton_energy_plus = 0.;
	 double proton_pz_minus = 999.;
	 double proton_px_minus = 999.;
	 double proton_py_minus = 999.;
	 double proton_energy_minus = 0.;

	 for(vector<MyGenPart>::iterator it_genpart = genPart->begin();
	       it_genpart != genPart->end(); ++it_genpart){

	    double eta_gen = it_genpart->Eta();
	    int status = it_genpart->status;
	    int id = it_genpart->pdgId;

	    if (status != 1) continue; // final state particles
	    double energy_gen = it_genpart->Energy();
	    double px_gen = it_genpart->Px();
	    double py_gen = it_genpart->Py();
	    double pz_gen = it_genpart->Pz();

	    genEPlusPz +=  (energy_gen + pz_gen);
	    genEMinusPz += (energy_gen - pz_gen);

	    if (id != 2212) continue; // select protons 

	    double proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
	    double pz_cut = 0.75*proton_pi;
	    if ( fabs(pz_gen) < pz_cut) continue;

	    if (pz_gen > proton_pz_plus)  { 
	       proton_pz_plus = pz_gen; proton_energy_plus = energy_gen; 
	       proton_px_plus = px_gen; proton_py_plus = py_gen;
	    } 
	    if (pz_gen < proton_pz_minus) { 
	       proton_pz_minus = pz_gen; proton_energy_minus = energy_gen; 
	       proton_px_minus = px_gen; proton_py_minus = py_gen;
	    } 
	 }
	 double xi_plus_gen = genEPlusPz/8000.;
	 double xi_minus_gen = genEMinusPz/8000.;
	 double xi_proton_plus = -1.;
	 double xi_proton_minus = -1.;
	 double t_proton_plus = 0.;
	 double t_proton_minus = 0.;
	 double thx_proton_plus = 0.; 
	 double thy_proton_plus = 0.; 
	 double thx_proton_minus = 0.; 
	 double thy_proton_minus = 0.; 

	 bool proton_minus_rp_accept_120 = false;
	 bool proton_minus_rp_accept_121 = false;
	 bool proton_minus_rp_accept_124 = false;
	 bool proton_minus_rp_accept_125 = false;
	 bool proton_minus_rp_accept_020 = false;

	 bool proton_plus_rp_accept_020 = false;
	 bool proton_plus_rp_accept_021 = false;
	 bool proton_plus_rp_accept_024 = false;
	 bool proton_plus_rp_accept_025 = false;
	 bool proton_plus_rp_accept_120 = false;

	 std::map<int,std::vector<double> > proton_plus_pars;
	 std::map<int,std::vector<double> > proton_minus_pars;

	 if(proton_pz_plus > 0.){
	    xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
	    t_proton_plus = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
	    thx_proton_plus = atan(-proton_px_plus/proton_pi);
	    thy_proton_plus = atan(proton_py_plus/proton_pi);

	    //FIXME
	    double out_x, out_thx, out_y, out_thy, out_xi;
	    proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
	    proton_plus_pars[20] = std::vector<double>(5,0.);
	    proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
	    proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
	    proton_plus_pars[20][4] = out_xi;

	    proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21);
	    //proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
	    //proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
	    proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
	    proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
	    proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);
	 }

	 if(proton_pz_minus < 0.){ 
	    xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
	    t_proton_minus = -2*( (proton_pi*proton_energy_minus) - (proton_pi*proton_pz_minus) );
	    thx_proton_minus = atan(-proton_px_minus/proton_pi);
	    thy_proton_minus = atan(proton_py_minus/proton_pi);

	    double out_x, out_thx, out_y, out_thy, out_xi;
	    proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
	    proton_minus_pars[120] = std::vector<double>(5,0.);
	    proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
	    proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
	    proton_minus_pars[120][4] = out_xi;

	    proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121);
	    proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124);
	    proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125);
	    proton_minus_rp_accept_020 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 20);
	 }

	 histosTH1F["xi_proton_plus"]->Fill( xi_proton_plus , event_weight );
	 histosTH1F["t_proton_plus"]->Fill( fabs(t_proton_plus) , event_weight );
	 histosTH1F["xi_proton_minus"]->Fill( xi_proton_minus , event_weight );
	 histosTH1F["t_proton_minus"]->Fill( fabs(t_proton_minus) , event_weight );

	 histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
	 histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );
	 histosTH1F["thx_proton_minus"]->Fill( thx_proton_minus , event_weight );
	 histosTH1F["thy_proton_minus"]->Fill( thy_proton_minus , event_weight );

	 histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
	 histosTH2F["proton_minus_xi_vs_t"]->Fill( fabs(t_proton_minus) , xi_proton_minus , event_weight );

	 bool proton_plus_rp_accept = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 ) ||
	    ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 );
	 bool proton_minus_rp_accept = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 ) ||
	    ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 );

	 if(proton_plus_rp_accept){
	    histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight ); 
	    histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
	 }
	 if(proton_minus_rp_accept){
	    histosTH1F["xi_proton_minus_accepted"]->Fill( xi_proton_minus , event_weight );
	    histosTH1F["t_proton_minus_accepted"]->Fill( fabs(t_proton_minus) , event_weight ); 
	    histosTH2F["proton_minus_xi_vs_t_accepted"]->Fill( fabs(t_proton_minus) , xi_proton_minus , event_weight );
	 }

	 if(proton_plus_rp_accept_020){
	    histosTH1F["xi_proton_plus_accepted_020"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted_020"]->Fill( fabs(t_proton_plus) , event_weight ); 
	    histosTH1F["posx_proton_plus_accepted_020"]->Fill( proton_plus_pars[20][0] , event_weight );
	    histosTH1F["posy_proton_plus_accepted_020"]->Fill( proton_plus_pars[20][1] , event_weight );
	    histosTH2F["pos_y_vs_x_proton_plus_accepted_020"]->Fill( proton_plus_pars[20][0], proton_plus_pars[20][1] , event_weight );
	 }
	 if(proton_plus_rp_accept_021){
	    histosTH1F["xi_proton_plus_accepted_021"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted_021"]->Fill( fabs(t_proton_plus) , event_weight ); 
	 }
	 if(proton_plus_rp_accept_024){
	    histosTH1F["xi_proton_plus_accepted_024"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted_024"]->Fill( fabs(t_proton_plus) , event_weight ); 
	 }
	 if(proton_plus_rp_accept_025){
	    histosTH1F["xi_proton_plus_accepted_025"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted_025"]->Fill( fabs(t_proton_plus) , event_weight ); 
	 }
	 if(proton_plus_rp_accept_120){
	    histosTH1F["xi_proton_plus_accepted_120"]->Fill( xi_proton_plus , event_weight );
	    histosTH1F["t_proton_plus_accepted_120"]->Fill( fabs(t_proton_plus) , event_weight ); 
	 }

	 if(proton_minus_rp_accept_120){
	    histosTH1F["xi_proton_minus_accepted_120"]->Fill( xi_proton_minus , event_weight );
	    histosTH1F["t_proton_minus_accepted_120"]->Fill( fabs(t_proton_minus) , event_weight ); 
	    histosTH1F["posx_proton_minus_accepted_120"]->Fill( proton_minus_pars[120][0] , event_weight );
	    histosTH1F["posy_proton_minus_accepted_120"]->Fill( proton_minus_pars[120][1] , event_weight );
	    histosTH2F["pos_y_vs_x_proton_minus_accepted_120"]->Fill( proton_minus_pars[120][0], proton_minus_pars[120][1] , event_weight );
	 }
	 if(proton_minus_rp_accept_020){
	    histosTH1F["xi_proton_minus_accepted_020"]->Fill( xi_proton_minus , event_weight );
	    histosTH1F["t_proton_minus_accepted_020"]->Fill( fabs(t_proton_minus) , event_weight ); 
	 }
      }

      //-------------------------------------------------------------------------------------------------
      // Vertices
      int n_vertices_selected = 0;
      /*int idx_vtx_max_sumpt = -1;
      double sumpt_max = 0.;*/
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
         int idx_vtx = it_vtx - vertex_coll->begin();
         //if( it_vtx->SumPtTracks > sumpt_max ){ idx_vtx_max_sumpt = idx_vtx; sumpt_max = it_vtx->SumPtTracks; }

         if( it_vtx->fake ) continue;
         if( !it_vtx->validity ) continue;
         ++n_vertices_selected;

         histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
         histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
         histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );
         histosTH1F["vtx_ndof"]->Fill( it_vtx->ndof, event_weight );
         histosTH1F["vtx_chi2"]->Fill( it_vtx->chi2, event_weight );
      }
      //histosTH1F["vtx_sumpt_max"]->Fill( idx_vtx_max_sumpt, event_weight );
      histosTH1F["vertex_multiplicity"]->Fill( n_vertices_selected, event_weight );

      //MyVertex const& primaryVertex = vertex_coll->at(0);
      MyVertex& primaryVertex = vertex_coll->at(0);
      histosTH1F["prim_vtx_zpos"]->Fill( primaryVertex.z, event_weight );
      histosTH1F["prim_vtx_xpos"]->Fill( primaryVertex.x, event_weight );
      histosTH1F["prim_vtx_ypos"]->Fill( primaryVertex.y, event_weight );
      
      histosTH1F["prim_vtx_ndof"]->Fill( primaryVertex.ndof, event_weight );
      histosTH1F["prim_vtx_chi2"]->Fill( primaryVertex.chi2, event_weight );
      histosTH1F["prim_vtx_chi2n"]->Fill( primaryVertex.chi2n(), event_weight );
      histosTH1F["prim_vtx_ntracks"]->Fill( primaryVertex.ntracks, event_weight );
      histosTH1F["prim_vtx_sumpt"]->Fill( primaryVertex.SumPtTracks, event_weight );

      double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
      bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
                              primaryVertex.ndof > 4 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      if(selectVertex && !select_Vertex) continue;

      histosTH1F["prim_vtx_zpos_after_vtx_sel"]->Fill( primaryVertex.z, event_weight );
      histosTH1F["prim_vtx_xpos_after_vtx_sel"]->Fill( primaryVertex.x, event_weight );
      histosTH1F["prim_vtx_ypos_after_vtx_sel"]->Fill( primaryVertex.y, event_weight );

      histosTH1F["prim_vtx_ndof_after_vtx_sel"]->Fill( primaryVertex.ndof, event_weight );
      histosTH1F["prim_vtx_chi2_after_vtx_sel"]->Fill( primaryVertex.chi2, event_weight );
      histosTH1F["prim_vtx_chi2n_after_vtx_sel"]->Fill( primaryVertex.chi2n(), event_weight );
      histosTH1F["prim_vtx_ntracks_after_vtx_sel"]->Fill( primaryVertex.ntracks, event_weight );
      histosTH1F["prim_vtx_sumpt_after_vtx_sel"]->Fill( primaryVertex.SumPtTracks, event_weight );

      histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

      // Tracks
      int n_tracks_selected = 0;
      for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
         histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
         histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
         histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

         if( it_trk->Pt() < 0.5 ) continue;
         if( fabs( it_trk->Eta() ) > 2.5 ) continue;
         if( ( it_trk->dz / it_trk->edz ) > 5. ) continue;
         if( ( it_trk->d0 / it_trk->ed0 ) > 5. ) continue;
    
         /*
         outtrack.quality[0] = intrack.quality(TrackBase::qualityByName("loose"));
         outtrack.quality[1] = intrack.quality(TrackBase::qualityByName("tight"));
         outtrack.quality[2] = intrack.quality(TrackBase::qualityByName("highPurity"));
         outtrack.quality[3] = intrack.quality(TrackBase::qualityByName("confirmed"));
         outtrack.quality[4] = intrack.quality(TrackBase::qualityByName("goodIterative"));
         */ 
         if( !it_trk->quality[2] ) continue;

         ++n_tracks_selected;
      }
      //histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );

      // Jets
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin() ; it_jet != pfJet_coll->end() ; ++it_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;

         MyBaseJet const& basejet = it_jet->mapjet[ jetCorrName ];
         histosTH1F["jet_pt"]->Fill( basejet.Pt(), event_weight );
         if(basejet.Pt() > 0.) histosTH1F["jet_eta"]->Fill( basejet.Eta(), event_weight );
         histosTH1F["jet_phi"]->Fill( basejet.Phi(), event_weight );
      }

      bool select_leadingJet = false;
      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
	 histosTH1F["leadingJet_pt"]->Fill( leadingJet.Pt(), event_weight );
	 if(leadingJet.Pt() > 0.) histosTH1F["leadingJet_eta"]->Fill( leadingJet.Eta(), event_weight );
	 histosTH1F["leadingJet_phi"]->Fill( leadingJet.Phi(), event_weight );

         histosTH1F["leadingJet_jec"]->Fill( leadingJet.jec, event_weight );

         //bool looseJetId = ( pfJet_coll->at(0) ).LooseJetId;
         bool looseJetId = loosePFJetID(pfJet_coll->at(0),jetCollName);
         histosTH1F["leadingJet_looseJetId"]->Fill( looseJetId ? 1 : 0, event_weight );
         //bool tightJetId = ( pfJet_coll->at(0) ).TightJetId;
         bool tightJetId = tightPFJetID(pfJet_coll->at(0),jetCollName);
         histosTH1F["leadingJet_tightJetId"]->Fill( tightJetId ? 1 : 0, event_weight );

         select_leadingJet = (leadingJet.Pt() >= ptJetMin && fabs(leadingJet.Eta()) <= etaJetMax && looseJetId);
      }

      bool select_secondJet = false;
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[ jetCorrName ];
	 histosTH1F["secondJet_pt"]->Fill( secondJet.Pt(), event_weight );
	 if(secondJet.Pt() > 0.) histosTH1F["secondJet_eta"]->Fill( secondJet.Eta(), event_weight );
	 histosTH1F["secondJet_phi"]->Fill( secondJet.Phi(), event_weight );

         histosTH1F["secondJet_jec"]->Fill( secondJet.jec, event_weight );

         //bool looseJetId = ( pfJet_coll->at(1) ).LooseJetId;
         bool looseJetId = loosePFJetID(pfJet_coll->at(1),jetCollName);
         histosTH1F["secondJet_looseJetId"]->Fill( looseJetId ? 1 : 0, event_weight );
         //bool tightJetId = ( pfJet_coll->at(1) ).TightJetId;
         bool tightJetId = tightPFJetID(pfJet_coll->at(1),jetCollName);
         histosTH1F["secondJet_tightJetId"]->Fill( tightJetId ? 1 : 0, event_weight );

         select_secondJet = (secondJet.Pt() >= ptJetMin && fabs(secondJet.Eta()) <= etaJetMax && looseJetId);
      }
      bool select_Jet = ( select_leadingJet && select_secondJet);
      if(selectJets && !select_Jet) continue;

      histosTH1F["EventSelection"]->Fill( "Jet", event_weight );

      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
	 histosTH1F["leadingJet_pt_after_jet_sel"]->Fill( leadingJet.Pt(), event_weight );
	 histosTH1F["leadingJet_eta_after_jet_sel"]->Fill( leadingJet.Eta(), event_weight );
	 histosTH1F["leadingJet_phi_after_jet_sel"]->Fill( leadingJet.Phi(), event_weight );
	 histosTH1F["leadingJet_jec_after_jet_sel"]->Fill( leadingJet.jec, event_weight );
      }
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[ jetCorrName ];
	 histosTH1F["secondJet_pt_after_jet_sel"]->Fill( secondJet.Pt(), event_weight );
	 histosTH1F["secondJet_eta_after_jet_sel"]->Fill( secondJet.Eta(), event_weight );
	 histosTH1F["secondJet_phi_after_jet_sel"]->Fill( secondJet.Phi(), event_weight );
	 histosTH1F["secondJet_jec_after_jet_sel"]->Fill( secondJet.jec, event_weight );
      }

      // Particle-flow
      vector<MyPFCand> particles_sorted;
      double etaEdgeLow = -999.0;
      double etaEdgeHigh = 999.0;

      double pfEPlusPz = 0.;
      double pfEMinusPz = 0.;
      double xiCorrFactor = 1.0;
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin();
	                             it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();

         // Apply thresholds
         if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;

         if( eta >= etaEdgeLow && eta <= etaEdgeHigh ) particles_sorted.push_back(*it_pfcand);
         pfEPlusPz  += energy + pz;
         pfEMinusPz += energy - pz;
      }
      // Xi (CMS)
      double pfXiPlusReco = xiCorrFactor*pfEPlusPz/8000.;
      double pfXiMinusReco = xiCorrFactor*pfEMinusPz/8000.;

      // Find eta_min & eta_max
      double pfEtaMin = -999.;
      double pfEtaMax = 999.;
      if( particles_sorted.size() > 0 ){
	 std::stable_sort(particles_sorted.begin(), particles_sorted.end(), sortByEta);
	 pfEtaMin = particles_sorted.at(0).Eta();
	 pfEtaMax = particles_sorted.at(particles_sorted.size() - 1).Eta();
      }

      //if( selectEtaMax && (pfEtaMax > 3.0) ) continue;
      if( selectEtaMax && (pfEtaMax > etaMaxThreshold) ) continue;
      histosTH1F["EventSelection"]->Fill( "EtaMax", event_weight );

      //if( selectEtaMin && (pfEtaMin < -3.0) ) continue;
      if( selectEtaMin && (pfEtaMin < -etaMaxThreshold) ) continue;
      histosTH1F["EventSelection"]->Fill( "EtaMin", event_weight );

      // TOTEM T2
      int n_t2_tracks_selected = 0;
      int n_t2_tracks_selected_zplus = 0;
      int n_t2_tracks_selected_zminus = 0;
      if(!isMC){
	 vector<double> const& t2_trk_entryX = t2_event->TrkEntryX;
	 vector<double> const& t2_trk_entryY = t2_event->TrkEntryY;
	 vector<double> const& t2_trk_entryZ =  t2_event->TrkEntryZ;
	 vector<double> const& t2_trk_chiProb =  t2_event->TrkChiProb;

	 size_t n_t2_tracks = t2_trk_chiProb.size();
	 for(size_t i_t2_trk = 0; i_t2_trk < n_t2_tracks; ++i_t2_trk){
	    double trk_entryZ = t2_trk_entryZ[i_t2_trk];
	    int zside = ( trk_entryZ >= 0. ) ? 1 : -1;
	    if( zside > 0 )
	       histosTH1F["t2_track_chi2Prob_zplus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );
	    else
	       histosTH1F["t2_track_chi2Prob_zminus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );

	    // Select tracks
	    if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

	    ++n_t2_tracks_selected;
	    if( zside > 0 ) ++n_t2_tracks_selected_zplus;
	    else            ++n_t2_tracks_selected_zminus;

	    if( zside > 0 ){
	       histosTH1F["t2_track_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
	       histosTH1F["t2_track_entryY_zplus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
	       histosTH2F["t2_track_entryY_vs_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
	    } else{
	       histosTH1F["t2_track_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
	       histosTH1F["t2_track_entryY_zminus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
	       histosTH2F["t2_track_entryY_vs_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
	    }
	 }

	 if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;
	 histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Plus", event_weight );

	 if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;
	 histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Minus", event_weight );
      }

      bool proton_right_valid = false;
      bool proton_left_valid = false;
      if(!isMC){
	 proton_right_valid = rec_proton_right->valid;
	 proton_left_valid = rec_proton_left->valid;
	 if( selectSingleArmRecProton && (proton_right_valid && proton_left_valid) ) continue;
	 histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

	 if( selectDoubleArmRecProton && !(proton_right_valid && proton_left_valid) ) continue;
	 histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

	 bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
	 bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
	 if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
	 histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

	 if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
	 histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );
      }
      //-------------------
      // After selection 
      //-------------------
      histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );

      histosTH1F["pf_etaMax"]->Fill( pfEtaMax, event_weight ); 
      histosTH1F["pf_etaMin"]->Fill( pfEtaMin, event_weight );

      double pfDeltaEta = pfEtaMax - pfEtaMin;
      histosTH1F["pf_deltaEta"]->Fill( pfDeltaEta, event_weight ); 

      histosTH1F["pf_EPlusPz"]->Fill( pfEPlusPz, event_weight );
      histosTH1F["pf_EMinusPz"]->Fill( pfEMinusPz, event_weight );
      histosTH1F["pf_xiPlus"]->Fill( pfXiPlusReco, event_weight );
      histosTH1F["pf_xiMinus"]->Fill( pfXiMinusReco, event_weight );
      histosTH1F["pf_logXiPlus"]->Fill( log10(pfXiPlusReco), event_weight );
      histosTH1F["pf_logXiMinus"]->Fill( log10(pfXiMinusReco), event_weight );

      histosTH1F["t2_track_multiplicity_zplus"]->Fill( n_t2_tracks_selected_zplus, event_weight );
      histosTH1F["t2_track_multiplicity_zminus"]->Fill( n_t2_tracks_selected_zminus, event_weight );
      histosTH2F["t2_track_multiplicity_vs_track_multiplicity"]->Fill( n_tracks_selected, n_t2_tracks_selected, event_weight );
      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
         histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), n_t2_tracks_selected, event_weight );
      }

      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin();
	                             it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );

	 if(partType == MyPFCand::X)
	    histosTH2F["energyVsEtaUndefined"]->Fill( eta, energy, event_weight );
	 else if(partType == MyPFCand::h)
	    histosTH2F["energyVsEtaChargedHadron"]->Fill( eta, energy, event_weight ); 
	 else if(partType == MyPFCand::e) 
	    histosTH2F["energyVsEtaElectron"]->Fill( eta, energy, event_weight );
	 else if(partType == MyPFCand::mu) 
	    histosTH2F["energyVsEtaMuon"]->Fill( eta, energy, event_weight );
	 else if(partType == MyPFCand::gamma) 
	    histosTH2F["energyVsEtaGamma"]->Fill( eta, energy, event_weight );
	 else if(partType == MyPFCand::h0) 
	    histosTH2F["energyVsEtaNeutralHadron"]->Fill( eta, energy, event_weight );
	 else if(partType == MyPFCand::h_HF){ 
	    histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );
	    /*if( part->ecalEnergy() > 0. ) histosTH2F["energyVsEtaHadronHFEcalEnergy"]->Fill( eta, energy, event_weight );
	    else                          histosTH2F["energyVsEtaHadronHFNoEcalEnergy"]->Fill( eta, energy, event_weight );*/
	 } else if(partType == MyPFCand::egamma_HF) 
	    histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );
      }

      /*for(vector<MyFSCHit>::iterator it_hit = fscHits_coll->begin() ; it_hit != fscHits_coll->end() ; ++it_hit){
         histosTH1F["fscHit_energy"]->Fill( it_hit->energy, event_weight );
         histosTH1F["fscHit_time"]->Fill( it_hit->time, event_weight );
      }*/

      // RP protons
      if(!isMC){
         RPRootDumpTrackInfo* rp_track_info_020 = rp_track_info[20];
         RPRootDumpTrackInfo* rp_track_info_120 = rp_track_info[120];

         bool rp_track_valid_020 = rp_track_info_020->valid;
         if( rp_track_valid_020 ){
            double rp_track_posx_020 = rp_track_info_020->x;
            double rp_track_posy_020 = rp_track_info_020->y;
            histosTH1F["rp_track_posx_020"]->Fill( rp_track_posx_020, event_weight );
            histosTH1F["rp_track_posy_020"]->Fill( rp_track_posy_020, event_weight );
            histosTH2F["rp_track_pos_y_vs_x_020"]->Fill( rp_track_posx_020, rp_track_posy_020, event_weight );
         }

         bool rp_track_valid_120 = rp_track_info_120->valid;
         if( rp_track_valid_120 ){
            double rp_track_posx_120 = rp_track_info_120->x;
            double rp_track_posy_120 = rp_track_info_120->y;
            histosTH1F["rp_track_posx_120"]->Fill( rp_track_posx_120, event_weight );
            histosTH1F["rp_track_posy_120"]->Fill( rp_track_posy_120, event_weight );
            histosTH2F["rp_track_pos_y_vs_x_120"]->Fill( rp_track_posx_120, rp_track_posy_120, event_weight );
         }

	 //bool proton_right_valid = rec_proton_right->valid;
	 double chi2_proton_right = rec_proton_right->chi2;
	 double chindf_proton_right = rec_proton_right->chindf;
	 double xi_proton_right = rec_proton_right->xi;
	 //bool good_proton_right = proton_right_valid && (chi2_proton_right/chindf_proton_right > 1);
	 bool good_proton_right = proton_right_valid && (xi_proton_right < 0.);
	 if( good_proton_right ){
	    //double xi_proton_right = rec_proton_right->xi;
	    double t_proton_right = rec_proton_right->t;
	    histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight );
	    histosTH1F["proton_right_xi"]->Fill( -xi_proton_right, event_weight );
	    histosTH1F["proton_right_t"]->Fill( -t_proton_right, event_weight );

	    if(-xi_proton_right > 0.){
	       xi_proton_right = -xi_proton_right;
	       histosTH1F["proton_right_logXi"]->Fill( log10(xi_proton_right), event_weight );
	       histosTH1F["pf_xiMinus_minus_proton_right_xi"]->Fill( (pfXiMinusReco - xi_proton_right), event_weight );
	       histosTH2F["proton_right_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_right), event_weight );
	       histosTH2F["proton_right_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_right), event_weight );
	       histosTH2F["proton_right_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, xi_proton_right, event_weight );
	       histosTH2F["proton_right_xi_vs_pf_xiMinus"]->Fill( pfXiMinusReco, xi_proton_right, event_weight );
	       histosTH2F["proton_right_logXi_vs_t"]->Fill( -t_proton_right, log10(xi_proton_right), event_weight );
	    }

	    if( pfJet_coll->size() > 0 ){
	       MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
	       histosTH2F["proton_right_t_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), -t_proton_right, event_weight );
	    }
	 }

	 //bool proton_left_valid = rec_proton_left->valid;
	 double chi2_proton_left = rec_proton_left->chi2;
	 double chindf_proton_left = rec_proton_left->chindf;
	 double xi_proton_left = rec_proton_left->xi;
	 //bool good_proton_left = proton_left_valid && (chi2_proton_left/chindf_proton_left > 1);
	 bool good_proton_left = proton_left_valid && (xi_proton_left < 0.);
	 if( good_proton_left ){
	    //double xi_proton_left = rec_proton_left->xi;
	    double t_proton_left = rec_proton_left->t;
	    histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
	    histosTH1F["proton_left_xi"]->Fill( -xi_proton_left, event_weight );
	    histosTH1F["proton_left_t"]->Fill( -t_proton_left, event_weight );

	    if(-xi_proton_left > 0.){
	       xi_proton_left = -xi_proton_left;
	       histosTH1F["proton_left_logXi"]->Fill( log10(xi_proton_left), event_weight );
	       histosTH1F["pf_xiPlus_minus_proton_left_xi"]->Fill( (pfXiPlusReco - xi_proton_left), event_weight );
	       histosTH2F["proton_left_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_left), event_weight );
	       histosTH2F["proton_left_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_left), event_weight );
	       histosTH2F["proton_left_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, xi_proton_left, event_weight );
	       histosTH2F["proton_left_xi_vs_pf_xiMinus"]->Fill( pfXiMinusReco, xi_proton_left, event_weight );
	       histosTH2F["proton_left_logXi_vs_t"]->Fill( -t_proton_left, log10(xi_proton_left), event_weight );
	    }

	    if( pfJet_coll->size() > 0 ){
	       MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
	       histosTH2F["proton_left_t_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), -t_proton_left, event_weight );
	    }
	 }
	 bool proton_pair_valid = rec_proton_pair->valid;
	 double chi2_proton_pair = rec_proton_pair->chi2;
	 double chindf_proton_pair = rec_proton_pair->chindf;
	 double xi_proton_pair_right = rec_proton_pair->xir;
	 double xi_proton_pair_left = rec_proton_pair->xil;
	 //bool good_proton_pair = proton_pair_valid && (chi2_proton_pair/chindf_proton_pair > 2);
	 bool good_proton_pair = proton_pair_valid && (xi_proton_pair_right < 0.) && (xi_proton_pair_left < 0.);
	 if( good_proton_pair ){
	    histosTH1F["proton_pair_chi2"]->Fill( chi2_proton_pair, event_weight );

	    //double xi_proton_pair_right = rec_proton_pair->xir;
	    double t_proton_pair_right = rec_proton_pair->tr;
	    histosTH1F["proton_pair_right_xi"]->Fill( -xi_proton_pair_right, event_weight );
	    histosTH1F["proton_pair_right_t"]->Fill( -t_proton_pair_right, event_weight );
	    if(-xi_proton_pair_right > 0.)
	       histosTH1F["proton_pair_right_logXi"]->Fill( log10(-xi_proton_pair_right), event_weight );
	    /*if(xi_proton_pair_right > 0.)
	      histosTH1F["proton_pair_right_logXi"]->Fill( log10(xi_proton_pair_right), event_weight );*/

	    //double xi_proton_pair_left = rec_proton_pair->xil;
	    double t_proton_pair_left = rec_proton_pair->tl;
	    histosTH1F["proton_pair_left_xi"]->Fill( -xi_proton_pair_left, event_weight );
	    histosTH1F["proton_pair_left_t"]->Fill( -t_proton_pair_left, event_weight );
	    if(-xi_proton_pair_left > 0.)
	       histosTH1F["proton_pair_left_logXi"]->Fill( log10(-xi_proton_pair_left), event_weight );
	    /*if(xi_proton_pair_left > 0.)
	      histosTH1F["proton_pair_left_logXi"]->Fill( log10(xi_proton_pair_left), event_weight );*/
	 }
      }
    } // End of loop over events in a file
    
    // Close current file
    file->Close();
    
  } // End of loop over files
  
  
  // Output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();
  
  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
                                  it_histo != histosTH1F.end(); ++it_histo)
     (*it_histo).second->Write();
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

  output->Close();
}
