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
#include <TSystemFile.h>
#include <TSystemDirectory.h>

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

#include "EventMetaData.h"
#include "NtupleInfo.h"
#include "RPEvent.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPattern.h"
#include "RPRootDumpPatternInfo.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "T1Event.h"
#include "T2Event.h"
#include "TriggerData.h"

#include "analysis_tools.h"
#include "rp_aperture_config.h"
#include "deltaPhi.h"
//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>

using namespace std;


void ana_psel_DataMC_Jpsi(vector<string> const& fileNames, string const& outputFileName = "output.root",const Double_t t_proton_down_=0.03, const Double_t t_proton_up_=1.0,const Int_t Bin_mass=100 ,const Int_t nevt_max = -1){

	bool isMC  = true;
	bool verbose = false;
	//string treeName ="evt";
	//string treeName ="cms_totem";
	string treeName = (!isMC) ? "cms_totem" : "evt";
	double etaMaxThreshold = 2.0;

	bool selectBunchCrossing = false;
	bool selectVertex = true;
	bool selectMuons = true;
	bool selectTrack = true;
	bool selectEtaMax = false;
	bool selectEtaMin = false;
	bool selectZeroHitsT2Plus = false;
	bool selectZeroHitsT2Minus = false;
	bool selectSingleArmRecProton = true;
	//bool selectSingleArmRecProton = false;
	bool selectDoubleArmRecProton = false;
	//bool selectDoubleArmRecProton = true;
	bool selectElastic = false;
	bool selectNonElastic = false;
	//bool selectNonElastic = true;
	bool selectRPPlusAccept = false;
	bool selectRPMinusAccept = true;

	// Declaration of histograms
	map<string,TH1F*> histosTH1F;

	vector<int> bunchCrossingList;
	bunchCrossingList.push_back(27);
	bunchCrossingList.push_back(649);
	bunchCrossingList.push_back(2991);


	vector<string> selectHLTPathNames;
	selectHLTPathNames.push_back("HLT_L1_DoubleMu0");


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
	//if(verbose)cout<<"pass MC 1"<<endl;
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



	vector<string> selections;
	selections.push_back("All");
	selections.push_back("BunchCrossing");
	selections.push_back("HLT");
	selections.push_back("Vertex");
	selections.push_back("Muons");
	selections.push_back("EtaMax");
	selections.push_back("EtaMin");
	selections.push_back("ZeroHitsT2Plus");
	selections.push_back("ZeroHitsT2Minus");
	selections.push_back("SingleArmRP");
	selections.push_back("DoubleArmRP");
	selections.push_back("Elastic");
	selections.push_back("NonElastic");
	selections.push_back("RPProton");
	selections.push_back("MCRPPlusAccept");
	selections.push_back("MCRPMinusAccept");
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

	histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
	histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
	histosTH1F["vtx_ndof"] = new TH1F("vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["vtx_chi2"] = new TH1F("vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);

	histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertices" , 30 , 0 , 30);
	histosTH1F["vertex_multiplicity_after_vtx_sel"] = new TH1F("vertex_multiplicity_after_vtx_sel", "n vertices after vtx sel" , 30 , 0 , 30);
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
	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 100 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_rapidity"] = new TH1F("track_rapidity", "#rapidity(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 100 , -M_PI , M_PI);
	histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);

	histosTH1F["muon_pt"] = new TH1F("muon_pt", "p_{T}(muon)" , 100 , 0. , 100.);
	histosTH1F["muon_eta"] = new TH1F("muon_eta", "#eta(muon)" , 100 , -5.2 , 5.2);
	histosTH1F["muon_phi"] = new TH1F("muon_phi", "#phi(muon)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muon_multiplicity"] = new TH1F("muon_multiplicity", "n muons" , 100 , 0 , 100);

	//Muon1 info

	histosTH1F["muon1_pt"] = new TH1F("muon1_pt", "p_{T}(mu1)" , 100 , 0. , 100.);
	histosTH1F["muon1_eta"] = new TH1F("muon1_eta", "#eta(mu1)" , 100 , -5.2 , 5.2);
	histosTH1F["muon1_phi"] = new TH1F("muon1_phi", "#phi(muon1)" , 100 , -1.2*M_PI ,1.2*M_PI);
	histosTH1F["muon1_rapidity"] = new TH1F("muon1_rapidity", "y(mu1)" , 100 , -15. , 15.);
	//Muon1 info

	histosTH1F["muon2_pt"] = new TH1F("muon2_pt", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["muon2_eta"] = new TH1F("muon2_eta", "#eta(mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["muon2_phi"] = new TH1F("muon2_phi", "#phi(muon2)" , 100 , -1.2*M_PI ,1.2*M_PI);
	histosTH1F["muon2_rapidity"] = new TH1F("muon2_rapidity", "y(mu2)" , 100 , -15. , 15.);

	//Deltas info
	histosTH1F["muonDeltaPt"] = new TH1F("muonDeltaPt", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta"] = new TH1F("muonDeltaEta", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi"] = new TH1F("muonDeltaPhi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY"] = new TH1F("muonDeltaY", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//Dimuon
	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2"] = new TH1F("dimuon_pt2", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt"] = new TH1F("dimuon_pt", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta"] = new TH1F("dimuon_eta", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity"] = new TH1F("dimuon_rapidity", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity"] = new TH1F("dimuon_multiplicity", "n dimuons" , 100 , 0 , 100); 
	/*histosTH1F["proton_right_xi_t_cut"] = new TH1F("proton_right_xi_t_selected", "#xi" , 29 , xi_bins);
	  histosTH1F["proton_right_t_t_cut"] = new TH1F("proton_right_t_t_selected", "|t|" , 11 , tbins);*/
	//dimuon_proton_right
	//histosTH1F["dimuon_proton_right_xi"] = new TH1F("dimuon_proton_right_xi", "#xi" , 29 , xi_bins);
	//histosTH1F["dimuon_proton_right_t"] = new TH1F("dimuon_proton_right_t", "-t" , 11 , tbins);     
	histosTH1F["dimuon_proton_right_xi"] = new TH1F("dimuon_proton_right_xi", "#xi" , 100,0. ,1. );
	histosTH1F["dimuon_proton_right_t"] = new TH1F("dimuon_proton_right_t", "-t" , 100 , 0.,5.);
	histosTH1F["muonDeltaPt_proton_right"] = new TH1F("muonDeltaPt_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_proton_right"] = new TH1F("muonDeltaEta_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_proton_right"] = new TH1F("muonDeltaPhi_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_proton_right"] = new TH1F("muonDeltaY_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	histosTH1F["dimuon_mass_proton_right"] = new TH1F("dimuon_mass_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2_proton_right"] = new TH1F("dimuon_pt2_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt_proton_right"] = new TH1F("dimuon_pt_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta_proton_right"] = new TH1F("dimuon_eta_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_proton_right"] = new TH1F("dimuon_rapidity_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_proton_right"] = new TH1F("dimuon_multiplicity_proton_right", "n dimuons" , 100 , 0 , 100);
	//_t_cut_proton_right
	Float_t tbins_matrix[34] = { 0.0144, 0.0256, 0.04, 0.0576, 0.0784, 0.1024, 0.1296, 0.16, 0.1936, 0.2304, 0.2704, 0.3136, 0.36, 0.4096, 0.4624, 0.5184, 0.5778, 0.64, 0.7056, 0.7744, 0.8464, 0.9216, 1., 1.0816, 1.1664, 1.2544, 1.3456, 1.44, 1.5376, 1.6384, 1.7424, 1.8496, 1.96, 2.0736};
	histosTH1F["dimuon_mass_t_cut_proton_right"] = new TH1F("dimuon_mass_t_cut_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt2_t_cut_proton_right"] = new TH1F("dimuon_pt2_t_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_pt_t_cut_proton_right"] = new TH1F("dimuon_pt_t_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_eta_t_cut_proton_right"] = new TH1F("dimuon_eta_t_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_right"] = new TH1F("dimuon_rapidity_t_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_right"] = new TH1F("dimuon_multiplicity_t_cut_proton_right", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_right_xi_t_cut"] = new TH1F("proton_right_xi_t_selected", "#xi" , 100, 0. ,5.0);
	histosTH1F["proton_right_t_t_cut"] = new TH1F("proton_right_t_t_selected", "-t" , 100, 0., 5.0);
	//Deltas info
	histosTH1F["muonDeltaPt_t_cut_proton_right"] = new TH1F("muonDeltaPt_t_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_t_cut_proton_right"] = new TH1F("muonDeltaEta_t_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_t_cut_proton_right"] = new TH1F("muonDeltaPhi_t_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 ,- 2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_t_cut_proton_right"] = new TH1F("muonDeltaY_t_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);



	//dimuon_proton_left
	histosTH1F["dimuon_proton_left_xi"] = new TH1F("dimuon_proton_left_xi", "#xi" , 100, 0.,1.);
	histosTH1F["dimuon_proton_left_t"] = new TH1F("dimuon_proton_left_t", "-t" , 100, 0., 5.);
	histosTH1F["dimuon_mass_proton_left"] = new TH1F("dimuon_mass_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_proton_left"] = new TH1F("dimuon_pt_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_proton_left"] = new TH1F("dimuon_pt2_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_proton_left"] = new TH1F("dimuon_eta_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_proton_left"] = new TH1F("dimuon_rapidity_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_proton_left"] = new TH1F("dimuon_multiplicity_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_proton_left"] = new TH1F("muonDeltaPt_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_proton_left"] = new TH1F("muonDeltaEta_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_proton_left"] = new TH1F("muonDeltaPhi_proton_left", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_proton_left"] = new TH1F("muonDeltaY_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//_t_cut_proton_left
	histosTH1F["dimuon_mass_t_cut_proton_left"] = new TH1F("dimuon_mass_t_cut_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_t_cut_proton_left"] = new TH1F("dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_t_cut_proton_left"] = new TH1F("dimuon_pt2_t_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_t_cut_proton_left"] = new TH1F("dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_left"] = new TH1F("dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_left"] = new TH1F("dimuon_multiplicity_t_cut_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_left_xi_t_cut"] = new TH1F("proton_left_xi_t_selected", "#xi" ,  100, 0.,1.);
	histosTH1F["proton_left_t_t_cut"] = new TH1F("proton_left_t_t_selected", "-t" , 100 , 0., 5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_t_cut_proton_left"] = new TH1F("muonDeltaPt_t_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_t_cut_proton_left"] = new TH1F("muonDeltaEta_t_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_t_cut_proton_left"] = new TH1F("muonDeltaPhi_t_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,- 2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_t_cut_proton_left"] = new TH1F("muonDeltaY_t_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);



	//jpsi mass selection
	histosTH1F["jpsi_mass"] = new TH1F("jpsi_mass", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt"] = new TH1F("jpsi_pt", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2"] = new TH1F("jpsi_pt2", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta"] = new TH1F("jpsi_eta", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity"] = new TH1F("jpsi_rapidity", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity"] = new TH1F("jpsi_multiplicity", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi"] = new TH1F("muonDeltaPt_jpsi", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi"] = new TH1F("muonDeltaEta_jpsi", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi"] = new TH1F("muonDeltaPhi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi"] = new TH1F("muonDeltaY_jpsi", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//Dimuon


	//proton_right
	histosTH1F["jpsi_proton_right_xi"] = new TH1F("jpsi_proton_right_xi", "#xi" , 100, 0., 1.);
	histosTH1F["jpsi_proton_right_t"] = new TH1F("jpsi_proton_right_t", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_mass_proton_right"] = new TH1F("jpsi_mass_proton_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_proton_right"] = new TH1F("jpsi_pt_proton_right", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_proton_right"] = new TH1F("jpsi_pt2_proton_right", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_proton_right"] = new TH1F("jpsi_eta_proton_right", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_proton_right"] = new TH1F("jpsi_rapidity_proton_right", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity_proton_right"] = new TH1F("jpsi_multiplicity_proton_right", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_proton_right"] = new TH1F("muonDeltaPt_jpsi_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_proton_right"] = new TH1F("muonDeltaEta_jpsi_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_proton_right"] = new TH1F("muonDeltaPhi_jpsi_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_proton_right"] = new TH1F("muonDeltaY_jpsi_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//proton_left
	histosTH1F["jpsi_proton_left_xi"] = new TH1F("jpsi_proton_left_xi", "#xi" , 100, 0.,1. );
	histosTH1F["jpsi_proton_left_t"] = new TH1F("jpsi_proton_left_t", "-t" , 100, 0.,5.);
	histosTH1F["jpsi_mass_proton_left"] = new TH1F("jpsi_mass_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_proton_left"] = new TH1F("jpsi_pt_proton_left", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_proton_left"] = new TH1F("jpsi_pt2_proton_left", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_proton_left"] = new TH1F("jpsi_eta_proton_left", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_proton_left"] = new TH1F("jpsi_rapidity_proton_left", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_multiplicity_proton_left"] = new TH1F("jpsi_multiplicity_proton_left", "n jpsi(dimuons)" , 100 , 0 , 100);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_proton_left"] = new TH1F("muonDeltaPt_jpsi_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_proton_left"] = new TH1F("muonDeltaEta_jpsi_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_proton_left"] = new TH1F("muonDeltaPhi_jpsi_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_proton_left"] = new TH1F("muonDeltaY_jpsi_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);        

	//_t_cut_proton_right
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_right"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_right"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_proton_right"] = new TH1F("jpsi_dimuon_pt2_t_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 2000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_right"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_right"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_proton_right", "n dimuons" , 100 , 0 , 20);

	histosTH1F["muonDeltaPt_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaPt_jpsi_t_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaEta_jpsi_t_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaPhi_jpsi_t_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_proton_right"] = new TH1F("muonDeltaY_jpsi_t_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	histosTH1F["jpsi_proton_right_xi_t_cut"] = new TH1F("jpsi_proton_right_xi_t_selected", "#xi" , 100, 0.,1.);
	histosTH1F["jpsi_proton_right_t_t_cut"] = new TH1F("jpsi_proton_right_t_t_selected", "-t" , 100, 0.,5.);



	//histosTH1F["xi_all_before"] = new TH1F("xi_all_before", "#xi" ,100, -1.,1. ); 

	//_t_cut_proton_left
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_left"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt2_t_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 2000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_left"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_proton_left", "n jpsi_dimuons" , 100 , 0 , 20);

	histosTH1F["jpsi_proton_left_xi_t_cut"] = new TH1F("jpsi_proton_left_xi_t_selected", "#xi" ,100, 0.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut"] = new TH1F("jpsi_proton_left_t_t_selected", "-t" , 100, 0., 5.);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaPt_jpsi_t_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaEta_jpsi_t_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaPhi_jpsi_t_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI ,2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_proton_left"] = new TH1F("muonDeltaY_jpsi_t_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//t_cut_xi_cut_proton_left
	histosTH1F["jpsi_dimuon_mass_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_mass_t_cut_xi_cut_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_pt_t_cut_xi_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_dimuon_pt2_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_pt2_t_cut_xi_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 2000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_eta_t_cut_xi_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_rapidity_t_cut_xi_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_left"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_left", "n jpsi_dimuons" , 100 , 0 , 20);
	histosTH1F["jpsi_proton_left_xi_t_cut_xi_cut"] = new TH1F("jpsi_proton_left_xi_selected_t_and_xi ", "#xi" , 30, -1.,1. );
	histosTH1F["jpsi_proton_left_t_t_cut_xi_cut"] = new TH1F("jpsi_proton_left_t_selected_t_and_xi", "-t" ,  30, 0., 5.);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPt_jpsi_t_cut_xi_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaEta_jpsi_t_cut_xi_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPhi_jpsi_t_cut_xi_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 ,-2.2*M_PI ,2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaY_jpsi_t_cut_xi_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//
	histosTH1F["pf_etaMax"] = new TH1F("pf_etaMax","#eta^{max}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_etaMin"] = new TH1F("pf_etaMin","#eta^{min}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_deltaEta"] = new TH1F("pf_deltaEta","#Delta#eta",100,0.,10.);
	histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",24,binningEPlusPz);
	histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",24,binningEPlusPz);
	histosTH1F["pf_xiPlus"] = new TH1F("pf_xiPlus","#xi^{+}",100,-1.,1.);
	histosTH1F["pf_xiMinus"] = new TH1F("pf_xiMinus","#xi^{-}",100,-1.,1.);
	histosTH1F["pf_logXiPlus"] = new TH1F("pf_logXiPlus","log(#xi^{+})",100,-5.,0.);
	histosTH1F["pf_logXiMinus"] = new TH1F("pf_logXiMinus","log(#xi^{-})",100,-5.,0.);

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
	histosTH1F["rp_track_posx_021"] = new TH1F("rp_track_posx_021", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_021"] = new TH1F("rp_track_posy_021", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_022"] = new TH1F("rp_track_posx_022", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_022"] = new TH1F("rp_track_posy_022", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_023"] = new TH1F("rp_track_posx_023", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_023"] = new TH1F("rp_track_posy_023", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_024"] = new TH1F("rp_track_posx_024", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_024"] = new TH1F("rp_track_posy_024", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_025"] = new TH1F("rp_track_posx_025", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_025"] = new TH1F("rp_track_posy_025", "y(RP track)" , 500, -50., 50.);

	histosTH1F["rp_track_posx_120"] = new TH1F("rp_track_posx_120", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_120"] = new TH1F("rp_track_posy_120", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_121"] = new TH1F("rp_track_posx_121", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_121"] = new TH1F("rp_track_posy_121", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_122"] = new TH1F("rp_track_posx_122", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_122"] = new TH1F("rp_track_posy_122", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_123"] = new TH1F("rp_track_posx_123", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_123"] = new TH1F("rp_track_posy_123", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_124"] = new TH1F("rp_track_posx_124", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_124"] = new TH1F("rp_track_posy_124", "y(RP track)" , 500, -50., 50.);
	histosTH1F["rp_track_posx_125"] = new TH1F("rp_track_posx_125", "x(RP track)" , 200, -10., 10.);
	histosTH1F["rp_track_posy_125"] = new TH1F("rp_track_posy_125", "y(RP track)" , 500, -50., 50.);

	histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 200 , 0. , 5.);
	histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 100 , 0. , 100.);
	histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 200 , -1. , 1.);
	histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 200 , 0. , 5.);
	histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 100 , 0. , 100.);

	/*histosTH1F["proton_pair_right_xi"] = new TH1F("proton_pair_right_xi", "#xi" , 200 , -1. , 1.);
	  histosTH1F["proton_pair_right_logXi"] = new TH1F("proton_pair_right_logXi","log(#xi)",200,-5.,0.);
	  histosTH1F["proton_pair_right_t"] = new TH1F("proton_pair_right_t", "-t" , 200 , 0. , 5.);
	  histosTH1F["proton_pair_left_xi"] = new TH1F("proton_pair_left_xi", "#xi" , 200 , -1. , 1.);
	  histosTH1F["proton_pair_left_logXi"] = new TH1F("proton_pair_left_logXi","log(#xi)",200,-5.,0.);
	  histosTH1F["proton_pair_left_t"] = new TH1F("proton_pair_left_t", "-t" , 200 , 0. , 5.);
	  histosTH1F["proton_pair_chi2"] = new TH1F("proton_pair_chi2", "#chi^{2}" , 100 , 0. , 100.);

	  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);*/


	histosTH1F["proton_pair_right_xi"] = new TH1F("proton_pair_right_xi", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_pair_right_logXi"] = new TH1F("proton_pair_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_right_t"] = new TH1F("proton_pair_right_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_left_xi"] = new TH1F("proton_pair_left_xi", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_pair_left_logXi"] = new TH1F("proton_pair_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_left_t"] = new TH1F("proton_pair_left_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_chi2"] = new TH1F("proton_pair_chi2", "#chi^{2}" , 100 , 0. , 100.);

	histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 100 , -1. , 1.);
	histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 100 , -1. , 1.);

	histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);

	histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" , 500, 0., 5.);
	histosTH1F["thx_proton_minus"] = new TH1F("thx_proton_minus", "thx_proton_minus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_minus"] = new TH1F("thy_proton_minus", "thy_proton_minus" , 200, -5e-4, 5e-4);

	histosTH1F["xi_proton_t_range_plus"] = new TH1F("xi_proton_t_range_plus", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_plus"] = new TH1F("t_proton_xi_range_plus", "t_proton_plus" , 500, 0., 5.);

	histosTH1F["xi_proton_t_range_minus"] = new TH1F("xi_proton_t_range_minus", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_minus"] = new TH1F("t_proton_xi_range_minus", "t_proton_minus" , 500, 0., 5.);

	//FIXME
	histosTH1F["xi_proton_plus_accepted"] = new TH1F("xi_proton_plus_accepted", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted"] = new TH1F("t_proton_plus_accepted", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_minus_accepted"] = new TH1F("xi_proton_minus_accepted", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted"] = new TH1F("t_proton_minus_accepted", "t_proton_minus" , 500, 0., 5.);

	histosTH1F["xi_proton_t_range_plus_accepted"] = new TH1F("xi_proton_t_range_plus_accepted", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_plus_accepted"] = new TH1F("t_proton_xi_range_plus_accepted", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_t_range_minus_accepted"] = new TH1F("xi_proton_t_range_minus_accepted", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_minus_accepted"] = new TH1F("t_proton_xi_range_minus_accepted", "t_proton_minus" , 500, 0., 5.);

	histosTH1F["xi_proton_plus_selected"] = new TH1F("xi_proton_plus_selected", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_selected"] = new TH1F("t_proton_plus_selected", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_minus_selected"] = new TH1F("xi_proton_minus_selected", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_selected"] = new TH1F("t_proton_minus_selected", "t_proton_minus" , 500, 0., 5.);

	histosTH1F["xi_proton_t_range_plus_selected"] = new TH1F("xi_proton_t_range_plus_selected", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_plus_selected"] = new TH1F("t_proton_xi_range_plus_selected", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_t_range_minus_selected"] = new TH1F("xi_proton_t_range_minus_selected", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_xi_range_minus_selected"] = new TH1F("t_proton_xi_range_minus_selected", "t_proton_minus" , 500, 0., 5.);

	// RP stations
	histosTH1F["xi_proton_plus_accepted_020"] = new TH1F("xi_proton_plus_accepted_020", "xi_proton_plus" , 200, 0., 1.);

	histosTH1F["t_proton_plus_accepted_020"] = new TH1F("t_proton_plus_accepted_020", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["posx_proton_plus_accepted_020"] = new TH1F("posx_proton_plus_accepted_020", "posx_proton_plus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_plus_accepted_020"] = new TH1F("posy_proton_plus_accepted_020", "posy_proton_plus" , 200, -0.05, 0.05);

	histosTH1F["xi_proton_plus_accepted_021"] = new TH1F("xi_proton_plus_accepted_021", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_021"] = new TH1F("t_proton_plus_accepted_021", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_022"] = new TH1F("xi_proton_plus_accepted_022", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_022"] = new TH1F("t_proton_plus_accepted_022", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_023"] = new TH1F("xi_proton_plus_accepted_023", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_023"] = new TH1F("t_proton_plus_accepted_023", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_024"] = new TH1F("xi_proton_plus_accepted_024", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_024"] = new TH1F("t_proton_plus_accepted_024", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_025"] = new TH1F("xi_proton_plus_accepted_025", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_025"] = new TH1F("t_proton_plus_accepted_025", "t_proton_plus" , 500, 0., 5.);
	histosTH1F["xi_proton_plus_accepted_120"] = new TH1F("xi_proton_plus_accepted_120", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_120"] = new TH1F("t_proton_plus_accepted_120", "t_proton_plus" , 500, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_120"] = new TH1F("xi_proton_minus_accepted_120", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_120"] = new TH1F("t_proton_minus_accepted_120", "t_proton_minus" , 500, 0., 5.);
	histosTH1F["posx_proton_minus_accepted_120"] = new TH1F("posx_proton_minus_accepted_120", "posx_proton_minus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_minus_accepted_120"] = new TH1F("posy_proton_minus_accepted_120", "posy_proton_minus" , 200, -0.05, 0.05);

	histosTH1F["xi_proton_minus_accepted_020"] = new TH1F("xi_proton_minus_accepted_020", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_020"] = new TH1F("t_proton_minus_accepted_020", "t_proton_minus" , 500, 0., 5.);


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
	histosTH2F["rp_track_pos_y_vs_x_021"] = new TH2F("rp_track_pos_y_vs_x_021", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_022"] = new TH2F("rp_track_pos_y_vs_x_022", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_023"] = new TH2F("rp_track_pos_y_vs_x_023", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_024"] = new TH2F("rp_track_pos_y_vs_x_024", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_025"] = new TH2F("rp_track_pos_y_vs_x_025", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_121"] = new TH2F("rp_track_pos_y_vs_x_121", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_122"] = new TH2F("rp_track_pos_y_vs_x_122", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_123"] = new TH2F("rp_track_pos_y_vs_x_123", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_124"] = new TH2F("rp_track_pos_y_vs_x_124", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_125"] = new TH2F("rp_track_pos_y_vs_x_125", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

	histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	histosTH2F["proton_plus_xi_vs_t_selected"] = new TH2F("proton_plus_xi_vs_t_selected","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t_selected"] = new TH2F("proton_minus_xi_vs_t_selected","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	histosTH2F["pos_y_vs_x_proton_plus_020"] = new TH2F("pos_y_vs_x_proton_plus_020", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_021"] = new TH2F("pos_y_vs_x_proton_plus_021", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_022"] = new TH2F("pos_y_vs_x_proton_plus_022", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_023"] = new TH2F("pos_y_vs_x_proton_plus_023", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_024"] = new TH2F("pos_y_vs_x_proton_plus_024", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_025"] = new TH2F("pos_y_vs_x_proton_plus_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_120"] = new TH2F("pos_y_vs_x_proton_minus_120", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_121"] = new TH2F("pos_y_vs_x_proton_minus_121", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_122"] = new TH2F("pos_y_vs_x_proton_minus_122", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_123"] = new TH2F("pos_y_vs_x_proton_minus_123", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_124"] = new TH2F("pos_y_vs_x_proton_minus_124", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_125"] = new TH2F("pos_y_vs_x_proton_minus_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);

	histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_021"] = new TH2F("pos_y_vs_x_proton_plus_accepted_021", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_022"] = new TH2F("pos_y_vs_x_proton_plus_accepted_022", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_023"] = new TH2F("pos_y_vs_x_proton_plus_accepted_023", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_024"] = new TH2F("pos_y_vs_x_proton_plus_accepted_024", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_plus_accepted_025"] = new TH2F("pos_y_vs_x_proton_plus_accepted_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_121"] = new TH2F("pos_y_vs_x_proton_minus_accepted_121", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_122"] = new TH2F("pos_y_vs_x_proton_minus_accepted_122", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_123"] = new TH2F("pos_y_vs_x_proton_minus_accepted_123", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_124"] = new TH2F("pos_y_vs_x_proton_minus_accepted_124", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_125"] = new TH2F("pos_y_vs_x_proton_minus_accepted_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
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
	//===========MC
	// Reco Muons rp_accept
	histosTH1F["dimuon_mass_rpplus_accept"] = new TH1F("dimuon_mass_rpplus_accept", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_rpplus_accept"] = new TH1F("dimuon_pt_rpplus_accept", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_rpplus_accept"] = new TH1F("dimuon_pt2_rpplus_accept", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_rpplus_accept"] = new TH1F("dimuon_eta_rpplus_accept", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_rpplus_accept"] = new TH1F("dimuon_rapidity_rpplus_accept", "y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["dimuon_multiplicity_rp_accept"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_rpplus_accept"] = new TH1F("muonDeltaPt_rpplus_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_rpplus_accept"] = new TH1F("muonDeltaEta_rpplus_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi_rpplus_accept"] = new TH1F("muonDeltaPhi_rpplus_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDphi_rpplus_accept"] = new TH1F("muonDphi_rpplus_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDeltaY_rpplus_accept"] = new TH1F("muonDeltaY_rpplus_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);


	// Reco Muons rp_accept && _tsel
	histosTH1F["dimuon_mass_rpplus_accept_tsel"] = new TH1F("dimuon_mass_rpplus_accept_tsel", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_rpplus_accept_tsel"] = new TH1F("dimuon_pt_rpplus_accept_tsel", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_rpplus_accept_tsel"] = new TH1F("dimuon_pt2_rpplus_accept_tsel", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_rpplus_accept_tsel"] = new TH1F("dimuon_eta_rpplus_accept_tsel", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_rpplus_accept_tsel"] = new TH1F("dimuon_rapidity_rpplus_accept_tsel", "y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["dimuon_multiplicity_rp_accept_tsel"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_rpplus_accept_tsel"] = new TH1F("muonDeltaPt_rpplus_accept_tsel", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_rpplus_accept_tsel"] = new TH1F("muonDeltaEta_rpplus_accept_tsel", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi_rpplus_accept_tsel"] = new TH1F("muonDeltaPhi_rpplus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDphi_rpplus_accept_tsel"] = new TH1F("muonDphi_rpplus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDeltaY_rpplus_accept_tsel"] = new TH1F("muonDeltaY_rpplus_accept_tsel", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["dimuon_xi_plus_rpplus_accept_tsel"]= new TH1F("dimuon_xi_plus_rpplus_accept_tsel", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["dimuon_t_plus_rpplus_accept_tsel"] = new TH1F("dimuon_t_plus_rpplus_accept_tsel", "|t|^{+}" , 100,0.,5.0);

	//jpsi mass selection
	histosTH1F["jpsi_mass_rpplus_accept"] = new TH1F("jpsi_mass_rpplus_accept", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_rpplus_accept"] = new TH1F("jpsi_pt_rpplus_accept", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_rpplus_accept"] = new TH1F("jpsi_pt2_rpplus_accept", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_rpplus_accept"] = new TH1F("jpsi_eta_rpplus_accept", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_rpplus_accept"] = new TH1F("jpsi_rapidity_rpplus_accept", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["jpsi_multiplicity_rpplus_accept"] = new TH1F("jpsi_multiplicity_rp_accept", "n jpsi(dimuons)" , 100 , 0 , 100);
	histosTH1F["jpsi_xi_plus_rpplus_accept"]= new TH1F("jpsi_xi_plus_rpplus_accept", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["jpsi_t_plus_rpplus_accept"] = new TH1F("jpsi_t_plus_rpplus_accept", "|t|^{+}" , 100,0.,5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_rpplus_accept"] = new TH1F("muonDeltaPt_jpsi_rpplus_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_rpplus_accept"] = new TH1F("muonDeltaEta_jpsi_rpplus_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_rpplus_accept"] = new TH1F("muonDeltaPhi_jpsi_rpplus_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDphi_jpsi_rpplus_accept"] = new TH1F("muonDphi_jpsi_rpplus_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_rpplus_accept"] = new TH1F("muonDeltaY_jpsi_rpplus_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//jpsi rp accept && _tsel
	//jpsi mass selection
	histosTH1F["jpsi_mass_rpplus_accept_tsel"] = new TH1F("jpsi_mass_rpplus_accept_tsel", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_rpplus_accept_tsel"] = new TH1F("jpsi_pt_rpplus_accept_tsel", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_rpplus_accept_tsel"] = new TH1F("jpsi_pt2_rpplus_accept_tsel", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_rpplus_accept_tsel"] = new TH1F("jpsi_eta_rpplus_accept_tsel", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_rpplus_accept_tsel"] = new TH1F("jpsi_rapidity_rpplus_accept_tsel", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["jpsi_multiplicity_rp_accept_tsel"] = new TH1F("jpsi_multiplicity_rp_accept_tsel", "n jpsi(dimuons)" , 100 , 0 , 100);
	histosTH1F["jpsi_xi_plus_rpplus_accept_tsel"]= new TH1F("jpsi_xi_plus_rpplus_accept_tsel", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["jpsi_t_plus_rpplus_accept_tsel"] = new TH1F("jpsi_t_plus_rpplus_accept_tsel", "|t|^{+}" , 100,0.,5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_rpplus_accept_tsel"] = new TH1F("muonDeltaPt_jpsi_rpplus_accept_tsel", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_rpplus_accept_tsel"] = new TH1F("muonDeltaEta_jpsi_rpplus_accept_tsel", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_rpplus_accept_tsel"] = new TH1F("muonDeltaPhi_jpsi_rpplus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDphi_jpsi_rpplus_accept_tsel"] = new TH1F("muonDphi_jpsi_rpplus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_rpplus_accept_tsel"] = new TH1F("muonDeltaY_jpsi_rpplus_accept_tsel", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//Minus
	// Reco Muons rp_accept
	histosTH1F["dimuon_mass_rpminus_accept"] = new TH1F("dimuon_mass_rpminus_accept", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_rpminus_accept"] = new TH1F("dimuon_pt_rpminus_accept", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_rpminus_accept"] = new TH1F("dimuon_pt2_rpminus_accept", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_rpminus_accept"] = new TH1F("dimuon_eta_rpminus_accept", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_rpminus_accept"] = new TH1F("dimuon_rapidity_rpminus_accept", "y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["dimuon_multiplicity_rp_accept"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_rpminus_accept"] = new TH1F("muonDeltaPt_rpminus_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_rpminus_accept"] = new TH1F("muonDeltaEta_rpminus_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi_rpminus_accept"] = new TH1F("muonDeltaPhi_rpminus_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDphi_rpminus_accept"] = new TH1F("muonDphi_rpminus_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDeltaY_rpminus_accept"] = new TH1F("muonDeltaY_rpminus_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);


	// Reco Muons rp_accept && _tsel
	histosTH1F["dimuon_mass_rpminus_accept_tsel"] = new TH1F("dimuon_mass_rpminus_accept_tsel", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt_rpminus_accept_tsel"] = new TH1F("dimuon_pt_rpminus_accept_tsel", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2_rpminus_accept_tsel"] = new TH1F("dimuon_pt2_rpminus_accept_tsel", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta_rpminus_accept_tsel"] = new TH1F("dimuon_eta_rpminus_accept_tsel", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_rpminus_accept_tsel"] = new TH1F("dimuon_rapidity_rpminus_accept_tsel", "y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["dimuon_multiplicity_rp_accept_tsel"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);
	histosTH1F["muonDeltaPt_rpminus_accept_tsel"] = new TH1F("muonDeltaPt_rpminus_accept_tsel", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_rpminus_accept_tsel"] = new TH1F("muonDeltaEta_rpminus_accept_tsel", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
	histosTH1F["muonDeltaPhi_rpminus_accept_tsel"] = new TH1F("muonDeltaPhi_rpminus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDphi_rpminus_accept_tsel"] = new TH1F("muonDphi_rpminus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muonDeltaY_rpminus_accept_tsel"] = new TH1F("muonDeltaY_rpminus_accept_tsel", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["dimuon_xi_minus_rpminus_accept_tsel"]= new TH1F("dimuon_xi_minus_rpminus_accept_tsel", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["dimuon_t_minus_rpminus_accept_tsel"] = new TH1F("dimuon_t_minus_rpminus_accept_tsel", "|t|^{+}" , 100,0.,5.0);

	//jpsi mass selection
	histosTH1F["jpsi_mass_rpminus_accept"] = new TH1F("jpsi_mass_rpminus_accept", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_rpminus_accept"] = new TH1F("jpsi_pt_rpminus_accept", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_rpminus_accept"] = new TH1F("jpsi_pt2_rpminus_accept", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_rpminus_accept"] = new TH1F("jpsi_eta_rpminus_accept", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_rpminus_accept"] = new TH1F("jpsi_rapidity_rpminus_accept", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["jpsi_multiplicity_rpplus_accept"] = new TH1F("jpsi_multiplicity_rp_accept", "n jpsi(dimuons)" , 100 , 0 , 100);
	histosTH1F["jpsi_xi_minus_rpminus_accept"]= new TH1F("jpsi_xi_minus_rpminus_accept", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["jpsi_t_minus_rpminus_accept"] = new TH1F("jpsi_t_minus_rpminus_accept", "|t|^{+}" , 100,0.,5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_rpminus_accept"] = new TH1F("muonDeltaPt_jpsi_rpminus_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_rpminus_accept"] = new TH1F("muonDeltaEta_jpsi_rpminus_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_rpminus_accept"] = new TH1F("muonDeltaPhi_jpsi_rpminus_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDphi_jpsi_rpminus_accept"] = new TH1F("muonDphi_jpsi_rpminus_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_rpminus_accept"] = new TH1F("muonDeltaY_jpsi_rpminus_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

	//jpsi rp accept && _tsel
	//jpsi mass selection
	histosTH1F["jpsi_mass_rpminus_accept_tsel"] = new TH1F("jpsi_mass_rpminus_accept_tsel", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_pt_rpminus_accept_tsel"] = new TH1F("jpsi_pt_rpminus_accept_tsel", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
	histosTH1F["jpsi_pt2_rpminus_accept_tsel"] = new TH1F("jpsi_pt2_rpminus_accept_tsel", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["jpsi_eta_rpminus_accept_tsel"] = new TH1F("jpsi_eta_rpminus_accept_tsel", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_rapidity_rpminus_accept_tsel"] = new TH1F("jpsi_rapidity_rpminus_accept_tsel", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
	//histosTH1F["jpsi_multiplicity_rp_accept_tsel"] = new TH1F("jpsi_multiplicity_rp_accept_tsel", "n jpsi(dimuons)" , 100 , 0 , 100);
	histosTH1F["jpsi_xi_minus_rpminus_accept_tsel"]= new TH1F("jpsi_xi_minus_rpminus_accept_tsel", "#xi^{+}" , 100 , 0.,1.);
	histosTH1F["jpsi_t_minus_rpminus_accept_tsel"] = new TH1F("jpsi_t_minus_rpminus_accept_tsel", "|t|^{+}" , 100,0.,5.0);

	//Deltas info
	histosTH1F["muonDeltaPt_jpsi_rpminus_accept_tsel"] = new TH1F("muonDeltaPt_jpsi_rpminus_accept_tsel", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["muonDeltaEta_jpsi_rpminus_accept_tsel"] = new TH1F("muonDeltaEta_jpsi_rpminus_accept_tsel", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
	histosTH1F["muonDeltaPhi_jpsi_rpminus_accept_tsel"] = new TH1F("muonDeltaPhi_jpsi_rpminus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDphi_jpsi_rpminus_accept_tsel"] = new TH1F("muonDphi_jpsi_rpminus_accept_tsel", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
	histosTH1F["muonDeltaY_jpsi_rpminus_accept_tsel"] = new TH1F("muonDeltaY_jpsi_rpminus_accept_tsel", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
	//============

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	//===================
	//int i_tot = 0 , nevt_tot = 0;
	vector<TString>* vfiles = new vector<TString>; 
	for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );

	//Declaration of tree and its branches variables
	//TTree* tree = new TTree(treeName.c_str(),"");
	TTree* tree = NULL;
	MyEvtId*           evtId        = NULL;
	MyL1TrigOld*       l1Trig       = NULL;  
	MyHLTrig*          hltTrig      = NULL;
	vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	//vector<MyPFJet>*   pfJet_coll   = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	vector<MyFSCHit>*  fscHits_coll = NULL;
	vector<MyFSCDigi>* fscDigis_coll = NULL;
	MyGenKin*          genKin        = NULL;
	//===================
	T1Event* t1_event = NULL;
	T2Event* t2_event = NULL; 
	RPRootDumpReconstructedProton* rec_proton_left  = NULL;
	RPRootDumpReconstructedProton* rec_proton_right = NULL;
	RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
	RPRootDumpReconstructedProton* sim_proton_right = NULL;
	RPRootDumpReconstructedProton* sim_proton_left = NULL;
	map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
	map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
	map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
	map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
	//===================
	//if(!isMC){
	TString outtxt_right = outputFileName;
	TString outtxt_left = outputFileName;

	outtxt_right.ReplaceAll("root","txt_right");
	outtxt_left.ReplaceAll("root","txt_left");  

	//const string left_ = "left_"+ outtxt_left ;
	//std::string right_ = "right_"+ outtxt_right ;
	ofstream outstring_left(outtxt_left); 
	ofstream outstring_right(outtxt_right);//}
	//===================
	//===================
	if(isMC) rp_aperture_config();

	//if(verbose)cout<<"rp_aperture_config"<<endl;
	int n_vertices_selected = 0;
	int n_select_Vertex_After_vtx_cut =0;
	int n_tracks_selected = 0;
	int n_muon_selected = 0;
	int n_dimuons_selected = 0;
	int n_dimuons_t_selected_proton_left = 0;
	//int n_dimuons_t_at_xi_selected_proton_left = 0;
	int n_dimuons_t_selected_proton_right = 0;
	//int n_dimuons_t_at_xi_selected_proton_right = 0;
	int n_jpsi_selected = 0;
	int n_jpsi_t_selected_proton_left = 0;
	//int n_jpsi_t_at_xi_selected_proton_left = 0;
	int n_jpsi_t_selected_proton_right = 0;
	//int n_jpsi_t_at_xi_selected_proton_right = 0;
	int n_jpsi_proton_left = 0;
	int n_jpsi_proton_right=0;
	int n_dimuons_proton_left = 0;
	int n_dimuons_proton_right=0;
	int n_dimuons_rp_selected_plus_side = 0;
	int n_dimuons_rp_selected_tsel_plus_side =0;
	int n_dimuons_rp_selected_minus_side = 0;
	int n_dimuons_rp_selected_tsel_minus_side =0;

	int n_jpsi_selected_rp_accept_plus_side = 0;
	int n_jpsi_selected_rp_accept_tsel_plus_side = 0;
	int n_jpsi_selected_rp_accept_minus_side = 0;
	int n_jpsi_selected_rp_accept_tsel_minus_side = 0;

	int i_tot = 0 , nevt_tot = 0;
	if(verbose)cout<<"i_tot"<<endl;
	//starting Loop over files, stops at end of list of files or when reached nevt_max
	for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){

		TFile* file = TFile::Open(*itfiles,"READ");

		//getting the tree form the current file
		tree = (TTree*) file->Get( treeName.c_str() );

		//Getting number of events
		int nev = int(tree->GetEntriesFast());

		nevt_tot += nev;
		cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

		// Add branches to TTree ----------------------------------------------------------------------
		if(verbose)cout<<"branches"<<endl;

		if(isMC){
			tree->SetBranchAddress("evtId",&evtId);
			tree->SetBranchAddress("generalTracks",&track_coll); 
			tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
			tree->SetBranchAddress("muons",&muon_coll);
			tree->SetBranchAddress("particleFlow",&pFlow_coll);
			tree->SetBranchAddress("genKin",&genKin);
			tree->SetBranchAddress("genPart",&genPart);
			//if(verbose)cout<<"add branches to ttree"<<endl;

		}else{ 
			tree->SetBranchAddress("cmsEvtUA",&evtId);
			tree->SetBranchAddress("cmsTrigUA",&l1Trig);
			tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
			tree->SetBranchAddress("cmsTracksUA",&track_coll); 
			tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
			//tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
			//tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
			tree->SetBranchAddress("cmsMuonsUA",&muon_coll);
			tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
			tree->SetBranchAddress("branchT1EV.",&t1_event);
			tree->SetBranchAddress("branchT2EV.",&t2_event);
			tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
			tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
			tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
			tree->SetBranchAddress("sim_prot_left.",&sim_proton_left);
			tree->SetBranchAddress("sim_prot_right.",&sim_proton_right);//}
			std::vector<unsigned int> rp_list;
			rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(22); rp_list.push_back(23); rp_list.push_back(24); rp_list.push_back(25);
			rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(122); rp_list.push_back(123); rp_list.push_back(124); rp_list.push_back(125);
			char br_name[200];
			char name[200];
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
	}


	//if(verbose)cout<<"pass MC 2"<<endl;




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
		bool passed_HLTMuon= false;
		string HLT_muon = "HLT_L1DoubleMu0_v1"; 
		//bool passed_vrt= false;      
		bool passed_muon = false;
		bool mu1_selected = false;
		bool mu2_selected = false;
		bool chmuons = false;
		bool select_Track = false;
		histosTH1F["EventSelection"]->Fill( "All", event_weight );
		//========================================================

		bool proton_plus_rp_accept = false;
		bool proton_minus_rp_accept = false;
		bool proton_plus_xi_range = false;
		bool proton_plus_t_range = false;
		bool proton_minus_xi_range = false;
		bool proton_minus_t_range = false;
		double xi_proton_plus = -1.;
		double xi_proton_minus = -1.;
		double t_proton_plus = 0.;
		double t_proton_minus = 0.;
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
				//if(verbose)cout<<"pass MC 3 GenPart"<<endl;
				//double eta_gen = it_genpart->Eta();
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

				//double proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
				//double pz_cut = 0.75*proton_pi;
				//if ( fabs(pz_gen) < pz_cut) continue;

				if (pz_gen > proton_pz_plus)  { 
					proton_pz_plus = pz_gen; proton_energy_plus = energy_gen; 
					proton_px_plus = px_gen; proton_py_plus = py_gen;
				} 
				if (pz_gen < proton_pz_minus) { 
					proton_pz_minus = pz_gen; proton_energy_minus = energy_gen; 
					proton_px_minus = px_gen; proton_py_minus = py_gen;
				} 
			}
			//double xi_plus_gen = genEPlusPz/8000.;
			//double xi_minus_gen = genEMinusPz/8000.;
			/*double xi_proton_plus = -1.;
			  double xi_proton_minus = -1.;
			  double t_proton_plus = 0.;
			  double t_proton_minus = 0.;*/
			double thx_proton_plus = 0.; 
			double thy_proton_plus = 0.; 
			double thx_proton_minus = 0.; 
			double thy_proton_minus = 0.; 

			bool proton_minus_rp_accept_120 = false;
			bool proton_minus_rp_accept_121 = false;
			bool proton_minus_rp_accept_122 = false;
			bool proton_minus_rp_accept_123 = false;
			bool proton_minus_rp_accept_124 = false;
			bool proton_minus_rp_accept_125 = false;
			bool proton_minus_rp_accept_020 = false;

			bool proton_plus_rp_accept_020 = false;
			bool proton_plus_rp_accept_021 = false;
			bool proton_plus_rp_accept_022 = false;
			bool proton_plus_rp_accept_023 = false;
			bool proton_plus_rp_accept_024 = false;
			bool proton_plus_rp_accept_025 = false;
			bool proton_plus_rp_accept_120 = false;

			std::map<int,std::vector<double> > proton_plus_pars;
			std::map<int,std::vector<double> > proton_minus_pars;

			if(proton_pz_plus > 0.){
				xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
				//t_proton_plus = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
				TLorentzVector vec_pi(0.,0.,proton_pi,proton_pi);
				TLorentzVector vec_pf(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus);
				TLorentzVector vec_t = (vec_pf - vec_pi);
				t_proton_plus = vec_t.Mag2();

				thx_proton_plus = atan(-proton_px_plus/proton_pi);
				thy_proton_plus = atan(proton_py_plus/proton_pi);

				//FIXME
				double out_x, out_thx, out_y, out_thy, out_xi;
				proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[20] = std::vector<double>(5,0.);
				proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
				proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
				proton_plus_pars[20][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_020"]->Fill( proton_plus_pars[20][0], proton_plus_pars[20][1] , event_weight );

				//proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21);
				proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[21] = std::vector<double>(5,0.);
				proton_plus_pars[21][0] = out_x; proton_plus_pars[21][1] = out_y;
				proton_plus_pars[21][2] = out_thx; proton_plus_pars[21][3] = out_thy;
				proton_plus_pars[21][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_021"]->Fill( proton_plus_pars[21][0], proton_plus_pars[21][1] , event_weight );

				//proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
				proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[22] = std::vector<double>(5,0.);
				proton_plus_pars[22][0] = out_x; proton_plus_pars[22][1] = out_y;
				proton_plus_pars[22][2] = out_thx; proton_plus_pars[22][3] = out_thy;
				proton_plus_pars[22][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_022"]->Fill( proton_plus_pars[22][0], proton_plus_pars[22][1] , event_weight );

				//proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
				proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[23] = std::vector<double>(5,0.);
				proton_plus_pars[23][0] = out_x; proton_plus_pars[23][1] = out_y;
				proton_plus_pars[23][2] = out_thx; proton_plus_pars[23][3] = out_thy;
				proton_plus_pars[23][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_023"]->Fill( proton_plus_pars[23][0], proton_plus_pars[23][1] , event_weight );

				//proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
				proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[24] = std::vector<double>(5,0.);
				proton_plus_pars[24][0] = out_x; proton_plus_pars[24][1] = out_y;
				proton_plus_pars[24][2] = out_thx; proton_plus_pars[24][3] = out_thy;
				proton_plus_pars[24][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_024"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );

				//proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
				proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25, out_x, out_thx, out_y, out_thy, out_xi);
				proton_plus_pars[25] = std::vector<double>(5,0.);
				proton_plus_pars[25][0] = out_x; proton_plus_pars[25][1] = out_y;
				proton_plus_pars[25][2] = out_thx; proton_plus_pars[25][3] = out_thy;
				proton_plus_pars[25][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_plus_025"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );

				proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);

				histosTH1F["xi_proton_plus"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus"]->Fill( fabs(t_proton_plus) , event_weight );
				//cout<<"t_proton_plus"<<fabs(t_proton_plus)<<endl;
				proton_plus_xi_range = ( xi_proton_plus >= 0.03 && xi_proton_plus < 0.2 );
				proton_plus_t_range = ( fabs(t_proton_plus) < 1.0 );
				if(proton_plus_t_range){
					histosTH1F["xi_proton_t_range_plus"]->Fill( xi_proton_plus , event_weight );
					if(proton_plus_xi_range){
						histosTH1F["t_proton_xi_range_plus"]->Fill( fabs(t_proton_plus) , event_weight );}
				}

				histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
				histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );

				histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
			}

			if(proton_pz_minus < 0.){ 
				xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
				//t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) );
				TLorentzVector vec_pi(0.,0.,-proton_pi,proton_pi);
				TLorentzVector vec_pf(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus);
				TLorentzVector vec_t = (vec_pf - vec_pi);
				t_proton_minus = vec_t.Mag2();

				thx_proton_minus = atan(-proton_px_minus/proton_pi);
				thy_proton_minus = atan(proton_py_minus/proton_pi);

				double out_x, out_thx, out_y, out_thy, out_xi;
				proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[120] = std::vector<double>(5,0.);
				proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
				proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
				proton_minus_pars[120][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_120"]->Fill( proton_minus_pars[120][0], proton_minus_pars[120][1] , event_weight );

				//proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121);
				proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[121] = std::vector<double>(5,0.);
				proton_minus_pars[121][0] = out_x; proton_minus_pars[121][1] = out_y;
				proton_minus_pars[121][2] = out_thx; proton_minus_pars[121][3] = out_thy;
				proton_minus_pars[121][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_121"]->Fill( proton_minus_pars[121][0], proton_minus_pars[121][1] , event_weight );

				//proton_minus_rp_accept_122 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 122);
				proton_minus_rp_accept_122 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 122, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[122] = std::vector<double>(5,0.);
				proton_minus_pars[122][0] = out_x; proton_minus_pars[122][1] = out_y;
				proton_minus_pars[122][2] = out_thx; proton_minus_pars[122][3] = out_thy;
				proton_minus_pars[122][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_122"]->Fill( proton_minus_pars[122][0], proton_minus_pars[122][1] , event_weight );

				//proton_minus_rp_accept_123 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 123);
				proton_minus_rp_accept_123 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 123, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[123] = std::vector<double>(5,0.);
				proton_minus_pars[123][0] = out_x; proton_minus_pars[123][1] = out_y;
				proton_minus_pars[123][2] = out_thx; proton_minus_pars[123][3] = out_thy;
				proton_minus_pars[123][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_123"]->Fill( proton_minus_pars[123][0], proton_minus_pars[123][1] , event_weight );

				//proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124);
				proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[124] = std::vector<double>(5,0.);
				proton_minus_pars[124][0] = out_x; proton_minus_pars[124][1] = out_y;
				proton_minus_pars[124][2] = out_thx; proton_minus_pars[124][3] = out_thy;
				proton_minus_pars[124][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_124"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );

				//proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125);
				proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125, out_x, out_thx, out_y, out_thy, out_xi);
				proton_minus_pars[125] = std::vector<double>(5,0.);
				proton_minus_pars[125][0] = out_x; proton_minus_pars[125][1] = out_y;
				proton_minus_pars[125][2] = out_thx; proton_minus_pars[125][3] = out_thy;
				proton_minus_pars[125][4] = out_xi;
				histosTH2F["pos_y_vs_x_proton_minus_125"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );

				proton_minus_rp_accept_020 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 20);

				histosTH1F["xi_proton_minus"]->Fill( xi_proton_minus , event_weight );
				histosTH1F["t_proton_minus"]->Fill( fabs(t_proton_minus) , event_weight );

				proton_minus_xi_range = ( xi_proton_minus >= 0.03 && xi_proton_minus < 0.2 );
				proton_minus_t_range = ( fabs(t_proton_minus)> 0.03 && fabs(t_proton_minus) < 1.0 );
				if(proton_minus_t_range){
					histosTH1F["xi_proton_t_range_minus"]->Fill( xi_proton_minus , event_weight );
					if(proton_minus_xi_range){
						histosTH1F["t_proton_xi_range_minus"]->Fill( fabs(t_proton_minus) , event_weight );}
				}

				histosTH1F["thx_proton_minus"]->Fill( thx_proton_minus , event_weight );
				histosTH1F["thy_proton_minus"]->Fill( thy_proton_minus , event_weight );

				histosTH2F["proton_minus_xi_vs_t"]->Fill( fabs(t_proton_minus) , xi_proton_minus , event_weight );
			}

			// Check RP combinations  
			proton_plus_rp_accept = (proton_pz_plus > 0.) && 
				( ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 )|| 
				  ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 )||
				  ( proton_plus_rp_accept_022 && proton_plus_rp_accept_023 ));
			proton_minus_rp_accept = (proton_pz_minus < 0.) && 
				( ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 )||
				  ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 )||
				  ( proton_minus_rp_accept_122 && proton_minus_rp_accept_123 ) );

			if( proton_plus_rp_accept ){
				histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight ); 
				histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
				if(proton_plus_t_range){
					histosTH1F["xi_proton_t_range_plus_accepted"]->Fill( xi_proton_plus , event_weight );
					if(proton_plus_xi_range){
						histosTH1F["t_proton_xi_range_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight );}
				}
			}

			if( proton_minus_rp_accept ){
				histosTH1F["xi_proton_minus_accepted"]->Fill( xi_proton_minus , event_weight );
				histosTH1F["t_proton_minus_accepted"]->Fill( fabs(t_proton_minus) , event_weight ); 
				histosTH2F["proton_minus_xi_vs_t_accepted"]->Fill( fabs(t_proton_minus) , xi_proton_minus , event_weight );

				if(proton_minus_t_range){
					histosTH1F["xi_proton_t_range_minus_accepted"]->Fill( xi_proton_minus , event_weight );
					if(proton_minus_xi_range){
						histosTH1F["t_proton_xi_range_minus_accepted"]->Fill( fabs(t_proton_minus) , event_weight );}
				}
			}

			// RP stations
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
				histosTH2F["pos_y_vs_x_proton_plus_accepted_021"]->Fill( proton_plus_pars[21][0], proton_plus_pars[21][1] , event_weight );
			}
			if(proton_plus_rp_accept_022){
				histosTH1F["xi_proton_plus_accepted_022"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted_022"]->Fill( fabs(t_proton_plus) , event_weight ); 
				histosTH2F["pos_y_vs_x_proton_plus_accepted_022"]->Fill( proton_plus_pars[22][0], proton_plus_pars[22][1] , event_weight );
			}
			if(proton_plus_rp_accept_023){
				histosTH1F["xi_proton_plus_accepted_023"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted_023"]->Fill( fabs(t_proton_plus) , event_weight ); 
				histosTH2F["pos_y_vs_x_proton_plus_accepted_023"]->Fill( proton_plus_pars[23][0], proton_plus_pars[23][1] , event_weight );
			}
			if(proton_plus_rp_accept_024){
				histosTH1F["xi_proton_plus_accepted_024"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted_024"]->Fill( fabs(t_proton_plus) , event_weight ); 
				histosTH2F["pos_y_vs_x_proton_plus_accepted_024"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
			}
			if(proton_plus_rp_accept_025){
				histosTH1F["xi_proton_plus_accepted_025"]->Fill( xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted_025"]->Fill( fabs(t_proton_plus) , event_weight ); 
				histosTH2F["pos_y_vs_x_proton_plus_accepted_025"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
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
			if(proton_minus_rp_accept_121){
				histosTH2F["pos_y_vs_x_proton_minus_accepted_121"]->Fill( proton_minus_pars[121][0], proton_minus_pars[121][1] , event_weight );
			}
			if(proton_minus_rp_accept_122){
				histosTH2F["pos_y_vs_x_proton_minus_accepted_122"]->Fill( proton_minus_pars[122][0], proton_minus_pars[122][1] , event_weight );
			}
			if(proton_minus_rp_accept_123){
				histosTH2F["pos_y_vs_x_proton_minus_accepted_123"]->Fill( proton_minus_pars[123][0], proton_minus_pars[123][1] , event_weight );
			}
			if(proton_minus_rp_accept_124){
				histosTH2F["pos_y_vs_x_proton_minus_accepted_124"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );
			}
			if(proton_minus_rp_accept_125){
				histosTH2F["pos_y_vs_x_proton_minus_accepted_125"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );
			}

			if(proton_minus_rp_accept_020){
				histosTH1F["xi_proton_minus_accepted_020"]->Fill( xi_proton_minus , event_weight );
				histosTH1F["t_proton_minus_accepted_020"]->Fill( fabs(t_proton_minus) , event_weight ); 
			}
		}
		//if(verbose)cout<<"fim genpart"<<endl;
		//-------------------------------------------------------------------------------------------------
		// Event selection
		//-------------------------------------------------------------------------------------------------
		// Bunch crossing & HLT
		if(!isMC){
			//==========================================================
			//int orbitNumber = evtId->Orbit;
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
					// if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );

					if( hltName == HLT_muon){ 
						passed_HLTMuon = true;

						if( it_hlt->second ) histosTH1F["hltTrigFired"]->Fill( hltName.c_str(), event_weight );
					}
				}


			}

			if(!passed_HLTMuon) continue;

			histosTH1F["EventSelection"]->Fill( "HLT", event_weight );
		}
		//-------------------------------------------------------------------------------------------------
		//-------------------------------------------------------------------------------------------------
		/*int idx_vtx_max_sumpt = -1;
		  double sumpt_max = 0.;*/
		for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
			//int idx_vtx = it_vtx - vertex_coll->begin();
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
		// Include Muon selection REf CMS AN AN 12-067 
		//bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
		//                        primaryVertex.ndof < 10 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0); 

		if(selectVertex && !select_Vertex) continue;

		++n_select_Vertex_After_vtx_cut;


		histosTH1F["vertex_multiplicity_after_vtx_sel"]->Fill( n_select_Vertex_After_vtx_cut, event_weight );  
		histosTH1F["prim_vtx_zpos_after_vtx_sel"]->Fill( primaryVertex.z, event_weight );
		histosTH1F["prim_vtx_xpos_after_vtx_sel"]->Fill( primaryVertex.x, event_weight );
		histosTH1F["prim_vtx_ypos_after_vtx_sel"]->Fill( primaryVertex.y, event_weight );

		histosTH1F["prim_vtx_ndof_after_vtx_sel"]->Fill( primaryVertex.ndof, event_weight );
		histosTH1F["prim_vtx_chi2_after_vtx_sel"]->Fill( primaryVertex.chi2, event_weight );
		histosTH1F["prim_vtx_chi2n_after_vtx_sel"]->Fill( primaryVertex.chi2n(), event_weight );
		histosTH1F["prim_vtx_ntracks_after_vtx_sel"]->Fill( primaryVertex.ntracks, event_weight );
		histosTH1F["prim_vtx_sumpt_after_vtx_sel"]->Fill( primaryVertex.SumPtTracks, event_weight );

		histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

		int prim_vtx_id = primaryVertex.id;

		// Tracks
		//int n_tracks_selected = 0;
		for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
			histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
			histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
			histosTH1F["track_rapidity"]->Fill( it_trk->Rapidity(), event_weight );
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

			select_Track = true;
			++n_tracks_selected;

		}
		histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );
		if(selectTrack && !select_Track) continue;  

		// Muons variables
		//int n_muon_selected = 0;
		double chmu1 = 0.; double chmu2 = 0.;
		double phimu1 = 0.; double ptmu1 = 0.;double etamu1 = 0.;double ymu1 = 0.;
		double phimu2 = 0.; double ptmu2 = 0.;double etamu2 = 0.;double ymu2 = 0.;
		double deltaphi = 0.; double deltaeta = 0.; double deltapt = 0.; double deltay = 0.; double Dphi = 0.;
		double dimuon_mass = 0.; double dimuon_pt=0.; double dimuon_pt2=0.; double dimuon_eta =0.;
		double dimuon_rapidity = 0.; double jpsi_mass = 0.; double jpsi_pt=0.;
		double jpsi_pt2=0.; double jpsi_eta =0.; double jpsi_rapidity = 0.;
		// Muons
		vector<MyMuon> muons_selected;
		for(vector<MyMuon>::iterator it_muon = muon_coll->begin() ; it_muon != muon_coll->end() ; ++it_muon){

			if( !(it_muon->IsTrackerMuon || it_muon->IsGlobalMuon) ) continue;

			MyTracks const& muon_innerTrack = it_muon->innerTrack;
			bool muon_id = it_muon->TMOneStationAngTight &&
				muon_innerTrack.chi2n < 1.8 &&
				muon_innerTrack.nValidPixelHits > 0 &&
				muon_innerTrack.vtxdxy[prim_vtx_id] < 3. &&
				muon_innerTrack.vtxdz[prim_vtx_id] < 30.;

			if( !muon_id ) continue;

			++n_muon_selected;

			histosTH1F["muon_pt"]->Fill( it_muon->Pt(), event_weight );
			histosTH1F["muon_eta"]->Fill( it_muon->Eta(), event_weight );
			histosTH1F["muon_phi"]->Fill( it_muon->Phi(), event_weight );

			muons_selected.push_back( *it_muon );

			histosTH1F["muon_multiplicity"]->Fill( n_muon_selected, event_weight );
		}
		bool select_Muons = ( muons_selected.size() >= 2 );
		if(selectMuons && !select_Muons) continue;
		histosTH1F["EventSelection"]->Fill( "Muons", event_weight );
		if(verbose)cout<<"muons selected size (>=2) :"<<muons_selected.size()<<endl;
		//----------------------------------------------------------------------------------------
		// Muon pairs
		for(vector<MyMuon>::iterator it_mu1 = muons_selected.begin() ;
				it_mu1 != muons_selected.end() ; ++it_mu1){
			histosTH1F["muon1_pt"]->Fill( it_mu1->Pt(), event_weight );
			histosTH1F["muon1_eta"]->Fill( it_mu1->Eta(), event_weight );
			histosTH1F["muon1_rapidity"]->Fill( it_mu1->Rapidity(), event_weight );
			histosTH1F["muon1_phi"]->Fill(it_mu1->Phi(), event_weight );
			phimu1 = it_mu1->Phi();  ptmu1 = it_mu1->Pt(); 
			etamu1 = it_mu1->Eta(); ymu1 = it_mu1->Rapidity();
			chmu1 = it_mu1->charge;
			mu1_selected = true;

			for(vector<MyMuon>::iterator it_mu2 = muons_selected.begin() ;
					it_mu2 != muons_selected.end() ; ++it_mu2){
				bool os_muons = ( it_mu1->charge*it_mu2->charge < 0. );
				if( !os_muons ) continue;
				++n_dimuons_selected;
				mu2_selected = true;
				//cout<<"mu2_selected :"<<mu2_selected<<endl;
				histosTH1F["muon2_pt"]->Fill( it_mu2->Pt(), event_weight );
				histosTH1F["muon2_eta"]->Fill( it_mu2->Eta(), event_weight );
				histosTH1F["muon2_rapidity"]->Fill( it_mu2->Rapidity(), event_weight );
				histosTH1F["muon2_phi"]->Fill(it_mu2->Phi(), event_weight );
				phimu2 = it_mu2->Phi(); ptmu2 = it_mu2->Pt();
				etamu2 = it_mu2->Eta(); ymu2 = it_mu2->Rapidity();
				chmu2 = it_mu2->charge;

				//...
				TLorentzVector& muon1_lorentz = *it_mu1;
				TLorentzVector& muon2_lorentz = *it_mu2;
				TLorentzVector dimuon_lorentz(0.,0.,0.,0.);
				dimuon_lorentz += muon1_lorentz;
				dimuon_lorentz += muon2_lorentz;
				//histosTH1F["dimuon_multiplicity"]->Fill(n_dimuons_selected, event_weight );
				histosTH1F["dimuon_mass"]->Fill( dimuon_lorentz.M(), event_weight );
				histosTH1F["dimuon_pt"]->Fill( dimuon_lorentz.Pt(), event_weight );
				histosTH1F["dimuon_pt2"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );
				histosTH1F["dimuon_eta"]->Fill( dimuon_lorentz.Eta(), event_weight );
				histosTH1F["dimuon_rapidity"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
				dimuon_mass = dimuon_lorentz.M();  
				dimuon_pt=dimuon_lorentz.Pt(); 
				dimuon_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
				dimuon_eta =dimuon_lorentz.Eta();
				dimuon_rapidity = dimuon_lorentz.Rapidity();
				//cout<<"Delta Definitions :"<<endl;		
				deltapt = fabs(muon1_lorentz.Pt() - muon2_lorentz.Pt());
				deltaeta = fabs(muon1_lorentz.Eta() - muon2_lorentz.Eta());
				deltaphi = fabs(phimu1 - phimu2);
				if(deltaphi > M_PI)deltaphi = (2*M_PI - deltaphi);
				deltay = fabs(muon1_lorentz.Rapidity() - muon2_lorentz.Rapidity());
				Dphi = std::fabs(deltaPhi(phimu1 ,phimu2));  
				//cout<<"Delta Phi :"<<deltaphi<<endl;
				histosTH1F["muonDeltaPt"]->Fill(deltapt, event_weight );
				histosTH1F["muonDeltaEta"]->Fill(deltaeta, event_weight );	
				histosTH1F["muonDeltaPhi"]->Fill(deltaphi, event_weight );
				histosTH1F["muonDeltaY"]->Fill(deltay, event_weight ); 
				//histosTH1F["muonDphi"]->Fill(Dphi, event_weight );
				//cout<<"mass_dimuon = "<< dimuon_lorentz.M() << endl;
				////////////////////////////////////////////////////////////////////////////////////////////
				//double Jpsi_mass = dimuon_lorentz.M();
				//cout<<"jpsi_mass_dimuon = "<< Jpsi_mass << endl;                              
				if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){
					//cout<<"jpsi_mass = "<< dimuon_mass << endl;
					++n_jpsi_selected;
					//double Dphi_jpsi = std::fabs(deltaPhi(phimu1 ,phimu2));
					//histosTH1F["jpsi_multiplicity"]->Fill(n_dimuons_selected, event_weight );
					histosTH1F["jpsi_mass"]->Fill( dimuon_lorentz.M(), event_weight );
					histosTH1F["jpsi_pt"]->Fill( dimuon_lorentz.Pt(), event_weight );
					histosTH1F["jpsi_pt2"]->Fill( (dimuon_lorentz.Pt()*dimuon_lorentz.Pt()), event_weight );
					histosTH1F["jpsi_eta"]->Fill( dimuon_lorentz.Eta(), event_weight );
					histosTH1F["jpsi_rapidity"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
					histosTH1F["muonDeltaPt_jpsi"]->Fill(deltapt, event_weight );
					histosTH1F["muonDeltaEta_jpsi"]->Fill(deltaeta, event_weight );	
					histosTH1F["muonDeltaPhi_jpsi"]->Fill(deltaphi, event_weight );
					//histosTH1F["muonDphi_jpsi"]->Fill(Dphi_jpsi, event_weight );
					histosTH1F["muonDeltaY_jpsi"]->Fill(deltay, event_weight ); 
					jpsi_mass = dimuon_lorentz.M();  
					jpsi_pt=dimuon_lorentz.Pt(); 
					jpsi_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
					jpsi_eta =dimuon_lorentz.Eta();
					jpsi_rapidity = dimuon_lorentz.Rapidity();
					//cout<<"finalizando jpsi mass cut"<<endl; 
				}   
				//cout<<"final do loop mass range"<<endl;
			}
		}
		if(!mu2_selected)continue;

		//-----------------------------------------------------------------------------------------
		// Particle-flow
		vector<MyPFCand> particles_sorted;
		double etaEdgeLow = -999.0;
		double etaEdgeHigh = 999.0;

		double pfEPlusPz = 0.;
		double pfEMinusPz = 0.;
		double xiCorrFactor = 1.0;
		for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin();
				it_pfcand != pFlow_coll->end(); ++it_pfcand){
			//int partType = it_pfcand->particleId;
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
		////if( selectEtaMax && (pfEtaMax > etaMaxThreshold) ) continue;
		histosTH1F["EventSelection"]->Fill( "EtaMax", event_weight );

		//if( selectEtaMin && (pfEtaMin < -3.0) ) continue;
		////if( selectEtaMin && (pfEtaMin < -etaMaxThreshold) ) continue;
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
		//=============================================================================
		if(!isMC){
			proton_right_valid = rec_proton_right->valid;
			proton_left_valid = rec_proton_left->valid;
			bool double_arm_rec_proton = (proton_right_valid && proton_left_valid);
			bool single_arm_rec_proton = (proton_right_valid || proton_left_valid) && !double_arm_rec_proton;

			if( selectSingleArmRecProton && !single_arm_rec_proton ) continue;
			histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

			if( selectDoubleArmRecProton && !double_arm_rec_proton ) continue;
			histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

			bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
			bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
			if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
			histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

			if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
			histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );
		} else if(isMC){
			if( selectRPPlusAccept && !proton_plus_rp_accept ) continue;
			histosTH1F["EventSelection"]->Fill( "MCRPPlusAccept", event_weight );

			if( selectRPMinusAccept && !proton_minus_rp_accept ) continue;
			histosTH1F["EventSelection"]->Fill( "MCRPMinusAccept", event_weight );
		}


		// RP protons
		if(!isMC){
			RPRootDumpTrackInfo* rp_track_info_020 = rp_track_info[20];
			RPRootDumpTrackInfo* rp_track_info_021 = rp_track_info[21];
			RPRootDumpTrackInfo* rp_track_info_022 = rp_track_info[22];
			RPRootDumpTrackInfo* rp_track_info_023 = rp_track_info[23];
			RPRootDumpTrackInfo* rp_track_info_024 = rp_track_info[24];
			RPRootDumpTrackInfo* rp_track_info_025 = rp_track_info[25];

			RPRootDumpTrackInfo* rp_track_info_120 = rp_track_info[120];
			RPRootDumpTrackInfo* rp_track_info_121 = rp_track_info[121];
			RPRootDumpTrackInfo* rp_track_info_122 = rp_track_info[122];
			RPRootDumpTrackInfo* rp_track_info_123 = rp_track_info[123];
			RPRootDumpTrackInfo* rp_track_info_124 = rp_track_info[124];
			RPRootDumpTrackInfo* rp_track_info_125 = rp_track_info[125];

			bool rp_track_valid_020 = rp_track_info_020->valid;
			if( rp_track_valid_020 ){
				double rp_track_posx_020 = rp_track_info_020->x;
				double rp_track_posy_020 = rp_track_info_020->y;
				histosTH1F["rp_track_posx_020"]->Fill( rp_track_posx_020, event_weight );
				histosTH1F["rp_track_posy_020"]->Fill( rp_track_posy_020, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_020"]->Fill( rp_track_posx_020, rp_track_posy_020, event_weight );
			}
			bool rp_track_valid_021 = rp_track_info_021->valid;
			if( rp_track_valid_021 ){
				double rp_track_posx_021 = rp_track_info_021->x;
				double rp_track_posy_021 = rp_track_info_021->y;
				histosTH1F["rp_track_posx_021"]->Fill( rp_track_posx_021, event_weight );
				histosTH1F["rp_track_posy_021"]->Fill( rp_track_posy_021, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_021"]->Fill( rp_track_posx_021, rp_track_posy_021, event_weight );
			}
			bool rp_track_valid_022 = rp_track_info_022->valid;
			if( rp_track_valid_022 ){
				double rp_track_posx_022 = rp_track_info_022->x;
				double rp_track_posy_022 = rp_track_info_022->y;
				histosTH1F["rp_track_posx_022"]->Fill( rp_track_posx_022, event_weight );
				histosTH1F["rp_track_posy_022"]->Fill( rp_track_posy_022, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_022"]->Fill( rp_track_posx_022, rp_track_posy_022, event_weight );
			}
			bool rp_track_valid_023 = rp_track_info_023->valid;
			if( rp_track_valid_023 ){
				double rp_track_posx_023 = rp_track_info_023->x;
				double rp_track_posy_023 = rp_track_info_023->y;
				histosTH1F["rp_track_posx_023"]->Fill( rp_track_posx_023, event_weight );
				histosTH1F["rp_track_posy_023"]->Fill( rp_track_posy_023, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_023"]->Fill( rp_track_posx_023, rp_track_posy_023, event_weight );
			}
			bool rp_track_valid_024 = rp_track_info_024->valid;
			if( rp_track_valid_024 ){
				double rp_track_posx_024 = rp_track_info_024->x;
				double rp_track_posy_024 = rp_track_info_024->y;
				histosTH1F["rp_track_posx_024"]->Fill( rp_track_posx_024, event_weight );
				histosTH1F["rp_track_posy_024"]->Fill( rp_track_posy_024, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_024"]->Fill( rp_track_posx_024, rp_track_posy_024, event_weight );
			}
			bool rp_track_valid_025 = rp_track_info_025->valid;
			if( rp_track_valid_025 ){
				double rp_track_posx_025 = rp_track_info_025->x;
				double rp_track_posy_025 = rp_track_info_025->y;
				histosTH1F["rp_track_posx_025"]->Fill( rp_track_posx_025, event_weight );
				histosTH1F["rp_track_posy_025"]->Fill( rp_track_posy_025, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_025"]->Fill( rp_track_posx_025, rp_track_posy_025, event_weight );
			}

			bool rp_track_valid_120 = rp_track_info_120->valid;
			if( rp_track_valid_120 ){
				double rp_track_posx_120 = rp_track_info_120->x;
				double rp_track_posy_120 = rp_track_info_120->y;
				histosTH1F["rp_track_posx_120"]->Fill( rp_track_posx_120, event_weight );
				histosTH1F["rp_track_posy_120"]->Fill( rp_track_posy_120, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_120"]->Fill( rp_track_posx_120, rp_track_posy_120, event_weight );
			}
			bool rp_track_valid_121 = rp_track_info_121->valid;
			if( rp_track_valid_121 ){
				double rp_track_posx_121 = rp_track_info_121->x;
				double rp_track_posy_121 = rp_track_info_121->y;
				histosTH1F["rp_track_posx_121"]->Fill( rp_track_posx_121, event_weight );
				histosTH1F["rp_track_posy_121"]->Fill( rp_track_posy_121, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_121"]->Fill( rp_track_posx_121, rp_track_posy_121, event_weight );
			}
			bool rp_track_valid_122 = rp_track_info_122->valid;
			if( rp_track_valid_122 ){
				double rp_track_posx_122 = rp_track_info_122->x;
				double rp_track_posy_122 = rp_track_info_122->y;
				histosTH1F["rp_track_posx_122"]->Fill( rp_track_posx_122, event_weight );
				histosTH1F["rp_track_posy_122"]->Fill( rp_track_posy_122, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_122"]->Fill( rp_track_posx_122, rp_track_posy_122, event_weight );
			}
			bool rp_track_valid_123 = rp_track_info_123->valid;
			if( rp_track_valid_123 ){
				double rp_track_posx_123 = rp_track_info_123->x;
				double rp_track_posy_123 = rp_track_info_123->y;
				histosTH1F["rp_track_posx_123"]->Fill( rp_track_posx_123, event_weight );
				histosTH1F["rp_track_posy_123"]->Fill( rp_track_posy_123, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_123"]->Fill( rp_track_posx_123, rp_track_posy_123, event_weight );
			}
			bool rp_track_valid_124 = rp_track_info_124->valid;
			if( rp_track_valid_124 ){
				double rp_track_posx_124 = rp_track_info_124->x;
				double rp_track_posy_124 = rp_track_info_124->y;
				histosTH1F["rp_track_posx_124"]->Fill( rp_track_posx_124, event_weight );
				histosTH1F["rp_track_posy_124"]->Fill( rp_track_posy_124, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_124"]->Fill( rp_track_posx_124, rp_track_posy_124, event_weight );
			}
			bool rp_track_valid_125 = rp_track_info_125->valid;
			if( rp_track_valid_125 ){
				double rp_track_posx_125 = rp_track_info_125->x;
				double rp_track_posy_125 = rp_track_info_125->y;
				histosTH1F["rp_track_posx_125"]->Fill( rp_track_posx_125, event_weight );
				histosTH1F["rp_track_posy_125"]->Fill( rp_track_posy_125, event_weight );
				histosTH2F["rp_track_pos_y_vs_x_125"]->Fill( rp_track_posx_125, rp_track_posy_125, event_weight );
			}

			double chi2_proton_right = rec_proton_right->chi2;
			double chindf_proton_right = rec_proton_right->chindf;
			double xi_proton_right = -rec_proton_right->xi;
			//bool good_proton_right = proton_right_valid && (chi2_proton_right/chindf_proton_right > 1);
			bool good_proton_right = proton_right_valid && (xi_proton_right >= 0.);
			if( good_proton_right ){
				double t_proton_right = rec_proton_right->t;
				histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight );
				histosTH1F["proton_right_xi"]->Fill( xi_proton_right, event_weight );
				histosTH1F["proton_right_t"]->Fill( -t_proton_right, event_weight );

				histosTH2F["proton_right_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, xi_proton_right, event_weight );
				histosTH2F["proton_right_xi_vs_pf_xiMinus"]->Fill( pfXiMinusReco, xi_proton_right, event_weight );
				if(xi_proton_right > 0.){
					histosTH1F["proton_right_logXi"]->Fill( log10(xi_proton_right), event_weight );
					histosTH2F["proton_right_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_right), event_weight );
					histosTH2F["proton_right_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_right), event_weight );
					histosTH2F["proton_right_logXi_vs_t"]->Fill( -t_proton_right, log10(xi_proton_right), event_weight );
				}

				histosTH1F["pf_xiMinus_minus_proton_right_xi"]->Fill( (pfXiMinusReco - xi_proton_right), event_weight );
				/*if( (pfXiMinusReco - xi_proton_right) < 0. ){
				  histosTH1F["proton_right_xi_selected"]->Fill( xi_proton_right, event_weight );
				  histosTH1F["proton_right_t_selected"]->Fill( -t_proton_right, event_weight );
				  }*/
			}

			double chi2_proton_left = rec_proton_left->chi2;
			double chindf_proton_left = rec_proton_left->chindf;
			double xi_proton_left = -rec_proton_left->xi;
			//bool good_proton_left = proton_left_valid && (chi2_proton_left/chindf_proton_left > 1);
			bool good_proton_left = proton_left_valid && (xi_proton_left >= 0.);
			if( good_proton_left ){
				double t_proton_left = rec_proton_left->t;
				histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
				histosTH1F["proton_left_xi"]->Fill( xi_proton_left, event_weight );
				histosTH1F["proton_left_t"]->Fill( -t_proton_left, event_weight );

				histosTH2F["proton_left_xi_vs_pf_xiPlus"]->Fill( pfXiPlusReco, xi_proton_left, event_weight );
				histosTH2F["proton_left_xi_vs_pf_xiMinus"]->Fill( pfXiMinusReco, xi_proton_left, event_weight );
				if(xi_proton_left > 0.){
					histosTH1F["proton_left_logXi"]->Fill( log10(xi_proton_left), event_weight );
					histosTH2F["proton_left_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_left), event_weight );
					histosTH2F["proton_left_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_left), event_weight );
					histosTH2F["proton_left_logXi_vs_t"]->Fill( -t_proton_left, log10(xi_proton_left), event_weight );
				}

				/*if( pfJet_coll->size() > 0 ){
				  MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
				  histosTH2F["proton_left_t_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), -t_proton_left, event_weight );
				  }*/

				histosTH1F["pf_xiPlus_minus_proton_left_xi"]->Fill( (pfXiPlusReco - xi_proton_left), event_weight );

				/*if( (pfXiPlusReco - xi_proton_left) < 0. ){
				  histosTH1F["proton_left_xi_selected"]->Fill( xi_proton_left, event_weight );
				  histosTH1F["proton_left_t_selected"]->Fill( -t_proton_left, event_weight );
				  }*/
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

			//-------------------
			// After selection 
			//-------------------
			if( good_proton_right ){
				n_dimuons_proton_right++;
				double t_proton_right = rec_proton_right->t;
				histosTH1F["dimuon_proton_right_xi"]->Fill( xi_proton_right, event_weight );
				histosTH1F["dimuon_proton_right_t"]->Fill( fabs(t_proton_right), event_weight );

				histosTH1F["dimuon_mass_proton_right"]->Fill( dimuon_mass, event_weight );
				histosTH1F["dimuon_pt_proton_right"]->Fill( dimuon_pt, event_weight );
				histosTH1F["dimuon_pt2_proton_right"]->Fill( dimuon_pt2, event_weight );
				histosTH1F["dimuon_eta_proton_right"]->Fill( dimuon_eta, event_weight );
				histosTH1F["dimuon_rapidity_proton_right"]->Fill( dimuon_rapidity, event_weight );			
				histosTH1F["muonDeltaPt_proton_right"]->Fill(deltapt, event_weight );
				histosTH1F["muonDeltaEta_proton_right"]->Fill(deltaeta, event_weight );
				histosTH1F["muonDeltaPhi_proton_right"]->Fill(deltaphi, event_weight );
				histosTH1F["muonDeltaY_proton_right"]->Fill(deltay, event_weight );
				if (fabs(t_proton_right)>t_proton_down_ && fabs(t_proton_right)<t_proton_up_){
					//if(proton_plus_t_range){
					n_dimuons_t_selected_proton_right++;
					if(verbose) cout<< "t proton right side (abs):"<< "["<< fabs(t_proton_right)<< "]"<< endl;

					histosTH1F["dimuon_mass_t_cut_proton_right"]->Fill( dimuon_mass, event_weight );
					histosTH1F["dimuon_pt_t_cut_proton_right"]->Fill(dimuon_pt , event_weight );

					histosTH1F["dimuon_pt2_t_cut_proton_right"]->Fill( dimuon_pt2, event_weight );
					histosTH1F["dimuon_eta_t_cut_proton_right"]->Fill( dimuon_eta, event_weight );
					histosTH1F["dimuon_rapidity_t_cut_proton_right"]->Fill( dimuon_rapidity, event_weight );
					histosTH1F["proton_right_xi_t_cut"]->Fill( xi_proton_right, event_weight );
					histosTH1F["proton_right_t_t_cut"]->Fill( fabs(t_proton_right), event_weight );
					histosTH1F["muonDeltaPt_t_cut_proton_right"]->Fill(deltapt, event_weight );
					histosTH1F["muonDeltaEta_t_cut_proton_right"]->Fill(deltaeta, event_weight );	
					histosTH1F["muonDeltaPhi_t_cut_proton_right"]->Fill(deltaphi, event_weight );
					histosTH1F["muonDeltaY_t_cut_proton_right"]->Fill(deltay, event_weight ); 



				}
				}
				if( good_proton_left ){
					n_dimuons_proton_left++;
					double t_proton_left = rec_proton_left->t;
					histosTH1F["dimuon_proton_left_xi"]->Fill( xi_proton_left, event_weight );
					histosTH1F["dimuon_proton_left_t"]->Fill( fabs(t_proton_left), event_weight ); 

					histosTH1F["dimuon_mass_proton_left"]->Fill( dimuon_mass, event_weight );
					histosTH1F["dimuon_pt_proton_left"]->Fill( dimuon_pt, event_weight );
					histosTH1F["dimuon_pt2_proton_left"]->Fill( dimuon_pt2, event_weight );
					histosTH1F["dimuon_eta_proton_left"]->Fill( dimuon_eta, event_weight );
					histosTH1F["dimuon_rapidity_proton_left"]->Fill( dimuon_rapidity, event_weight );
					histosTH1F["muonDeltaPt_proton_left"]->Fill(deltapt, event_weight );
					histosTH1F["muonDeltaEta_proton_left"]->Fill(deltaeta, event_weight );
					histosTH1F["muonDeltaPhi_proton_left"]->Fill(deltaphi, event_weight );
					histosTH1F["muonDeltaY_proton_left"]->Fill(deltay, event_weight );
					if (fabs(t_proton_left)>t_proton_down_ && fabs(t_proton_left)<t_proton_up_){
						//if(proton_minus_t_range){
						n_dimuons_t_selected_proton_left++;	
						if(verbose) cout<< "t proton left side (abs):"<< "["<< fabs(t_proton_left)<< "]"<< endl;


						histosTH1F["dimuon_mass_t_cut_proton_left"]->Fill( dimuon_mass, event_weight );
						histosTH1F["dimuon_pt_t_cut_proton_left"]->Fill( dimuon_pt, event_weight );
						histosTH1F["dimuon_pt2_t_cut_proton_left"]->Fill(dimuon_pt2, event_weight );

						histosTH1F["dimuon_eta_t_cut_proton_left"]->Fill( dimuon_eta, event_weight );
						histosTH1F["dimuon_rapidity_t_cut_proton_left"]->Fill( dimuon_rapidity, event_weight );
						histosTH1F["proton_left_xi_t_cut"]->Fill( xi_proton_left, event_weight );
						histosTH1F["proton_left_t_t_cut"]->Fill( fabs(t_proton_left), event_weight );
						histosTH1F["muonDeltaPt_t_cut_proton_left"]->Fill(deltapt, event_weight );
						histosTH1F["muonDeltaEta_t_cut_proton_left"]->Fill(deltaeta, event_weight );	
						histosTH1F["muonDeltaPhi_t_cut_proton_left"]->Fill(deltaphi, event_weight );
						histosTH1F["muonDeltaY_t_cut_proton_left"]->Fill(deltay, event_weight ); 



					}
					}
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					if( good_proton_right ){
						if((dimuon_mass > 3.0) && (dimuon_mass < 3.2)){

							n_jpsi_proton_right++;
							double t_proton_right = rec_proton_right->t;
							histosTH1F["jpsi_proton_right_xi"]->Fill( xi_proton_right, event_weight );
							histosTH1F["jpsi_proton_right_t"]->Fill( fabs(t_proton_right), event_weight );

							//histosTH1F["jpsi_multiplicity_proton_right"]->Fill(n_dimuons_selected, event_weight );
							histosTH1F["jpsi_mass_proton_right"]->Fill( dimuon_mass, event_weight );
							histosTH1F["jpsi_pt_proton_right"]->Fill( dimuon_pt, event_weight );
							histosTH1F["jpsi_pt2_proton_right"]->Fill( dimuon_pt2, event_weight );
							histosTH1F["jpsi_eta_proton_right"]->Fill( dimuon_eta, event_weight );
							histosTH1F["jpsi_rapidity_proton_right"]->Fill( dimuon_rapidity, event_weight );
							histosTH1F["muonDeltaPt_jpsi_proton_right"]->Fill(deltapt, event_weight );
							histosTH1F["muonDeltaEta_jpsi_proton_right"]->Fill(deltaeta, event_weight );
							histosTH1F["muonDeltaPhi_jpsi_proton_right"]->Fill(deltaphi, event_weight );
							histosTH1F["muonDeltaY_jpsi_proton_right"]->Fill(deltay, event_weight );

							if (fabs(t_proton_right)>t_proton_down_ && fabs(t_proton_right)<t_proton_up_){
								//if(proton_plus_t_range){
								n_jpsi_t_selected_proton_right++;
								if(verbose) cout<< "t proton right side (abs):"<< "["<< fabs(t_proton_right)<< "]"<< endl;
								outstring_right << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;

								histosTH1F["jpsi_dimuon_mass_t_cut_proton_right"]->Fill( dimuon_mass, event_weight );
								histosTH1F["jpsi_dimuon_pt_t_cut_proton_right"]->Fill( dimuon_pt, event_weight );
								histosTH1F["jpsi_dimuon_pt2_t_cut_proton_right"]->Fill( dimuon_pt2, event_weight );
								histosTH1F["jpsi_dimuon_eta_t_cut_proton_right"]->Fill( dimuon_eta, event_weight );
								histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_right"]->Fill( dimuon_rapidity, event_weight );
								histosTH1F["jpsi_proton_right_xi_t_cut"]->Fill( xi_proton_right, event_weight );
								histosTH1F["jpsi_proton_right_t_t_cut"]->Fill( fabs(t_proton_right), event_weight );
								histosTH1F["muonDeltaPt_jpsi_t_cut_proton_right"]->Fill(deltapt, event_weight );
								histosTH1F["muonDeltaEta_jpsi_t_cut_proton_right"]->Fill(deltaeta, event_weight );
								histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_right"]->Fill(deltaphi, event_weight );
								histosTH1F["muonDeltaY_jpsi_t_cut_proton_right"]->Fill(deltay, event_weight );

							} // t right proton
							} // right proton
						} //jpsi mass
						if( good_proton_left ){
							if((dimuon_mass > 3.0) && (dimuon_mass < 3.2)){
								n_jpsi_proton_left++;
								double t_proton_left = rec_proton_left->t;

								histosTH1F["jpsi_proton_left_xi"]->Fill(xi_proton_left, event_weight );
								histosTH1F["jpsi_proton_left_t"]->Fill( fabs(t_proton_left), event_weight );

								histosTH1F["jpsi_mass_proton_left"]->Fill( dimuon_mass, event_weight );
								histosTH1F["jpsi_pt_proton_left"]->Fill( dimuon_pt, event_weight );
								histosTH1F["jpsi_pt2_proton_left"]->Fill( dimuon_pt2, event_weight );
								histosTH1F["jpsi_eta_proton_left"]->Fill( dimuon_eta, event_weight );
								histosTH1F["jpsi_rapidity_proton_left"]->Fill( dimuon_rapidity, event_weight );
								histosTH1F["muonDeltaPt_jpsi_proton_left"]->Fill(deltapt, event_weight );
								histosTH1F["muonDeltaEta_jpsi_proton_left"]->Fill(deltaeta, event_weight );
								histosTH1F["muonDeltaPhi_jpsi_proton_left"]->Fill(deltaphi, event_weight );
								histosTH1F["muonDeltaY_jpsi_proton_left"]->Fill(deltay, event_weight );                       			

								if (fabs(t_proton_left)>t_proton_down_ && fabs(t_proton_left)<t_proton_up_){
									//if(proton_minus_t_range){
									n_jpsi_t_selected_proton_left++;
									outstring_left << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;

									if(verbose) cout<< "t proton left side (abs):"<< "["<< fabs(t_proton_left)<< "]"<< endl;


									histosTH1F["jpsi_dimuon_mass_t_cut_proton_left"]->Fill( dimuon_mass, event_weight );
									histosTH1F["jpsi_dimuon_pt_t_cut_proton_left"]->Fill( dimuon_pt, event_weight );
									histosTH1F["jpsi_dimuon_pt2_t_cut_proton_left"]->Fill( dimuon_pt2, event_weight );
									histosTH1F["jpsi_dimuon_eta_t_cut_proton_left"]->Fill( dimuon_eta, event_weight );
									histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left"]->Fill( dimuon_rapidity, event_weight );
									histosTH1F["jpsi_proton_left_xi_t_cut"]->Fill( xi_proton_left, event_weight );
									histosTH1F["jpsi_proton_left_t_t_cut"]->Fill( fabs(t_proton_left), event_weight );

									histosTH1F["muonDeltaPt_jpsi_t_cut_proton_left"]->Fill(deltapt, event_weight );
									histosTH1F["muonDeltaEta_jpsi_t_cut_proton_left"]->Fill(deltaeta, event_weight );
									histosTH1F["muonDeltaPhi_jpsi_t_cut_proton_left"]->Fill(deltaphi, event_weight );
									histosTH1F["muonDeltaY_jpsi_t_cut_proton_left"]->Fill(deltay, event_weight );


								} // t, proton left
								} // good proton left
							}//jpsi mass


						} //!isMC
						//===========================================================================
						//-------------------
						// After selection 
						//-------------------

						//-------------------
						// Generator-level proton distributions
						//-------------------
						if(isMC){
							if( proton_plus_rp_accept ){
								histosTH1F["xi_proton_plus_selected"]->Fill( xi_proton_plus , event_weight );
								histosTH1F["t_proton_plus_selected"]->Fill( fabs(t_proton_plus) , event_weight ); 
								histosTH2F["proton_plus_xi_vs_t_selected"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
								//dimuons              
								++n_dimuons_rp_selected_plus_side;	
								//histosTH1F["dimuon_multiplicity_rp_accept"]->Fill(n_dimuons_rp_selected_plus_side, event_weight );
								histosTH1F["dimuon_mass_rpplus_accept"]->Fill( dimuon_mass, event_weight );
								histosTH1F["dimuon_pt_rpplus_accept"]->Fill( dimuon_pt, event_weight );
								histosTH1F["dimuon_pt2_rpplus_accept"]->Fill(dimuon_pt2, event_weight );
								histosTH1F["dimuon_eta_rpplus_accept"]->Fill( dimuon_eta, event_weight );
								histosTH1F["dimuon_rapidity_rpplus_accept"]->Fill( dimuon_rapidity, event_weight );
								histosTH1F["muonDeltaPt_rpplus_accept"]->Fill(deltapt, event_weight );
								histosTH1F["muonDeltaEta_rpplus_accept"]->Fill(deltaeta, event_weight );	
								histosTH1F["muonDphi_rpplus_accept"]->Fill(Dphi, event_weight );
								histosTH1F["muonDeltaPhi_rpplus_accept"]->Fill(deltaphi, event_weight );
								histosTH1F["muonDeltaY_rpplus_accept"]->Fill(deltay, event_weight ); 

								//jpsi
								if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){  

									// cout<<"jpsi_mass = "<< jpsi_mass << endl;
									//++n_jpsi_rp_accept_selected;
									++n_jpsi_selected_rp_accept_plus_side;
									//histosTH1F["jpsi_multiplicity_rp_accept"]->Fill(n_dimuons_selected, event_weight );
									histosTH1F["jpsi_mass_rpplus_accept"]->Fill(jpsi_mass, event_weight );
									histosTH1F["jpsi_pt_rpplus_accept"]->Fill( jpsi_pt, event_weight );
									histosTH1F["jpsi_pt2_rpplus_accept"]->Fill( jpsi_pt2, event_weight );
									histosTH1F["jpsi_eta_rpplus_accept"]->Fill( jpsi_eta, event_weight );
									histosTH1F["jpsi_rapidity_rpplus_accept"]->Fill( jpsi_rapidity, event_weight );
									histosTH1F["muonDeltaPt_jpsi_rpplus_accept"]->Fill(deltapt, event_weight );
									histosTH1F["muonDeltaEta_jpsi_rpplus_accept"]->Fill(deltaeta, event_weight );	
									histosTH1F["muonDeltaPhi_jpsi_rpplus_accept"]->Fill(deltaphi, event_weight );
									histosTH1F["muonDeltaY_jpsi_rpplus_accept"]->Fill(deltay, event_weight );  

									histosTH1F["jpsi_t_plus_rpplus_accept"]->Fill( fabs(t_proton_plus), event_weight );
									histosTH1F["jpsi_xi_plus_rpplus_accept"]->Fill( xi_proton_plus, event_weight );
									histosTH1F["muonDphi_jpsi_rpplus_accept"]->Fill(Dphi, event_weight );
								}

								if(proton_plus_t_range){
									histosTH1F["xi_proton_t_range_plus_selected"]->Fill( xi_proton_plus , event_weight );
									//----------------------------------------------
									//dimuons              
									++n_dimuons_rp_selected_tsel_plus_side;	
									//histosTH1F["dimuon_multiplicity_rp_accept"]->Fill(n_dimuons_rp_selected_plus_side, event_weight );
									histosTH1F["dimuon_mass_rpplus_accept_tsel"]->Fill( dimuon_mass, event_weight );
									histosTH1F["dimuon_pt_rpplus_accept_tsel"]->Fill( dimuon_pt, event_weight );
									histosTH1F["dimuon_pt2_rpplus_accept_tsel"]->Fill(dimuon_pt2, event_weight );
									histosTH1F["dimuon_eta_rpplus_accept_tsel"]->Fill( dimuon_eta, event_weight );
									histosTH1F["dimuon_rapidity_rpplus_accept_tsel"]->Fill( dimuon_rapidity, event_weight );
									histosTH1F["muonDeltaPt_rpplus_accept_tsel"]->Fill(deltapt, event_weight );
									histosTH1F["muonDeltaEta_rpplus_accept_tsel"]->Fill(deltaeta, event_weight );	
									histosTH1F["muonDphi_rpplus_accept_tsel"]->Fill(Dphi, event_weight );
									histosTH1F["muonDeltaPhi_rpplus_accept_tsel"]->Fill(deltaphi, event_weight );
									histosTH1F["muonDeltaY_rpplus_accept_tsel"]->Fill(deltay, event_weight ); 

									//jpsi
									if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){  
										//cout<< " rp accept after selection |t| && mass range: "<< " |t| proton plus "<< fabs(t_proton_plus)<<endl;
										//cout<<"jpsi_mass = "<< jpsi_mass << endl;
										++n_jpsi_selected_rp_accept_tsel_plus_side;
										//histosTH1F["jpsi_multiplicity_rp_accept"]->Fill(n_dimuons_selected, event_weight );
										histosTH1F["jpsi_mass_rpplus_accept_tsel"]->Fill(jpsi_mass, event_weight );
										histosTH1F["jpsi_pt_rpplus_accept_tsel"]->Fill( jpsi_pt, event_weight );
										histosTH1F["jpsi_pt2_rpplus_accept_tsel"]->Fill( jpsi_pt2, event_weight );
										histosTH1F["jpsi_eta_rpplus_accept_tsel"]->Fill( jpsi_eta, event_weight );
										histosTH1F["jpsi_rapidity_rpplus_accept_tsel"]->Fill( jpsi_rapidity, event_weight );
										histosTH1F["muonDeltaPt_jpsi_rpplus_accept_tsel"]->Fill(deltapt, event_weight );
										histosTH1F["muonDeltaEta_jpsi_rpplus_accept_tsel"]->Fill(deltaeta, event_weight );	
										histosTH1F["muonDeltaPhi_jpsi_rpplus_accept_tsel"]->Fill(deltaphi, event_weight );
										histosTH1F["muonDeltaY_jpsi_rpplus_accept_tsel"]->Fill(deltay, event_weight );  
										//cout<<"jpsi_t_plus_rpplus_accept_tsel"<<fabs(t_proton_plus)<<endl;
										histosTH1F["jpsi_t_plus_rpplus_accept_tsel"]->Fill( fabs(t_proton_plus), event_weight );
										histosTH1F["jpsi_xi_plus_rpplus_accept_tsel"]->Fill( xi_proton_plus, event_weight );
										histosTH1F["muonDphi_jpsi_rpplus_accept_tsel"]->Fill(Dphi, event_weight );
									}

									if(proton_plus_xi_range){
										histosTH1F["t_proton_xi_range_plus_selected"]->Fill( fabs(t_proton_plus) , event_weight );}
								}
							}
							if( proton_minus_rp_accept ){
								histosTH1F["xi_proton_minus_selected"]->Fill( xi_proton_minus , event_weight );
								histosTH1F["t_proton_minus_selected"]->Fill( fabs(t_proton_minus) , event_weight ); 
								histosTH2F["proton_minus_xi_vs_t_selected"]->Fill( fabs(t_proton_minus) , xi_proton_minus , event_weight );
								//dimuons              
								++n_dimuons_rp_selected_minus_side;	

								histosTH1F["dimuon_mass_rpminus_accept"]->Fill( dimuon_mass, event_weight );
								histosTH1F["dimuon_pt_rpminus_accept"]->Fill( dimuon_pt, event_weight );
								histosTH1F["dimuon_pt2_rpminus_accept"]->Fill(dimuon_pt2, event_weight );
								histosTH1F["dimuon_eta_rpminus_accept"]->Fill( dimuon_eta, event_weight );
								histosTH1F["dimuon_rapidity_rpminus_accept"]->Fill( dimuon_rapidity, event_weight );
								histosTH1F["muonDeltaPt_rpminus_accept"]->Fill(deltapt, event_weight );
								histosTH1F["muonDeltaEta_rpminus_accept"]->Fill(deltaeta, event_weight );	
								histosTH1F["muonDphi_rpminus_accept"]->Fill(Dphi, event_weight );
								histosTH1F["muonDeltaPhi_rpminus_accept"]->Fill(deltaphi, event_weight );
								histosTH1F["muonDeltaY_rpminus_accept"]->Fill(deltay, event_weight ); 

								//jpsi
								if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){  

									//cout<<"jpsi_mass = "<< jpsi_mass << endl;
									++n_jpsi_selected_rp_accept_minus_side;

									histosTH1F["jpsi_mass_rpminus_accept"]->Fill(jpsi_mass, event_weight );
									histosTH1F["jpsi_pt_rpminus_accept"]->Fill( jpsi_pt, event_weight );
									histosTH1F["jpsi_pt2_rpminus_accept"]->Fill( jpsi_pt2, event_weight );
									histosTH1F["jpsi_eta_rpminus_accept"]->Fill( jpsi_eta, event_weight );
									histosTH1F["jpsi_rapidity_rpminus_accept"]->Fill( jpsi_rapidity, event_weight );
									histosTH1F["muonDeltaPt_jpsi_rpminus_accept"]->Fill(deltapt, event_weight );
									histosTH1F["muonDeltaEta_jpsi_rpminus_accept"]->Fill(deltaeta, event_weight );	
									histosTH1F["muonDeltaPhi_jpsi_rpminus_accept"]->Fill(deltaphi, event_weight );
									histosTH1F["muonDeltaY_jpsi_rpminus_accept"]->Fill(deltay, event_weight );  

									histosTH1F["jpsi_t_minus_rpminus_accept"]->Fill( fabs(t_proton_minus), event_weight );
									histosTH1F["jpsi_xi_minus_rpminus_accept"]->Fill( xi_proton_minus, event_weight );
									histosTH1F["muonDphi_jpsi_rpminus_accept"]->Fill(Dphi, event_weight );
								}




								if(proton_minus_t_range){
									histosTH1F["xi_proton_t_range_minus_selected"]->Fill( xi_proton_minus , event_weight );
									//----------------------------------------------
									//dimuons              
									++n_dimuons_rp_selected_tsel_minus_side;	

									histosTH1F["dimuon_mass_rpminus_accept_tsel"]->Fill( dimuon_mass, event_weight );
									histosTH1F["dimuon_pt_rpminus_accept_tsel"]->Fill( dimuon_pt, event_weight );
									histosTH1F["dimuon_pt2_rpminus_accept_tsel"]->Fill(dimuon_pt2, event_weight );
									histosTH1F["dimuon_eta_rpminus_accept_tsel"]->Fill( dimuon_eta, event_weight );
									histosTH1F["dimuon_rapidity_rpminus_accept_tsel"]->Fill( dimuon_rapidity, event_weight );
									histosTH1F["muonDeltaPt_rpminus_accept_tsel"]->Fill(deltapt, event_weight );
									histosTH1F["muonDeltaEta_rpminus_accept_tsel"]->Fill(deltaeta, event_weight );	
									histosTH1F["muonDphi_rpminus_accept_tsel"]->Fill(Dphi, event_weight );
									histosTH1F["muonDeltaPhi_rpminus_accept_tsel"]->Fill(deltaphi, event_weight );
									histosTH1F["muonDeltaY_rpminus_accept_tsel"]->Fill(deltay, event_weight ); 

									//jpsi
									if(((dimuon_mass > 3.0) && (dimuon_mass < 3.2))){  
										//cout<< " rp accept after selection |t| && mass range: "<< " |t| proton minus "<< fabs(t_proton_minus)<<endl;
										//cout<<"jpsi_mass = "<< jpsi_mass << endl;
										++n_jpsi_selected_rp_accept_tsel_minus_side;

										histosTH1F["jpsi_mass_rpminus_accept_tsel"]->Fill(jpsi_mass, event_weight );
										histosTH1F["jpsi_pt_rpminus_accept_tsel"]->Fill( jpsi_pt, event_weight );
										histosTH1F["jpsi_pt2_rpminus_accept_tsel"]->Fill( jpsi_pt2, event_weight );
										histosTH1F["jpsi_eta_rpminus_accept_tsel"]->Fill( jpsi_eta, event_weight );
										histosTH1F["jpsi_rapidity_rpminus_accept_tsel"]->Fill( jpsi_rapidity, event_weight );
										histosTH1F["muonDeltaPt_jpsi_rpminus_accept_tsel"]->Fill(deltapt, event_weight );
										histosTH1F["muonDeltaEta_jpsi_rpminus_accept_tsel"]->Fill(deltaeta, event_weight );	
										histosTH1F["muonDeltaPhi_jpsi_rpminus_accept_tsel"]->Fill(deltaphi, event_weight );
										histosTH1F["muonDeltaY_jpsi_rpminus_accept_tsel"]->Fill(deltay, event_weight );  
										//cout<<"jpsi_t_plus_rpminus_accept_tsel"<<fabs(t_proton_minus)<<endl;
										histosTH1F["jpsi_t_minus_rpminus_accept_tsel"]->Fill( fabs(t_proton_minus), event_weight );
										histosTH1F["jpsi_xi_minus_rpminus_accept_tsel"]->Fill( xi_proton_minus, event_weight );
										histosTH1F["muonDphi_jpsi_rpminus_accept_tsel"]->Fill(Dphi, event_weight );
									}



									if(proton_minus_xi_range){
										histosTH1F["t_proton_xi_range_minus_selected"]->Fill( fabs(t_proton_minus) , event_weight );}
								}
							}
						}

						//-------------------
						// Detector-level distributions
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
						/*if( pfJet_coll->size() > 0 ){
						  MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
						  histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), n_t2_tracks_selected, event_weight );
						  }*/

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

						///////////////////////////////////////////////////////////////////////////////

					} // End of loop over events in a file

					cout<<"Total of evts="<< nev << endl << *itfiles << endl;
					cout<<"n_vertices_selected ="<< n_vertices_selected << endl;
					cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
					cout<<"n_tracks_selected="<< n_tracks_selected<<endl;
					cout<<"n_muon_selected="<< n_muon_selected<<endl;
					cout<<"n_dimuons_selected="<< n_dimuons_selected<<endl;
					cout<<"n_dimuons_proton_left="<<n_dimuons_proton_left<<endl;
					cout<<"n_dimuons_proton_right="<<n_dimuons_proton_right<<endl;
					cout<<"n_dimuons_t_selected_proton_left="<< n_dimuons_t_selected_proton_left<<endl;
					//cout<<"n_dimuons_t_at_xi_selected_proton_left="<< n_dimuons_t_at_xi_selected_proton_left<<endl;
					cout<<"n_dimuons_t_selected_proton_right="<< n_dimuons_t_selected_proton_right<<endl;
					//cout<<"n_dimuons_t_at_xi_selected_proton_right="<< n_dimuons_t_at_xi_selected_proton_right<<endl;	
					cout<<"n_jpsi_selected="<< n_jpsi_selected<<endl;
					cout<<"n_jpsi_proton_left="<<n_jpsi_proton_left<<endl;
					cout<<"n_jpsi_proton_right="<<n_jpsi_proton_right<<endl;
					cout<<"n_jpsi_t_selected_proton_left="<< n_jpsi_t_selected_proton_left<<endl;
					//cout<<"n_jpsi_t_at_xi_selected_proton_left="<< n_jpsi_t_at_xi_selected_proton_left<<endl;
					cout<<"n_jpsi_t_selected_proton_right="<< n_jpsi_t_selected_proton_right<<endl;
					//cout<<"n_jpsi_t_at_xi_selected_proton_right="<< n_jpsi_t_at_xi_selected_proton_right<<endl;	
					cout<<"n_dimuons_rp_accept_plus_side="<< n_dimuons_rp_selected_plus_side <<endl;
					cout<<"n_dimuons_rp_accept_tsel_plus_side="<< n_dimuons_rp_selected_tsel_plus_side <<endl;
					cout<<"n_jpsi_selected_rp_accept_plus_side ="<< n_jpsi_selected_rp_accept_plus_side<<endl;
					cout<<"n_jpsi_selected_rp_accept_tsel_plus_side ="<< n_jpsi_selected_rp_accept_tsel_plus_side<<endl;

					cout<<"n_dimuons_rp_accept_minus_side="<< n_dimuons_rp_selected_minus_side <<endl;
					cout<<"n_dimuons_rp_accept_tsel_minus_side="<< n_dimuons_rp_selected_tsel_minus_side <<endl;
					cout<<"n_jpsi_selected_rp_accept_minus_side ="<< n_jpsi_selected_rp_accept_minus_side<<endl;
					cout<<"n_jpsi_selected_rp_accept_tsel_minus_side ="<< n_jpsi_selected_rp_accept_tsel_minus_side<<endl;
					//  }	
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
