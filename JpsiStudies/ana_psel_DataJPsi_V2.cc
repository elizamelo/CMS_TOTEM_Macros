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

void ana_psel_DataJPsi_V2(vector<string> const& fileNames, string const& outputFileName = "output.root", const Double_t t_proton_down_=0.0, const Double_t t_proton_up_=1.0, const Double_t xi_proton_down_=0.0,const Double_t xi_proton_up_=0.2,const Int_t Bin_mass=100 ,const Int_t nevt_max = -1){

	//bool isMC  = false;
	bool verbose = false;
	string treeName = "cms_totem";
	double etaMaxThreshold = 2.0;


	bool selectBunchCrossing = false;
	bool selectVertex = true;
	bool selectMuons = true;
        bool selectTrack = true;
	bool selectEtaMax = false;
	bool selectEtaMin = false;
	bool selectSingleArmRecProton = true;
        //bool selectSingleArmRecProton = false;
	bool selectDoubleArmRecProton = false;
	//bool selectDoubleArmRecProton = true;
        bool selectElastic = false;
	bool selectNonElastic = false;
        //bool selectNonElastic = true;

	vector<int> bunchCrossingList;
	bunchCrossingList.push_back(27);
	bunchCrossingList.push_back(649);
	bunchCrossingList.push_back(2991);

	// bool selectHLTORorAND = true; // OR=true - AND=false
	vector<string> selectHLTPathNames;
	/*selectHLTPathNames.push_back("HLT_L1DoubleJet20part1_v1");
	  selectHLTPathNames.push_back("HLT_L1DoubleJet20part2_v1");*/
	selectHLTPathNames.push_back("HLT_L1_DoubleMu0");
	//selectHLTPathNames.push_back("HLT_L1Tech53_MB_2_v1");
	//selectHLTPathNames.push_back("HLT_L1Tech53_MB_3_v1");
	//selectHLTPathNames.push_back("HLT_ZeroBias_v7");
	//selectHLTPathNames.push_back("HLT_RomanPots_Tech52_v1");

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
	selections.push_back("Muons");
	selections.push_back("EtaMax");
	selections.push_back("EtaMin");
	selections.push_back("SingleArmRP");
	selections.push_back("DoubleArmRP");
	selections.push_back("Elastic");
	selections.push_back("NonElastic");
	int nBinsEventSelection = selections.size();
	histosTH1F["EventSelection"] = new TH1F("EventSelection","EventSelection",nBinsEventSelection,0,nBinsEventSelection);

        //Float_t tbins[12] = {0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.24, 0.30, 0.36, 0.45, 0.65, 1.};
	//Float_t xi_bins[30] = {0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15,0.20};
	
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


        //t_cut_xi_cut_proton_right
        /*histosTH1F["dimuon_mass_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_mass_t_cut_xi_cut_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["dimuon_pt2_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_pt2_t_cut_xi_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
        histosTH1F["dimuon_pt_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_pt_t_cut_xi_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["dimuon_eta_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_eta_t_cut_xi_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_rapidity_t_cut_xi_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_multiplicity_t_cut_xi_cut_proton_right", "n dimuons" , 100 , 0 , 100);
        histosTH1F["proton_right_xi_t_cut_xi_cut"] = new TH1F("proton_right_xi_selected_t_and_xi", "#xi", 100 ,0.,1.);
        histosTH1F["proton_right_t_t_cut_xi_cut"] = new TH1F("proton_right_t_selected_t_and_xi", "-t" , 100, 0.,5.);


        //Muon1 info
        //histosTH1F["muon1_mass"] = new TH1F("muon1_mass", "mass_{#mu 1}" , Bin_mass , 0. , 5.);
        histosTH1F["muon1_pt_t_cut_xi_cut_proton_right"] = new TH1F("muon1_pt_t_cut_xi_cut_proton_right", "p_{T}(mu1)" , 100 , 0. , 100.);
        histosTH1F["muon1_eta_t_cut_xi_cut_proton_right"] = new TH1F("muon1_eta_t_cut_xi_cut_proton_right", "#eta(mu1)" , 100 , -5.2 , 5.2);
        histosTH1F["muon1_phi_t_cut_xi_cut_proton_right"] = new TH1F("muon1_phi_t_cut_xi_cut_proton_right", "#phi(muon1)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muon1_rapidity_t_cut_xi_cut_proton_right"] = new TH1F("muon1_rapidity_t_cut_xi_cut_proton_right", "y(mu1)" , 100 , -15. , 15.);
        //Muon1 info
        //histosTH1F["muon2_mass"] = new TH1F("muon2_mass", "mass_{#mu 2}" , Bin_mass , 0. , 5.);
        histosTH1F["muon2_pt_t_cut_xi_cut_proton_right"] = new TH1F("muon2_pt_t_cut_xi_cut_proton_right", "p_{T}(mu2)" , 100 , 0. , 100.);
        histosTH1F["muon2_eta_t_cut_xi_cut_proton_right"] = new TH1F("muon2_eta_t_cut_xi_cut_proton_right", "#eta(mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["muon2_phi_t_cut_xi_cut_proton_right"] = new TH1F("muon2_phi_t_cut_xi_cut_proton_right", "#phi(muon2)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muon2_rapidity_t_cut_xi_cut_proton_right"] = new TH1F("muon2_rapidity_t_cut_xi_cut_proton_right", "y(mu2)" , 100 , -15. , 15.);
        //Deltas info
        histosTH1F["muonDeltaPt_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaPt_t_cut_xi_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaEta_t_cut_xi_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["muonDeltaPhi_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaPhi_t_cut_xi_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 ,- 2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaY_t_cut_xi_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);*/

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

        //t_cut_xi_cut_proton_left
       /* histosTH1F["dimuon_mass_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_mass_t_cut_xi_cut_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["dimuon_pt_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_pt_t_cut_xi_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["dimuon_pt2_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_pt2_t_cut_xi_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 1000.);

        histosTH1F["dimuon_eta_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_eta_t_cut_xi_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_rapidity_t_cut_xi_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_multiplicity_t_cut_xi_cut_proton_left", "n dimuons" , 100 , 0 , 100);
        histosTH1F["proton_left_xi_t_cut_xi_cut"] = new TH1F("proton_left_xi_selected_t_and_xi ", "#xi" ,  100, 0.,1.);
        histosTH1F["proton_left_t_t_cut_xi_cut"] = new TH1F("proton_left_t_selected_t_and_xi", "-t" , 100, 0., 5.);

         //Muon1 info
        //histosTH1F["muon1_mass"] = new TH1F("muon1_mass", "mass_{#mu 1}" , Bin_mass , 0. , 5.);
        histosTH1F["muon1_pt_t_cut_xi_cut_proton_left"] = new TH1F("muon1_pt_t_cut_xi_cut_proton_left", "p_{T}(mu1)" , 100 , 0. , 100.);
        histosTH1F["muon1_eta_t_cut_xi_cut_proton_left"] = new TH1F("muon1_eta_t_cut_xi_cut_proton_left", "#eta(mu1)" , 100 , -5.2 , 5.2);
        histosTH1F["muon1_phi_t_cut_xi_cut_proton_left"] = new TH1F("muon1_phi_t_cut_xi_cut_proton_left", "#phi(muon1)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muon1_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("muon1_rapidity_t_cut_xi_cut_proton_left", "y(mu1)" , 100 , -15. , 15.);
        //Muon1 info
        //histosTH1F["muon2_mass"] = new TH1F("muon2_mass", "mass_{#mu 2}" , Bin_mass , 0. , 5.);
        histosTH1F["muon2_pt_t_cut_xi_cut_proton_left"] = new TH1F("muon2_pt_t_cut_xi_cut_proton_left", "p_{T}(mu2)" , 100 , 0. , 100.);
        histosTH1F["muon2_eta_t_cut_xi_cut_proton_left"] = new TH1F("muon2_eta_t_cut_xi_cut_proton_left", "#eta(mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["muon2_phi_t_cut_xi_cut_proton_left"] = new TH1F("muon2_phi_t_cut_xi_cut_proton_left", "#phi(muon2)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muon2_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("muon2_rapidity_t_cut_xi_cut_proton_left", "y(mu2)" , 100 , -15. , 15.);
        //Deltas info
        histosTH1F["muonDeltaPt_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPt_t_cut_xi_cut_proton_left", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaEta_t_cut_xi_cut_proton_left", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["muonDeltaPhi_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaPhi_t_cut_xi_cut_proton_left", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_t_cut_xi_cut_proton_left"] = new TH1F("muonDeltaY_t_cut_xi_cut_proton_left", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
        */




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

       //t_cut_xi_cut_proton_right
       /* histosTH1F["jpsi_dimuon_mass_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_mass_t_cut_xi_cut_proton_right", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["jpsi_dimuon_pt2_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_pt2_t_cut_xi_cut_proton_right", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 2000.);
        histosTH1F["jpsi_dimuon_pt_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_pt_t_cut_xi_cut_proton_right", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
        histosTH1F["jpsi_dimuon_eta_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_eta_t_cut_xi_cut_proton_right", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["jpsi_dimuon_rapidity_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_rapidity_t_cut_xi_cut_proton_right", "y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_right"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_right", "n jpsi_dimuons" , 100 , 0 , 20);
        histosTH1F["jpsi_proton_right_xi_t_cut_xi_cut"] = new TH1F("jpsi_proton_right_xi_selected_t_and_xi", "#xi" ,100, 0.,1. );
        histosTH1F["jpsi_proton_right_t_t_cut_xi_cut"] = new TH1F("jpsi_proton_right_t_selected_t_and_xi", "-t" , 100, 0.,5.);*/

       //Deltas info
       /* histosTH1F["muonDeltaPt_jpsi_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaPt_jpsi_t_cut_xi_cut_proton_right", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_jpsi_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaEta_jpsi_t_cut_xi_cut_proton_right", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["muonDeltaPhi_jpsi_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaPhi_jpsi_t_cut_xi_cut_proton_right", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_jpsi_t_cut_xi_cut_proton_right"] = new TH1F("muonDeltaY_jpsi_t_cut_xi_cut_proton_right", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);*/



	//_t_cut_proton_left
	histosTH1F["jpsi_dimuon_mass_t_cut_proton_left"] = new TH1F("jpsi_dimuon_mass_t_cut_proton_left", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["jpsi_dimuon_pt_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 100 , 0. , 20.);
        histosTH1F["jpsi_dimuon_pt2_t_cut_proton_left"] = new TH1F("jpsi_dimuon_pt2_t_cut_proton_left", "p_{T}^{2}(mu1,mu2)" , 100 , 0. , 2000.);
	histosTH1F["jpsi_dimuon_eta_t_cut_proton_left"] = new TH1F("jpsi_dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["jpsi_dimuon_rapidity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_left"] = new TH1F("jpsi_dimuon_multiplicity_t_cut_proton_left", "n jpsi_dimuons" , 100 , 0 , 20);

        histosTH1F["jpsi_proton_left_xi_t_cut"] = new TH1F("jpsi_proton_left_xi_t_selected", "#xi" ,30, 0.,1. );
        histosTH1F["jpsi_proton_left_t_t_cut"] = new TH1F("jpsi_proton_left_t_t_selected", "-t" , 30, 0., 5.);

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

	/*histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
	  histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);*/

	histosTH1F["t2_track_chi2Prob_zplus"] = new TH1F("t2_track_chi2Prob_zplus", "#chi^{2}" , 100 , 0. , 1.);
	histosTH1F["t2_track_entryX_zplus"] = new TH1F("t2_track_entryX_zplus", "x_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_entryY_zplus"] = new TH1F("t2_track_entryY_zplus", "y_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_multiplicity_zplus"] = new TH1F("t2_track_multiplicity_zplus", "n tracks" , 100 , 0 , 100);
	histosTH1F["t2_track_chi2Prob_zminus"] = new TH1F("t2_track_chi2Prob_zminus", "#chi^{2}" , 100 , 0. , 1.);
	histosTH1F["t2_track_entryX_zminus"] = new TH1F("t2_track_entryX_zminus", "x_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_entryY_zminus"] = new TH1F("t2_track_entryY_zminus", "y_{trk}" , 160 , -160. , 160.);
	histosTH1F["t2_track_multiplicity_zminus"] = new TH1F("t2_track_multiplicity_zminus", "n tracks" , 100 , 0 , 100);

	histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 100 , 0. , 1.);
        histosTH1F["proton_right_xi_total"] = new TH1F("proton_right_xi_all", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",100,-5.,0.);
	histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 100 , 0. , 100.);
	histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 100 , 0. , 1.);
        histosTH1F["proton_left_xi_total"] = new TH1F("proton_left_xi_all", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 100 , 0. , 100.);

        histosTH1F["proton_left_thx"] = new TH1F("thx_proton_left", "thx_proton_left" , 100, -5e-4, 5e-4);
        histosTH1F["proton_left_thy"] = new TH1F("thx_proton_left", "thx_proton_left" , 100, -5e-4, 5e-4);
	histosTH1F["tx_proton_left"] = new TH1F("tx_proton_left", "tx_proton_left" , 100, -10., 10.);
	histosTH1F["ty_proton_left "] = new TH1F("ty_proton_left", "ty_proton_left" , 100, -50., 50.);
	histosTH1F["x0_proton_left"] = new TH1F("x0_proton_left", "x0_proton_left" , 100, -10., 10.);
        histosTH1F["y0_proton_left "] = new TH1F("y0_proton_left", "y0_proton_left" , 100, -50., 50.);

        histosTH1F["proton_right_thx"] = new TH1F("thx_proton_right", "thx_proton_right" , 100, -5e-4, 5e-4);
        histosTH1F["proton_right_thy"] = new TH1F("thx_proton_right", "thx_proton_right" , 100, -5e-4, 5e-4);
	histosTH1F["tx_proton_right"] = new TH1F("tx_proton_right", "tx_proton_right" , 100, -10., 10.);
	histosTH1F["ty_proton_right"] = new TH1F("ty_proton_right", "ty_proton_right" , 100, -50., 50.);
	histosTH1F["x0_proton_right"] = new TH1F("x0_proton_right", "x0_proton_right" , 100, -10., 10.);
        histosTH1F["y0_proton_right"] = new TH1F("y0_proton_right", "y0_proton_right" , 100, -50., 50.);

        //  histosTH1F["rp_track_posx_020"] = new TH1F("rp_track_posx_020", "x(RP track)" , 200, -10., 10.);
        //  histosTH1F["rp_track_posy_020"] = new TH1F("rp_track_posy_020", "y(RP track)" , 500, -50., 50.);


	histosTH1F["proton_pair_right_xi"] = new TH1F("proton_pair_right_xi", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_pair_right_logXi"] = new TH1F("proton_pair_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_right_t"] = new TH1F("proton_pair_right_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_left_xi"] = new TH1F("proton_pair_left_xi", "#xi" , 100 , -1. , 1.);
	histosTH1F["proton_pair_left_logXi"] = new TH1F("proton_pair_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_pair_left_t"] = new TH1F("proton_pair_left_t", "-t" , 100 , 0. , 5.);
	histosTH1F["proton_pair_chi2"] = new TH1F("proton_pair_chi2", "#chi^{2}" , 100 , 0. , 100.);

	histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 100 , -1. , 1.);
	histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 100 , -1. , 1.);

	map<string,TH2F*> histosTH2F;  
	histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
	histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
	histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

	histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 100, -5., 0., 100, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 100, -5., 0., 100, -5., 0.);
	histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 100, -5., 0., 100, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 100, -5., 0., 100, -5., 0.);
	histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 100, 0., 5., 100, -5., 0.);
	histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 100, 0., 5., 100, -5., 0.);

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	//===================

	vector<TString>* vfiles = new vector<TString>; 
	for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );

	//Declaration of tree and its branches variables
	//TTree* tree = new TTree(treeName.c_str(),"");
	TTree* tree = NULL;
	MyEvtId*           evtId        = NULL;
	MyL1TrigOld*       l1Trig       = NULL;  
	MyHLTrig*          hltTrig      = NULL;
	//vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	//vector<MyPFJet>*   pfJet_coll   = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	//vector<MyFSCHit>*  fscHits_coll = NULL;
	//vector<MyFSCDigi>* fscDigis_coll = NULL;
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
        TString outtxt_right = outputFileName;
	TString outtxt_left = outputFileName;

	outtxt_right.ReplaceAll("root","txt_right");
	outtxt_left.ReplaceAll("root","txt_left");  

	//const string left_ = "left_"+ outtxt_left ;
	//std::string right_ = "right_"+ outtxt_right ;
	ofstream outstring_left(outtxt_left); 
	ofstream outstring_right(outtxt_right);

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
		tree->SetBranchAddress("cmsEvtUA",&evtId);
		tree->SetBranchAddress("cmsTrigUA",&l1Trig);
		tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
		tree->SetBranchAddress("cmsTracksUA",&track_coll); 
		tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
		//tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
		//tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
		tree->SetBranchAddress("cmsMuonsUA",&muon_coll);
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

		//if(isMC) tree->SetBranchAddress("genPart",&genPart);

		/*//Getting number of events
		  int nev = int(tree->GetEntriesFast());
		  nevt_tot += nev;
		  cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/


                        int n_vertices_selected = 0;
                        int n_select_Vertex_After_vtx_cut =0;
                        int n_tracks_selected = 0;
                        int n_muon_selected = 0;
                        int n_dimuons_selected = 0;
                        int n_dimuons_t_selected_proton_left = 0;
                        int n_dimuons_t_at_xi_selected_proton_left = 0;
                        int n_dimuons_t_selected_proton_right = 0;
                        int n_dimuons_t_at_xi_selected_proton_right = 0;
                        int n_jpsi_selected = 0;
                        int n_jpsi_t_selected_proton_left = 0;
                        int n_jpsi_t_at_xi_selected_proton_left = 0;
                        int n_jpsi_t_selected_proton_right = 0;
                        int n_jpsi_t_at_xi_selected_proton_right = 0;
                        int n_jpsi_proton_left = 0;
                        int n_jpsi_proton_right=0;
                        int n_dimuons_proton_left = 0;
                        int n_dimuons_proton_right=0;

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




			//-------------------------------------------------------------------------------------------------
			//filling pt distribution for the generated particles
			//ie those from pythia generator, without reconstruction
			/*if(isMC){
			  for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
			  pt_gen->Fill(p->Pt());
			  }*/

			//-------------------------------------------------------------------------------------------------
			// Vertices
/*			int n_vertices_selected = 0;
			int n_select_Vertex_After_vtx_cut =0; 
			int n_tracks_selected = 0;
			int n_muon_selected = 0;
			int n_dimuons_selected = 0;
        		int n_dimuons_t_selected_proton_left = 0;
        		int n_dimuons_t_at_xi_selected_proton_left = 0;
        		int n_dimuons_t_selected_proton_right = 0;
        		int n_dimuons_t_at_xi_selected_proton_right = 0;
*/
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
                        double deltaphi = 0.; double deltaeta = 0.; double deltapt = 0.; double deltay = 0.; //double Dphi = 0.;
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
                                        double Dphi = std::fabs(deltaPhi(phimu1 ,phimu2));  
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
				       // int partType = it_pfcand->particleId;
				double eta = it_pfcand->Eta();
				double energy = it_pfcand->Energy();
				double pz = it_pfcand->Pz();

				// Apply thresholds
				//if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;

				//if( eta >= etaEdgeLow && eta <= etaEdgeHigh ) particles_sorted.push_back(*it_pfcand);
				//pfEPlusPz  += energy + pz;
				//pfEMinusPz += energy - pz;
			}

			// Find eta_min & eta_max
			/*double pfEtaMin = -999.;
			double pfEtaMax = 999.;
			if( particles_sorted.size() > 0 ){
				std::stable_sort(particles_sorted.begin(), particles_sorted.end(), sortByEta);
				pfEtaMin = particles_sorted.at(0).Eta();
				pfEtaMax = particles_sorted.at(particles_sorted.size() - 1).Eta();
			}
			histosTH1F["pf_etaMax"]->Fill( pfEtaMax, event_weight ); 
			histosTH1F["pf_etaMin"]->Fill( pfEtaMin, event_weight );

			//if( selectEtaMax && (pfEtaMax > 3.0) ) continue;
			if( selectEtaMax && (pfEtaMax > etaMaxThreshold) ) continue;
			histosTH1F["EventSelection"]->Fill( "EtaMax", event_weight );
                         */
			//if( selectEtaMin && (pfEtaMin < -3.0) ) continue;
			//if( selectEtaMin && (pfEtaMin < -etaMaxThreshold) ) continue;
			//histosTH1F["EventSelection"]->Fill( "EtaMin", event_weight );

			bool proton_right_valid = rec_proton_right->valid;
			bool proton_left_valid = rec_proton_left->valid;
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

			/*double pfDeltaEta = pfEtaMax - pfEtaMin;
			histosTH1F["pf_deltaEta"]->Fill( pfDeltaEta, event_weight ); 
                        cout<<"selectSingleArmRecProton: "<<selectSingleArmRecProton<<endl;
			double pfXiPlusReco = xiCorrFactor*pfEPlusPz/8000.;
			double pfXiMinusReco = xiCorrFactor*pfEMinusPz/8000.;
			histosTH1F["pf_EPlusPz"]->Fill( pfEPlusPz, event_weight );
			histosTH1F["pf_EMinusPz"]->Fill( pfEMinusPz, event_weight );
			histosTH1F["pf_xiPlus"]->Fill( pfXiPlusReco, event_weight );
			histosTH1F["pf_xiMinus"]->Fill( pfXiMinusReco, event_weight );
			histosTH1F["pf_logXiPlus"]->Fill( log10(pfXiPlusReco), event_weight );
			histosTH1F["pf_logXiMinus"]->Fill( log10(pfXiMinusReco), event_weight );
                        */
			/*for(vector<MyFSCHit>::iterator it_hit = fscHits_coll->begin() ; it_hit != fscHits_coll->end() ; ++it_hit){
			  histosTH1F["fscHit_energy"]->Fill( it_hit->energy, event_weight );
			  histosTH1F["fscHit_time"]->Fill( it_hit->time, event_weight );
			  }*/

			vector<double> const& t2_trk_entryX = t2_event->TrkEntryX;
			vector<double> const& t2_trk_entryY = t2_event->TrkEntryY;
			vector<double> const& t2_trk_entryZ =  t2_event->TrkEntryZ;
			vector<double> const& t2_trk_chiProb =  t2_event->TrkChiProb;

			int n_t2_tracks_selected = 0;
			int n_t2_tracks_selected_zplus = 0;
			int n_t2_tracks_selected_zminus = 0;
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
			histosTH1F["t2_track_multiplicity_zplus"]->Fill( n_t2_tracks_selected_zplus, event_weight );
			histosTH1F["t2_track_multiplicity_zminus"]->Fill( n_t2_tracks_selected_zminus, event_weight );
			histosTH2F["t2_track_multiplicity_vs_track_multiplicity"]->Fill( n_tracks_selected, n_t2_tracks_selected, event_weight );



			// RP protons 
                        //double thx, thy, phi, t, tx, ty, xi, x0, y0, chi2, chindf;
			//bool proton_right_valid = rec_proton_right->valid;
			double chi2_proton_right = rec_proton_right->chi2;
                        double thx_proton_right = rec_proton_right->thx;
                        double thy_proton_right = rec_proton_right->thy;
                        double tx_proton_right = rec_proton_right->tx;
                        double ty_proton_right = rec_proton_right->ty;
                        double x0_proton_right = rec_proton_right->x0;
                        double y0_proton_right = rec_proton_right->y0;         
          
			//double chindf_proton_right = rec_proton_right->chindf;
			double xi_proton_right = rec_proton_right->xi;
			//bool good_proton_right = proton_right_valid && (chi2_proton_right/chindf_proton_right > 1);
			bool good_proton_right = proton_right_valid && (xi_proton_right < 0.);
			if( good_proton_right ){
				//double xi_proton_right = rec_proton_right->xi;
				double t_proton_right = rec_proton_right->t;
				histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight );
				histosTH1F["proton_right_xi"]->Fill(abs( xi_proton_right), event_weight );
				histosTH1F["proton_right_xi_total"]->Fill( xi_proton_right, event_weight );
				histosTH1F["proton_right_t"]->Fill( -t_proton_right, event_weight );
				histosTH1F["proton_right_thx"]->Fill(proton_right_thx, event_weight );
				histosTH1F["proton_right_thy"]->Fill(proton_right_thy, event_weight );
				histosTH1F["tx_proton_right"]->Fill(tx_proton_right , event_weight );
				histosTH1F["ty_proton_right "]->Fill(ty_proton_right , event_weight );
				histosTH1F["x0_proton_right"]->Fill(x0_proton_right , event_weight );
				histosTH1F["y0_proton_right "]->Fill(y0_proton_right , event_weight );
                  
                 		if(-xi_proton_right > 0.){
					xi_proton_right = -xi_proton_right;
					histosTH1F["proton_right_logXi"]->Fill( log10(xi_proton_right), event_weight );
					//histosTH1F["pf_xiMinus_minus_proton_right_xi"]->Fill( (pfXiMinusReco - xi_proton_right), event_weight );
					//histosTH2F["proton_right_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_right), event_weight );
					//histosTH2F["proton_right_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_right), event_weight );
					//histosTH2F["proton_right_logXi_vs_t"]->Fill( -t_proton_right, log10(xi_proton_right), event_weight );
				}
			}

			//bool proton_left_valid = rec_proton_left->valid;
			double chi2_proton_left = rec_proton_left->chi2;
			//double chindf_proton_left = rec_proton_left->chindf;
			double xi_proton_left = rec_proton_left->xi;
                        double thx_proton_left = rec_proton_left->thx;
                        double thy_proton_left = rec_proton_left->thy;
                        double tx_proton_left = rec_proton_left->tx;
                        double ty_proton_left = rec_proton_left->ty;
                        double x0_proton_left = rec_proton_left->x0;
                        double y0_proton_left = rec_proton_left->y0;  
			//bool good_proton_left = proton_left_valid && (chi2_proton_left/chindf_proton_left > 1);
			bool good_proton_left = proton_left_valid && (xi_proton_left < 0.);
			if( good_proton_left ){
				//double xi_proton_left = rec_proton_left->xi;
				double t_proton_left = rec_proton_left->t;
				histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
				histosTH1F["proton_left_xi"]->Fill( abs(xi_proton_left), event_weight );
				histosTH1F["proton_left_xi_total"]->Fill( xi_proton_left, event_weight );
				histosTH1F["proton_left_t"]->Fill( -t_proton_left, event_weight );
				histosTH1F["proton_left_thx"]->Fill(proton_left_thx, event_weight );
				histosTH1F["proton_left_thy"]->Fill(proton_left_thy, event_weight );
				histosTH1F["tx_proton_left"]->Fill(tx_proton_left , event_weight );
				histosTH1F["ty_proton_left "]->Fill(ty_proton_left , event_weight );
				histosTH1F["x0_proton_left"]->Fill(x0_proton_left , event_weight );
				histosTH1F["y0_proton_left "]->Fill(y0_proton_left , event_weight );

				if(-xi_proton_left > 0.){
					xi_proton_left = -xi_proton_left;
					histosTH1F["proton_left_logXi"]->Fill( log10(xi_proton_left), event_weight );
					//histosTH1F["pf_xiPlus_minus_proton_left_xi"]->Fill( (pfXiPlusReco - xi_proton_left), event_weight );
					//histosTH2F["proton_left_logXi_vs_pf_logXiPlus"]->Fill( log10(pfXiPlusReco),log10(xi_proton_left), event_weight );
					//histosTH2F["proton_left_logXi_vs_pf_logXiMinus"]->Fill( log10(pfXiMinusReco),log10(xi_proton_left), event_weight );
					//histosTH2F["proton_left_logXi_vs_t"]->Fill( -t_proton_left, log10(xi_proton_left), event_weight );
				}
			}

			bool proton_pair_valid = rec_proton_pair->valid;
			double chi2_proton_pair = rec_proton_pair->chi2;
			//double chindf_proton_pair = rec_proton_pair->chindf;
			double xi_proton_pair_right = rec_proton_pair->xir;
			double xi_proton_pair_left = rec_proton_pair->xil;
			//bool good_proton_pair = proton_pair_valid && (chi2_proton_pair/chindf_proton_pair > 2);
			bool good_proton_pair = proton_pair_valid && (xi_proton_pair_right < 0.) && (xi_proton_pair_left < 0.);
			if( good_proton_pair ){
				histosTH1F["proton_pair_chi2"]->Fill( chi2_proton_pair, event_weight );

				//double xi_proton_pair_right = rec_proton_pair->xir;
				double t_proton_pair_right = rec_proton_pair->tr;
				histosTH1F["proton_pair_right_xi"]->Fill( xi_proton_pair_right, event_weight );
				histosTH1F["proton_pair_right_t"]->Fill( -t_proton_pair_right, event_weight );
				if(-xi_proton_pair_right > 0.)
					histosTH1F["proton_pair_right_logXi"]->Fill( log10(-xi_proton_pair_right), event_weight );
				/*if(xi_proton_pair_right > 0.)
				  histosTH1F["proton_pair_right_logXi"]->Fill( log10(xi_proton_pair_right), event_weight );*/

				//double xi_proton_pair_left = rec_proton_pair->xil;
				double t_proton_pair_left = rec_proton_pair->tl;
				histosTH1F["proton_pair_left_xi"]->Fill( xi_proton_pair_left, event_weight );
				histosTH1F["proton_pair_left_t"]->Fill( -t_proton_pair_left, event_weight );
				if(-xi_proton_pair_left > 0.)
					histosTH1F["proton_pair_left_logXi"]->Fill( log10(-xi_proton_pair_left), event_weight );
				/*if(xi_proton_pair_left > 0.)
				  histosTH1F["proton_pair_left_logXi"]->Fill( log10(xi_proton_pair_left), event_weight );*/
			}



			//------------------
			
					
					if( good_proton_right ){
                                           n_dimuons_proton_right++;
					   double t_proton_right = rec_proton_right->t;
				           histosTH1F["dimuon_proton_right_xi"]->Fill( xi_proton_right, event_weight );
                                           histosTH1F["dimuon_proton_right_t"]->Fill( fabs(t_proton_right), event_weight );
                                           //histosTH1F["dimuon_multiplicity_proton_right"]->Fill(n_dimuons_selected, event_weight );
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
                                                         n_dimuons_t_selected_proton_right++;
			                                if(verbose) cout<< "t proton right side (abs):"<< "["<< t_proton_down_  << ","<< t_proton_up_ << "]"<< endl;
							//histosTH1F["dimuon_multiplicity_t_cut_proton_right"]->Fill(n_dimuons_selected, event_weight );
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


							/*if (-xi_proton_right<xi_proton_down_ || -xi_proton_right>xi_proton_up_){
							     n_dimuons_t_at_xi_selected_proton_right++;
                                                              if(verbose) cout<< "xi proton right side (abs):"<< "["<< xi_proton_down_  << ","<< xi_proton_up_ << "]"<< endl;
                                                                
                                                                histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_right"]->Fill(n_dimuons_selected, event_weight );
								histosTH1F["dimuon_mass_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.M(), event_weight );
								histosTH1F["dimuon_pt_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Pt(), event_weight );
                                                                histosTH1F["dimuon_pt2_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );

								histosTH1F["dimuon_eta_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Eta(), event_weight );
								histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
								histosTH1F["proton_right_xi_t_cut_xi_cut"]->Fill( xi_proton_right, event_weight );
                                                                histosTH1F["proton_right_t_t_cut_xi_cut"]->Fill( fabs(t_proton_right), event_weight );  
                                                                ///
                                                                histosTH1F["muon1_pt_t_cut_xi_cut_proton_right"]->Fill( muon1_lorentz.Pt(), event_weight );
								histosTH1F["muon1_eta_t_cut_xi_cut_proton_right"]->Fill( muon1_lorentz.Eta(), event_weight );
								histosTH1F["muon1_phi_t_cut_xi_cut_proton_right"]->Fill( muon1_lorentz.Phi(), event_weight );
								histosTH1F["muon1_rapidity_t_cut_xi_cut_proton_right"]->Fill( muon1_lorentz.Rapidity(),event_weight );
								///
                                                                histosTH1F["muon2_pt_t_cut_xi_cut_proton_right"]->Fill( muon2_lorentz.Pt(), event_weight );
								histosTH1F["muon2_eta_t_cut_xi_cut_proton_right"]->Fill( muon2_lorentz.Eta(), event_weight );
								histosTH1F["muon2_phi_t_cut_xi_cut_proton_right"]->Fill( muon2_lorentz.Phi(), event_weight );
								histosTH1F["muon2_rapidity_t_cut_xi_cut_proton_right"]->Fill( muon2_lorentz.Rapidity(),event_weight );

								histosTH1F["muonDeltaPt_t_cut_xi_cut_proton_right"]->Fill(deltapt, event_weight );
                                                                histosTH1F["muonDeltaEta_t_cut_xi_cut_proton_right"]->Fill(deltaeta, event_weight );	
                                                                histosTH1F["muonDeltaPhi_t_cut_xi_cut_proton_right"]->Fill(deltaphi, event_weight );
                                                                histosTH1F["muonDeltaY_t_cut_xi_cut_proton_right"]->Fill(deltay, event_weight ); 
							}*/
						}
					}
					if( good_proton_left ){
                                           n_dimuons_proton_left++;
	              			   double t_proton_left = rec_proton_left->t;
                                           histosTH1F["dimuon_proton_left_xi"]->Fill( xi_proton_left, event_weight );
                                           histosTH1F["dimuon_proton_left_t"]->Fill( fabs(t_proton_left), event_weight ); 
                                           //histosTH1F["dimuon_multiplicity_proton_left"]->Fill(n_dimuons_selected, event_weight );
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
							n_dimuons_t_selected_proton_left++;	
                                                        if(verbose) cout<< "t proton left side (abs) l1242:"<< "["<< t_proton_down_  << ","<< t_proton_up_ << "]"<< endl;
 
							//histosTH1F["dimuon_multiplicity_t_cut_proton_left"]->Fill(n_dimuons_selected, event_weight );
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

							/*if (-xi_proton_left<xi_proton_down_ || -xi_proton_left>xi_proton_up_){
                                                                n_dimuons_t_at_xi_selected_proton_left++;
                                                                if(verbose) cout<< "xi proton left side (abs):"<<"["<< xi_proton_down_  << ","<< xi_proton_up_ << "]"<< endl; 
                                                                 
								histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_left"]->Fill(n_dimuons_selected, event_weight );
								histosTH1F["dimuon_mass_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.M(), event_weight );
                                                                histosTH1F["dimuon_pt2_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );

								histosTH1F["dimuon_pt_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Pt(), event_weight );
								histosTH1F["dimuon_eta_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Eta(), event_weight );
								histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
 								histosTH1F["proton_left_xi_t_cut_xi_cut"]->Fill( xi_proton_left, event_weight );
                                                                histosTH1F["proton_left_t_t_cut_xi_cut"]->Fill( fabs(t_proton_left), event_weight );
                                                                ///
                                                                histosTH1F["muon1_pt_t_cut_xi_cut_proton_left"]->Fill( muon1_lorentz.Pt(), event_weight );
								histosTH1F["muon1_eta_t_cut_xi_cut_proton_left"]->Fill( muon1_lorentz.Eta(), event_weight );
								histosTH1F["muon1_phi_t_cut_xi_cut_proton_left"]->Fill( muon1_lorentz.Phi(), event_weight );
								histosTH1F["muon1_rapidity_t_cut_xi_cut_proton_left"]->Fill( muon1_lorentz.Rapidity(),event_weight );
								///
                                                                histosTH1F["muon2_pt_t_cut_xi_cut_proton_left"]->Fill( muon2_lorentz.Pt(), event_weight );
								histosTH1F["muon2_eta_t_cut_xi_cut_proton_left"]->Fill( muon2_lorentz.Eta(), event_weight );
								histosTH1F["muon2_phi_t_cut_xi_cut_proton_left"]->Fill( muon2_lorentz.Phi(), event_weight );
								histosTH1F["muon2_rapidity_t_cut_xi_cut_proton_left"]->Fill( muon2_lorentz.Rapidity(),event_weight );

								histosTH1F["muonDeltaPt_t_cut_xi_cut_proton_left"]->Fill(deltapt, event_weight );
                                                                histosTH1F["muonDeltaEta_t_cut_xi_cut_proton_left"]->Fill(deltaeta, event_weight );	
                                                                histosTH1F["muonDeltaPhi_t_cut_xi_cut_proton_left"]->Fill(deltaphi, event_weight );
                                                                histosTH1F["muonDeltaY_t_cut_xi_cut_proton_left"]->Fill(deltay, event_weight ); 
							}*/

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
                                               n_jpsi_t_selected_proton_right++;
			                       if(verbose) cout<< "t proton right side (abs):"<< "["<< t_proton_down_  << ","<< t_proton_up_ << "]"<< endl;
                                               outstring_right << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;
					       //histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_right"]->Fill(n_dimuons_selected, event_weight );
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

					       /*if (-xi_proton_right<xi_proton_down_ || -xi_proton_right>xi_proton_up_){
					          n_jpsi_t_at_xi_selected_proton_right++;
                                                  if(verbose) cout<< "xi proton right side (abs):"<< "["<< xi_proton_down_  << ","<< xi_proton_up_ << "]"<< endl;
                                                  outstring_right << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;
              
                                                  histosTH1F["jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_right"]->Fill(n_dimuons_selected, event_weight );
						  histosTH1F["jpsi_dimuon_mass_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.M(), event_weight );
						  histosTH1F["jpsi_dimuon_pt_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Pt(), event_weight );
                                                  histosTH1F["jpsi_dimuon_pt2_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );

						  histosTH1F["jpsi_dimuon_eta_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Eta(), event_weight );
						  histosTH1F["jpsi_dimuon_rapidity_t_cut_xi_cut_proton_right"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
						  histosTH1F["jpsi_proton_right_xi_t_cut_xi_cut"]->Fill( xi_proton_right, event_weight );
                                                  histosTH1F["jpsi_proton_right_t_t_cut_xi_cut"]->Fill( fabs(t_proton_right), event_weight );  

					          histosTH1F["muonDeltaPt_jpsi_t_cut_xi_cut_proton_right"]->Fill(deltapt, event_weight );
                                                  histosTH1F["muonDeltaEta_jpsi_t_cut_xi_cut_proton_right"]->Fill(deltaeta, event_weight );	
                                                  histosTH1F["muonDeltaPhi_jpsi_t_cut_xi_cut_proton_right"]->Fill(deltaphi, event_weight );
                                                  histosTH1F["muonDeltaY_jpsi_t_cut_xi_cut_proton_right"]->Fill(deltay, event_weight ); 
                                                  
					       } // xi, t, right proton */
					    } // t right proton
                                          } // right proton
                                        } //jpsi mass
                                          if( good_proton_left ){
                                            if((dimuon_mass > 3.0) && (dimuon_mass < 3.2)){
                                            n_jpsi_proton_left++;
					    double t_proton_left = rec_proton_left->t;
					    
                                            histosTH1F["jpsi_proton_left_xi"]->Fill(xi_proton_left, event_weight );
                                            histosTH1F["jpsi_proton_left_t"]->Fill( fabs(t_proton_left), event_weight );
                                            //histosTH1F["jpsi_multiplicity_proton_left"]->Fill(n_dimuons_selected, event_weight );
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
						  n_jpsi_t_selected_proton_left++;
                                                  outstring_left << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;
						  	
                                                  if(verbose) cout<< "t proton left side (abs):"<< "["<< t_proton_down_  << ","<< t_proton_up_ << "]"<< endl;
 
						  histosTH1F["jpsi_dimuon_multiplicity_t_cut_proton_left"]->Fill(n_dimuons_selected, event_weight );
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

					       /*if (-xi_proton_left<xi_proton_down_ || -xi_proton_left>xi_proton_up_){
                                                  n_jpsi_t_at_xi_selected_proton_left++;
                                                  
                                                  if(verbose) cout<< "xi proton left side (abs):"<<"["<< xi_proton_down_  << ","<< xi_proton_up_ << "]"<< endl; 
                                                  outstring_left << evtId->Run << ":"<< evtId->LumiSect << ":"<< evtId->Evt << endl;
               
						  histosTH1F["jpsi_dimuon_multiplicity_t_cut_xi_cut_proton_left"]->Fill(n_dimuons_selected, event_weight );
						  histosTH1F["jpsi_dimuon_mass_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.M(), event_weight );
						  histosTH1F["jpsi_dimuon_pt2_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Pt()*dimuon_lorentz.Pt(), event_weight );

                                                  histosTH1F["jpsi_dimuon_pt_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Pt(), event_weight );
						  histosTH1F["jpsi_dimuon_eta_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Eta(), event_weight );
						  histosTH1F["jpsi_dimuon_rapidity_t_cut_xi_cut_proton_left"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
 						  histosTH1F["jpsi_proton_left_xi_t_cut_xi_cut"]->Fill( xi_proton_left, event_weight );
                                                  histosTH1F["jpsi_proton_left_t_t_cut_xi_cut"]->Fill( fabs(t_proton_left), event_weight );
                                                  histosTH1F["muonDeltaPt_jpsi_t_cut_xi_cut_proton_left"]->Fill(deltapt, event_weight );
                                                  histosTH1F["muonDeltaEta_jpsi_t_cut_xi_cut_proton_left"]->Fill(deltaeta, event_weight );	
                                                  histosTH1F["muonDeltaPhi_jpsi_t_cut_xi_cut_proton_left"]->Fill(deltaphi, event_weight );
                                                  histosTH1F["muonDeltaY_jpsi_t_cut_xi_cut_proton_left"]->Fill(deltay, event_weight ); 

					       } //xi, t, proton left*/

					    } // t, proton left
					  } // good proton left
                                    }//jpsi mass
				//} //muon 1
			//} //muon 2

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
                cout<<"n_dimuons_t_at_xi_selected_proton_left="<< n_dimuons_t_at_xi_selected_proton_left<<endl;
                cout<<"n_dimuons_t_selected_proton_right="<< n_dimuons_t_selected_proton_right<<endl;
                cout<<"n_dimuons_t_at_xi_selected_proton_right="<< n_dimuons_t_at_xi_selected_proton_right<<endl;	
                cout<<"n_jpsi_selected="<< n_jpsi_selected<<endl;
                cout<<"n_jpsi_proton_left="<<n_jpsi_proton_left<<endl;
                cout<<"n_jpsi_proton_right="<<n_jpsi_proton_right<<endl;
                cout<<"n_jpsi_t_selected_proton_left="<< n_jpsi_t_selected_proton_left<<endl;
                cout<<"n_jpsi_t_at_xi_selected_proton_left="<< n_jpsi_t_at_xi_selected_proton_left<<endl;
                cout<<"n_jpsi_t_selected_proton_right="<< n_jpsi_t_selected_proton_right<<endl;
                cout<<"n_jpsi_t_at_xi_selected_proton_right="<< n_jpsi_t_at_xi_selected_proton_right<<endl;	
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
