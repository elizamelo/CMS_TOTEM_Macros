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
#include <TString.h>

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

//#include "analysis_tools.h"
#include "rp_aperture_config.h"


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

double etaBinsHCALBoundaries[] = {-5.205, -4.903, -4.730,
	-4.552, -4.377, -4.204, -4.027, -3.853, -3.677, -3.503, -3.327, -3.152,
	-3.000, -2.868, -2.650, -2.500,
	-2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479,
	-1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.870, -0.783,
	-0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
	0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696,
	0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392,
	1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322,
	2.500, 2.650, 2.868, 3.000,
	3.152, 3.327, 3.503, 3.677, 3.853, 4.027, 4.204, 4.377, 4.552,
	4.730, 4.903, 5.205}; // 41 + 41 bins



void analysis_UATree_DYmuons(vector<string> const& fileNames, string const& outputFileName = "output.root", const Double_t t_proton_down_=0.003, const Double_t t_proton_up_=1.0, const Double_t xi_proton_down_=0.01,const Double_t xi_proton_up_=0.15,const Int_t Bin_mass=200 ,const Int_t nevt_max = 57754){

	bool isMC  = true;
	bool verbose = false;
	string treeName = "evt";//"cms_totem";
	double etaMaxThreshold = 2.0;


	bool selectElastic = false;
	//bool selectNonElastic = false;
	bool selectNonElastic = true;

	//==============================
	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


	// Declaration of histograms
	map<string,TH1F*> histosTH1F;

	/*	vector<string> selections;
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
	 */
	Float_t tbins[12] = {0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.24, 0.30, 0.36, 0.45, 0.65, 1.};
	Float_t xi_bins[29] = {0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15};

	/*	for(size_t k = 0; k < selections.size(); ++k)
		histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

		histosTH1F["bunchCrossingNumber"] = new TH1F("bunchCrossingNumber", "bunchCrossingNumber" , 3900 , 0 , 3900);

		histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
		histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

		int nBinsHLT = hltPathNames.size(); 
		histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);
		for(size_t k = 0; k < nBinsHLT; ++k) 
		histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1), hltPathNames[k].c_str() );
	 */
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
	// MC gen histograms

	if(isMC){
		histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 100);

	}


	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
	histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);

	histosTH1F["muon_pt"] = new TH1F("muon_pt", "p_{T}(muon)" , 150 , 0. , 150.);
	histosTH1F["muon_eta"] = new TH1F("muon_eta", "#eta(muon)" , 200 , -5.2 , 5.2);
	histosTH1F["muon_phi"] = new TH1F("muon_phi", "#phi(muon)" , 200 , -M_PI , M_PI);
	histosTH1F["muon_multiplicity"] = new TH1F("muon_multiplicity", "n muons" , 100 , 0 , 100);

	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt"] = new TH1F("dimuon_pt", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta"] = new TH1F("dimuon_eta", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity"] = new TH1F("dimuon_rapidity", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity"] = new TH1F("dimuon_multiplicity", "n dimuons" , 100 , 0 , 100); 

	//_t_cut_proton_right
	histosTH1F["dimuon_mass_t_cut_proton_right"] = new TH1F("dimuon_mass_t_cut_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt_t_cut_proton_right"] = new TH1F("dimuon_pt_t_cut_proton_right", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta_t_cut_proton_right"] = new TH1F("dimuon_eta_t_cut_proton_right", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_right"] = new TH1F("dimuon_rapidity_t_cut_proton_right", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_right"] = new TH1F("dimuon_multiplicity_t_cut_proton_right", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_right_xi_t_cut"] = new TH1F("proton_right_xi_t_selected", "#xi" , 28 , xi_bins);
	histosTH1F["proton_right_t_t_cut"] = new TH1F("proton_right_t_t_selected", "|t|" , 11 , tbins);

	//t_cut_xi_cut_proton_right
	histosTH1F["dimuon_mass_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_mass_t_cut_xi_cut_proton_right", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_pt_t_cut_xi_cut_proton_right", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_eta_t_cut_xi_cut_proton_right", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_rapidity_t_cut_xi_cut_proton_right", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_right"] = new TH1F("dimuon_multiplicity_t_cut_xi_cut_proton_right", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_right_xi_t_cut_xi_cut"] = new TH1F("proton_right_xi_selected_t_and_xi", "#xi" , 28 , xi_bins);
	histosTH1F["proton_right_t_t_cut_xi_cut"] = new TH1F("proton_right_t_selected_t_and_xi", "|t|" , 11 , tbins);

	//_t_cut_proton_left
	histosTH1F["dimuon_mass_t_cut_proton_left"] = new TH1F("dimuon_mass_t_cut_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt_t_cut_proton_left"] = new TH1F("dimuon_pt_t_cut_proton_left", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta_t_cut_proton_left"] = new TH1F("dimuon_eta_t_cut_proton_left", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_proton_left"] = new TH1F("dimuon_rapidity_t_cut_proton_left", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_proton_left"] = new TH1F("dimuon_multiplicity_t_cut_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_left_xi_t_cut"] = new TH1F("proton_left_xi_t_selected", "#xi" , 28 , xi_bins);
	histosTH1F["proton_left_t_t_cut"] = new TH1F("proton_left_t_t_selected", "|t|" , 11 , tbins);
	//t_cut_xi_cut_proton_left
	histosTH1F["dimuon_mass_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_mass_t_cut_xi_cut_proton_left", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_pt_t_cut_xi_cut_proton_left", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_eta_t_cut_xi_cut_proton_left", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_rapidity_t_cut_xi_cut_proton_left", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_t_cut_xi_cut_proton_left"] = new TH1F("dimuon_multiplicity_t_cut_xi_cut_proton_left", "n dimuons" , 100 , 0 , 100);
	histosTH1F["proton_left_xi_t_cut_xi_cut"] = new TH1F("proton_left_xi_selected_t_and_xi ", "#xi" , 28 ,xi_bins);
	histosTH1F["proton_left_t_t_cut_xi_cut"] = new TH1F("proton_left_t_selected_t_and_xi", "|t|" , 11 , tbins);


	histosTH1F["pf_etaMax"] = new TH1F("pf_etaMax","#eta^{max}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_etaMin"] = new TH1F("pf_etaMin","#eta^{min}",82,etaBinsHCALBoundaries);
	histosTH1F["pf_deltaEta"] = new TH1F("pf_deltaEta","#Delta#eta",100,0.,10.);
	histosTH1F["pf_EPlusPz"] = new TH1F("pf_EPlusPz","sum(E + pz)",24,binningEPlusPz);
	histosTH1F["pf_EMinusPz"] = new TH1F("pf_EMinusPz","sum(E - pz)",24,binningEPlusPz);
	histosTH1F["pf_xiPlus"] = new TH1F("pf_xiPlus","#xi^{+}",200,-1.,1.);
	histosTH1F["pf_xiMinus"] = new TH1F("pf_xiMinus","#xi^{-}",200,-1.,1.);
	histosTH1F["pf_logXiPlus"] = new TH1F("pf_logXiPlus","log(#xi^{+})",200,-5.,0.);
	histosTH1F["pf_logXiMinus"] = new TH1F("pf_logXiMinus","log(#xi^{-})",200,-5.,0.);

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

	histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 28 , xi_bins);
	histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 11 , tbins);
	histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 100 , 0. , 100.);
	histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 28 ,xi_bins);
	histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",200,-5.,0.);
	histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 11 , tbins);
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

	map<string,TH2F*> histosTH2F;  
	histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
	histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
	histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

	histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
	histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
	histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 200, 0., 5., 200, -5., 0.);

// MC 
  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 50 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi+", "#xi^{+}" , 82 , 0,4);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi-", "#xi^{-}" , 82 , 0,4);
  histosTH1F["logxi_plus"] = new TH1F("logxi+", "Log #xi^{+}" , 82 , -3,0.5);
  histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 82 , -3,0.5);
  histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 82 , -3,0.5);
  histosTH1F["correction"] = new TH1F("correction", "Correction factor" , 82 , 0,2);
  histosTH1F["resolution_after"] = new TH1F("resolution_after", "Resolution" , 82 , -2,2);
  histosTH1F["resolution_before"] = new TH1F("resolution_before", "Resolution" , 82 , -2,2);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);

  histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , 200, 0., 1.);
  histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["t_proton_plus_tbin1"] = new TH1F("t_proton_plus_tbin1", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["t_proton_plus_tbin2"] = new TH1F("t_proton_plus_tbin2", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["t_proton_plus_tbin3"] = new TH1F("t_proton_plus_tbin3", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["t_proton_plus_tbin4"] = new TH1F("t_proton_plus_tbin4", "t_proton_plus" , 200, 0., 5.);
  histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
  histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);




//MC proton info
//map<string,TH2F*> histosTH2F;
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
  histosTH2F["xi_plus_reco_gen"] = new TH2F("xi+","xi+",82,0,0.5,82,0,0.5);
  histosTH2F["xi_minus_reco_gen"] = new TH2F("xi-","xi-",82,0,0.5,82,0,0.5);
  histosTH2F["logxi_plus_reco_gen"] = new TH2F("logxi+","xi+",82,-3,0.5,82,-3,0.5);
  histosTH2F["logxi_minus_reco_gen"] = new TH2F("logxi-","xi-",82,-3,0.5,82,-3,0.5);
//	ofstream outstring_right(outtxt_right);

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
		tree->SetBranchAddress("evtId",&evtId);
//		tree->SetBranchAddress("cmsEvtUA",&evtId);
//		tree->SetBranchAddress("cmsTrigUA",&l1Trig);
//		tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
		tree->SetBranchAddress("cmsTracksUA",&track_coll); 
		tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
		//tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
		//tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
		tree->SetBranchAddress("cmsMuonsUA",&muon_coll);
		tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
		tree->SetBranchAddress("genKin",&genKin);
    		tree->SetBranchAddress("genPart",&genPart);
		//tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
		//tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
/*		tree->SetBranchAddress("branchT2EV.",&t2_event);
		tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
		tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
		tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
*/		
			/*//rp_track_info[20] = NULL;
		  tree->SetBranchAddress("track_rp_20.", &rp_track_info[20]);*/
		
/*
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
*/

		/*readRPBranches(tree,
		  rec_proton_left,
		  rec_proton_right,
		  rec_proton_pair,
		  rp_track_info,
		  rp_digi_info,
		  rp_par_patterns_info,
		  rp_nonpar_patterns_info,
		  rp_multi_track_info);*/

		if(isMC) tree->SetBranchAddress("genPart",&genPart);

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
//			bool passed_HLTMuon= false;
//			string HLT_muon = "HLT_L1DoubleMu0_v1"; 
			bool passed_vrt= false;      
			bool passed_muon = false;

/*
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
*/



			//-------------------------------------------------------------------------------------------------
			//filling pt distribution for the generated particles
			//ie those from pythia generator, without reconstruction
			if(isMC){
			  for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
			  pt_gen->Fill(p->Pt());
			  }

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
			histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );


		
      //GenPart
      double genEPlusPz = 0;
      double genEMinusPz = 0;
      double cm = 8000;
      double proton_pi = 4000;
      double proton_pz_plus=-999;
      double proton_px_plus = -999.;
      double proton_py_plus = -999.;
      double proton_energy_plus = 0.;
      double proton_pz_minus=999;
      double proton_px_minus = 999.;
      double proton_py_minus = 999.;
      double proton_energy_minus = 0.;
      double px_gen;
      double py_gen;
      double pz_gen;
      double energy_gen;
      double proton_pf;

      for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){

         double eta_gen = it_genpart->Eta();
         int status = it_genpart->status;
         int id = it_genpart->pdgId;

         if (status != 1) continue; //final state for the particles
         if (id != 2212) continue;
//       if (eta_gen<4.9 && eta_gen>-4.9){ xi_posit_gen = true;
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
            proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
            double pz_cut = 0.7*proton_pi;
            if (fabs(pz_gen) < pz_cut) continue;
		            genEPlusPz += (energy_gen + pz_gen);
            genEMinusPz += (energy_gen - pz_gen);

         if (pz_gen > proton_pz_plus) {
             proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
             proton_px_plus = px_gen; proton_py_plus = py_gen;
         }
         if (pz_gen < proton_pz_minus) {
             proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
             proton_px_minus = px_gen; proton_py_minus = py_gen;
         }

      }

      double xi_plus_gen = genEPlusPz/cm; //cout<<xi1_gen<<endl;
      double xi_minus_gen = genEMinusPz/cm;
      double xi_proton_plus = -1.;
      double xi_proton_minus = -1.;
      double t_proton_plus = 0.;
      double t_proton_minus = 0.;
      double thx_proton_plus = 0.;
      double thy_proton_plus = 0.;
      double thx_proton_minus = 0.;
      double thy_proton_minus = 0.;

      ++nevents_gen ;
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
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
         proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);

      }
   if(proton_pz_minus < 0.){
         xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) );

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

	  // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999.;
/*      double cm = 8000;*/

      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
         int partType = it_pfcand->particleId;
         double eta = it_pfcand->Eta();
         double energy = it_pfcand->Energy();
         double pz = it_pfcand->Pz();

         // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;

         soma1 += (energy + pz);
         soma2 += (energy - pz);

         if (eta > eta_max) {eta_max = eta; PF_eta_max = true;}
         if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}


         histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );
/*         if(partType == MyPFCand::h || partType == MyPFCand::h0 || partType == MyPFCand::h_HF){
           if(eta<-5 && eta>5) continue;
           Double_t mx = genPart->
           Double_t xigen = mx*mx/cm;
         }*/

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
            histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );}
            /*if( part->ecalEnergy() > 0. ) histosTH2F["energyVsEtaHadronHFEcalEnergy"]->Fill( eta, energy, event_weight );
            else                          histosTH2F["energyVsEtaHadronHFNoEcalEnergy"]->Fill( eta, energy, event_weight );*/
         else if(partType == MyPFCand::egamma_HF)
            histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );

       }
       if(!PF_eta_max) continue;
       if(!PF_eta_min) continue;
       weight_total_PF_selected += event_weight;

       ++nevents_pf;

       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;

       double correction = xi_plus_Reco/xi_plus_gen;
       double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
       double xi_reconst = xi_plus_Reco/0.8;
       double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;

    //rp_accept

      if (fabs(t_proton_plus)<0.03 || fabs(t_proton_plus)>1)continue;
      if (xi_proton_plus<0.03 || xi_proton_plus>0.1)continue;
      histosTH1F["xi_proton_plus"]->Fill( xi_proton_plus , event_weight );
      histosTH1F["t_proton_plus"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["t_proton_plus_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
      histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );
      histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );

      bool proton_plus_rp_accept = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 );


      if(!proton_plus_rp_accept)continue;
      histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
      histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
      histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["t_proton_plus_accepted_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );


       histosTH1F["vtx_zpos_after"]->Fill( zpos, event_weight );
       histosTH1F["vtx_xpos_after"]->Fill( xpos, event_weight );
       histosTH1F["vtx_ypos_after"]->Fill( ypos, event_weight );

       histosTH2F["xi_plus_reco_gen"]->Fill( xi_plus_gen, xi_plus_Reco, event_weight );
       histosTH2F["xi_minus_reco_gen"]->Fill( xi_minus_gen, xi_minus_Reco, event_weight );
       histosTH2F["logxi_plus_reco_gen"]->Fill( log10(xi_plus_gen), log10(xi_plus_Reco), event_weight );
       histosTH2F["logxi_minus_reco_gen"]->Fill( log10(xi_minus_gen), log10(xi_minus_Reco), event_weight );

       histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
       histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
       histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
       histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );
       histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
       histosTH1F["logxi_plus"]->Fill( log10(xi_plus_Reco), event_weight  );
	
			// Muons
			//int n_muon_selected = 0;
/*			vector<MyMuon> muons_selected;
			for(vector<MyMuon>::iterator it_muon = muon_coll->begin() ; it_muon != muon_coll->end() ; ++it_muon){

				if( !(it_muon->IsTrackerMuon || it_muon->IsGlobalMuon) ) continue;

				MyTracks const& muon_innerTrack = it_muon->innerTrack;
				bool muon_id = it_muon->TMOneStationAngTight &&
					muon_innerTrack.chi2n < 1.8 &&
					muon_innerTrack.nValidPixelHits > 0 &&
					muon_innerTrack.vtxdxy[prim_vtx_id] < 3. &&
					muon_innerTrack.vtxdz[prim_vtx_id] < 30.;

				if( !muon_id ) continue;

					primaryVertex.ndof > 4 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
			// Include Muon selection REf CMS AN AN 12-067 
			//bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity &&
			//                        primaryVertex.ndof < 10 && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0); 

*/			if(selectVertex && !select_Vertex) continue;

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
			histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );


		
      //GenPart
      double genEPlusPz = 0;
      double genEMinusPz = 0;
      double cm = 8000;
      double proton_pi = 4000;
      double proton_pz_plus=-999;
      double proton_px_plus = -999.;
      double proton_py_plus = -999.;
      double proton_energy_plus = 0.;
      double proton_pz_minus=999;
      double proton_px_minus = 999.;
      double proton_py_minus = 999.;
      double proton_energy_minus = 0.;
      double px_gen;
      double py_gen;
      double pz_gen;
      double energy_gen;
      double proton_pf;

      for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){

         double eta_gen = it_genpart->Eta();
         int status = it_genpart->status;
         int id = it_genpart->pdgId;

         if (status != 1) continue; //final state for the particles
         if (id != 2212) continue;
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
            proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
            double pz_cut = 0.7*proton_pi;
            if (fabs(pz_gen) < pz_cut) continue;
		            genEPlusPz += (energy_gen + pz_gen);
            genEMinusPz += (energy_gen - pz_gen);

         if (pz_gen > proton_pz_plus) {
             proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
             proton_px_plus = px_gen; proton_py_plus = py_gen;
         }
         if (pz_gen < proton_pz_minus) {
             proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
             proton_px_minus = px_gen; proton_py_minus = py_gen;
         }

      }

      double xi_plus_gen = genEPlusPz/cm; //cout<<xi1_gen<<endl;
      double xi_minus_gen = genEMinusPz/cm;
      double xi_proton_plus = -1.;
      double xi_proton_minus = -1.;
      double t_proton_plus = 0.;
      double t_proton_minus = 0.;
      double thx_proton_plus = 0.;
      double thy_proton_plus = 0.;
      double thx_proton_minus = 0.;
      double thy_proton_minus = 0.;

      ++nevents_gen ;
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
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
         proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);

      }
   if(proton_pz_minus < 0.){
         xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) );

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

	  // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999.;
/*      double cm = 8000;*/

      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
         int partType = it_pfcand->particleId;
         double eta = it_pfcand->Eta();
         double energy = it_pfcand->Energy();
         double pz = it_pfcand->Pz();

         // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;

         soma1 += (energy + pz);
         soma2 += (energy - pz);

         if (eta > eta_max) {eta_max = eta; PF_eta_max = true;}
         if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}


         histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );
/*         if(partType == MyPFCand::h || partType == MyPFCand::h0 || partType == MyPFCand::h_HF){
           if(eta<-5 && eta>5) continue;
           Double_t mx = genPart->
           Double_t xigen = mx*mx/cm;
         }*/

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
            histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );}
         else if(partType == MyPFCand::egamma_HF)
            histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );

       }
       if(!PF_eta_max) continue;
       if(!PF_eta_min) continue;
       weight_total_PF_selected += event_weight;

       ++nevents_pf;

       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;

       double correction = xi_plus_Reco/xi_plus_gen;
       double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
       double xi_reconst = xi_plus_Reco/0.8;
       double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;

    //rp_accept

      if (fabs(t_proton_plus)<0.03 || fabs(t_proton_plus)>1)continue;
      if (xi_proton_plus<0.03 || xi_proton_plus>0.1)continue;
      histosTH1F["xi_proton_plus"]->Fill( xi_proton_plus , event_weight );
      histosTH1F["t_proton_plus"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["t_proton_plus_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
      histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );
      histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );

      bool proton_plus_rp_accept = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 );


      if(!proton_plus_rp_accept)continue;
      histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
      histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
      histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight );
      histosTH1F["t_proton_plus_accepted_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );


       histosTH1F["vtx_zpos_after"]->Fill( zpos, event_weight );
       histosTH1F["vtx_xpos_after"]->Fill( xpos, event_weight );
       histosTH1F["vtx_ypos_after"]->Fill( ypos, event_weight );

       histosTH2F["xi_plus_reco_gen"]->Fill( xi_plus_gen, xi_plus_Reco, event_weight );
       histosTH2F["xi_minus_reco_gen"]->Fill( xi_minus_gen, xi_minus_Reco, event_weight );
       histosTH2F["logxi_plus_reco_gen"]->Fill( log10(xi_plus_gen), log10(xi_plus_Reco), event_weight );
       histosTH2F["logxi_minus_reco_gen"]->Fill( log10(xi_minus_gen), log10(xi_minus_Reco), event_weight );

       histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
       histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
       histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
       histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );
       histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
       histosTH1F["logxi_plus"]->Fill( log10(xi_plus_Reco), event_weight  );
	

		} // End of loop over events in a file
/////////////////////////////////////////////////////////////////////////////////////
		cout<<"Total of evts="<< nev << endl << *itfiles << endl;
		cout<<"n_vertices_selected ="<< n_vertices_selected << endl;
		cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
		cout<<"n_tracks_selected="<< n_tracks_selected<<endl;
		cout<<"n_muon_selected="<< n_muon_selected<<endl;
		cout<<"n_dimuons_selected="<< n_dimuons_selected<<endl;
		cout<<"n_dimuons_t_selected_proton_left="<< n_dimuons_t_selected_proton_left<<endl;
		cout<<"n_dimuons_t_at_xi_selected_proton_left="<< n_dimuons_t_at_xi_selected_proton_left<<endl;
		cout<<"n_dimuons_t_selected_proton_right="<< n_dimuons_t_selected_proton_right<<endl;
		cout<<"n_dimuons_t_at_xi_selected_proton_right="<< n_dimuons_t_at_xi_selected_proton_right<<endl;	

		// Close current file
		file->Close();
//		outstring_left.close();
//		outstring_right.close();
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
