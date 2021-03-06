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
#include <TRandom.h>
#include <TSystemDirectory.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

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

#include "rp_aperture_config.h"

#define PI 3.141592653589793
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

/*
 */
void CMS_MC_POMWIG_DY_plus_v01(vector<string> const& fileNames, string const& outputFileName = "output.root", const Double_t t_proton_down_=0.003, const Double_t t_proton_up_=1.0, const Double_t xi_proton_down_=0.01,const Double_t xi_proton_up_=0.15,const Int_t Bin_mass=200 ,const Int_t nevt_max = -1){

	//(string const& outputFileName = "mc_pomwig_plus_RP.root", const Int_t nevt_max = -1){
	bool isMC  = false;
	bool verbose = false;
	bool selectVertex = true;	
	bool selectMuons = true;
	string treeName = "evt";//"cms_totem";

	//string jetCollName = "ak5PFJets";
	//string jetCorrName = "ak5PFL2L3Residual";
	//  string jetCorrName = "ak5PFL2L3"; 

	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


	// Declaration of histograms
	map<string,TH1F*> histosTH1F;


	//Vertex plots
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

	//	histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
	//	histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
	//	histosTH1F["vtx_xpos_random"] = new TH1F("vtx_xpos_random", "x(vtx)" , 150 , -1.5 , 1.5);
	//	histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
	//	histosTH1F["vtx_ypos_random"] = new TH1F("vtx_ypos_random", "y(vtx)" , 150 , -1.5 , 1.5);
	//	histosTH1F["vtx_zpos_after"] = new TH1F("vtx_zpos_after", "z(vtx)" , 150 , -30. , 30.);
	//	histosTH1F["vtx_xpos_after"] = new TH1F("vtx_xpos_after", "x(vtx)" , 150 , -1.5 , 1.5);
	//	histosTH1F["vtx_ypos_after"] = new TH1F("vtx_ypos_after", "y(vtx)" , 150 , -1.5 , 1.5);

	//MC Gen 
	histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
	//
	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
	histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);

	// Reco Muons plots
	histosTH1F["muon_pt"] = new TH1F("muon_pt", "p_{T}(muon)" , 150 , 0. , 150.);
	histosTH1F["muon_eta"] = new TH1F("muon_eta", "#eta(muon)" , 200 , -5.2 , 5.2);
	histosTH1F["muon_phi"] = new TH1F("muon_phi", "#phi(muon)" , 200 , -M_PI , M_PI);
	histosTH1F["muon_multiplicity"] = new TH1F("muon_multiplicity", "n muons" , 100 , 0 , 100);

	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt"] = new TH1F("dimuon_pt", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta"] = new TH1F("dimuon_eta", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity"] = new TH1F("dimuon_rapidity", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity"] = new TH1F("dimuon_multiplicity", "n dimuons" , 100 , 0 , 100);

	// Reco Muons rp_accept
	histosTH1F["dimuon_mass_rp_accept"] = new TH1F("dimuon_mass_rp_accept", "mass(mu1,mu2)" , Bin_mass , 0. , 120.);
	histosTH1F["dimuon_pt_rp_accept"] = new TH1F("dimuon_pt_rp_accept", "p_{T}(mu1,mu2)" , 150 , 0. , 150.);
	histosTH1F["dimuon_eta_rp_accept"] = new TH1F("dimuon_eta_rp_accept", "#eta(mu1,mu2)" , 200 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity_rp_accept"] = new TH1F("dimuon_rapidity_rp_accept", "y(mu1,mu2)" , 200 , -15. , 15.);
	histosTH1F["dimuon_multiplicity_rp_accept"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);


	/*
	   histosTH1F["jet_pt"] = new TH1F("jet_pt", "p_{T}(jet)" , 150 , 0. , 150.);
	   histosTH1F["jet_eta"] = new TH1F("jet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
	   histosTH1F["jet_phi"] = new TH1F("jet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);

	   histosTH1F["leadingJet_pt"] = new TH1F("leadingJet_pt", "p_{T}(jet)" , 150 , 0. , 230.);
	   histosTH1F["leadingJet_eta"] = new TH1F("leadingJet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
	   histosTH1F["leadingJet_phi"] = new TH1F("leadingJet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);
	   histosTH1F["leadingJet_pt_selected"] = new TH1F("leadingJet_pt_selected", "p_{T}(jet)" , 150 , 0. , 230.);
	   histosTH1F["leadingJet_eta_selected"] = new TH1F("leadingJet_eta_selected", "#eta(jet)" , 200 , -5.2 , 5.2);

	   histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , 150 , 0. , 230.);
	   histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "#eta(jet)" , 200 , -5.2 , 5.2);
	   histosTH1F["secondJet_phi"] = new TH1F("secondJet_phi", "#phi(jet)" , 200 , -M_PI , M_PI);
	   histosTH1F["secondJet_pt_selected"] = new TH1F("secondJet_pt_selected", "p_{T}(jet)" , 150 , 0. , 230.);
	   histosTH1F["secondJet_eta_selected"] = new TH1F("secondJet_eta_selected", "#eta(jet)" , 200 , -5.2 , 5.2);

	   histosTH1F["DeltaPtJet"] = new TH1F("Delta_pt_Jet", "#Delta p_{T}(jet)" , 150 , 0. , 150.);
	   histosTH1F["DeltaEtaJet"] = new TH1F("Delta_eta_Jet", "#Delta#eta(jet)" , 200 , 0 , 5.2);
	   histosTH1F["DeltaPhiJet"] = new TH1F("Delta_phi_Jet", "#Delta#phi(jet)" , 200 , -M_PI , M_PI);
	   histosTH1F["Mass_Jet"] = new TH1F("mass_Jet", "Mass(jet)" , 200 , 0 , 450);
	 */
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

	histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" ,   30,-1.0,1.0);
	histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" ,   20,-1.0,1.0);
	histosTH1F["t_proton_plus_tbin1"] = new TH1F("t_proton_plus_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin2"] = new TH1F("t_proton_plus_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin3"] = new TH1F("t_proton_plus_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin4"] = new TH1F("t_proton_plus_tbin4", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);

	Float_t tbins[12] = {0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.24, 0.30, 0.36, 0.45, 0.65, 1.};
	Float_t xi_bins[29] = {0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15};

	//Float_t tbins[12] = {0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.24, 0.30, 0.36, 0.45, 0.65, 1.};
	//Float_t xi_bins[51] = {0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25};

	Float_t tbins_matrix[34] = { 0.0144, 0.0256, 0.04, 0.0576, 0.0784, 0.1024, 0.1296, 0.16, 0.1936, 0.2304, 0.2704, 0.3136, 0.36, 0.4096, 0.4624, 0.5184, 0.5778, 0.64, 0.7056, 0.7744, 0.8464, 0.9216, 1., 1.0816, 1.1664, 1.2544, 1.3456, 1.44, 1.5376, 1.6384, 1.7424, 1.8496, 1.96, 2.0736};
	histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" ,  30,-1.0,1.0);// 28, xi_bins);
	histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" ,   20,-1.0,1.0);//11, tbins);
	histosTH1F["t_proton_minus_accepted"] = new TH1F("t_proton_minus_accepted", "t_proton_minus",  11,-1.0,1.0);//, 11, tbins);
	histosTH1F["thx_proton_minus"] = new TH1F("thx_proton_minus", "thx_proton_minus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_minus"] = new TH1F("thy_proton_minus", "thy_proton_minus" , 200, -5e-4, 5e-4);
	histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" ,   30,-1.0,1.0);//28, xi_bins);
	histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" ,   11,-1.0,1.0);//11, tbins);
	histosTH1F["t_proton_plus_newbinning"] = new TH1F("t_proton_plus_newbinning", "t_proton_plus" , 33, tbins_matrix);
	histosTH1F["t_proton_plus_accepted"] = new TH1F("t_proton_plus_accepted", "t_proton_plus" ,   11,-1.0,1.0);//11, tbins);
	histosTH1F["t_proton_plus_accepted_newbinning"] = new TH1F("t_proton_plus_accepted_newbinning", "t_proton_plus" , 33, tbins_matrix);
	histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
	histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);


	//FIXME
	histosTH1F["xi_proton_minus_accepted"] = new TH1F("xi_proton_minus_accepted", "xi_proton_minus" ,  30,-1.0,1.0);// 28, xi_bins);
	histosTH1F["xi_proton_plus_accepted"] = new TH1F("xi_proton_plus_accepted", "xi_proton_plus" ,   30,-1.0,1.0);//28, xi_bins);

	histosTH1F["xi_proton_plus_accepted_020"] = new TH1F("xi_proton_plus_accepted_020", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_020"] = new TH1F("t_proton_plus_accepted_020", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["posx_proton_plus_accepted_020"] = new TH1F("posx_proton_plus_accepted_020", "posx_proton_plus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_plus_accepted_020"] = new TH1F("posy_proton_plus_accepted_020", "posy_proton_plus" , 200, -0.05, 0.05);
	histosTH1F["t_proton_plus_accepted_020_tbin1"] = new TH1F("t_proton_plus_accepted_020_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_020_tbin2"] = new TH1F("t_proton_plus_accepted_020_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_020_tbin3"] = new TH1F("t_proton_plus_accepted_020_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_020_tbin4"] = new TH1F("t_proton_plus_accepted_020_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_021"] = new TH1F("xi_proton_plus_accepted_021", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_021"] = new TH1F("t_proton_plus_accepted_021", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_021_tbin1"] = new TH1F("t_proton_plus_accepted_021_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_021_tbin2"] = new TH1F("t_proton_plus_accepted_021_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_021_tbin3"] = new TH1F("t_proton_plus_accepted_021_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_021_tbin4"] = new TH1F("t_proton_plus_accepted_021_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_022"] = new TH1F("xi_proton_plus_accepted_022", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_022"] = new TH1F("t_proton_plus_accepted_022", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_022_tbin1"] = new TH1F("t_proton_plus_accepted_022_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_022_tbin2"] = new TH1F("t_proton_plus_accepted_022_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_022_tbin3"] = new TH1F("t_proton_plus_accepted_022_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_022_tbin4"] = new TH1F("t_proton_plus_accepted_022_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_023"] = new TH1F("xi_proton_plus_accepted_023", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_023"] = new TH1F("t_proton_plus_accepted_023", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_023_tbin1"] = new TH1F("t_proton_plus_accepted_023_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_023_tbin2"] = new TH1F("t_proton_plus_accepted_023_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_023_tbin3"] = new TH1F("t_proton_plus_accepted_023_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_023_tbin4"] = new TH1F("t_proton_plus_accepted_023_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_024"] = new TH1F("xi_proton_plus_accepted_024", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_024"] = new TH1F("t_proton_plus_accepted_024", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_024_tbin1"] = new TH1F("t_proton_plus_accepted_024_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_024_tbin2"] = new TH1F("t_proton_plus_accepted_024_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_024_tbin3"] = new TH1F("t_proton_plus_accepted_024_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_024_tbin4"] = new TH1F("t_proton_plus_accepted_024_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_025"] = new TH1F("xi_proton_plus_accepted_025", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_025"] = new TH1F("t_proton_plus_accepted_025", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_025_tbin1"] = new TH1F("t_proton_plus_accepted_025_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_025_tbin2"] = new TH1F("t_proton_plus_accepted_025_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_025_tbin3"] = new TH1F("t_proton_plus_accepted_025_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_025_tbin4"] = new TH1F("t_proton_plus_accepted_025_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_plus_accepted_120"] = new TH1F("xi_proton_plus_accepted_120", "xi_proton_plus" , 200, 0., 1.);
	histosTH1F["t_proton_plus_accepted_120"] = new TH1F("t_proton_plus_accepted_120", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_120_tbin1"] = new TH1F("t_proton_plus_accepted_120_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_120_tbin2"] = new TH1F("t_proton_plus_accepted_120_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_120_tbin3"] = new TH1F("t_proton_plus_accepted_120_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_accepted_120_tbin4"] = new TH1F("t_proton_plus_accepted_120_tbin4", "t_proton_plus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_120"] = new TH1F("xi_proton_minus_accepted_120", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_120"] = new TH1F("t_proton_minus_accepted_120", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["posx_proton_minus_accepted_120"] = new TH1F("posx_proton_minus_accepted_120", "posx_proton_minus" , 200, -0.05, 0.05);
	histosTH1F["posy_proton_minus_accepted_120"] = new TH1F("posy_proton_minus_accepted_120", "posy_proton_minus" , 200, -0.05, 0.05);
	histosTH1F["t_proton_minus_accepted_120_tbin1"] = new TH1F("t_proton_minus_accepted_120_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_120_tbin2"] = new TH1F("t_proton_minus_accepted_120_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_120_tbin3"] = new TH1F("t_proton_minus_accepted_120_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_120_tbin4"] = new TH1F("t_proton_minus_accepted_120_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_121"] = new TH1F("xi_proton_minus_accepted_121", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_121"] = new TH1F("t_proton_minus_accepted_121", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_121_tbin1"] = new TH1F("t_proton_minus_accepted_121_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_121_tbin2"] = new TH1F("t_proton_minus_accepted_121_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_121_tbin3"] = new TH1F("t_proton_minus_accepted_121_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_121_tbin4"] = new TH1F("t_proton_minus_accepted_121_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_122"] = new TH1F("xi_proton_minus_accepted_122", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_122"] = new TH1F("t_proton_minus_accepted_122", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_122_tbin1"] = new TH1F("t_proton_minus_accepted_122_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_122_tbin2"] = new TH1F("t_proton_minus_accepted_122_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_122_tbin3"] = new TH1F("t_proton_minus_accepted_122_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_122_tbin4"] = new TH1F("t_proton_minus_accepted_122_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_123"] = new TH1F("xi_proton_minus_accepted_123", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_123"] = new TH1F("t_proton_minus_accepted_123", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_123_tbin1"] = new TH1F("t_proton_minus_accepted_123_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_123_tbin2"] = new TH1F("t_proton_minus_accepted_123_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_123_tbin3"] = new TH1F("t_proton_minus_accepted_123_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_123_tbin4"] = new TH1F("t_proton_minus_accepted_123_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_124"] = new TH1F("xi_proton_minus_accepted_124", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_124"] = new TH1F("t_proton_minus_accepted_124", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_124_tbin1"] = new TH1F("t_proton_minus_accepted_124_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_124_tbin2"] = new TH1F("t_proton_minus_accepted_124_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_124_tbin3"] = new TH1F("t_proton_minus_accepted_124_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_124_tbin4"] = new TH1F("t_proton_minus_accepted_124_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_125"] = new TH1F("xi_proton_minus_accepted_125", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_125"] = new TH1F("t_proton_minus_accepted_125", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_125_tbin1"] = new TH1F("t_proton_minus_accepted_125_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_125_tbin2"] = new TH1F("t_proton_minus_accepted_125_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_125_tbin3"] = new TH1F("t_proton_minus_accepted_125_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_125_tbin4"] = new TH1F("t_proton_minus_accepted_125_tbin4", "t_proton_minus" , 200, 0., 5.);

	histosTH1F["xi_proton_minus_accepted_020"] = new TH1F("xi_proton_minus_accepted_020", "xi_proton_minus" , 200, 0., 1.);
	histosTH1F["t_proton_minus_accepted_020"] = new TH1F("t_proton_minus_accepted_020", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_020_tbin1"] = new TH1F("t_proton_minus_accepted_020_tbin1", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_020_tbin2"] = new TH1F("t_proton_minus_accepted_020_tbin2", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_020_tbin3"] = new TH1F("t_proton_minus_accepted_020_tbin3", "t_proton_minus" , 200, 0., 5.);
	histosTH1F["t_proton_minus_accepted_020_tbin4"] = new TH1F("t_proton_minus_accepted_020_tbin4", "t_proton_minus" , 200, 0., 5.);


	map<string,TH2F*> histosTH2F;
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

	histosTH2F["rp_track_pos_y_vs_x_020"] = new TH2F("rp_track_pos_y_vs_x_020", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
	histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

	histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
	histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

	histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 200, -0.05, 0.05, 200, -0.05, 0.05);
	histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 200, -0.05, 0.05, 200, -0.05, 0.05);

	for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
		it->second->Sumw2();
	for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
		it->second->Sumw2();
	gStyle->SetPalette(1);
	//===================
	int i_tot = 0 , nevt_tot = 0;
	/*
	//   vector<TString>* vfiles = new vector<TString>(1,"merged_reduced_8372_198903_LP_Jets1_1_test_v1.root"); 
	//   vector<TString>* vfiles = new vector<TString>; 
	//   vfiles->push_back("UABaseTree_Pomwig_SD_8TeV_9_1_99h.root" ); 

	//const char *dirname="/storage1/lhuertas/Analysis/CMSTotem/MC/Pomwig_SDPlus_8TeV/";
	const char *dirname="/storage1/lhuertas/Analysis/CMSTotem/MC/Pomwig_SDPlus_8TeV/";
	//   const char *dirname="dcap://se-dcache.hepgrid.uerj.br:22138/pnfs/hepgrid.uerj.br/data/cms/store/user/lhuertas/MC/Pomwig_SD_8TeV/";
	const char *ext=".root";
	vector<TString>* vfiles = new vector<TString>; 

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) {
	TSystemFile *file;
	TString fname;
	TIter next(files);
	while ((file=(TSystemFile*)next())) {
	fname = file->GetName();
	if (!file->IsDirectory() && fname.EndsWith(ext)) {
	TString root_file = dirname + string(fname.Data());
	vfiles->push_back(root_file); cout<<root_file<<endl;      
	}
	}
	}
	 */


	vector<TString>* vfiles = new vector<TString>;
	for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );

	//Declaration of tree and its branches variables
	//  TTree* tree = new TTree(treeName.c_str(),"");
	TTree* tree = NULL;
	MyEvtId*           evtId        = NULL;
	vector<MyGenPart>* genPart      = NULL;
	vector<MyTracks>*  track_coll   = NULL;
	vector<MyVertex>*  vertex_coll  = NULL;
	vector<MyPFJet>*   pfJet_coll   = NULL;
	vector<MyMuon>*    muon_coll   = NULL;
	vector<MyGenJet>*   genJet_coll   = NULL;
	vector<MyPFCand>*  pFlow_coll   = NULL;
	MyGenKin*  genKin   = NULL;
	//===================

	rp_aperture_config();

	int n_vertices_selected = 0;
	int n_select_Vertex_After_vtx_cut =0;
	int n_tracks_selected = 0;
	int n_muon_selected = 0;
	int n_dimuons_selected = 0;
	int n_dimuons_rp_selected_plus_side = 0;
	int n_dimuons_rp_accept_plus_side = 0;
	int n_dimuons_rp_selected_minus_side = 0;
	int n_dimuons_rp_accept_minus_side = 0;
	//int n_dimuons_t_selected_proton_left = 0;
	//int n_dimuons_t_at_xi_selected_proton_left = 0;
	//int n_dimuons_t_selected_proton_right = 0;
	//int n_dimuons_t_at_xi_selected_proton_right = 0;
	//double nevents_jets = 0; 
	double nevents_pf = 0; 
	double nevents_gen = 0; 
	double nevents_total = 0; 
	//double events_jets = 0;
	//double events_pf = 0;
	double events_gen = 0;

	double nweight_total = 0; 
	double weight_total_leadingJet = 0; 
	double weight_total_secondJet = 0; 
	double weight_total_leadingJet_selected = 0; 
	double weight_total_secondJet_selected = 0; 
	double weight_total_Jet_selected = 0; 
	double weight_total_PF_selected = 0; 
	//starting Loop over files, stops at end of list of files or when reached nevt_max
	for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){

		TFile* file = TFile::Open(*itfiles,"READ");

		//getting the tree form the current file
		tree = (TTree*) file->Get( treeName.c_str() );

		//Getting number of events
		int nev = int(tree->GetEntriesFast());
		nevt_tot += nev;
		cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

		//adding branches to the tree ----------------------------------------------------------------------
		tree->SetBranchAddress("evtId",&evtId);
		tree->SetBranchAddress("generalTracks",&track_coll); 
		tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
		tree->SetBranchAddress("muons",&muon_coll);
		tree->SetBranchAddress("particleFlow",&pFlow_coll);
		tree->SetBranchAddress("genKin",&genKin);
		tree->SetBranchAddress("genPart",&genPart);


		//starting loop over events, stops when reached end of file or nevt_max
		for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

			//printing the % of events done every 10k evts
			if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

			//Filling the variables defined setting branches
			tree->GetEntry(i_evt);

			//bool passedHLT = false;
			bool passedvtx = false;
			//bool jet1_selected = false;
			//bool jet2_selected = false;
			bool pz_proton_max = false;
			bool PF_eta_max = false;
			bool PF_eta_min = false;
			bool xi_negat_gen = false;
			bool xi_posit_gen = false;

			//AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
			// double event_weight = genKin->genWeight; 
			double event_weight = 1.0;
			nweight_total += event_weight; 
			++nevents_total;

			//-------------------------------------------------------------------------------------------------
			//filling pt distribution for the generated particles
			//ie those from pythia generator, without reconstruction
			//	if(isMC){
			//	  for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
			//	  pt_gen->Fill(p->Pt());
			//	  }

			//-------------------------------------------------------------------------------------------------

			// Vertices
			for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
				int idx_vtx = it_vtx - vertex_coll->begin();

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

			//histosTH1F["EventSelection"]->Fill( "Vertex", event_weight );

			int prim_vtx_id = primaryVertex.id;


			// Tracks
			//int n_tracks_selected = 0;
			for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
				histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
				histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
				histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

				++n_tracks_selected;
			}
			histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );

			//include Muon and Dimuons
			//first muon	
			double  muon1_mass  = 0.;
			double muons1_pt = 0.;
			double muons1_eta = 0.;
			double muons1_rapidity = 0.; 
			//second  muon    
			double  muon2_mass  = 0.;
			double muons2_pt = 0.;
			double muons2_eta = 0.;
			double muons2_rapidity = 0.;   
			// dimuons
			double  dimuons_mass  = 0.;
			double dimuons_pt = 0.;
			double dimuons_eta = 0.;
			double dimuons_rapidity = 0.;	
			bool dimuons_true;




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
			// histosTH1F["EventSelection"]->Fill( "Muons", event_weight );
			//------------------
			// Muon pairs
			for(vector<MyMuon>::iterator it_mu1 = muons_selected.begin() ;
					it_mu1 != muons_selected.end() ; ++it_mu1){
				for(vector<MyMuon>::iterator it_mu2 = muons_selected.begin() ;
						it_mu2 != muons_selected.end() ; ++it_mu2){
					bool os_muons = ( it_mu1->charge*it_mu2->charge < 0. );

					if( !os_muons ) continue;
					++n_dimuons_selected;
					dimuons_true = true;
					//...
					TLorentzVector& muon1_lorentz = *it_mu1;
					TLorentzVector& muon2_lorentz = *it_mu2;
					TLorentzVector dimuon_lorentz(0.,0.,0.,0.);
					dimuon_lorentz += muon1_lorentz;
					dimuon_lorentz += muon2_lorentz;
					dimuons_mass  = dimuon_lorentz.M();
					dimuons_pt = dimuon_lorentz.Pt();
					dimuons_eta = dimuon_lorentz.Eta();
					dimuons_rapidity = dimuon_lorentz.Rapidity();			
					//                                        histosTH1F["dimuon_multiplicity"]->Fill(n_dimuons_selected, event_weight );
					histosTH1F["dimuon_mass"]->Fill( dimuons_mass, event_weight );
					histosTH1F["dimuon_pt"]->Fill( dimuons_pt, event_weight );
					histosTH1F["dimuon_eta"]->Fill( dimuons_eta, event_weight );
					histosTH1F["dimuon_rapidity"]->Fill( dimuons_rapidity, event_weight );



				} //Second Muons 
			}// First Muon


			if( !dimuons_true ) continue;
//			cout << dimuons_true<<endl;





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
				if(verbose) cout<<"Particle ID: " << id << " energy: "<< energy_gen <<" px_gen: " << px_gen << " py_gen: "<< py_gen << " pz_gen: "<< pz_gen<<" proton_pf: "<< proton_pf << endl; 
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

				if(verbose)  cout<< " xi_proton_plus: "<< xi_proton_plus<< " t_proton_plus:"<< t_proton_plus<< endl;

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

				if(verbose)  cout<< " xi_proton_minus: "<< xi_proton_minus << " t_proton_minus: "<< t_proton_minus << endl;		
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


			histosTH1F["correction"]->Fill( correction, event_weight ); 
			histosTH1F["resolution_before"]->Fill( resolution_before, event_weight ); 
			histosTH1F["resolution_after"]->Fill( resolution_after, event_weight ); 
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
			//histosTH2F["xi+Vseta_max"]->Fill( eta_max, xi1, event_weight ); 

			//rp_accept
			if (verbose) cout<< " |t| selection  : "<< " | "<< t_proton_down_ << " - " << t_proton_up_<< " |"<< endl;

			if (fabs(t_proton_plus)<t_proton_down_ || fabs(t_proton_plus)>t_proton_up_)continue;
			if(verbose) cout<< " rp accept after selection |t| : "<< " |t| proton plus "<< fabs(t_proton_plus)<< " xi proton plus "<< xi_proton_plus << endl;
			if (xi_proton_plus<xi_proton_down_ || xi_proton_plus>xi_proton_up_)continue; 

			if (verbose)  cout<< " rp accept after selection xi: "<< " |t| proton plus "<< fabs(t_proton_plus)<< " xi proton plus "<< xi_proton_plus << endl;

			histosTH1F["xi_proton_plus"]->Fill( xi_proton_plus , event_weight );
			histosTH1F["t_proton_plus"]->Fill( fabs(t_proton_plus) , event_weight );
			histosTH1F["t_proton_plus_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );
			histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
			histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );
			histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );

			bool proton_plus_rp_accept = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 );


			if(!proton_plus_rp_accept){// continue;
				// xi rp accepted 
				histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
				histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
				histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight );
				histosTH1F["t_proton_plus_accepted_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );

				// xi reco  plots
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

				histosTH1F["dimuon_mass_rp_accept"]->Fill( dimuons_mass, event_weight );
				histosTH1F["dimuon_pt_rp_accept"]->Fill( dimuons_pt, event_weight );
				histosTH1F["dimuon_eta_rp_accept"]->Fill( dimuons_eta, event_weight );
				histosTH1F["dimuon_rapidity_rp_accept"]->Fill( dimuons_rapidity, event_weight );

				//----------------------------------------------
			}
			//------------------


		}//end loop for events

		cout<<"Total of evts="<< nev << endl << *itfiles << endl;
		cout<<"n_vertices_selected ="<< n_vertices_selected << endl;
		cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
		cout<<"n_tracks_selected="<< n_tracks_selected<<endl;
		cout<<"n_muon_selected="<< n_muon_selected<<endl;
		cout<<"n_dimuons_selected="<< n_dimuons_selected<<endl;
		cout<<"n_dimuons_rp_accept_plus_side="<< n_dimuons_rp_accept_plus_side <<endl;
		//cout<<"n_dimuons_t_at_xi_selected_proton_right="<< n_dimuons_t_at_xi_selected_proton_right<<endl;



		file->Close();

	}//end of loop over files

	//output file
	TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
	output->cd();
	histosTH1F["Eta_max"]->GetXaxis()->SetTitle("#eta^{max}");
	histosTH1F["Eta_min"]->GetXaxis()->SetTitle("#eta^{min}");
	histosTH1F["Delta_eta_maxmin"]->GetXaxis()->SetTitle("#eta^{max}-#eta^{min}");
	histosTH1F["xi_plus_Reco"]->GetXaxis()->SetTitle("#xi^{+}");
	histosTH1F["xi_minus_Reco"]->GetXaxis()->SetTitle("#xi^{-}");
	histosTH2F["xi_plus_reco_gen"]->GetXaxis()->SetTitle("#xi^{+}_{gen}");
	histosTH2F["xi_plus_reco_gen"]->GetYaxis()->SetTitle("#xi^{+}_{reconst}");
	histosTH2F["xi_minus_reco_gen"]->GetXaxis()->SetTitle("#xi^{-}_{gen}");
	histosTH2F["xi_minus_reco_gen"]->GetYaxis()->SetTitle("#xi^{-}_{reconst}");
	histosTH2F["logxi_plus_reco_gen"]->GetXaxis()->SetTitle("Log_{10}#xi^{+}_{gen}");
	histosTH2F["logxi_plus_reco_gen"]->GetYaxis()->SetTitle("Log_{10}#xi^{+}_{reconst}");
	histosTH2F["logxi_minus_reco_gen"]->GetXaxis()->SetTitle("Log_{10}#xi^{-}_{gen}");
	histosTH2F["logxi_minus_reco_gen"]->GetYaxis()->SetTitle("Log_{10}#xi^{-}_{reconst}");
	//	histosTH1F["leadingJet_pt"]->GetXaxis()->SetTitle("p_T^{jet1}(GeV)");
	//	histosTH1F["secondJet_pt"]->GetXaxis()->SetTitle("#eta^{jet1}");
	//	histosTH1F["secondJet_pt"]->GetXaxis()->SetTitle("p_T^{jet2}(GeV)");
	//	histosTH1F["secondJet_eta"]->GetXaxis()->SetTitle("#eta^{jet2}");
	//	histosTH1F["DeltaPtJet"]->GetXaxis()->SetTitle("#Delta p_T^{jets}(GeV)");
	//	histosTH1F["DeltaEtaJet"]->GetXaxis()->SetTitle("#Delta#eta^{jets}");
	//	histosTH1F["DeltaPhiJet"]->GetXaxis()->SetTitle("#Delta#phi^{jets}");
	/* 
	   float cross_section = 2.1189e7; //pb cross section for pomwig
	   float luminity_HLT_L1Jet1_198902 = 0.015879;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
	   float luminity_HLT_L1Jet2_198902 = 0.015879;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
	   float luminity_HLT_L1Jet1_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
	   float luminity_HLT_L1Jet2_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
	   float luminity = luminity_HLT_L1Jet1_198902 + luminity_HLT_L1Jet2_198902 + luminity_HLT_L1Jet1_198903 + luminity_HLT_L1Jet2_198903;

	   float n_events = luminity*cross_section;

	   float f1 = (float) nevents_total;
	//   float f1 = (float) nweight_total;
	Double_t scale = n_events/f1;

	//   float f3 = (float) weight_total_PF_selected ; 
	//   Double_t scale_PF = 1.0/f3;

	histosTH1F["leadingJet_pt"]->Scale(scale);
	histosTH1F["leadingJet_phi"]->Scale(scale);
	histosTH1F["leadingJet_eta"]->Scale(scale);
	histosTH1F["leadingJet_pt_selected"]->Scale(scale);
	histosTH1F["leadingJet_eta_selected"]->Scale(scale);
	histosTH1F["secondJet_pt"]->Scale(scale);
	histosTH1F["secondJet_phi"]->Scale(scale);
	histosTH1F["secondJet_eta"]->Scale(scale);
	histosTH1F["secondJet_pt_selected"]->Scale(scale);
	histosTH1F["secondJet_eta_selected"]->Scale(scale);
	histosTH1F["Eta_max"]->Scale(scale);		 
	histosTH1F["Eta_min"]->Scale(scale);
	histosTH1F["Delta_eta_maxmin"]->Scale(scale);
	histosTH1F["xi_plus_Reco"]->Scale(scale);
	histosTH1F["xi_minus_Reco"]->Scale(scale);
	histosTH1F["logxi_plus"]->Scale(scale);
	histosTH1F["DeltaPtJet"]->Scale(scale);
	histosTH1F["DeltaEtaJet"]->Scale(scale);
	histosTH1F["DeltaPhiJet"]->Scale(scale);
	histosTH1F["Mass_Jet"]->Scale(scale);

	histosTH1F["xi_proton_plus_accepted"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted"]->Scale(scale);
	histosTH2F["proton_plus_xi_vs_t_accepted"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin1"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin2"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin3"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin4"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin5"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin6"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin7"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin8"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin9"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin10"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin11"]->Scale(scale);
	histosTH1F["t_proton_plus_accepted_tbin12"]->Scale(scale);
	histosTH1F["xi_proton_minus_accepted"]->Scale(scale);
	histosTH1F["t_proton_minus_accepted"]->Scale(scale);
	histosTH2F["proton_minus_xi_vs_t_accepted"]->Scale(scale);
	histosTH1F["t_proton_minus"]->Scale(scale);
	 */
	histosTH2F["xi_plus_reco_gen"]->SetOption("colz");
	histosTH2F["xi_minus_reco_gen"]->SetOption("colz");
	histosTH2F["logxi_plus_reco_gen"]->SetOption("colz");
	histosTH2F["logxi_minus_reco_gen"]->SetOption("colz");

	for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
			it_histo != histosTH1F.end(); ++it_histo){
		//   (*it_histo).second->Scale(scale);
		(*it_histo).second->Write();}
	for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
			it_histo != histosTH2F.end(); ++it_histo)
		(*it_histo).second->Write();

	output->Close();
}
