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

#include "deltaPhi.h"
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
void CMS_MC_POMPYT_JPSI_plus_v4(vector<string> const& fileNames, string const& outputFileName = "output.root", const Double_t t_proton_down_=0.0, const Double_t t_proton_up_=1.0, const Double_t xi_proton_down_=0.0,const Double_t xi_proton_up_=0.2,const Int_t Bin_mass=100 ,const Int_t nevt_max = -1){


	//bool isMC  = false;
	bool verbose = false;
	bool selectVertex = true;	
	bool selectMuons = true;
	string treeName = "evt";//"cms_totem";

	const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


	// Declaration of histograms
	map<string,TH1F*> histosTH1F;


	//Vertex plots
	histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 100 , -30. , 30.);
	histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 100 , -1.5 , 1.5);
	histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 100 , -1.5 , 1.5);
	histosTH1F["vtx_ndof"] = new TH1F("vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["vtx_chi2"] = new TH1F("vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);

	histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertices" , 30 , 0 , 30);
	histosTH1F["vertex_multiplicity_after_vtx_sel"] = new TH1F("vertex_multiplicity_after_vtx_sel", "n vertices after vtx sel" , 30 , 0 , 30);
	histosTH1F["prim_vtx_zpos"] = new TH1F("prim_vtx_zpos", "z(vtx)" , 100 , -30. , 30.);
	histosTH1F["prim_vtx_xpos"] = new TH1F("prim_vtx_xpos", "x(vtx)" , 100 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos"] = new TH1F("prim_vtx_ypos", "y(vtx)" , 100 , -1.5 , 1.5);

	histosTH1F["prim_vtx_ndof"] = new TH1F("prim_vtx_ndof", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["prim_vtx_chi2"] = new TH1F("prim_vtx_chi2", "chi2(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_chi2n"] = new TH1F("prim_vtx_chi2n", "chi2n(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_ntracks"] = new TH1F("prim_vtx_ntracks", "n_{trk}(vtx)" , 30 , 0 , 30);
	histosTH1F["prim_vtx_sumpt"] = new TH1F("prim_vtx_sumpt", "sum(p_{T})(vtx)" , 100 , 0. , 100.);

	histosTH1F["prim_vtx_zpos_after_vtx_sel"] = new TH1F("prim_vtx_zpos_after_vtx_sel", "z(vtx)" , 100 , -30. , 30.);
	histosTH1F["prim_vtx_xpos_after_vtx_sel"] = new TH1F("prim_vtx_xpos_after_vtx_sel", "x(vtx)" , 100 , -1.5 , 1.5);
	histosTH1F["prim_vtx_ypos_after_vtx_sel"] = new TH1F("prim_vtx_ypos_after_vtx_sel", "y(vtx)" , 100 , -1.5 , 1.5);

	histosTH1F["prim_vtx_ndof_after_vtx_sel"] = new TH1F("prim_vtx_ndof_after_vtx_sel", "ndof(vtx)" , 100 , 0. , 15.);
	histosTH1F["prim_vtx_chi2_after_vtx_sel"] = new TH1F("prim_vtx_chi2_after_vtx_sel", "chi2(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_chi2n_after_vtx_sel"] = new TH1F("prim_vtx_chi2n_after_vtx_sel", "chi2n(vtx)" , 100 , 0. , 10.);
	histosTH1F["prim_vtx_ntracks_after_vtx_sel"] = new TH1F("prim_vtx_ntracks_after_vtx_sel", "n_{trk}(vtx)" , 30 , 0 , 30);
	histosTH1F["prim_vtx_sumpt_after_vtx_sel"] = new TH1F("prim_vtx_sumpt_after_vtx_sel", "sum(p_{T})(vtx)" , 100 , 0. , 100.);	

	//MC Gen 

        histosTH1F["t_proton_minus_before"] = new TH1F("t_proton_minus_before", "|t|^{-}" , 100 , 0.,4.);
        histosTH1F["xi_proton_minus_before"] = new TH1F("xi_proton_minus_before", "#xi^{-}" , 100 , 0.,1.);
        histosTH1F["xi_proton_plus_before"] = new TH1F("xi_proton_plus_before", "xi_proton_plus" ,   100,0.,1.0);//28, xi_bins);
        histosTH1F["t_proton_plus_before"] = new TH1F("t_proton_plus_before", "t_proton_plus" ,   100,0.,4.0);//11, tbins);
        histosTH1F["xi_proton_plus_after"] = new TH1F("xi_proton_plus_after", "xi_proton_plus" ,   100,0,1.0);//28, xi_bins);
        histosTH1F["t_proton_plus_after"] = new TH1F("t_proton_plus_after", "t_proton_plus" ,   100,0.,4.0);//11, tbins);
	//
	histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 100 , 0. , 15.);
	histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 100 , -5.2 , 5.2);
	histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 100 , -M_PI , M_PI);
	histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);

	// Reco Muons plots
	histosTH1F["muon_pt"] = new TH1F("muon_pt", "p_{T}(muon)" , 100 , 0. , 100.);
	histosTH1F["muon_eta"] = new TH1F("muon_eta", "#eta(muon)" , 100 , -5.2 , 5.2);
	histosTH1F["muon_phi"] = new TH1F("muon_phi", "#phi(muon)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muon_multiplicity"] = new TH1F("muon_multiplicity", "n muons" , 100 , 0 , 100);

	histosTH1F["dimuon_mass"] = new TH1F("dimuon_mass", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
	histosTH1F["dimuon_pt"] = new TH1F("dimuon_pt", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
	histosTH1F["dimuon_pt2"] = new TH1F("dimuon_pt2", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
	histosTH1F["dimuon_eta"] = new TH1F("dimuon_eta", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["dimuon_rapidity"] = new TH1F("dimuon_rapidity", "y(mu1,mu2)" , 100 , -15. , 15.);
	histosTH1F["dimuon_multiplicity"] = new TH1F("dimuon_multiplicity", "n dimuons" , 100 , 0 , 100);
        histosTH1F["dimuon_xi_plus"]= new TH1F("dimuon_xi_plus", "#xi^{+}" , 100 , 0.,1.);
        histosTH1F["dimuon_t_plus"] = new TH1F("dimuon_t_plus", "|t|^{+}" , 100,0.,4.0);

        // Reco Muons rp_accept
        histosTH1F["dimuon_mass_rp_accept"] = new TH1F("dimuon_mass_rp_accept", "mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["dimuon_pt_rp_accept"] = new TH1F("dimuon_pt_rp_accept", "p_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["dimuon_pt2_rp_accept"] = new TH1F("dimuon_pt2_rp_accept", "p_{T}2(mu1,mu2)" , 100 , 0. , 1000.);
        histosTH1F["dimuon_eta_rp_accept"] = new TH1F("dimuon_eta_rp_accept", "#eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["dimuon_rapidity_rp_accept"] = new TH1F("dimuon_rapidity_rp_accept", "y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["dimuon_multiplicity_rp_accept"] = new TH1F("dimuon_multiplicity_rp_accept", "n dimuons" , 100 , 0 , 100);
        histosTH1F["muonDeltaPt_rp_accept"] = new TH1F("muonDeltaPt_rp_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_rp_accept"] = new TH1F("muonDeltaEta_rp_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
        histosTH1F["muonDeltaPhi_rp_accept"] = new TH1F("muonDeltaPhi_rp_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muonDphi_rp_accept"] = new TH1F("muonDphi_rp_accept", "#Delta#phi(mu1,mu2)" , 100 , -1.2*M_PI , 1.2*M_PI);
        histosTH1F["muonDeltaY_rp_accept"] = new TH1F("muonDeltaY_rp_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["dimuon_xi_plus_rp_accept"]= new TH1F("dimuon_xi_plus_rp_accept", "#xi^{+}" , 100 , 0.,1.);
        histosTH1F["dimuon_t_plus_rp_accept"] = new TH1F("dimuon_t_plus_rp_accept", "|t|^{+}" , 100,0.,4.0);

        //Muon1 info
        histosTH1F["muon1_pt"] = new TH1F("muon1_pt", "p_{T}(mu1)" , 100 , 0. , 100.);
	histosTH1F["muon1_eta"] = new TH1F("muon1_eta", "#eta(mu1)" , 100 , -5.2 , 5.2);
	histosTH1F["muon1_phi"] = new TH1F("muon1_phi", "#phi(muon1)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muon1_rapidity"] = new TH1F("muon1_rapidity", "y(mu1)" , 100 , -15. , 15.);
        //Muon2 info
	histosTH1F["muon2_pt"] = new TH1F("muon2_pt", "p_{T}(mu2)" , 100 , 0. , 100.);
	histosTH1F["muon2_eta"] = new TH1F("muon2_eta", "#eta(mu2)" , 100 , -5.2 , 5.2);
	histosTH1F["muon2_phi"] = new TH1F("muon2_phi", "#phi(muon2)" , 100 , -1.2*M_PI , 1.2*M_PI);
	histosTH1F["muon2_rapidity"] = new TH1F("muon2_rapidity", "y(mu2)" , 100 , -15. , 15.);
	//Deltas info TEST
        histosTH1F["muonDeltaPt_test"] = new TH1F("muonDeltaPt_test", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_test"] = new TH1F("muonDeltaEta_test", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
        histosTH1F["muonDeltaPhi_test"] = new TH1F("muonDeltaPhi_test", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_test"] = new TH1F("muonDeltaY_test", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
        //Deltas info
        histosTH1F["muonDeltaPt"] = new TH1F("muonDeltaPt", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta"] = new TH1F("muonDeltaEta", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);	
        histosTH1F["muonDeltaPhi"] = new TH1F("muonDeltaPhi", "#Delta#phi(mu1,mu2)" , 100 , -2*M_PI , 2*M_PI);
        histosTH1F["muonDphi"] = new TH1F("muonDphi", "#Delta#phi(mu1,mu2)" , 100 , -2*M_PI , 2*M_PI);
        histosTH1F["muonDeltaY"] = new TH1F("muonDeltaY", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);
        //jpsi mass selection
        histosTH1F["jpsi_mass"] = new TH1F("jpsi_mass", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["jpsi_pt"] = new TH1F("jpsi_pt", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
        histosTH1F["jpsi_pt2"] = new TH1F("jpsi_pt2", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
        histosTH1F["jpsi_eta"] = new TH1F("jpsi_eta", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["jpsi_rapidity"] = new TH1F("jpsi_rapidity", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["jpsi_multiplicity"] = new TH1F("jpsi_multiplicity", "n jpsi(dimuons)" , 100 , 0 , 100);
        histosTH1F["jpsi_xi_plus"]= new TH1F("jpsi_xi_plus", "#xi^{+}" , 100 , 0.,1.);
        histosTH1F["jpsi_t_plus"] = new TH1F("jpsi_t_plus", "|t|^{+}" , 100,0.,4.0);

        //Deltas info
        histosTH1F["muonDeltaPt_jpsi"] = new TH1F("muonDeltaPt_jpsi", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_jpsi"] = new TH1F("muonDeltaEta_jpsi", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["muonDeltaPhi_jpsi"] = new TH1F("muonDeltaPhi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDphi_jpsi"] = new TH1F("muonDphi_jpsi", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_jpsi"] = new TH1F("muonDeltaY_jpsi", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);

        //jpsi mass selection
        histosTH1F["jpsi_mass_rp_accept"] = new TH1F("jpsi_mass_rp_accept", "jpsi_mass(mu1,mu2)" , Bin_mass , 0. , 10.);
        histosTH1F["jpsi_pt_rp_accept"] = new TH1F("jpsi_pt_rp_accept", "jpsi p_{T}(mu1,mu2)" , 100 , 0. , 20.);
        histosTH1F["jpsi_pt2_rp_accept"] = new TH1F("jpsi_pt2_rp_accept", "jpsi p_{T}^{2}(mu1,mu2)" , 100 , 0. , 1000.);
        histosTH1F["jpsi_eta_rp_accept"] = new TH1F("jpsi_eta_rp_accept", "jpsi #eta(mu1,mu2)" , 100 , -5.2 , 5.2);
        histosTH1F["jpsi_rapidity_rp_accept"] = new TH1F("jpsi_rapidity_rp_accept", "jpsi y(mu1,mu2)" , 100 , -15. , 15.);
        histosTH1F["jpsi_multiplicity_rp_accept"] = new TH1F("jpsi_multiplicity_rp_accept", "n jpsi(dimuons)" , 100 , 0 , 100);
        histosTH1F["jpsi_xi_plus_rp_accept"]= new TH1F("jpsi_xi_plus_rp_accept", "#xi^{+}" , 100 , 0.,1.);
        histosTH1F["jpsi_t_plus_rp_accept"] = new TH1F("jpsi_t_plus_rp_accept", "|t|^{+}" , 100,0.,4.0);

        //Deltas info
        histosTH1F["muonDeltaPt_jpsi_rp_accept"] = new TH1F("muonDeltaPt_jpsi_rp_accept", "#Deltap_{T}(mu1,mu2)" , 100 , 0. , 100.);
        histosTH1F["muonDeltaEta_jpsi_rp_accept"] = new TH1F("muonDeltaEta_jpsi_rp_accept", "#Delta#eta(mu1,mu2)" , 100 , 0. , 10.);
        histosTH1F["muonDeltaPhi_jpsi_rp_accept"] = new TH1F("muonDeltaPhi_jpsi_rp_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDphi_jpsi_rp_accept"] = new TH1F("muonDphi_jpsi_rp_accept", "#Delta#phi(mu1,mu2)" , 100 , -2.2*M_PI , 2.2*M_PI);
        histosTH1F["muonDeltaY_jpsi_rp_accept"] = new TH1F("muonDeltaY_jpsi_rp_accept", "#Deltay(mu1,mu2)" , 100 , 0. , 10.);


	histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
	histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 50 , 0,11);
	histosTH1F["xi_plus_Reco"] = new TH1F("xi+_Reco", "#xi^{+}" , 82 , 0,4);
	histosTH1F["xi_minus_Reco"] = new TH1F("xi-_Reco", "#xi^{-}" , 82 , 0,4);
	histosTH1F["logxi_plus"] = new TH1F("logxi+", "Log #xi^{+}" , 82 , -3,0.5);
	histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 82 , -3,0.5);
	histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 82 , -3,0.5);
	histosTH1F["correction"] = new TH1F("correction", "Correction factor" , 82 , 0,2);
	histosTH1F["resolution_after"] = new TH1F("resolution_after", "Resolution" , 82 , -2,2);
	histosTH1F["resolution_before"] = new TH1F("resolution_before", "Resolution" , 82 , -2,2);

	histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
	histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);

	//histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" ,   30,-1.0,1.0);
	//histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" ,   20,-1.0,1.0);
	histosTH1F["t_proton_plus_tbin1"] = new TH1F("t_proton_plus_tbin1", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin2"] = new TH1F("t_proton_plus_tbin2", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin3"] = new TH1F("t_proton_plus_tbin3", "t_proton_plus" , 200, 0., 5.);
	histosTH1F["t_proton_plus_tbin4"] = new TH1F("t_proton_plus_tbin4", "t_proton_plus" , 200, 0., 5.);
	//histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
	//histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);


Float_t tbins_matrix[34] = { 0.0144, 0.0256, 0.04, 0.0576, 0.0784, 0.1024, 0.1296, 0.16, 0.1936, 0.2304, 0.2704, 0.3136, 0.36, 0.4096, 0.4624, 0.5184, 0.5778, 0.64, 0.7056, 0.7744, 0.8464, 0.9216, 1., 1.0816, 1.1664, 1.2544, 1.3456, 1.44, 1.5376, 1.6384, 1.7424, 1.8496, 1.96, 2.0736};
//histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" ,  30,-1.0,1.0);// 28, xi_bins);
//histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" ,   20,-1.0,1.0);//11, tbins);
histosTH1F["t_proton_minus_accepted"] = new TH1F("t_proton_minus_accepted", "t_proton_minus",  11,-1.0,1.0);//, 11, tbins);
histosTH1F["thx_proton_minus"] = new TH1F("thx_proton_minus", "thx_proton_minus" , 200, -5e-4, 5e-4);
histosTH1F["thy_proton_minus"] = new TH1F("thy_proton_minus", "thy_proton_minus" , 200, -5e-4, 5e-4);
//histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" ,   30,-1.0,1.0);//28, xi_bins);
//histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" ,   11,-1.0,1.0);//11, tbins);
histosTH1F["t_proton_plus_newbinning"] = new TH1F("t_proton_plus_newbinning", "t_proton_plus" , 33, tbins_matrix);
histosTH1F["t_proton_plus_accepted"] = new TH1F("t_proton_plus_accepted", "t_proton_plus" , 30,0.,1.0);//11, tbins);
histosTH1F["t_proton_plus_accepted_newbinning"] = new TH1F("t_proton_plus_accepted_newbinning", "t_proton_plus" , 33, tbins_matrix);
histosTH1F["thx_proton_plus"] = new TH1F("thx_proton_plus", "thx_proton_plus" , 200, -5e-4, 5e-4);
histosTH1F["thy_proton_plus"] = new TH1F("thy_proton_plus", "thy_proton_plus" , 200, -5e-4, 5e-4);
//histosTH1F["xi_plus_gen"] = new TH1F("xi+_gen", " #xi_{+}^{gen}" , 82 , 0.,3.0);
//histosTH1F["xi_minus_gen"] = new TH1F("xi-_gen", " #xi_{-}^{gen}" , 82 , 0.,3.0);


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


vector<TString>* vfiles = new vector<TString>;
for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );

//Declaration of tree and its branches variables
//  TTree* tree = new TTree(treeName.c_str(),"");
TTree* tree = NULL;
MyEvtId*           evtId        = NULL;
vector<MyGenPart>* genPart      = NULL;
vector<MyTracks>*  track_coll   = NULL;
vector<MyVertex>*  vertex_coll  = NULL;
//vector<MyPFJet>*   pfJet_coll   = NULL;
vector<MyMuon>*    muon_coll   = NULL;
//vector<MyGenJet>*   genJet_coll   = NULL;
vector<MyPFCand>*  pFlow_coll   = NULL;
MyGenKin*  genKin   = NULL;
//===================
/*TTree t1("t1","a simple Tree with simple variables");
Float_t eta_max, eta_min;
t1.Branch("eta_max",&eta_max,"eta_max/F");*/
rp_aperture_config();

int n_vertices_selected = 0;
int n_select_Vertex_After_vtx_cut =0;
int n_tracks_selected = 0;
int n_muon_selected = 0;
int n_dimuons_selected = 0;
int n_dimuons_rp_selected_plus_side = 0;
//int n_dimuons_rp_accept_plus_side = 0;
int n_jpsi_selected = 0;
int n_jpsi_selected_rp_accept_plus_side = 0;
double nevents_pf = 0; 
double nevents_gen = 0; 
double nevents_total = 0; 
double n_muon_selected_size = 0;
//double events_gen = 0;
double nweight_total = 0; 
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
		//bool passedvtx = false;
		bool mu1_selected = false;
		bool mu2_selected = false;
                bool chmuons = false;
		//bool pz_proton_max = false;
		bool PF_eta_max = false;
		bool PF_eta_min = false;
		//bool xi_negat_gen = false;
		//bool xi_posit_gen = false;

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
			//int idx_vtx = it_vtx - vertex_coll->begin();

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
	 	  //histosTH1F["EventSelection"]->Fill( "Muons", event_weight );
                  ++n_muon_selected_size;
                  if(verbose)cout<<"muons selected size (>=2) :"<<muons_selected.size()<<endl;
                //===============================================================
                //Muon pairs
                  double chmu1 = 0.; double chmu2 = 0.;
                  double phimu1 = 0.;double ptmu1 = 0.;double etamu1 = 0.;double ymu1 = 0.;
                  double phimu2 = 0.;double ptmu2 = 0.;double etamu2 = 0.;double ymu2 = 0.;
                  double deltaphi = 0.;double deltaeta = 0.;double deltapt = 0.;double deltay = 0.;double Dphi = 0.;
                  double dimuon_mass = 0.; double dimuon_pt=0.; double dimuon_pt2=0.;double dimuon_eta =0.;double dimuon_rapidity = 0.;
                  double jpsi_mass = 0.; double jpsi_pt=0.; double jpsi_pt2=0.;double jpsi_eta =0.;double jpsi_rapidity = 0.;
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
                                        histosTH1F["muon2_pt"]->Fill( it_mu2->Pt(), event_weight );
					histosTH1F["muon2_eta"]->Fill( it_mu2->Eta(), event_weight );
					histosTH1F["muon2_rapidity"]->Fill( it_mu2->Rapidity(), event_weight );
					histosTH1F["muon2_phi"]->Fill(it_mu2->Phi(), event_weight );
                                        phimu2 = it_mu2->Phi();ptmu2 = it_mu2->Pt();
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
                                   					
                                        deltapt = fabs(muon1_lorentz.Pt() - muon2_lorentz.Pt());
                                        deltaeta = fabs(muon1_lorentz.Eta() - muon2_lorentz.Eta());
                                        deltaphi = fabs(phimu1 - phimu2);
                                        if(deltaphi > M_PI)deltaphi = (2*M_PI - deltaphi);
                                        deltay = fabs(muon1_lorentz.Rapidity() - muon2_lorentz.Rapidity());
                                        Dphi = std::fabs(deltaPhi(phimu1 ,phimu2));  
                                      
                                        histosTH1F["muonDeltaPt"]->Fill(deltapt, event_weight );
                                        histosTH1F["muonDeltaEta"]->Fill(deltaeta, event_weight );	
                                        histosTH1F["muonDeltaPhi"]->Fill(deltaphi, event_weight );
                                        histosTH1F["muonDeltaY"]->Fill(deltay, event_weight ); 
                                        histosTH1F["muonDphi"]->Fill(Dphi, event_weight );
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        double Jpsi_mass = dimuon_lorentz.M();
                                        //cout<<"jpsi_mass_dimuon = "<< Jpsi_mass << endl;                              
                                        if(((Jpsi_mass > 3.0) && (Jpsi_mass < 3.2))){
                                           //cout<<"jpsi_mass = "<< jpsi_mass << endl;
                                           ++n_jpsi_selected;
                                        
                                          //histosTH1F["jpsi_multiplicity"]->Fill(n_dimuons_selected, event_weight );
                                          histosTH1F["jpsi_mass"]->Fill( dimuon_lorentz.M(), event_weight );
                                          histosTH1F["jpsi_pt"]->Fill( dimuon_lorentz.Pt(), event_weight );
                                          histosTH1F["jpsi_pt2"]->Fill( (dimuon_lorentz.Pt()*dimuon_lorentz.Pt()), event_weight );
                                          histosTH1F["jpsi_eta"]->Fill( dimuon_lorentz.Eta(), event_weight );
                                          histosTH1F["jpsi_rapidity"]->Fill( dimuon_lorentz.Rapidity(), event_weight );
					  histosTH1F["muonDeltaPt_jpsi"]->Fill(deltapt, event_weight );
                                          histosTH1F["muonDeltaEta_jpsi"]->Fill(deltaeta, event_weight );	
                                          histosTH1F["muonDeltaPhi_jpsi"]->Fill(deltaphi, event_weight );
                                          histosTH1F["muonDphi_jpsi"]->Fill(Dphi, event_weight );
                                          histosTH1F["muonDeltaY_jpsi"]->Fill(deltay, event_weight );  
                                          //histosTH1F["jpsi_t_plus"]->Fill( fabs(t_proton_plus), event_weight );
                                          //histosTH1F["jpsi_xi_plus"]->Fill( xi_proton_plus, event_weight );
                                          jpsi_mass = dimuon_lorentz.M();  
                                          jpsi_pt=dimuon_lorentz.Pt(); 
                                          jpsi_pt2=dimuon_lorentz.Pt()*dimuon_lorentz.Pt(); 
                                          jpsi_eta =dimuon_lorentz.Eta();
                                          jpsi_rapidity = dimuon_lorentz.Rapidity();
                                        }   
                        }
                    }

                //===========
                chmuons = ( chmu1*chmu2 < 0. );
                cout<<"chmuons before"<<chmuons<<endl;
                if( !chmuons ) continue;
                cout<<"chmuons after"<<chmuons<<endl;
                if( !mu2_selected ) continue;
                 cout<<"mu2_selected"<<mu2_selected<<endl;
                //===============================================================
		//GenPart
		double genEPlusPz = 0.;
		double genEMinusPz = 0.;
		double cm = 8000.;
		double proton_pi = 4000;
		double proton_pz_plus=-999;
		double proton_px_plus = -999.;
		double proton_py_plus = -999.;
		double proton_energy_plus = 0;
		double proton_pz_minus= 999;
		double proton_px_minus = 999.;
		double proton_py_minus = 999.;
		double proton_energy_minus = 0;
		double px_gen;
		double py_gen;
		double pz_gen;
		double energy_gen;
		double proton_pf=0.;
                double proton_in =0.;

		for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){

			//double eta_gen = it_genpart->Eta();
			int status = it_genpart->status;
			int id = it_genpart->pdgId;

			if (status != 1) continue; //final state for the particles
			if (id != 2212) continue;

			energy_gen = it_genpart->Energy();
			px_gen = it_genpart->Px();
			py_gen = it_genpart->Py();
			pz_gen = it_genpart->Pz();
			//proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
                        proton_pf = sqrt(energy_gen*energy_gen+py_gen*py_gen+pz_gen*pz_gen);
                        proton_in = sqrt(2*(proton_pi*proton_pi));
			double pz_cut = 0.7*proton_pi;

			if(verbose) cout<<"Particle ID: " << id << " energy: "<< energy_gen <<" px_gen: " << px_gen << " py_gen: "<< py_gen << " pz_gen: "<< pz_gen<<" proton_pf: "<< proton_pf << endl; 

			//if (fabs(pz_gen) < pz_cut) continue;

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
		/*double t_plus = 0.;
		double fatora = 0.;
		double fatorb = 0.;*/
		double t_proton_plus = 0.;
		//double t_dir_plus = 0.;
		double t_proton_minus = 0.;
		//double t_dir_minus = 0.;
		double thx_proton_plus = 0.;
		double thy_proton_plus = 0.;
		double thx_proton_minus = 0.;
		double thy_proton_minus = 0.;

		++nevents_gen ;

		bool proton_minus_rp_accept_120 = false;
		bool proton_minus_rp_accept_121 = false;
		//bool proton_minus_rp_accept_122 = false;
		//bool proton_minus_rp_accept_123 = false;
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
			
                        TLorentzVector p2(0,0,4000,4000);
                        TLorentzVector p3(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus); // 4 momentum of p3
                        t_proton_plus =(p2-p3).M2();
			thx_proton_plus = atan(-proton_px_plus/proton_pi);
			thy_proton_plus = atan(proton_py_plus/proton_pi);
			cout<<"t_proton_plus :"<< t_proton_plus <<endl;
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

		histosTH1F["xi_proton_plus_before"]->Fill( xi_proton_plus , event_weight );
		histosTH1F["t_proton_plus_before"]->Fill( fabs(t_proton_plus) , event_weight );
                cout<<"t_proton_plus_before ="<<fabs(t_proton_plus)<<endl;
            
                if(((jpsi_mass > 3.0) && (jpsi_mass < 3.2))){
                cout<<"jpsi_t_plus"<<fabs(t_proton_plus)<<endl;
                histosTH1F["jpsi_t_plus"]->Fill( fabs(t_proton_plus), event_weight );
                histosTH1F["jpsi_xi_plus"]->Fill( xi_proton_plus, event_weight );}
                

		if(proton_pz_minus < 0.){
			xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;

                        TLorentzVector p1(0,0,-4000,4000);
                        TLorentzVector p4(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus); // 4 momentum of p3
                        t_proton_plus=(p1-p4).M2();
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

		histosTH1F["xi_proton_minus_before"]->Fill( xi_proton_minus , event_weight );
		histosTH1F["t_proton_minus_before"]->Fill( fabs(t_proton_minus) , event_weight );                                               
                //histosTH1F["t_dir_minus_before"]->Fill( t_dir_minus , event_weight );
	
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





			
                                  
//////////////////////////////////////////////////////
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

		histosTH1F["thx_proton_plus"]->Fill( thx_proton_plus , event_weight );
		histosTH1F["thy_proton_plus"]->Fill( thy_proton_plus , event_weight );
		histosTH2F["proton_plus_xi_vs_t"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );	 
		//histosTH2F["xi+Vseta_max"]->Fill( eta_max, xi1, event_weight ); 

		//rp_accept
		if (verbose) cout<< " |t| selection  : "<< " | "<< t_proton_down_ << " - " << t_proton_up_<< " |"<< endl;

		if (fabs(t_proton_plus)<t_proton_down_ || fabs(t_proton_plus)>t_proton_up_)continue;
		if(verbose) cout<< " rp accept after selection |t| : "<< " |t| proton plus "<< fabs(t_proton_plus)<< " xi proton plus "<< xi_proton_plus << endl;
		if (xi_proton_plus<xi_proton_down_ || xi_proton_plus>xi_proton_up_)continue; 

		if (verbose)  cout<< " rp accept after selection xi: "<< " |t| proton plus "<< fabs(t_proton_plus)<< " xi proton plus "<< xi_proton_plus << endl;

		histosTH1F["xi_proton_plus_after"]->Fill( xi_proton_plus , event_weight );
		histosTH1F["t_proton_plus_after"]->Fill( fabs(t_proton_plus) , event_weight );
                cout<<" t_proton_plus_after ="<< fabs(t_proton_plus)<<endl;
		histosTH1F["t_proton_plus_newbinning"]->Fill( fabs(t_proton_plus) , event_weight );


                if (verbose)cout<< " 020: " << proton_plus_rp_accept_020 <<endl;
                if (verbose)cout<< " 021: " << proton_plus_rp_accept_021 <<endl;
                if (verbose)cout<< " 022: " << proton_plus_rp_accept_022 <<endl;
                if (verbose)cout<< " 023: " << proton_plus_rp_accept_023 <<endl;
                if (verbose)cout<< " 024: " << proton_plus_rp_accept_024 <<endl;
                if (verbose)cout<< " 025: " << proton_plus_rp_accept_025 <<endl;


		bool proton_plus_rp_accept = (( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 ) ||( proton_plus_rp_accept_022 && proton_plus_rp_accept_023 ) );


		if(!proton_plus_rp_accept) continue; 
		// xi rp accepted 
		histosTH1F["xi_proton_plus_accepted"]->Fill( xi_proton_plus , event_weight );
		histosTH2F["proton_plus_xi_vs_t_accepted"]->Fill( fabs(t_proton_plus) , xi_proton_plus , event_weight );
		histosTH1F["t_proton_plus_accepted"]->Fill( fabs(t_proton_plus) , event_weight );
                cout<<"t_proton_plus_accepted ="<<fabs(t_proton_plus)<<endl;
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
		//----------------------------------------------
                //dimuons              
		++n_dimuons_rp_selected_plus_side;	
                //histosTH1F["dimuon_multiplicity_rp_accept"]->Fill(n_dimuons_rp_selected_plus_side, event_weight );
                histosTH1F["dimuon_mass_rp_accept"]->Fill( dimuon_mass, event_weight );
                histosTH1F["dimuon_pt_rp_accept"]->Fill( dimuon_pt, event_weight );
                histosTH1F["dimuon_pt2_rp_accept"]->Fill(dimuon_pt2, event_weight );
                histosTH1F["dimuon_eta_rp_accept"]->Fill( dimuon_eta, event_weight );
                histosTH1F["dimuon_rapidity_rp_accept"]->Fill( dimuon_rapidity, event_weight );
                histosTH1F["muonDeltaPt_rp_accept"]->Fill(deltapt, event_weight );
                histosTH1F["muonDeltaEta_rp_accept"]->Fill(deltaeta, event_weight );	
                histosTH1F["muonDphi_rp_accept"]->Fill(Dphi, event_weight );
                histosTH1F["muonDeltaPhi_rp_accept"]->Fill(deltaphi, event_weight );
                histosTH1F["muonDeltaY_rp_accept"]->Fill(deltay, event_weight ); 
                histosTH1F["dimuon_t_plus_rp_accept"]->Fill( fabs(t_proton_plus), event_weight );
		cout<<"dimuon_t_plus_rp_accept= "<< fabs(t_proton_plus)<<endl;
                histosTH1F["dimuon_xi_plus_rp_accept"]->Fill( xi_proton_plus, event_weight );
                //jpsi  
                //if(!((jpsi_mass > 3.0) && (jpsi_mass < 3.2)))continue;
                if(((jpsi_mass > 3.0) && (jpsi_mass < 3.2))){
                   cout<<"jpsi_mass = "<< jpsi_mass << endl;
                  //++n_jpsi_rp_accept_selected;
                  ++n_jpsi_selected_rp_accept_plus_side;
                  //histosTH1F["jpsi_multiplicity_rp_accept"]->Fill(n_dimuons_selected, event_weight );
                  histosTH1F["jpsi_mass_rp_accept"]->Fill(jpsi_mass, event_weight );
                  histosTH1F["jpsi_pt_rp_accept"]->Fill( jpsi_pt, event_weight );
                  histosTH1F["jpsi_pt2_rp_accept"]->Fill( jpsi_pt2, event_weight );
                  histosTH1F["jpsi_eta_rp_accept"]->Fill( jpsi_eta, event_weight );
                  histosTH1F["jpsi_rapidity_rp_accept"]->Fill( jpsi_rapidity, event_weight );
                  histosTH1F["muonDeltaPt_jpsi_rp_accept"]->Fill(deltapt, event_weight );
                  histosTH1F["muonDeltaEta_jpsi_rp_accept"]->Fill(deltaeta, event_weight );	
                  histosTH1F["muonDeltaPhi_jpsi_rp_accept"]->Fill(deltaphi, event_weight );
                  histosTH1F["muonDeltaY_jpsi_rp_accept"]->Fill(deltay, event_weight );  
                  
                  histosTH1F["jpsi_t_plus_rp_accept"]->Fill( fabs(t_proton_plus), event_weight );
                  histosTH1F["jpsi_xi_plus_rp_accept"]->Fill( xi_proton_plus, event_weight );
                  histosTH1F["muonDphi_jpsi_rp_accept"]->Fill(Dphi, event_weight );
                 }
                                //} //Second Muons 
                       // }// First Muon

	
			


			//------------------



		}//end loop for events

		cout<<"Total of evts="<< nev << endl << *itfiles << endl;
                cout<<"n_vertices_selected ="<< n_vertices_selected << endl;
                cout<<"n_select_Vertex_After_vtx_cut="<< n_select_Vertex_After_vtx_cut<<endl;
                cout<<"n_tracks_selected="<< n_tracks_selected<<endl;
                cout<<"n_muon_selected="<< n_muon_selected<<endl;
                cout<<"n_muon_selected_size ="<<n_muon_selected_size<<endl;
                cout<<"n_dimuons_selected="<< n_dimuons_selected<<endl;
		cout<<"n_dimuons_rp_accept_plus_side="<< n_dimuons_rp_selected_plus_side <<endl;
                cout<<"n_jpsi_selected ="<< n_jpsi_selected<<endl;
		cout<<"n_jpsi_selected_rp_accept_plus_side ="<< n_jpsi_selected_rp_accept_plus_side<<endl;
                cout<<"nevents_total for weight ="<< nevents_total<<endl;

 

		file->Close();

	}//end of loop over files

	//output file
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

