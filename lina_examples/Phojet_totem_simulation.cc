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

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

//TOTEM CLASSES
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

//#include "rp_aperture_config.h"

#define PI 3.141592653589793
using namespace std;


void Phojet_totem_simulation(string const& outputFileName = "phojet_totem_simulation.root", const Int_t nevt_max = -1){
  
  bool verbose = false;
  string treeName = "TotemNtuple";
  
  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


  map<string,TH1F*> histosTH1F;

  Float_t tbins[12] = {0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.24, 0.30, 0.36, 0.45, 0.65, 1.};
  histosTH1F["proton_rec_right_xi"] = new TH1F("proton_rec_right_xi", "proton_rec_right_xi" , 200, -1, 1);
  histosTH1F["proton_rec_right_t"] = new TH1F("proton_rec_right_t", "proton_rec_right_t" , 11, tbins);
  histosTH1F["proton_rec_right_chi2"] = new TH1F("proton_rec_right_chi2", "proton_rec_right_chi2" , 200, 0, 100);
  histosTH1F["proton_rec_left_xi"] = new TH1F("proton_rec_left_xi", "proton_rec_left_xi" , 200, -1, 1);
  histosTH1F["proton_rec_left_t"] = new TH1F("proton_rec_left_t", "proton_rec_left_t" , 11, tbins);
  histosTH1F["proton_rec_left_chi2"] = new TH1F("proton_rec_left_chi2", "proton_rec_left_chi2" , 200, 0, 100);
  histosTH1F["proton_sim_right_xi"] = new TH1F("proton_sim_right_xi", "proton_sim_right_xi" , 200, -1, 1);
  histosTH1F["proton_sim_right_t"] = new TH1F("proton_sim_right_t", "proton_sim_right_t" , 11, tbins);
  histosTH1F["proton_sim_right_chi2"] = new TH1F("proton_sim_right_chi2", "proton_sim_right_chi2" , 200, 0, 100);
  histosTH1F["proton_sim_left_xi"] = new TH1F("proton_sim_left_xi", "proton_sim_left_xi" , 200, -1, 1);
  histosTH1F["proton_sim_left_t"] = new TH1F("proton_sim_left_t", "proton_sim_left_t" , 11, tbins);
  histosTH1F["proton_sim_left_chi2"] = new TH1F("proton_sim_left_chi2", "proton_sim_left_chi2" , 200, 0, 100);

  map<string,TH2F*> histosTH2F;
  histosTH2F["proton_right_t_sim_vs_rec"] = new TH2F("proton_right_t_sim_vs_rec","proton_right_t_sim_vs_rec", 11 , tbins, 11, tbins);
  histosTH2F["proton_left_t_sim_vs_rec"] = new TH2F("proton_left_t_sim_vs_rec","proton_left_t_sim_vs_rec", 11 , tbins, 11, tbins);


  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  gStyle->SetPalette(1);
  //===================
  int i_tot = 0 , nevt_tot = 0;


  vector<TString>* vfiles = new vector<TString>; 
  vfiles->push_back("/storage1/lhuertas/Analysis/CMSTotem/MC/TotemSimulation/July2012_90m/RPphojetDPEbeta90energy4p0TeV_July2012_totem_ntuple.root");
   

  //Declaration of tree and its branches variables
  TTree* tree = new TTree(treeName.c_str(),"");
  T1Event* t1_event = NULL;
  T2Event* t2_event = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  RPRootDumpReconstructedProton* sim_proton_right = NULL;
  RPRootDumpReconstructedProton* sim_proton_left = NULL;
  RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
  map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
  //rp_aperture_config();

  
  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get( treeName.c_str() );


    //Getting number of events
    int nev = int(tree->GetEntries());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("branchT1EV.",&t1_event);
    tree->SetBranchAddress("branchT2EV.",&t2_event);
    tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
    tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
    tree->SetBranchAddress("sim_prot_left.",&sim_proton_left);
    tree->SetBranchAddress("sim_prot_right.",&sim_proton_right);
    tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
    std::vector<unsigned int> rp_list;
    rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(22); rp_list.push_back(23); rp_list.push_back(24); rp_list.push_back(25);
    rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(122); rp_list.push_back(123); rp_list.push_back(124); rp_list.push_back(125);
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
  
    /*//Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/
  
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){
    
    //printing the % of events done every 10k evts
    if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
  //    double event_weight = genKin->genWeight; 
      double event_weight = 1.0;
      
      // RP reconstructed protons
      bool proton_rec_right_valid = rec_proton_right->valid;
      double chi2_proton_rec_right = rec_proton_right->chi2;
      double chindf_proton_rec_right = rec_proton_right->chindf;
      double xi_proton_rec_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_rec_right = rec_proton_right->t;
      bool good_proton_rec_right = proton_rec_right_valid; // && (xi_proton_right < 0.);

      bool proton_rec_left_valid = rec_proton_left->valid;
      double chi2_proton_rec_left = rec_proton_left->chi2;
      double chindf_proton_rec_left = rec_proton_left->chindf;
      double xi_proton_rec_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_rec_left = rec_proton_left->t;
      bool good_proton_rec_left = proton_rec_left_valid; // && (xi_proton_left < 0.);
 

      // RP simulated protons
      bool proton_sim_right_valid = sim_proton_right->valid;
      double chi2_proton_sim_right = sim_proton_right->chi2;
      double chindf_proton_sim_right = sim_proton_right->chindf;
      double xi_proton_sim_right = sim_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_sim_right = sim_proton_right->t;
      bool good_proton_sim_right = proton_sim_right_valid; // && (xi_proton_right < 0.);

      bool proton_sim_left_valid = sim_proton_left->valid;
      double chi2_proton_sim_left = sim_proton_left->chi2;
      double chindf_proton_sim_left = sim_proton_left->chindf;
      double xi_proton_sim_left = sim_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_sim_left = sim_proton_left->t;
      bool good_proton_sim_left = proton_sim_left_valid; // && (xi_proton_left < 0.);


      if( good_proton_rec_right || good_proton_sim_right ){
        if(-t_proton_rec_right>=0.03 && -t_proton_rec_right<=1. && -xi_proton_rec_right>=0.03 && -xi_proton_rec_right<=0.1 && -t_proton_sim_right>=0.03 && -t_proton_sim_right<=1. && -xi_proton_sim_right>=0.03 && -xi_proton_sim_right<=0.1) {
          histosTH1F["proton_rec_right_chi2"]->Fill( chi2_proton_rec_right, event_weight );
          histosTH1F["proton_rec_right_xi"]->Fill(-xi_proton_rec_right, event_weight );
          histosTH1F["proton_rec_right_t"]->Fill( -t_proton_rec_right, event_weight );
          histosTH1F["proton_sim_right_chi2"]->Fill( chi2_proton_sim_right, event_weight );
          histosTH1F["proton_sim_right_xi"]->Fill(-xi_proton_sim_right, event_weight );
          histosTH1F["proton_sim_right_t"]->Fill( -t_proton_sim_right, event_weight );
          histosTH2F["proton_right_t_sim_vs_rec"]->Fill(-t_proton_rec_right,-t_proton_sim_right, event_weight ); 

        }
      }
/*      if( good_proton_sim_right ){
        if(-t_proton_sim_right>=0.03 && -t_proton_sim_right<=1. && -xi_proton_sim_right>=0.03 && -xi_proton_sim_right<=0.1) {
          histosTH1F["proton_sim_right_chi2"]->Fill( chi2_proton_sim_right, event_weight );
          histosTH1F["proton_sim_right_xi"]->Fill(-xi_proton_sim_right, event_weight );
          histosTH1F["proton_sim_right_t"]->Fill( -t_proton_sim_right, event_weight );
        }
      }
*/


      if( good_proton_rec_left ){
        if(-t_proton_rec_left>=0.03 && -t_proton_rec_left<=1. && -xi_proton_rec_left>=0.03 && -xi_proton_rec_left<=0.1) {
          histosTH1F["proton_rec_left_chi2"]->Fill( chi2_proton_rec_left, event_weight );
          histosTH1F["proton_rec_left_xi"]->Fill(-xi_proton_rec_left, event_weight );
          histosTH1F["proton_rec_left_t"]->Fill( -t_proton_rec_left, event_weight );
        }
      }

      
      if( good_proton_sim_left ){
        if(-t_proton_sim_left>=0.03 && -t_proton_sim_left<=1. && -xi_proton_sim_left>=0.03 && -xi_proton_sim_left<=0.1) {
          histosTH1F["proton_sim_left_chi2"]->Fill( chi2_proton_sim_left, event_weight );
          histosTH1F["proton_sim_left_xi"]->Fill(-xi_proton_sim_left, event_weight );
          histosTH1F["proton_sim_left_t"]->Fill( -t_proton_sim_left, event_weight );
        }
      }

    }//end loop for events

   file->Close();
 
  }//end of loop over files
     
  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();
  histosTH1F["proton_rec_right_chi2"]->GetXaxis()->SetTitle("#chi^{2}");
  histosTH1F["proton_rec_right_xi"]->GetXaxis()->SetTitle("#xi_{TOTEM}");
  histosTH1F["proton_rec_right_t"]->GetXaxis()->SetTitle("|t|(GeV^{2})");
  histosTH1F["proton_sim_right_chi2"]->GetXaxis()->SetTitle("#chi^{2}");
  histosTH1F["proton_sim_right_xi"]->GetXaxis()->SetTitle("#xi_{TOTEM}");
  histosTH1F["proton_sim_right_t"]->GetXaxis()->SetTitle("|t|(GeV^{2})");
  histosTH1F["proton_rec_left_chi2"]->GetXaxis()->SetTitle("#chi^{2}");
  histosTH1F["proton_rec_left_xi"]->GetXaxis()->SetTitle("#xi_{TOTEM}");
  histosTH1F["proton_rec_left_t"]->GetXaxis()->SetTitle("|t|(GeV^{2})");
  histosTH1F["proton_sim_left_chi2"]->GetXaxis()->SetTitle("#chi^{2}");
  histosTH1F["proton_sim_left_xi"]->GetXaxis()->SetTitle("#xi_{TOTEM}");
  histosTH1F["proton_sim_left_t"]->GetXaxis()->SetTitle("|t|(GeV^{2})");
  histosTH2F["proton_right_t_sim_vs_rec"]->GetXaxis()->SetTitle("|t|^{rec}(GeV^{2})");
  histosTH2F["proton_right_t_sim_vs_rec"]->GetYaxis()->SetTitle("|t|^{gen}(GeV^{2})");
 
  
  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin(); it_histo != histosTH1F.end(); ++it_histo){
     (*it_histo).second->GetYaxis()->SetTitle("Events");
     (*it_histo).second->Write();
  }
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

  output->Close();
}
