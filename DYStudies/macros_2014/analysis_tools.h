#include "TTree.h"

#include "MyPFCand.h"
#include "MyPFJet.h"

#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpReconstructedProtonPair.h"

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>

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
float minEVarBin = 0.; 
float binningEPlusPz[]={minEVarBin,5.,10.,15.,20.,25.,30.,40.,50.,60.,70.,80.,90.,
                                   100.,120.,140.,160.,180.,200.,225.,250.,
                                   275.,300.,350.,400.}; // 24 bins

bool sortByEta( const TLorentzVector& a, const TLorentzVector& b){ 
   return a.Eta() < b.Eta();
}

bool loosePFJetID(MyPFJet& pfJet, string const& coll){
  return ( pfJet.nconstituent > 1 && pfJet.fhad_ne < 0.99 && pfJet.fem_ne < 0.99 &&
           ( ( fabs((pfJet.mapjet[coll]).Eta()) <= 2.4 && pfJet.fhad_ch > 0 && pfJet.multi_ch > 0 && pfJet.fem_ch < 0.99 ) || fabs((pfJet.mapjet[coll]).Eta()) > 2.4 ) );
}
bool tightPFJetID(MyPFJet& pfJet, string const& coll){
  return ( pfJet.nconstituent > 1 && pfJet.fhad_ne < 0.90 && pfJet.fem_ne < 0.90 &&
           ( ( fabs((pfJet.mapjet[coll]).Eta()) <= 2.4 && pfJet.fhad_ch > 0 && pfJet.multi_ch > 0 && pfJet.fem_ch < 0.99 ) || fabs((pfJet.mapjet[coll]).Eta()) > 2.4 ) );
}

//--------
TBranch* checkAndGetBranch(TTree* tree, std::string branchName){
   TBranch *branch = tree->GetBranch(branchName.c_str());
   if (!branch) {
      string dotBranchName = branchName + ".";
      branch = tree->GetBranch(dotBranchName.c_str());
   }
   if (!branch) {
      tree->Print();
      std::cout << " No data branch " << branchName << " found in input file!" << std::endl;
      exit(-1);
   }
   return branch;
}

void readRPBranches(TTree *tree_totem,
                    RPRootDumpReconstructedProton *recProtLeft,
                    RPRootDumpReconstructedProton *recProtRight,
                    RPRootDumpReconstructedProtonPair *recProtonPair, 
                    std::map<unsigned int, RPRootDumpTrackInfo *>& track_info,
                    std::map<unsigned int, RPRootDumpDigiInfo *>& digi_info,
                    std::map<unsigned int, RPRootDumpPatternInfo *>& par_patterns_info,
                    std::map<unsigned int, RPRootDumpPatternInfo *>& nonpar_patterns_info, 
                    std::map<unsigned int, std::vector<RPRootDumpTrackInfo> *>& multi_track_info,
                    bool includeDigi = false, bool includePatterns = false, bool includeMultiTracks = false) {

   checkAndGetBranch(tree_totem, "rec_prot_left")->SetAddress( &recProtLeft );
   checkAndGetBranch(tree_totem, "rec_prot_right")->SetAddress( &recProtRight );
   checkAndGetBranch(tree_totem, "rec_prot_pair")->SetAddress( &recProtonPair );

   //the digis, patterns and tracks are stored separately
   //for every RP, this loops iterate over RP ids
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
         //track_info[id] = NULL;
	 checkAndGetBranch(tree_totem, br_name)->SetAddress( &track_info[id] );

	 if (includeDigi) {
	    sprintf(br_name, "digi_rp_%u.", id);
	    checkAndGetBranch(tree_totem, br_name)->SetAddress( &digi_info[id] );
	 }
	 if (includePatterns) {
	    sprintf(br_name, "par_patterns_rp_%u.", id);
	    checkAndGetBranch(tree_totem, br_name)->SetAddress( &par_patterns_info[id] );

	    sprintf(br_name, "nonpar_patterns_rp_%u.", id);
	    checkAndGetBranch(tree_totem, br_name)->SetAddress( &nonpar_patterns_info[id] );
	 }
	 if (includeMultiTracks) {
	    sprintf(br_name, "multi_track_rp_%u", id);
	    checkAndGetBranch(tree_totem, br_name)->SetAddress( &multi_track_info[id] );
	 }
      }
   }
}

//--------
enum calo_region_t {Barrel,Endcap,Transition,Forward};

typedef std::map<int,std::pair<double,double> > ThresholdsPerType;
typedef std::map<int,ThresholdsPerType>    ThresholdsPerRegion;
 
void resetPFThresholds(ThresholdsPerType& thresholdsPFlow){
  thresholdsPFlow[MyPFCand::X] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::h] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::e] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::mu] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::gamma] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::h0] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::h_HF] = std::make_pair(-1.,-1.);
  thresholdsPFlow[MyPFCand::egamma_HF] = std::make_pair(-1.,-1.);
}

bool pflowThreshold(MyPFCand const& part, ThresholdsPerRegion const& thresholdMap){

   bool accept = true;

   //FIXME
   double eta = part.Eta();
   // HF eta rings 29, 40, 41
   /*if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 2.976) ) || 
         (fabs(eta) >= 4.730) ) return false;*/
   // HF eta rings 29, 30, 40, 41
   if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || 
         (fabs(eta) >= 4.730) ) return false;

   int region = -1;
   if( (fabs(eta) >= 0.) && (fabs(eta) < 1.4) ) region = Barrel;
   else if( (fabs(eta) >= 1.4) && (fabs(eta) < 2.6) ) region = Endcap;
   else if( (fabs(eta) >= 2.6) && (fabs(eta) < 3.2) ) region = Transition;
   else if( (fabs(eta) >= 3.2) ) region = Forward;
   ThresholdsPerType const& thresholds = thresholdMap.find(region)->second;
   
   double ptThreshold = -1.0;
   double eThreshold = -1.0;
   int partType = part.particleId;
   ThresholdsPerType::const_iterator it_threshold = thresholds.find(partType);
   if(it_threshold != thresholds.end()) {
      ptThreshold = it_threshold->second.first;
      eThreshold = it_threshold->second.second;
   }

   if(part.Pt() < ptThreshold) accept = false;
   if(part.Energy() < eThreshold) accept = false;

   return accept;
}

//---------------------
bool elastic_top45_bot56(std::map<unsigned int, RPRootDumpTrackInfo*>& track_info_){
   bool diag_top45_bot56_topol = track_info_[121]->valid && track_info_[125]->valid && 
                                 track_info_[20]->valid && track_info_[24]->valid;
   bool low_xi_56_top45_bot_56 = fabs(track_info_[125]->x - 1.60185e-03) < 5 * 2.88043e-01;
   bool low_xi_45_top45_bot_56 = fabs(track_info_[24]->x - 3.46543e-04) < 5 * 2.81866e-01;

   bool thxthx_divergence_cut_top45_bot_56 = fabs( 
      (track_info_[24]->x - track_info_[20]->x) / (-track_info_[24]->z + track_info_[20]->z) 
      - track_info_[24]->x / 1e3 / -1.865760089 * 0.055513006 
      + (track_info_[125]->x - track_info_[121]->x) / (track_info_[125]->z - track_info_[121]->z) 
      - track_info_[125]->x / 1e3 / -1.865760089 * 0.055513006 - -6.60482e-07 ) < 5 * 4.89e-6;

   bool div_y_cut_top45_bot_56 = fabs(
      (track_info_[24]->y / 1e3 / 263.143819469 + track_info_[20]->y / 1e3 / 237.668241862) / 2 
      + (track_info_[125]->y / 1e3 / 263.143819469 + track_info_[121]->y / 1e3 / 237.668241862) / 2 
      - 5.25e-7 ) < 5 * 3.2199e-6;

   /*bool not_elastic_top45_bot_56 = diag_top45_bot56_topol &&
                                   !low_xi_56_top45_bot_56 || 
                                   !low_xi_45_top45_bot_56 || 
                                   !thxthx_divergence_cut_top45_bot_56 || 
                                   !div_y_cut_top45_bot_56;*/
   bool elastic_top45_bot56 = diag_top45_bot56_topol &&
                              low_xi_56_top45_bot_56 && 
                              low_xi_45_top45_bot_56 && 
                              thxthx_divergence_cut_top45_bot_56 && 
                              div_y_cut_top45_bot_56;

   
   return elastic_top45_bot56;
}

bool elastic_bot45_top56(std::map<unsigned int, RPRootDumpTrackInfo*>& track_info_){
   bool diag_bot45_top56_topol = track_info_[120]->valid && track_info_[124]->valid &&
                                 track_info_[21]->valid && track_info_[25]->valid;

   bool low_xi_56_bot45_top56 = fabs(track_info_[124]->x - 0) < 5 * 2.99388e-01;
   bool low_xi_45_bot45_top56 = fabs(track_info_[25]->x - 0) < 5 * 2.90251e-01;

   bool x_div_bot45_top56 = fabs(
      (track_info_[124]->x - track_info_[120]->x) / (track_info_[124]->z - track_info_[120]->z) 
      - track_info_[124]->x / 1e3 / -1.865760089 * 0.055513006 
      + (track_info_[25]->x - track_info_[21]->x) / (-track_info_[25]->z + track_info_[21]->z) 
      - track_info_[25]->x / 1e3 / -1.865760089 * 0.055513006 - -3.31182e-07) < 5 * 5.69131e-06;

   bool div_y_cut_bot45_top56 = fabs(
      (track_info_[25]->y / 1e3 / 263.143819469 + track_info_[21]->y / 1e3 / 237.668241862) / 2 
      + (track_info_[124]->y / 1e3 / 263.143819469 + track_info_[120]->y / 1e3 / 237.668241862) / 2 
      - -4.797e-7) < 5 * 3.221e-6;

   /*bool not_elastic_bot45_top56 = diag_bot45_top56_topol && 
                                  !low_xi_56_bot45_top56 || 
                                  !low_xi_45_bot45_top56 || 
                                  !x_div_bot45_top56 || 
                                  !div_y_cut_bot45_top56;*/
   bool elastic_bot45_top56 = diag_bot45_top56_topol && 
                              low_xi_56_bot45_top56 && 
                              low_xi_45_bot45_top56 && 
                              x_div_bot45_top56 && 
                              div_y_cut_bot45_top56;

   return elastic_bot45_top56;
}
