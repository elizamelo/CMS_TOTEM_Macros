#include <map>
#include <TFile.h>
#include <TROOT.h>
#include <TDirectory.h>

//input: 
//x* [m]
//thx* [rad] (p_x/p_nominal, p_nominal=4TeV) !!!
//y* [m]
//thy* [rad] (p_y/p_nominal, p_nominal=4TeV) !!!
//xi* [delta_p/p_nominal, should be negative by definition]
//
//output:
//proton passed successfully lhc apertures and entered the RP (only vertical position verified)
//
//comment:
//For the time being the optics used is symmetric for both arms. The optics studies are ongoing which will 
//result in increased reconstruction resolution. Once completed we'll update you with the outcome.

#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"

bool protonRPFiducialIdeal(double st_x, double st_y, double rp_position);
bool protonRPDetected(double x, double thx, double y, double thy, double xi, int rp_id, bool verbose = false);
bool protonRPDetected(double x, double thx, double y, double thy, double xi, int rp_id,
                      double& out_x, double& out_thx, double& out_y, double& out_thy, double& out_xi,
                      bool verbose = false);

std::map<int, double> rp_positions;
std::map<int, LHCOpticsApproximator *> rp_aproximators;

void rp_aperture_config()
{
    TFile* file = TFile::Open("parametrization/parametrization_4000GeV_90p0_reco.root");

    rp_positions[120] = 7.15;
    rp_positions[124] = 8;
    rp_positions[121] = -7.4;
    rp_positions[125] = -8.25;
    //FIXME
    rp_positions[122] = 7.;
    rp_positions[123] = 6.4;

    rp_positions[20] = 7.25;
    rp_positions[24] = 8.1;
    rp_positions[21] = -7.5;
    rp_positions[25] = -8.35;
    //FIXME
    rp_positions[22] = 6.6;
    rp_positions[23] = 6.25;
    
    rp_aproximators[120] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_1_lhcb1");
    rp_aproximators[124] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_2_lhcb1");
    rp_aproximators[121] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_1_lhcb1");
    rp_aproximators[125] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_2_lhcb1");
    rp_aproximators[122] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_h_1_lhcb1");
    rp_aproximators[123] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_h_2_lhcb1");

    rp_aproximators[20] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_1_lhcb2");
    rp_aproximators[24] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_2_lhcb2");
    rp_aproximators[21] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_1_lhcb2");
    rp_aproximators[25] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_v_2_lhcb2");
    rp_aproximators[22] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_h_1_lhcb2");
    rp_aproximators[23] = (LHCOpticsApproximator *) file->Get("ip5_to_station_220_h_2_lhcb2");
}

bool protonRPFiducialIdeal(double st_x, double st_y, double rp_position){ //mm

   // Diamond dimensions
   double b = 19.53;
   double l1 = 18.75;
   double l2 = 32.56;
   double a1 = l1/( TMath::Sqrt( double(2.) ) );
   double a2 = l2/( TMath::Sqrt( double(2.) ) );

   bool border1 = ( st_y > rp_position);

   bool border2 = ( st_y >  st_x + ( rp_position - b/2 ) );
   bool border3 = ( st_y < -st_x + ( rp_position + b/2 + 2*a1 ) );

   bool border4 = ( st_y <  st_x + ( rp_position + b/2 + 2*a1 ) );
   bool border5 = ( st_y > -st_x + ( rp_position - b/2 ) );

   bool accept = border1 && border2 && border3 && border4 && border5;
   return accept;
}

bool protonRPDetected(double x, double thx, double y, double thy, double xi, int rp_id, bool verbose)
{
    double out_x, out_thx, out_y, out_thy, out_xi;
    return protonRPDetected(x, thx, y, thy, xi, rp_id, out_x, out_thx, out_y, out_thy, out_xi, verbose);
}

bool protonRPDetected(double x, double thx, double y, double thy, double xi, int rp_id,
                      double& out_x, double& out_thx, double& out_y, double& out_thy, double& out_xi,
                      bool verbose)
{
    double in[5];
    double out[5];

    in[0] = x; //x* [m]
    in[1] = thx; //thx* [rad]
    in[2] = y; //y* [m]
    in[3] = thy; //thy* [rad]
    in[4] = xi; //xi* [deltap/p]
    
    bool check_apertures = true;
    bool within_apertures = rp_aproximators[rp_id]->Transport(in, out, check_apertures);
    bool detected = false;

    if(within_apertures)
    {
        std::vector<int> rp_list_top_bot;
        rp_list_top_bot.push_back(20); rp_list_top_bot.push_back(21);
        rp_list_top_bot.push_back(24); rp_list_top_bot.push_back(25);
        rp_list_top_bot.push_back(120); rp_list_top_bot.push_back(121);
        rp_list_top_bot.push_back(124); rp_list_top_bot.push_back(125);
        std::vector<int> rp_list_hor;
        rp_list_hor.push_back(22); rp_list_hor.push_back(23);
        rp_list_hor.push_back(122); rp_list_hor.push_back(123);
	if( std::find(rp_list_top_bot.begin(), rp_list_top_bot.end(), rp_id) != rp_list_top_bot.end() ){
	   double sign = ( rp_id%2 == 0 ) ? (1) : (-1); //top or bottom pot
	   //detected = out[2]*1000*sign > rp_positions[rp_id]*sign;

	   double st_x = out[0]*1000*sign; 
	   double st_y = out[2]*1000*sign; 
	   double rp_pos = rp_positions[rp_id]*sign;
	   detected = protonRPFiducialIdeal(st_x,st_y,rp_pos);
	} else if( std::find(rp_list_hor.begin(), rp_list_hor.end(), rp_id) != rp_list_hor.end() ){
	   double st_x = -out[2]*1000; 
	   double st_y =  out[0]*1000; 
	   double rp_pos = rp_positions[rp_id];
	   detected = protonRPFiducialIdeal(st_x,st_y,rp_pos);
        } else cout << "WARNING: RP ID INVALID " << rp_id << endl;
    }

    if(verbose){
       cout<<"RP: "<<rp_id<<endl;
       cout<<"x="<<out[0]<<endl;
       cout<<"thx="<<out[1]<<endl;
       cout<<"y="<<out[2]<<endl;
       cout<<"thy="<<out[3]<<endl;
       cout<<"xi="<<out[4]<<endl;
       cout<<"within LHC aperture: "<<within_apertures<<endl;
       cout<<"detected: "<<detected<<endl;
    }

    out_x = out[0]; out_thx = out[1]; 
    out_y = out[2]; out_thy = out[3];
    out_xi = out[4];
 
    return detected;
}


