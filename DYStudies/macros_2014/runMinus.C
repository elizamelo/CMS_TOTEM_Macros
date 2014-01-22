{
	std::cout << "Loading ana_psel_Data_MC_DDY_v1.cc" << std::endl;
	// gROOT->ProcessLine(".L tdrstyle.C");
	gROOT->ProcessLine(".L ana_psel_Data_MC_DDY_v1.cc++");
	/////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Loading Data_Mll_20_60_RPMinusAccept" << std::endl;
	string input = "/afs/cern.ch/work/s/sfonseca/8369_8371_8272_CMS_TOTEM_NTuples_13JUN.root";
	string output = "/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/Workspace/DYStudies/SingleArmRecProton/new_hist/Data_Mll_20_60_RPMinusAccept.root"; 
	gROOT->ProcessLine("ana_psel_DataMC_DDY_v1(vector<string>(1,input),output,0.03,1.0,20.0,60.0,200,-1);");
	///////////////////////////////////////////////////////////////////////////////
	std::cout << "Loading Data_Mll_10_20_RPMinusAccept" << std::endl;
	//string input = "/storage2/sfonseca/CMSSW/TOTEM/CMSSW_5_3_2_patch4/src/test/MC/POMPTY/DYToMuMu_M_10To20/DYToMuMu_M_10To20.root";
	string output = "/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/Workspace/DYStudies/SingleArmRecProton/new_hist/Data_Mll_10_20_RPMinusAccept.root";
	gROOT->ProcessLine("ana_psel_DataMC_DDY_v1(vector<string>(1,input),output,0.03,1.0,10.0,20.0,200,-1);");
	///////////////////////////////////////////////////////////////////////////////////
	std::cout << "Loading Data_Mll_6_10_RPMinusAccept" << std::endl;
	//string input = "/storage2/sfonseca/CMSSW/TOTEM/CMSSW_5_3_2_patch4/src/test/MC/POMPTY/DYToMuMu_M_6To10/DYToMuMu_M_6To10.root";
	string output = "/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/Workspace/DYStudies/SingleArmRecProton/new_hist/Data_Mll_6_10_RPMinusAccept.root";				
	gROOT->ProcessLine("ana_psel_DataMC_DDY_v1(vector<string>(1,input),output,0.03,1.0,6.0,10.0,200,-1);");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




}

