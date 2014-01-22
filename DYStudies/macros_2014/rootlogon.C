{
   std::cout << "Loading rootlogon.C" << std::endl;
   gROOT->ProcessLine(".L tdrstyle.C");
   setTDRStyle();
   gStyle->SetOptStat(1111111);
   gStyle->SetHistLineWidth(2);
   //gStyle->SetHistFillStyle(1);
   gStyle->SetMarkerStyle(23);
   gStyle->SetMarkerSize(0.3);
   gStyle->SetErrorX(0.5);
   gStyle->SetPalette(1);

//Lina
   gSystem->Load("/storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotemEnv/TOTEMdataFormat/lib/libTOTEMdataFormat.so");
   gSystem->Load("/storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotemEnv/CMSdataFormat/lib/libCMSdataFormat.so");

     gSystem->AddIncludePath(" -I/storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotemEnv/CMSdataFormat/src ");
     gSystem->AddIncludePath(" -I/storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotemEnv/TOTEMdataFormat/src ");

     gSystem->AddIncludePath(" -I/afs/cern.ch/exp/totem/scratch/hubert/cmssw_code/CMSSW_4_2_4_SVN_6/CMSSW_4_2_4/src ");

     gSystem->Load("parametrization/libSimG4CMSTotemRPProtTranspPar");
 
 //Sandro Config ( it does not workinhg with MC samples) 05/09/2013
   //gSystem->Load("/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/TOTEMdataFormat/lib/libTOTEMdataFormat.so");
   //gSystem->Load("/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/CMSdataFormat/lib/libCMSdataFormat.so");

   //gSystem->AddIncludePath(" -I/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/CMSdataFormat/src ");
   //gSystem->AddIncludePath(" -I/storage2/sfonseca/CMSSW/TOTEM/CMSTotem/TOTEMdataFormat/src ");


   //gSystem->AddIncludePath(" -I/afs/cern.ch/exp/totem/scratch/hubert/cmssw_code/CMSSW_4_2_4_SVN_4/CMSSW_4_2_4/src ");

   //gSystem->Load("parametrization/libSimG4CMSTotemRPProtTranspPar");
   
   //gSystem->Load("libFWCoreFWLite");
   //AutoLibraryLoader::enable();
}
