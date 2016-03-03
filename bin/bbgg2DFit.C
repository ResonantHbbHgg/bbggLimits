// C++ headers
#include "HiggsAnalysis/bbggLimits/src/BrazilianFlag.cc"
#include "HiggsAnalysis/bbggLimits/interface/bbgg2DFitter.h"
#include "HiggsAnalysis/bbggLimits/interface/Colors.h"
#include "HiggsAnalysis/bbggLimits/interface/RunCombine.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
// RooFit headers
#include <RooWorkspace.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooArgSet.h>
#include "RooStats/HLFactory.h"
#include <RooDataSet.h>
#include <RooFormulaVar.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooBernstein.h>
#include <RooExtendPdf.h>
#include <RooMinimizer.h>
#include "RooStats/RooStatsUtils.h"
#include <RooProdPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>

// namespaces
using namespace std;
using namespace RooFit;
using namespace RooStats;

int main(int argc, const char* argv[])
{
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL) ;
	RooMsgService::instance().setSilentMode(true);
	std::vector<std::string>higgstrue;
  //Parameters to set
	Float_t mass = 1;
	Float_t lumi = 1;
	Bool_t doBands = 1;
  bool docombine=false;
	bool doBrazilianFlag=false;
	int version = 44;
	string analysisType = "a";
	string nonresFile = "b";
	Bool_t useSigTheoryUnc = 1;
	std::vector<std::pair<int,std::string>> sigMass;
	Bool_t doblinding = 1;
	Int_t NCAT = 1;
	bool addHiggs =1;
	string signalType = "";
	string signalDir = "";
	string energy = "13TeV";
	string cardName = "";
	string datadir="";
	string dataname="";
	std::string dirhiggs="";
	std::string path_dir="./bbggToolsResults";
  std::vector<string>folder_to_create;
	//name files
	std::string fileBkgName = "hgg.inputbkg_8TeV";
	std::string part="LT_GluGluTo";
	std::string part2="ToHHTo2B2G_M-";
	std::string end="_narrow_13TeV-madgraph.root";
  if(argc < 2) 
	{
		cout <<red<< "Please, parse json configuration file!"<<normal << endl;
		return 0;
  }
  if(argc < 3) 
	{
		cout <<yellow<< "Please, specify the folder which will hold the files!"<<normal << endl;
		cout <<yellow<<"Default folder : "<<path_dir<<" !"<<normal<<std::endl;
	}
  //Read config file
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( argv[1], pt );
  cout << "Reading input configuration file..." << endl;
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "" ) )
  {
  	cout << "Reading " << rowPair.first << " options..." << endl;
    if (rowPair.first == "signal") 
    {
      signalType = rowPair.second.get<std::string>("type");
      signalDir = rowPair.second.get<std::string>("dir");
      cardName = rowPair.second.get<std::string>("signalModelCard");
      BOOST_FOREACH (boost::property_tree::ptree::value_type const& v, rowPair.second.get_child("mass"))
      {  
								static int color=1;
                folder_to_create.push_back("/"+signalType+"_M"+v.second.data());
                std::string filename=signalDir+part+signalType+part2+v.second.data()+end;
                std::pair<int,std::string>a(stoi(v.second.data()),filename);
   							sigMass.push_back(a);
								if(color%2==0)cout<<grey;
								cout << "\t Signal type: " << signalType << endl;
     	 					cout << "\t Signal samples location: " << signalDir << endl;
								cout << "\t Signal model card: " << cardName << endl;
      					cout << "\t Signal mass: " << v.second.data() << endl;
								cout<<normal;
                color++;
      }
      if (mass == 0 ) 
			{
      	cout << "Mass == 0 means non-resonant analysis, therefore:" << endl;
        nonresFile = rowPair.second.get<std::string>("nonresFile");
        cout << "\t Non resonant file: " << nonresFile <<endl;
      }
     
    }
    if (rowPair.first == "other") 
		{
      lumi = rowPair.second.get<float>("integratedLumi");
      energy = rowPair.second.get<std::string>("energy");
      mass = rowPair.second.get<float>("higgsMass");
      addHiggs = rowPair.second.get<bool>("addHiggs");
      doblinding = rowPair.second.get<bool>("doBlinding");
      doBands = rowPair.second.get<bool>("doBands");
      NCAT = rowPair.second.get<int>("ncat");
			doBrazilianFlag=rowPair.second.get<bool>("doBrazilianFlag");
      if(NCAT>3)
      {
				std::cout<<red<<"Error NCAT>3 !!!"<<normal<<std::endl;
				std::exit(0);
      }
      version=rowPair.second.get<int>("version");
			docombine=rowPair.second.get<bool>("runCombine");
      useSigTheoryUnc = rowPair.second.get<bool>("useSigTheoryUnc");
      analysisType = rowPair.second.get<string>("analysisType");
      cout << "Running options: " << endl;
      cout << "\t Integrated luminosity: " << lumi << endl;
      cout << "\t Center of mass energy: " << energy << endl;
      cout << "\t Higgs mass: " << mass << endl;
      cout << "\t Add Higgs: " << addHiggs <<endl;
      cout << "\t Do blinded analysis: " << doblinding << endl;
      cout << "\t Calculate and show 1 and 2 sigma bands on background fit: " << doBands << endl;
      cout << "\t Number of categories to fit: " << NCAT << endl;
      cout << "\t Analysis type: " << analysisType << endl;
   		cout << "\t Add an uncertainty to the datacard for the SM diHiggs theory uncertainty: " << useSigTheoryUnc << endl;
    }
   	if (rowPair.first == "data") 
	 	{
    	datadir = rowPair.second.get<string>("dir");
			cout << "\t Data location: " << datadir<< endl;
      dataname = rowPair.second.get<string>("name");
			cout << "\t Data type: " << dataname << endl;
    }
    if (rowPair.first == "higgs") 
		{
			dirhiggs=rowPair.second.get<string>("dir");
			
			std::string type=rowPair.second.get<string>("type");
      BOOST_FOREACH (boost::property_tree::ptree::value_type const& v, rowPair.second.get_child("type"))
      { 
						static int  color=1;
      			higgstrue.push_back(v.second.data());
						if(color%2==0)cout<<grey;
            cout << "\t Higgs signal location: " << dirhiggs << endl;
						cout << "\t Higgs signal type: " << v.second.data() << endl;
						cout<<normal;
						color++;
			}

    }
	}
if(argc==3)
{
	path_dir=argv[2];
	path_dir+="_v"+to_string(version);
}
else path_dir+="v_"+to_string(version);

std::string fileBaseName = "hgg.mH"+to_string(mass)+"_8TeV";
boost::filesystem::path dire(path_dir);
boost::filesystem::create_directory(dire);
std::map<std::string,std::string>higgsfilename
{
		{"ggh_m125_powheg_8TeV","hgg.hig.mH"+to_string(mass)+"_8TeV.ggh"},
		{"tth_m125_8TeV","hgg.hig.mH"+to_string(mass)+"_8TeV.tth"},
		{"vbf_m125_8TeV","hgg.hig.mH"+to_string(mass)+"_8TeV.vbf"},
		{"wzh_m125_8TeV_zh","hgg.hig.mH"+to_string(mass)+"_8TeV.vh"},
		{"bbh_m125_8TeV","hgg.hig.mH"+to_string(mass)+"_8TeV.bbh"}
};
std::map<std::string,int>higgsNumber
{
		{"ggh_m125_powheg_8TeV",0},
		{"tth_m125_8TeV",1},
		{"vbf_m125_8TeV",2},
		{"wzh_m125_8TeV_zh",3},
		{"bbh_m125_8TeV",4}
};
  for(unsigned int i=0;i!=sigMass.size();++i)
  {
    int sigMas=sigMass[i].first;
    std::string signalDir2="";
    if(energy=="13TeV") signalDir2 = sigMass[i].second;
    else signalDir2=signalDir+"/v"+to_string(version)+"/v"+to_string(version)+"_"+analysisType+"/"+signalType+"_m"+to_string(sigMas)+"_8TeV_m"+to_string(sigMas)+".root";
    std::string folder_name=path_dir+"/"+signalType+"_M"+to_string(sigMas);
    std::string HLFactoryname=signalType+"_M"+to_string(sigMas);
		std::string ddata;
    if(energy=="13TeV")ddata=datadir+"LT_"+dataname+".root";
		else ddata=datadir+"/v"+to_string(version)+"/v"+to_string(version)+"_"+analysisType+"/"+dataname+"_m"+to_string(sigMas)+".root";
    cout<<"Data: "<<ddata<<endl;
    TString card_name(cardName); // put the model parameters here!
    HLFactory hlf(HLFactoryname.c_str(), card_name, false);
    RooWorkspace* w = hlf.GetWs();
    //Object
  	bbgg2DFitter TheFitter = bbgg2DFitter( w, sigMas, lumi,folder_name,energy,doblinding, NCAT, addHiggs);
  	TheFitter.style();
    
  	int opened=TheFitter.AddSigData( mass,signalDir2);
	if(opened==-1)
	{
		std::cout<<yellow<<" File "<<signalDir2<<" not found or not openeable !!"<<normal<<std::endl;
		std::cout<<yellow<<" File skipped !!"<<normal<<std::endl;
		continue;
	}
        boost::filesystem::path diree(path_dir+folder_to_create[i]);
 	if(boost::filesystem::create_directory(diree)) 
 	{
  	boost::filesystem::path dirwork(path_dir+folder_to_create[i]+"/workspaces");
		boost::filesystem::create_directory(dirwork);
		boost::filesystem::path dirdata(path_dir+folder_to_create[i]+"/datacards");
		boost::filesystem::create_directory(dirdata);
  	}
 	cout<<green<<"SIGNAL ADDED"<<normal<<endl;
  TheFitter.SigModelFit( mass); // constructing signal pdf
	cout<<green<<"SIGNAL FITTED"<<normal<<endl;
  	TheFitter.MakeSigWS( fileBaseName);
        cout<<green<<"SIGNAL'S WORKSPACE DONE"<<normal<<endl;
  	TheFitter.MakePlots( mass);
	cout<<green<<"SIGNAL'S PLOT DONE"<<normal<<endl;
	//
        for(unsigned int J=0;J!=higgstrue.size();++J)
	{
      std::string direc=dirhiggs+"/v"+to_string(version)+"/v"+to_string(version)+"_"+analysisType+"/"+higgstrue[J]+"_m"+to_string(sigMas)+".root";
  		TheFitter.AddHigData( mass,direc,higgsNumber[higgstrue[J]]);
		cout<<green<<"HIGGS "<<higgstrue[J]<<" ADDED"<<normal<<endl;
  		TheFitter.HigModelFit( mass,higgsNumber[higgstrue[J]]); // constructing higgs pdf
		cout<<green<<"HIGGS "<<higgstrue[J]<<" FITTED"<<normal<<endl;
  		TheFitter.MakeHigWS( higgsfilename[higgstrue[J]],higgsNumber[higgstrue[J]]);
		cout<<green<<"HIGGS' "<<higgstrue[J]<<"WORKSPACE DONE"<<normal<<endl;
	}		
        TheFitter.MakePlotsHiggs(mass,higgstrue,higgsNumber);
	cout<<green<<"HIGGS' PLOTS FITTED"<<normal<<endl;
	TheFitter.AddBkgData(ddata);
    	cout<<green<<"BKG ADDED"<<normal<<endl;
  	//TheFitter.PrintWorkspace();
        RooFitResult* fitresults = TheFitter.BkgModelFit( doBands,addHiggs,higgstrue,higgsNumber); // this is berestein 3
	cout<<green<<"BKG FITTED"<<normal<<endl;
        if(fitresults==nullptr)
	{
		std::cout<<"Problem with fitresults !!"<<std::endl;
		std::exit(3);
	}
  	TheFitter.MakeBkgWS( fileBkgName);
	cout<<green<<"BKG'S WORKSPACE DONE"<<normal<<endl;
  	// construct the models to fit
  	//
  	TheFitter.MakeDataCardonecatnohiggs( fileBaseName, fileBkgName, useSigTheoryUnc);
	cout<<green<<"DATACARD_ONE_CAT DONE"<<normal<<endl;
  	TheFitter.MakeDataCard( fileBaseName, fileBkgName,higgsfilename,useSigTheoryUnc,higgstrue,higgsNumber);
	cout<<green<<"DATACARD DONE"<<normal<<endl;
  	fitresults->Print();
  }
  if(docombine==true)
  {
	RunCombine(path_dir,doblinding);
  }
  if(doBrazilianFlag==true)
  {
	BrazilianFlag(path_dir);
  }
  return 0;		
} // close runfits
