// C++ headers
#include "HiggsAnalysis/bbggLimits/interface/BrazilianFlag.h"
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
#include <boost/optional/optional.hpp>
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
        std::map<int,std::vector<float>> ParamsForFits;
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
	std::string part="LT_output_GluGluTo";
	std::string part2="ToHHTo2B2G_M-";
	std::string end="_narrow_13TeV-madgraph.root";
	float minMggMassFit=100;
	float maxMggMassFit=180;
	float minMjjMassFit=60;
	float maxMjjMassFit=180;
        float minSigFitMgg=115;
	float maxSigFitMgg=135;
  	float minSigFitMjj=60;
	float maxSigFitMjj=180;
        float minHigMggFit=115;
	float maxHigMggFit=135;
  	float minHigMjjFit=60;
	float maxHigMjjFit=180;
	bool HH=false;
	bool base=true;
	bool low=false;
	bool obs=false;
	bool twotag=false;
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
		if(color%2==0)cout<<grey;
      		cout << "\t Signal type: " << signalType << endl;
     	 	cout << "\t Signal samples location: " << signalDir << endl;
      		cout << "\t Signal model card: " << cardName << endl;
                cout << "\t Signal mass: " << v.second.data() << endl;
		TString fdtc = TString(signalType);
		fdtc.ReplaceAll("Bulk", "");
                folder_to_create.push_back("/"+std::string(fdtc.Data())+"_M"+v.second.data());
                std::string filename=TString(TString(signalDir).ReplaceAll( std::string("MASS"), std::string(v.second.data())) ).Data()+part+signalType+part2+v.second.data()+end;
                std::pair<int,std::string>a(stoi(v.second.data()),filename);
   		sigMass.push_back(a);
                std::string nameparam="param_"+v.second.data();
                boost::optional< const boost::property_tree::ptree& > child = rowPair.second.get_child_optional( nameparam.c_str() );
		if( child )
		{
			std::vector<float>Params;
  			BOOST_FOREACH (boost::property_tree::ptree::value_type const& v2, rowPair.second.get_child( nameparam.c_str() ))
      			{
				Params.push_back(stof(v2.second.data()));
			}
			if(Params.size()!=12)
			{
				std::cout<<yellow<<"\t Some parameters are missing for M"+v.second.data()+", I will use the default ones "<<normal<<std::endl;
			}
			else ParamsForFits[stoi(v.second.data())]=Params;
		}


		
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
      minMggMassFit=rowPair.second.get<float>("minMggMassFit");
      maxMggMassFit=rowPair.second.get<float>("maxMggMassFit");
      minMjjMassFit=rowPair.second.get<float>("minMjjMassFit");
      maxMjjMassFit=rowPair.second.get<float>("maxMjjMassFit");
      minSigFitMgg=rowPair.second.get<float>("minSigFitMgg");
      maxSigFitMgg=rowPair.second.get<float>("maxSigFitMgg");
      minSigFitMjj=rowPair.second.get<float>("minSigFitMjj");
      maxSigFitMjj=rowPair.second.get<float>("maxSigFitMjj");
      minHigMggFit=rowPair.second.get<float>("minHigMggFit");
      maxHigMggFit=rowPair.second.get<float>("maxHigMggFit");
      minHigMjjFit=rowPair.second.get<float>("minHigMjjFit");
      maxHigMjjFit=rowPair.second.get<float>("maxHigMjjFit");
      useSigTheoryUnc = rowPair.second.get<bool>("useSigTheoryUnc");
      analysisType = rowPair.second.get<string>("analysisType");
      HH=rowPair.second.get<bool>("HH");
      base=rowPair.second.get<bool>("base");
      low=rowPair.second.get<bool>("low");
      obs=rowPair.second.get<bool>("obs");
      twotag=rowPair.second.get<bool>("twotag");
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
      		if(addHiggs) higgstrue.push_back(v.second.data());
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
if(boost::filesystem::is_directory(path_dir)==true)
{
	std::cout<<red<<"Folder already exists. Please change the folder name or delete it "<<normal<<std::endl;
	std::exit(1);
}
std::string masswthout0 = std::to_string (mass);
masswthout0.erase ( masswthout0.find_last_not_of('0') + 1, std::string::npos );
if ((masswthout0.size () > 0)& (masswthout0.back()=='.'))  masswthout0.resize (masswthout0.size () - 1);
std::string fileBaseName = "hgg.mH"+masswthout0+"_8TeV";
boost::filesystem::path dire(path_dir);
boost::filesystem::create_directory(dire);
std::system(("cp "+std::string(argv[1])+" "+path_dir).c_str());
std::map<std::string,std::string>higgsfilename
{
		{"ggh_m125_powheg_8TeV","hgg.hig.mH"+masswthout0+"_8TeV.ggh"},
		{"tth_m125_8TeV","hgg.hig.mH"+masswthout0+"_8TeV.tth"},
		{"vbf_m125_8TeV","hgg.hig.mH"+masswthout0+"_8TeV.vbf"},
		{"wzh_m125_8TeV_zh","hgg.hig.mH"+masswthout0+"_8TeV.vh"},
		{"bbh_m125_8TeV","hgg.hig.mH"+masswthout0+"_8TeV.bbh"}
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
    TString newType = TString(signalType).ReplaceAll("Bulk", "");
    signalType = newType;
    if(energy=="13TeV") signalDir2 = sigMass[i].second;
    else signalDir2=signalDir+"/v"+to_string(version)+"/v"+to_string(version)+"_"+analysisType+"/"+signalType+"_m"+to_string(sigMas)+"_8TeV_m"+to_string(sigMas)+".root";
    std::string folder_name=path_dir+"/"+signalType+"_M"+to_string(sigMas);
    std::string HLFactoryname=signalType+"_M"+to_string(sigMas);
		std::string ddata;
    if(energy=="13TeV")ddata=std::string(TString( TString(datadir).ReplaceAll("MASS", TString::Itoa(sigMass[i].first, 10))).Data())+"LT_"+dataname+".root";
    else ddata=datadir+"/v"+to_string(version)+"/v"+to_string(version)+"_"+analysisType+"/"+dataname+"_m"+to_string(sigMas)+".root";
    cout<<"Data: "<<ddata<<endl;
    TString card_name(cardName); // put the model parameters here!
    HLFactory hlf(HLFactoryname.c_str(), card_name, false);
    RooWorkspace* w = hlf.GetWs();
    //Object
    bbgg2DFitter TheFitter;
    if(ParamsForFits.find(sigMas)!=ParamsForFits.end())
    {
	TheFitter.Initialize( w, sigMas, lumi,folder_name,energy,doblinding, NCAT, addHiggs,ParamsForFits[sigMas][0],ParamsForFits[sigMas][1],ParamsForFits[sigMas][2],ParamsForFits[sigMas][3],ParamsForFits[sigMas][4],ParamsForFits[sigMas][5],ParamsForFits[sigMas][6],ParamsForFits[sigMas][7],ParamsForFits[sigMas][8],ParamsForFits[sigMas][9],ParamsForFits[sigMas][10],ParamsForFits[sigMas][11]);
    }
    else TheFitter.Initialize( w, sigMas, lumi,folder_name,energy,doblinding, NCAT, addHiggs,minMggMassFit,maxMggMassFit,minMjjMassFit,maxMjjMassFit,minSigFitMgg,maxSigFitMgg,minSigFitMjj,maxSigFitMjj,minHigMggFit,maxHigMggFit,minHigMjjFit,maxHigMjjFit);
    TheFitter.style();
    TheFitter.SetType(signalType);
    
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
	if (addHiggs) {
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
        }
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
  	TheFitter.MakeDataCard( fileBaseName, fileBkgName,higgsfilename,useSigTheoryUnc,higgstrue,higgsNumber);
	cout<<green<<"DATACARD DONE"<<normal<<endl;
  	fitresults->Print();

	//MakeFancyDatacard
	//GetNumber of Observed Events:
	float sigExp[NCAT];
        float bkgObs[NCAT];
	for (int cc = 0; cc < NCAT; cc++){
		sigExp[cc] = -1;
		bkgObs[cc] = -1;
	}
	TString sigExpStr = TString("");
	TString bkgObsStr = TString("");
        for (int cc = 0; cc < NCAT; cc++) {
		sigExp[cc] = TheFitter.GetSigExpectedCats(cc);
		bkgObs[cc] = TheFitter.GetObservedCats(cc);
		sigExpStr += TString::Format("%f", sigExp[cc]);
		bkgObsStr += TString::Format("%f", bkgObs[cc]);
		if(cc < NCAT - 1) {
			sigExpStr += TString(",");
			bkgObsStr += TString(",");
		}
	}
	TString DCcommand = TString("python scripts/DataCardMaker.py -r -f ") + TString(folder_name) + TString::Format(" -n %d ", NCAT) + TString("-s ") + sigExpStr + TString(" -o ") + bkgObsStr;
	std::cout << DCcommand << std::endl;
	system(DCcommand);
  }
  if(docombine==true)
  {
	RunCombine(path_dir,doblinding);
  }
  if(doBrazilianFlag==true)
  {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(1) << lumi;
	system(TString("python scripts/ResPlotter.py -f" + path_dir+ " -l " + ss.str() ));
	BrazilianFlag(path_dir,HH,base,low,obs,twotag,energy,lumi);
  }
  return 0;		
} // close runfits
