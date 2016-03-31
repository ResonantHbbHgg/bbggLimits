#include "HiggsAnalysis/bbggLimits/interface/bbggHighMassFitter.h"
#include "HiggsAnalysis/bbggLimits/interface/Colors.h"
#include <string>
#include "TString.h"
#include <vector>
#include "RooStats/HLFactory.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include<iostream>                                                                 
using namespace RooFit;
using namespace RooStats ;

int main(int argc, char *argv[]) 
{
	std::cout<<green<<"HIGH MASS FITTER "<<normal<<std::endl;
	std::exit(1);
	bbggHighMassFitter TheFitter;
        float mass=0;
	TString signalname;
        int signalsample=1;
	unsigned int NCAT=3;
    	if (signalsample==0)
      	{ 
		signalname="HH";
        }
        TString fileBaseName("CMS_jj_"+signalname+TString::Format("_%.0f_13TeV", mass));
        std::vector<std::string> cat_names;
        cat_names.push_back("CMS_jj_4btag_cat0");
        cat_names.push_back("CMS_jj_3btag_HPHP_cat1");
        cat_names.push_back("CMS_jj_3btag_HPLP_cat2");
        TString fileBkgName("CMS_jj_bkg_HH_13TeV");
        TString card_name("Xvv_reduced_models_Bkg_HH_13TeV.rs");
        HLFactory hlf("HLFactory", card_name, false);
        RooWorkspace* w = hlf.GetWs();
        RooFitResult* fitresults[NCAT];
	double MMIN=0;
	double MMAX=0;
        bool dobands=false;
        w->var("mgg")->setMin(MMIN);
        w->var("mgg")->setMax(MMAX);
	std::string filePOSTfix="";
	TString mainCut("1");
	std::string folder_name="";
        TheFitter.Initialize(w,MMIN,MMAX,filePOSTfix,0.0,0.0,0,mainCut,0.0,0.0,0.0,folder_name,3);
	// Add data to the workspace
        TheFitter.AddSigData(mass, signalsample,cat_names);
        TheFitter.AddBkgData(cat_names);
        TheFitter.SigModelFit(mass, signalname, cat_names);
        TheFitter.BkgModelFit(dobands,cat_names, fitresults);
        TheFitter.MakeSigWS(fileBaseName, signalname,cat_names);
        TheFitter.MakeBkgWS(fileBkgName,cat_names);
        for(unsigned int i=0;i!=NCAT;++i)
	{
		TheFitter.MakeDataCard_1Channel(fileBaseName, fileBkgName,i, signalname, signalsample, cat_names, mass);
        }
        TheFitter.MakePlots(mass, fitresults, signalname,cat_names);
	return 0;
}
