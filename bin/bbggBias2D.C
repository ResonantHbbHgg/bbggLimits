//////////////////////////////////////////////////////////////////////////
//
// 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #505
// 
// Reading and writing ASCII configuration files
//
//
//
// 07/2008 - Wouter Verkerke 
// 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TString.h"
#include "TChain.h"
#include "TLegend.h"
#include "RooMinuit.h"
#include "TH1F.h"
#include "RooChi2Var.h"
#include "RooFitResult.h"
#include "RooCategory.h"
//#include "RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "TFile.h"

#include "HiggsAnalysis/bbggLimits/interface/bbggFittingTools.h"

//Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

using namespace RooFit;
namespace po = boost::program_options;


int main(int argc, const char* argv[])
{
  // C r e a t e  w o r k s p a c e 
  // -------------------------------
 
  if(argc < 2) {
    std::cout << "Please, parse JSON configuration file!" << std::endl;
    return 0;
  }

  //Parameters
  std::string rootName;
  std::string doBlinding;
  std::string blindingCut;
  std::string prefix;
  std::string plotsDir;
  std::string error;
  std::string selection;
  std::vector<std::string> vars;
  std::vector<std::string> varNames;
  std::vector<std::string> nbins;
  std::vector<std::string> functions;
  std::vector<std::string> functionsToFit;
  std::vector<std::string> functionsToPlot;
  std::vector<std::string> categories;
  std::vector<std::string> categoriesCuts;
  std::vector<std::string> biasFunctions;
  std::vector<TH1F> funcHistograms;
  //signal parameters
  std::string signal_File;
  std::string signalFit;
  std::vector<std::string> signal_Functions;
  std::vector<std::string> signal_functionsToFit;
  std::vector<std::string> signal_functionsToPlot;
  double signal_normalization_lumi_fb;
  double signal_normalization_totEvs;
  double signal_normalization_xsec;
  double bkgNorm_up, bkgNorm_down, bkgNorm;


  //Read json
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( argv[1], pt );
  rootName = pt.get_child("file").data();
  prefix = pt.get_child("prefix").data();
  plotsDir = pt.get_child("plotsDir").data();
  error = pt.get_child("error").data();
  selection = pt.get_child("selection").data();
  signal_File = pt.get_child("signal_File").data();
  signal_normalization_lumi_fb = TString(pt.get_child("signal_norm_lumi_fb").data().c_str()).Atof();
  signal_normalization_totEvs = TString(pt.get_child("signal_norm_totEvs").data().c_str()).Atof();
  signal_normalization_xsec = TString(pt.get_child("signal_norm_xsec").data().c_str()).Atof();
  bkgNorm_up = TString( pt.get_child("bkgNorm_up").data().c_str()).Atof();
  bkgNorm_down = TString( pt.get_child("bkgNorm_down").data().c_str()).Atof();
  bkgNorm = TString( pt.get_child("bkgNorm").data().c_str()).Atof();

  signalFit = pt.get_child("signalFit").data();
  
  std::cout << signal_normalization_lumi_fb << " " << signal_normalization_totEvs << " " << signal_normalization_xsec << std::endl;

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "vars" ) ){
    std::cout << rowPair.first << "\t" << rowPair.second.data() << std::endl;
    vars.push_back(rowPair.first);
    varNames.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "nbins" ) ){
    nbins.push_back(rowPair.second.data());
  }

  if(nbins.size() != vars.size()){
    std::cout << "Size of nbins vector must be the same as variables vector..." << std::endl;
    return 0;
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functions" ) ){
    functions.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functionsToFit" ) ){
    functionsToFit.push_back(rowPair.second.data());
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functionsToPlot" ) ){
    functionsToPlot.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "biasFunctions" ) ){
    biasFunctions.push_back(rowPair.second.data());
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "categories" ) ){
    categories.push_back( rowPair.first );
    categoriesCuts.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "blinding" ) ){
    blindingCut = rowPair.second.data();
    doBlinding = rowPair.first;
    break;
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "signal_Functions" ) ){
      signal_Functions.push_back(rowPair.second.data());
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "signal_functionsToFit" ) ){
      signal_functionsToFit.push_back(rowPair.second.data());
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "signal_functionsToPlot" ) ){
      signal_functionsToPlot.push_back(rowPair.second.data());
  }
  
  RooArgSet* TreeVars = new RooArgSet();
  
  for ( unsigned int cat = 0; cat < categories.size(); cat++) {
      RooWorkspace* w = new RooWorkspace("w");
      
      for( unsigned int i = 0; i < vars.size(); i++){
          // Initialize fit variable
          w->factory( vars[i].c_str() );
          std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
          std::cout << "Variable to be fit: " << sVar << std::endl;
          w->var( sVar.c_str() )->SetTitle(varNames[i].c_str());
          w->var( sVar.c_str() )->setUnit("GeV");
          TreeVars->add( *w->var( sVar.c_str() ) );
      }
      //Add extra things for TreeVars
//      RooRealVar* evWeight = new RooRealVar("evWeight","Total event weight",0,100,"");
//      TreeVars->add(*evWeight);
//      std::string cut = "1";
      
      //Get distribution from file
      TChain* treeTemp = new TChain("TCVARS");
      treeTemp->AddFile(TString(rootName));
      TTree* tree = (TTree*) treeTemp->CopyTree(std::string(selection + " && " + categoriesCuts[cat]).c_str() );
      std::cout << "Reading file: " << rootName << "\n\t with " << tree->GetEntries() << " entries \n\t Making the following cut: " << categoriesCuts[cat] << std::endl;
//      RooDataSet* data = new RooDataSet( "data", "data", tree, *TreeVars, cut.c_str(), "evWeight" );
      RooDataSet* data = new RooDataSet( "data", "data", tree, *TreeVars);
          
      //Get signal distribution from file
      TChain* treeTempSig = new TChain("TCVARS");
      treeTempSig->AddFile(TString(signal_File));
      std::string toCut = selection + " && " + categoriesCuts[cat];
      if(TString(categoriesCuts[cat]).Contains("< 0")) toCut = selection;
      TTree* treeSig = (TTree*) treeTempSig->CopyTree( toCut.c_str() );
      std::cout << "Reading signal file: " << signal_File << "\n\t with " << treeSig->GetEntries() << " entries \n\t Making the following cut: " << categoriesCuts[cat] << std::endl;
//      RooDataSet* dataSig = new RooDataSet( "dataSig", "dataSig", treeSig, *TreeVars, cut.c_str(), "evWeight" );
      RooDataSet* dataSig = new RooDataSet( "dataSig", "dataSig", treeSig, *TreeVars);
        
      w->import(*data);
      w->import(*dataSig);
      
      //Initialize all needed PDFs for signal
      for ( unsigned int f = 0; f < signal_Functions.size(); f++){
          TString thisFunction(signal_Functions[f]);
          std::cout << "Signal function added: " << thisFunction << std::endl;
          w->factory( thisFunction.Data() );
      }
      
      //Initialize all needed PDFs for background
      for ( unsigned int f = 0; f < functions.size(); f++){
        TString thisFunction(functions[f]);
        std::cout << "Function added: " << thisFunction << std::endl;
        w->factory( thisFunction.Data() );
      }
      
      std::vector<bbggFittingTools::FitRes> fitresults = bbggFittingTools::FitFunctions(w, signal_functionsToFit, dataSig, 1);
      
      for( unsigned int i = 0; i < vars.size(); i++){
          std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
          bbggFittingTools::PlotCurves(w, signal_functionsToPlot, fitresults, dataSig, sVar, nbins[i], plotsDir+"/signal"+sVar, 1, 0, -1);
          
          //Make signal parameters constant
          for( unsigned int ff = 0; ff < signal_functionsToFit.size(); ff++){
              TObjArray* funcNames = (TObjArray*) TString(signal_functionsToFit[ff]).Tokenize(":");
              const char* modelName = ((TObjString*) funcNames->At(0) )->String().Data();
          
              RooArgSet* params = (RooArgSet*) w->pdf( modelName )->getParameters( *w->var(sVar.c_str()) );
              TIterator* iter(params->createIterator());
              for (TObject *par = iter->Next(); par != 0; par = iter->Next()) {
                  RooRealVar *rrv = dynamic_cast<RooRealVar *>(par);
		  if( TString(rrv->GetName()).EqualTo("mjj") || TString(rrv->GetName()).EqualTo("mgg") ) continue;
		  std::cout << "##### Setting variable constant: " << rrv->GetName() << std::endl;
                  if(rrv) rrv->setConstant(kTRUE);
              }
          }
      }
      
      std::vector<bbggFittingTools::FitRes> bkgresults = bbggFittingTools::FitFunctions(w, functionsToFit, data);
      for( unsigned int i = 0; i < vars.size(); i++){
          std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
          bbggFittingTools::PlotCurves(w, functionsToPlot, bkgresults, data, sVar, nbins[i], plotsDir+"/bkg"+sVar, 0, 1);
      }
      
      //Do bias study business: create multipdf, etc
      RooCategory pdf_index("pdf_index","Index of Pdf which is active");
      RooArgList mypdfs;
      for( unsigned int bias = 0; bias < biasFunctions.size(); bias++){
	  std::cout << "Adding functions to multipdf! " << biasFunctions[bias].c_str() << std::endl;
	  const char* modelName = biasFunctions[bias].c_str();
	  mypdfs.add(* w->pdf( modelName ) );
      }
      RooMultiPdf multipdf("roomultipdf", "All Pdfs", pdf_index, mypdfs);
      RooRealVar norm("roomultipdf_norm", "Number of background events", bkgNorm, bkgNorm_down, bkgNorm_up);
//      RooRealVar norm("roomultipdf_norm", "Number of background events", 0.000001, 5000000);
      
      //Create bias workspace
      TFile * fout = new TFile("bias2d_pdfs.root", "RECREATE");
      RooWorkspace bW("Analysis2D", "Analysis2D");
      bW.import(pdf_index);
      bW.import(norm);
      bW.import(multipdf);
      bW.import( *w->pdf( signalFit.c_str() ) );
      bW.import(*data);
      bW.Print();
      bW.Write();
      fout->Close();
  }
  
  return 0;
}
