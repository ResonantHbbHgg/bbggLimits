//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 10:49:49 2015 by ROOT version 6.02/05
// from TTree bbggSelectionTree/Flat tree for HH->bbgg analyses (after pre selection)
// found on file: bbgg2DFitter_DoubleEG_37.root
//////////////////////////////////////////////////////////

#ifndef bbgg2DFitter_h
#define bbgg2DFitter_h
// C++ headers
//#include <iostream>
//#include <sstream>
#include <string>
//#include <cmath>
// ROOT headers
#include <TROOT.h>
//#include <TSystem.h>
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
//#include <RooWorkspace.h>
//#include <RooFitResult.h>
//#include <RooRealVar.h>
//#include <RooCategory.h>
//#include <RooArgSet.h>
//#include <RooStats/HLFactory.h>
//#include <RooDataSet.h>
//#include <RooFormulaVar.h>
//#include <RooGenericPdf.h>
//#include <RooPlot.h>
//#include <RooAbsPdf.h>
//#include <RooBernstein.h>
//#include <RooExtendPdf.h>
//#include <RooMinimizer.h>
//#include <RooStats/RooStatsUtils.h>
//#include <RooProdPdf.h>
//#include <RooExponential.h>
//#include <RooPolynomial.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
//#include <cmath>
//#include <Math/LorentzVector.h>
//#include <algorithm>
#include <string>
//#include <utility>

// namespaces
//using namespace std;
//using namespace RooFit;
//using namespace RooStats;

class RooWorkspace;
struct RooFitResult;
struct RooArgSet;

class bbgg2DFitter {
public :
//   Parameters
   Bool_t _doblinding;
   Int_t _NCAT;
   Int_t _sigMass;
   bool _addHiggs;
   float _lumi;
   TString _cut;
   std::string _energy;
   float _minMggMassFit;
   float _maxMggMassFit;
   float _minMjjMassFit;
   float _maxMjjMassFit;
   float _minSigFitMgg;
   float _maxSigFitMgg;
  	float _minSigFitMjj;
	float _maxSigFitMjj;
 float _minHigMggFit;
	float _maxHigMggFit;
  float _minHigMjjFit;
	float _maxHigMjjFit;
   //Workspace
   RooWorkspace* _w;
   std::string _folder_name;
   bbgg2DFitter() {}
   void Initialize(RooWorkspace* workspace, Int_t SigMass, float Lumi,std::string folder_name,std::string energy, Bool_t doBlinding, Int_t nCat, bool AddHiggs,float minMggMassFit,float maxMggMassFit,float minMjjMassFit,float maxMjjMassFit,float minSigFitMgg,float maxSigFitMgg,float minSigFitMjj,float maxSigFitMjj,float minHigMggFit,float maxHigMggFit,float minHigMjjFit,float maxHigMjjFit);
   virtual ~bbgg2DFitter() { }
   void SetCut(TString cut) {_cut = cut;}
   RooArgSet* defineVariables(); //DONE
   int AddSigData(float mass, TString signalfile); //DONE
   void AddHigData(float mass, TString signalfile, int higgschannel); //DONE
   void AddBkgData(TString datafile); //DONE
   void SigModelFit(float mass); //DONE
   void HigModelFit(float mass, int higgschannel); //DONE
   RooFitResult* BkgModelFit(Bool_t,bool,std::vector<std::string>higgstrue,std::map<std::string,int>higgsNumber); //DONE
   void MakePlots(float mass); //DONE
   void MakePlotsHiggs(float mass,std::vector<std::string>higgstrue,std::map<std::string,int>higgsNumber); //DONE
   void MakeSigWS(std::string filename); //DONE
   void MakeHigWS(std::string filename, int higgschannel); //DONE
   void MakeBkgWS(std::string filename); //DONE
   // const char* filenameh0, const char* filenameh1, const char* filenameh2, const char* filenameh4);
   void MakeDataCard(std::string filename, std::string filename1,std::map<std::string,std::string>higgsfilename, Bool_t,std::vector<std::string>,std::map<std::string,int>higgsNumber); //DONE
   void SetConstantParams(const RooArgSet* params); //DONE
   void PrintWorkspace();// {_w->Print("v");}
   TStyle * style(); //DONE
   
   ClassDef(bbgg2DFitter,0);
};

#endif

/*
#ifdef bbgg2DFitter_cxx
bbgg2DFitter::bbgg2DFitter(RooWorkspace* workspace, Int_t SigMass, float Lumi,std::string folder_name,std::string energy, Bool_t doBlinding = false, Int_t nCat = 0, bool AddHiggs = true)
{
    _doblinding = doBlinding;
    _NCAT = nCat;
    _sigMass = SigMass;
    _addHiggs = AddHiggs;
    _w = new RooWorkspace(*workspace);
    _lumi = Lumi;
    _cut = "1";
    _folder_name=folder_name;
    _energy=energy;
    
}

#endif // #ifdef bbgg2DFitter_cxx
*/
