//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 10:49:49 2015 by ROOT version 6.02/05
// from TTree bbggSelectionTree/Flat tree for HH->bbgg analyses (after pre selection)
// found on file: bbgg2DFitter_DoubleEG_37.root
//////////////////////////////////////////////////////////

#ifndef bbgg2DFitter_h
#define bbgg2DFitter_h
// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <unordered_map>
#include <map>
#include <algorithm>

// ROOT headers
#include <TROOT.h>
//#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGaxis.h>

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
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooMoment.h>

// namespaces
//using namespace std;
using namespace RooFit;
using namespace RooStats;


class RooWorkspace;
class RooRealVar;
class RooAbsPdf;
struct RooFitResult;
struct RooArgSet;

class bbgg2DFitter {

 private:

  std::vector<std::string> _singleHiggsNames;
  std::map<std::string,int> _singleHiggsMap;
  std::map<std::string,std::string> _singleHiggsWSfileNames;

  //   Parameters

  Int_t _verbLvl;

  Int_t _nonResWeightIndex;
  std::string _wName;
  Int_t _NR_MassRegion;

  Bool_t _doblinding;
  Int_t _NCAT;
  Int_t _sigMass;
  bool _addHiggs;
  float _lumi;
  TString _cut;
  std::string _energy;
  std::string _signalType;
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
  int _fitStrategy = 2;
  std::map<int,float> sigExpec;
  std::map<int,float> bkgExpec;
  std::map<int,float> dataObs;
  //Workspace
  RooWorkspace* _w;
  std::string _folder_name;

  TCanvas *_c1, *_c2;

 public :
   bbgg2DFitter() {}
   virtual ~bbgg2DFitter() { if (_c1) _c1->Close(); if (_c2) _c2->Close(); }
   void Initialize(RooWorkspace* workspace, Int_t SigMass, float Lumi,std::string folder_name,
		   std::string energy, Bool_t doBlinding, Int_t nCat, bool AddHiggs,
		   float minMggMassFit,float maxMggMassFit,float minMjjMassFit,float maxMjjMassFit,
		   float minSigFitMgg,float maxSigFitMgg,float minSigFitMjj,float maxSigFitMjj,
		   float minHigMggFit,float maxHigMggFit,float minHigMjjFit,float maxHigMjjFit,
		   Int_t doNRW=-2, std::string logFileName="");
   void SetVerbosityLevel(Int_t v) {_verbLvl=v;}
   void SetCut(TString cut) {_cut = cut;}
   void SetType(std::string tp) { _signalType = tp; }
   RooArgSet* defineVariables(); //DONE
   int AddSigData(float mass, TString signalfile); //DONE
   void AddHigData(float mass, TString signalfile, int higgschannel); //DONE
   void AddBkgData(TString datafile); //DONE
   void SigModelFit(float mass); //DONE
   void HigModelFit(float mass, int higgschannel); //DONE
   RooFitResult* BkgModelFit(Bool_t,bool); //DONE
   RooFitResult* BkgModelFit(Bool_t m,bool h,std::vector<std::string>higgstrue,std::map<std::string,int>higgsNumber) {
     return BkgModelFit(m, h);} //DONE
   void MakePlots(float mass); //DONE
   void SetFitStrategy( int st) { _fitStrategy = st; }
   void MakePlotsHiggs(float mass);
   void MakePlotsHiggs(float mass,std::vector<std::string>higgstrue,std::map<std::string,int>higgsNumber) {
     MakePlotsHiggs(mass);}
   void MakeSigWS(std::string filename); //DONE
   void MakeHigWS(std::string filename, int higgschannel); //DONE
   void MakeBkgWS(std::string filename); //DONE
   void MakeFitsForBias(std::string biasConfig, std::string outputFile);
   // const char* filenameh0, const char* filenameh1, const char* filenameh2, const char* filenameh4);
   void MakeDataCard(std::string filename, std::string filename1, Bool_t ); //DONE
   void MakeDataCard(std::string filename, std::string filename1,std::map<std::string,std::string>higgsfilename, Bool_t bb,std::vector<std::string>,std::map<std::string,int>higgsNumber) {
     MakeDataCard(filename, filename1, bb);
   }

   void SetConstantParams(const RooArgSet* params); //DONE
   void PrintWorkspace();// {_w->Print("v");}
   TStyle * style(); //DONE
   void SetSigExpectedCats(int cat, float expec) {
     if(sigExpec.find(cat) != sigExpec.end() ){std::cout << "[SetSigExpectedCats] Cat already set!" << std::endl;} else { sigExpec[cat] = expec; }}
   void SetBkgExpectedCats(int cat, float expec) {
     if(bkgExpec.find(cat) != bkgExpec.end() ){std::cout << "[SetBkgExpectedCats] Cat already set!" << std::endl;} else { bkgExpec[cat] = expec; }}
   void SetObservedCats(int cat, float observ) {
     if(dataObs.find(cat) != dataObs.end() ){std::cout << "[DataObservedCats] Cat already set!" << std::endl;} else { dataObs[cat] = observ; }}

   float GetSigExpectedCats(int cat) {
     if(sigExpec.find(cat) == sigExpec.end() ){std::cout << "[GetSigExpectedCats] Cat not found! Cat=" <<cat<< std::endl; return -1;}
     else { return sigExpec[cat]; }}
   float GetBkgExpectedCats(int cat) {
     if(bkgExpec.find(cat) == bkgExpec.end() ){std::cout << "[GetBkgExpectedCats] Cat not found! Cat="<<cat<< std::endl; return -1;}
     else { return bkgExpec[cat]; }}
   float GetObservedCats(int cat) {
     if(dataObs.find(cat) == dataObs.end() ){std::cout << "[GetObservedCats] Cat not found! Cat ="<<cat << std::endl; return -1;}
     else { return dataObs[cat]; }}

   std::vector<float> EffectiveSigma(RooRealVar* mass, RooAbsPdf* binned_pdf, float wmin, float wmax, float step, float epsilon);

   ClassDef(bbgg2DFitter,0);
};

#endif
