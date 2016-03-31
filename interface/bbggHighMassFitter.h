#ifndef bbggHighMassFitter_h
#define bbggHighMassFitter_h

// C++ headers
//#include <iostream>
//#include <sstream>
#include <string>
//#include <cmath>
// ROOT headers
#include <RooHist.h>
#include <TROOT.h>
//#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <RooFitResult.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <vector>
#include <RooWorkspace.h>
#include <RooHistPdf.h>
#include <RooCurve.h>
#include <RooPlot.h>
#include <TGraph.h>
class bbggHighMassFitter {
public :
   //parameters
   Int_t _NCAT;
   Double_t _MMIN /*= 1000.*/;
   Double_t _MMAX /*= 3000*/;
   std::string _filePOSTfix/*=""*/;
   double _analysisLumi /*= 2.1977*/; // Luminosity you use in your analysis
   double _nEventsInSignalMC /*= 0.*/; //number of events in Signal MC sample
   int _iGraviton /*= 0*/;
   TString _mainCut/*("1")*/;
   double _signalScaler/*=0*/;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb
   double _scaleFactorHP/*=1*/;// already done on 1 GeV Histo level
   double _scaleFactorLP/*=1*/;// already done on 1 GeV Histo level
   RooWorkspace* _w;
   std::string _folder_name;
   void Initialize(RooWorkspace* w,Double_t MMIN,Double_t MMAX,std::string filePOSTfix,double analysisLumi,double nEventsInSignalMC,int iGraviton,TString mainCut,double signalScaler,double scaleFactorHP, double scaleFactorLP,std::string folder_name,Int_t NCAT);
   bbggHighMassFitter() {}
   virtual ~bbggHighMassFitter();
   RooArgSet* defineVariables();
   void AddSigData(Float_t, int, std::vector<std::string>);
   void AddBkgData(std::vector<std::string>);
   void SigModelFit(Float_t, TString signalname, std::vector<std::string>);
   void BkgModelFit(Bool_t, std::vector<std::string>, RooFitResult** fitresults);
   void MakePlots(Float_t, RooFitResult** , TString signalname, std::vector<std::string>);
   void MakeSigWS(const char* filename, TString signalname, std::vector<std::string>);
   void MakeBkgWS(const char* filename, std::vector<std::string>);
   void SetConstantParams(const RooArgSet* params);
   void MakeDataCard_1Channel(const char* fileBaseName, const char* fileBkgName, int iChan, TString signalname, int signalsample, std::vector<std::string> cat_names, double mass);
   //Double_t effSigma(TH1 *hist);
   TStyle * style(); //DONE
   //ClassDef(bbggHighMassFitter,0);
};

#endif

