#ifndef bbggLTMaker_h
#define bbggLTMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH3F.h>
#include <iostream>
#include <fstream>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <cmath>
#include <Math/LorentzVector.h>
#include <algorithm>
#include <string>
#include <utility>

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


using namespace std;

// Total points for the Non-resonant re-weighting
//1507 (of initial grid) + 12 (benchmarks) = 1519
#define NRWTOT 1519

// Fixed size dimensions of array or collections stored in the TTree if any.

class bbggLTMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

//   Output file and tree
   TTree *outTree;
   TFile *outFile;
   
   ofstream f_event_dump;
   

   ULong64_t o_evt;
   UInt_t o_run;

   Float_t o_NRWeights[NRWTOT];
   TFile *NRwFile, *NRwFile2;
   TH2F *NR_Wei_Hists[NRWTOT];
   
   Int_t           o_category;
   Int_t           o_isSignal;
   Double_t        o_normalization;
   Double_t	   o_preweight;
   Double_t        o_btagweight;
   Double_t        o_weight;
   Double_t        o_bbMass;
   Double_t        o_ggMass;
   Double_t        o_bbggMass;
   Double_t        o_met;
   Double_t        o_phoevWeight;
   Double_t        o_ljet_bdis;
   Double_t        o_sjet_bdis;
   Double_t        o_HHTagger;
   Double_t        o_ttHTagger;
   Double_t        o_jt1diffweight;
   Double_t        o_jt2diffweight;
   Double_t        o_diffweight;
   Double_t        jet1PT;
   Double_t	   jet2PT;
   Double_t	   jet1ETA;
   Double_t	   jet2ETA;
   std::string outFileName;
   double mtotMin;
   double mtotMax;
   double normalization;
   double normalizationNR;
   double btagWP_loose;
   double btagWP_medium;
   double btagWP_tight;
   double cosThetaStarCutLow;
   double cosThetaStarCutHigh;
   double mvaCat0_lm;
   double mvaCat1_lm;
   double mvaCat0_hm;
   double mvaCat1_hm;
   double LowMassLeadingJetBtagCut;
   double HighMassLeadingJetBtagCut;
   double LowMassSubLeadingJetBtagCut;
   double HighMassSubLeadingJetBtagCut;
   int photonCR;
   int doKinFit;
   int doMX;
   int doNoCat;
   int doCatNonRes;
   int doCatLowMass;
   int doCatHighMass;
   int doCatMVA;
   int bVariation;
   int phoVariation;
   int trigVariation;
   int doNonResWeights;
   int photonCRNormToSig;
   int GenDiPhotonFilter;
   double massThreshold;
   bool isCustMVA;
   bool isRes;
   bool isETH;
   typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

   // Declaration of leaf types
   vector<double>   *genWeights;
   Double_t         genTotalWeight;
   LorentzVector    *leadingPhoton;
   vector<int>      *leadingPhotonID;
   vector<int>      *leadingPhotonISO;
   Int_t            leadingPhotonEVeto;
   LorentzVector    *subleadingPhoton;
   vector<int>      *subleadingPhotonID;
   vector<int>      *subleadingPhotonISO;
   Int_t            subleadingPhotonEVeto;
   LorentzVector    *diphotonCandidate;
   Int_t            nPromptInDiPhoton;
   LorentzVector    *leadingJet;
   LorentzVector    *leadingJet_KF;
   Float_t          leadingJet_bDis;
   LorentzVector    *subleadingJet;
   LorentzVector    *subleadingJet_KF;
   Float_t          subleadingJet_bDis;
   LorentzVector    *dijetCandidate;
   LorentzVector    *dijetCandidate_KF;
   LorentzVector    *diHiggsCandidate;
   LorentzVector    *diHiggsCandidate_KF;
   LorentzVector    *MET;
   Float_t	    CosThetaStar;
   Float_t	    CosThetaStar_CS;
   Float_t          ttHTagger;
   Float_t          HHTagger;
   Float_t          HHTagger_HM;
   Float_t          HHTagger_LM;
   Float_t	    leadingPhotonR9full5x5;
   Float_t	    subleadingPhotonR9full5x5;
   Int_t	isSignal;
   Int_t	isPhotonCR;
   Int_t	leadingJet_flavour;
   Int_t	subleadingJet_flavour;
   Int_t	leadingJet_hadFlavour;
   Int_t	subleadingJet_hadFlavour;
   Int_t        njets;
   //std::vector<std::pair<double,float>> btmap;
   TString myDiffOpt;

   Double_t gen_mHH, gen_cosTheta;
   ULong64_t event;
   UInt_t run;

   
   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!

   TBranch        *b_gen_mHH;   //!
   TBranch        *b_gen_cosTheta;   //!

   TBranch        *b_genWeights; //!
   TBranch        *b_genTotalWeight;   //!
   TBranch        *b_leadingPhoton;   //!
   TBranch        *b_leadingPhotonID;   //!
   TBranch        *b_leadingPhotonISO;   //!
   TBranch        *b_leadingPhotonEVeto;   //!
   TBranch        *b_subleadingPhoton;   //!
   TBranch        *b_subleadingPhotonID;   //!
   TBranch        *b_subleadingPhotonISO;   //!
   TBranch        *b_subleadingPhotonEVeto;   //!
   TBranch        *b_diphotonCandidate;   //!
   TBranch        *b_nPromptInDiPhoton;   //!
   TBranch        *b_leadingJet;   //!
   TBranch        *b_leadingJet_KF;   //!
   TBranch        *b_leadingJet_bDis;   //!
   TBranch        *b_subleadingJet;   //!
   TBranch        *b_subleadingJet_KF;   //!
   TBranch        *b_subleadingJet_bDis;   //!
   TBranch        *b_dijetCandidate;   //!
   TBranch        *b_dijetCandidate_KF;   //!
   TBranch        *b_diHiggsCandidate;   //!
   TBranch        *b_diHiggsCandidate_KF;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_njets;
   TBranch	  *b_isSignal;	//!
   TBranch	  *b_isPhotonCR;  //!
   TBranch	  *b_leadingJet_flavour;
   TBranch	  *b_subleadingJet_flavour;
   TBranch	  *b_leadingJet_hadFlavour;
   TBranch	  *b_subleadingJet_hadFlavour;
   TBranch	  *b_CosThetaStar;
   TBranch	  *b_CosThetaStar_CS;
   TBranch	  *b_leadingPhotonR9full5x5;
   TBranch	  *b_subleadingPhotonR9full5x5;
   TBranch        *b_ttHTagger;
   TBranch        *b_HHTagger;
   TBranch        *b_HHTagger_LM;
   TBranch        *b_HHTagger_HM;

   //Photon ID SF stuff
   TFile* photonidFile;
   TH2F* photonIDhist;
   TFile* csevFile;
   TH2F* csevhist;

   //Trigger stuff
   TFile* triggerFile;
   TH3F* ltriggerHist;
   TH3F* striggerHist;

   //BTaggin weights stuff
   BTagCalibration* calib;
   TFile* effsFile;
   TH2F* b_eff_medium;
   TH2F* b_eff_tight;
   TH2F* b_eff_loose;
   TH2F* c_eff_medium;
   TH2F* c_eff_tight;
   TH2F* c_eff_loose;
   TH2F* l_eff_medium;
   TH2F* l_eff_tight;
   TH2F* l_eff_loose;
   BTagCalibrationReader* b_reader_tight;
   BTagCalibrationReader* b_reader_tight_up;
   BTagCalibrationReader* b_reader_tight_down;
   BTagCalibrationReader* b_reader_medium;
   BTagCalibrationReader* b_reader_medium_up;
   BTagCalibrationReader* b_reader_medium_down;
   BTagCalibrationReader* b_reader_loose;
   BTagCalibrationReader* b_reader_loose_up;
   BTagCalibrationReader* b_reader_loose_down;
   BTagCalibrationReader* c_reader_tight;
   BTagCalibrationReader* c_reader_tight_up;
   BTagCalibrationReader* c_reader_tight_down;
   BTagCalibrationReader* c_reader_medium;
   BTagCalibrationReader* c_reader_medium_up;
   BTagCalibrationReader* c_reader_medium_down;
   BTagCalibrationReader* c_reader_loose;
   BTagCalibrationReader* c_reader_loose_up;
   BTagCalibrationReader* c_reader_loose_down;
   BTagCalibrationReader* l_reader_tight;
   BTagCalibrationReader* l_reader_tight_up;
   BTagCalibrationReader* l_reader_tight_down;
   BTagCalibrationReader* l_reader_medium;
   BTagCalibrationReader* l_reader_medium_up;
   BTagCalibrationReader* l_reader_medium_down;
   BTagCalibrationReader* l_reader_loose;
   BTagCalibrationReader* l_reader_loose_up;
   BTagCalibrationReader* l_reader_loose_down;

   //BTagCalibrationReader* b_diffreader_tight;
   //BTagCalibrationReader* c_diffreader_tight;
   //BTagCalibrationReader* l_diffreader_tight;

   BTagCalibrationReader* btag_reader;

   bbggLTMaker(TTree *tree=0, bool IsRes=0);
   virtual ~bbggLTMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

//mx usage
   void SetMax( double par ) { mtotMax = par; }
   void SetMin( double par ) { mtotMin = par; }
   void IsMX( int par ) { doMX = par; }
   void IsKinFit( int par ) { doKinFit = par; }
//normalization
   void SetNormalization(double norm, double norm2=1) { normalization = norm; normalizationNR = norm2;}
//categorization
   void DoNoCat( int cat ) { doNoCat = cat; }
   void DoCatNonRes( int cat ) { doCatNonRes = cat; }
   void DoCatLowMass( int cat ) { doCatLowMass = cat; }
   void DoCatHighMass( int cat ) { doCatHighMass = cat; }
   void DoCatMVA( int cat , float cat0_lm, float cat1_lm, float cat0_hm, float cat1_hm) { doCatMVA = cat; mvaCat0_lm = cat0_lm; mvaCat1_lm = cat1_lm; mvaCat0_hm = cat0_hm; mvaCat1_hm = cat1_hm;}
   void SetBTagWP_Tight( double par ) { btagWP_tight = par; }
   void SetBTagWP_Medium( double par ) { btagWP_medium = par; }
   void SetBTagWP_Loose( double par ) { btagWP_loose = par; }
//corrections
   void DoBVariation( int tt) { bVariation = tt; }
   void DoPhoVariation(int tt) { phoVariation = tt;}
   void DoTrigVariation(int tt) { trigVariation = tt;}
//other
   void IsPhotonCR( int pcr ) { photonCR = pcr; }
   void IsPhotonCRNormToSig( int pcr ) { photonCRNormToSig = pcr; }
   void SetCosThetaStarLow(float cut) { cosThetaStarCutLow = cut; }
   void SetCosThetaStarHigh(float cut) { cosThetaStarCutHigh = cut; }

//setup
   void BTagDiffWeightOpt( TString thisOpt) { myDiffOpt = thisOpt; }
   void SetOutFileName( std::string fname ) { outFileName = fname; }
   void BTagSetup(TString btagfile, TString effsfile);
   std::vector<std::pair<double,float>> BTagWeight(bbggLTMaker::LorentzVector jet1, int flavour1, bbggLTMaker::LorentzVector jet2, int flavour2, int variation=0);
   void SetupPhotonSF(TString idfile, TString evfile);
   float PhotonSF(LorentzVector pho, int phovar = 0);
   void SetMassThreshold(float par){ massThreshold = par;}
   void SetupTriggerSF(TString trig_file);
   float TriggerSF(LorentzVector lpho, float lr9, LorentzVector spho, float sr9, int var);

   void SetLowMassLeadingJetBtagCut ( double LMLJBTC ) { LowMassLeadingJetBtagCut = LMLJBTC;}
   void SetHighMassLeadingJetBtagCut ( double HMLJBTC ) { HighMassLeadingJetBtagCut = HMLJBTC;}
   void SetLowMassSubLeadingJetBtagCut ( double LMSJBTC ) { LowMassSubLeadingJetBtagCut = LMSJBTC;}
   void SetHighMassSubLeadingJetBtagCut ( double HMSJBTC ) { HighMassSubLeadingJetBtagCut = HMSJBTC;}

   void FilterGenDiPhotons() { GenDiPhotonFilter = 1;}
   
   void DoNRWeights(int doNRW) { doNonResWeights = doNRW; }

   void IsRes() {isRes=1;}
   void IsETH() {isETH=1;}

   void isCustomCatMVA() {isCustMVA = 1;}
};

#endif

#ifdef bbggLTMaker_cxx
bbggLTMaker::bbggLTMaker(TTree *tree, bool IsRes) : fChain(0)
{
   GenDiPhotonFilter = 0;
   isRes = IsRes;
   mtotMax = 12000.;
   mtotMin = 200.;
   normalization = 1.;
   normalizationNR = 1.;
   photonCR = 0;
   doMX = 1;
   doKinFit = 0;
   outFileName = "LT_output.root";
   doNoCat = 0;
   doCatNonRes = 0;
   doCatLowMass = 0;
   doCatHighMass = 0;
   doCatMVA = 0;
   mvaCat0_lm = -10;
   mvaCat1_lm = -10;
   mvaCat0_hm = -10;
   mvaCat1_hm = -10;
   btagWP_loose = 0.46;
   btagWP_medium = 0.8;
   btagWP_tight = 0.935;
   cosThetaStarCutLow = -100;
   cosThetaStarCutHigh = 100;
   bVariation = -999;
   phoVariation = -999;
   trigVariation = -999;
   doNonResWeights = 0;
   photonCRNormToSig = 0;
   massThreshold = 350;
   myDiffOpt = "central";
   LowMassLeadingJetBtagCut = -10;
   HighMassLeadingJetBtagCut = -10;
   LowMassSubLeadingJetBtagCut = -10;
   HighMassSubLeadingJetBtagCut = -10;
   njets = 0;
   Init(tree);
}

bbggLTMaker::~bbggLTMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   f_event_dump.close();
   //   photonidFile->Close();
   //   csevFile->Close();
   //   effsFile->Close();
}

Int_t bbggLTMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bbggLTMaker::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bbggLTMaker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   leadingPhoton = 0;
   subleadingPhoton = 0;
   diphotonCandidate = 0;
   leadingJet = 0;
   subleadingJet = 0;
   dijetCandidate = 0;
   diHiggsCandidate = 0;
   leadingJet_KF = 0;
   subleadingJet_KF = 0;
   dijetCandidate_KF = 0;
   diHiggsCandidate_KF = 0;
   MET = 0;
   genWeights = 0;
   leadingPhotonID = 0;
   leadingPhotonISO = 0;
   subleadingPhotonID = 0;
   subleadingPhotonISO = 0;
   CosThetaStar = 0;
   CosThetaStar_CS = 0;
   ttHTagger = 0;
   HHTagger = 0;
   HHTagger_LM = 0;
   HHTagger_HM = 0;
   isCustMVA=0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   //   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);

   fChain->SetBranchAddress("gen_mHH", &gen_mHH, &b_gen_mHH);
   fChain->SetBranchAddress("gen_cosTheta", &gen_cosTheta, &b_gen_cosTheta);

   fChain->SetBranchAddress("genWeights", &genWeights, &b_genWeights);
   fChain->SetBranchAddress("genTotalWeight", &genTotalWeight, &b_genTotalWeight);
   fChain->SetBranchAddress("leadingPhoton", &leadingPhoton, &b_leadingPhoton);
   fChain->SetBranchAddress("leadingPhotonID", &leadingPhotonID, &b_leadingPhotonID);
   fChain->SetBranchAddress("leadingPhotonISO", &leadingPhotonISO, &b_leadingPhotonISO);
   fChain->SetBranchAddress("leadingPhotonEVeto", &leadingPhotonEVeto, &b_leadingPhotonEVeto);
   fChain->SetBranchAddress("subleadingPhoton", &subleadingPhoton, &b_subleadingPhoton);
   fChain->SetBranchAddress("subleadingPhotonID", &subleadingPhotonID, &b_subleadingPhotonID);
   fChain->SetBranchAddress("subleadingPhotonISO", &subleadingPhotonISO, &b_subleadingPhotonISO);
   fChain->SetBranchAddress("subleadingPhotonEVeto", &subleadingPhotonEVeto, &b_subleadingPhotonEVeto);
   fChain->SetBranchAddress("diphotonCandidate", &diphotonCandidate, &b_diphotonCandidate);
   fChain->SetBranchAddress("nPromptInDiPhoton", &nPromptInDiPhoton, &b_nPromptInDiPhoton);
   fChain->SetBranchAddress("leadingJet", &leadingJet, &b_leadingJet);
   fChain->SetBranchAddress("leadingJet_KF", &leadingJet_KF, &b_leadingJet_KF);
   fChain->SetBranchAddress("leadingJet_bDis", &leadingJet_bDis, &b_leadingJet_bDis);
   fChain->SetBranchAddress("subleadingJet", &subleadingJet, &b_subleadingJet);
   fChain->SetBranchAddress("subleadingJet_KF", &subleadingJet_KF, &b_subleadingJet_KF);
   fChain->SetBranchAddress("subleadingJet_bDis", &subleadingJet_bDis, &b_subleadingJet_bDis);
   fChain->SetBranchAddress("dijetCandidate", &dijetCandidate, &b_dijetCandidate);
   fChain->SetBranchAddress("dijetCandidate_KF", &dijetCandidate_KF, &b_dijetCandidate_KF);
   fChain->SetBranchAddress("diHiggsCandidate", &diHiggsCandidate, &b_diHiggsCandidate);
   fChain->SetBranchAddress("diHiggsCandidate_KF", &diHiggsCandidate_KF, &b_diHiggsCandidate_KF);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("isSignal", &isSignal, &b_isSignal);
   fChain->SetBranchAddress("isPhotonCR", &isPhotonCR, &b_isPhotonCR);
   fChain->SetBranchAddress("leadingJet_flavour", &leadingJet_flavour, &b_leadingJet_flavour);
   fChain->SetBranchAddress("subleadingJet_flavour", &subleadingJet_flavour, &b_subleadingJet_flavour);
   fChain->SetBranchAddress("leadingJet_hadFlavour", &leadingJet_hadFlavour, &b_leadingJet_hadFlavour);
   fChain->SetBranchAddress("subleadingJet_hadFlavour", &subleadingJet_hadFlavour, &b_subleadingJet_hadFlavour);
   fChain->SetBranchAddress("CosThetaStar", &CosThetaStar, &b_CosThetaStar);
   fChain->SetBranchAddress("CosThetaStar_CS", &CosThetaStar_CS, &b_CosThetaStar_CS);
   fChain->SetBranchAddress("leadingPhotonR9full5x5", &leadingPhotonR9full5x5, &b_leadingPhotonR9full5x5);
   fChain->SetBranchAddress("subleadingPhotonR9full5x5", &subleadingPhotonR9full5x5, &b_subleadingPhotonR9full5x5);
   fChain->SetBranchAddress("ttHTagger", &ttHTagger, &b_ttHTagger);
   TString WhichTagger = "HHTagger";
   if(isETH) {
     WhichTagger = "HHTagger2017_transform";
     fChain->SetBranchAddress(TString(WhichTagger), &HHTagger, &b_HHTagger);
     //fChain->SetBranchAddress(TString(WhichTagger), &HHTagger_LM, &b_HHTagger_LM);
     //fChain->SetBranchAddress(TString(WhichTagger), &HHTagger_HM, &b_HHTagger_HM);
   }
   else {
     if(isRes) WhichTagger = "ResHHTagger";
     fChain->SetBranchAddress(WhichTagger, &HHTagger, &b_HHTagger);
     fChain->SetBranchAddress(TString(WhichTagger+"_LM"), &HHTagger_LM, &b_HHTagger_LM);
     fChain->SetBranchAddress(TString(WhichTagger+"_HM"), &HHTagger_HM, &b_HHTagger_HM);
   }
   Notify();
}

Bool_t bbggLTMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bbggLTMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bbggLTMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bbggLTMaker_cxx
