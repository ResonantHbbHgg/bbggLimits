//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 14 14:41:41 2016 by ROOT version 5.34/18
// from TTree bbggLTMaker/Flat tree for HH->bbgg analyses (after pre selection)
// found on file: /afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_sig/Hadd/output_GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph.root
//////////////////////////////////////////////////////////

#ifndef bbggLTMaker_h
#define bbggLTMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>

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


// Fixed size dimensions of array or collections stored in the TTree if any.

class bbggLTMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

//   Output file and tree
   TTree *outTree;
   TFile *outFile;
   Int_t           o_category;
   Double_t        o_normalization;
   Double_t	   o_preweight;
   Double_t        o_btagweight;
   Double_t        o_weight;
   Double_t        o_bbMass;
   Double_t        o_ggMass;
   Double_t        o_bbggMass;
   Double_t        o_phoevWeight;
   Double_t        jet1PT;
   Double_t	   jet2PT;
   Double_t	   jet1ETA;
   Double_t	   jet2ETA;
   std::string outFileName;
   double mtotMin;
   double mtotMax;
   double normalization;
   double btagWP;
   double btagWP_low;
   double btagWP_high;
   double cosThetaStarCut;
   int cosThetaStarCutCats;
   int photonCR;
   int doKinFit;
   int doMX;
   int tilt;
   int doNoCat;
   int doCatMixed;
   int doSingleCat;
   int bVariation;
   int phoVariation;
   double tiltWindow;
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
   Float_t	    CosThetaStar;
   Int_t	isSignal;
   Int_t	isPhotonCR;
   Int_t	leadingJet_flavour;
   Int_t	subleadingJet_flavour;
   Int_t	leadingJet_hadFlavour;
   Int_t	subleadingJet_hadFlavour;
   std::vector<std::pair<float,float>> btmap;

   // List of branches
   TBranch        *b_genWeights;   //!
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
   TBranch	  *b_isSignal;	//!
   TBranch	  *b_isPhotonCR;  //!
   TBranch	  *b_leadingJet_flavour;
   TBranch	  *b_subleadingJet_flavour;
   TBranch	  *b_leadingJet_hadFlavour;
   TBranch	  *b_subleadingJet_hadFlavour;
   TBranch	  *b_CosThetaStar;

   //Photon ID SF stuff
   TFile* photonidFile;
   TH2F* photonIDhist;
   TFile* csevFile;
   TH2F* csevhist;

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


   bbggLTMaker(TTree *tree=0);
   virtual ~bbggLTMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void SetMax( double max ){ mtotMax = max; }
   void SetMin( double min ){ mtotMin = min; }
   void SetNormalization(double norm) { normalization = norm; }
   void IsPhotonCR( int pcr ) { photonCR = pcr; }
   void IsMX( int mx ) { doMX = mx; }
   void IsKinFit( int kf ) { doKinFit = kf; }
   void SetOutFileName( std::string fname ) { outFileName = fname; }
   void SetBTagWP( double wp ) { btagWP = wp; }
   void SetBTagWP_High( double wp ) { btagWP_high= wp; }
   void SetBTagWP_Low( double wp ) { btagWP_low = wp; }
   void DoNoCat( int cat ) { doNoCat = cat; }
   void DoCatMixed( int cat ) { doCatMixed = cat; }
   void DoSingleCat( int cat ) { doSingleCat = cat; }
//   void SetTilt( int tt, double ttWind) { tilt = tt; tiltWindow = ttWind; }
   void SetTilt( int tt) { tilt = tt;}
   void DoBVariation( int tt) { bVariation = tt; }
   void BTagSetup(TString btagfile, TString effsfile);
   std::vector<std::pair<float,float>> BTagWeight(bbggLTMaker::LorentzVector jet1, int flavour1, bbggLTMaker::LorentzVector jet2, int flavour2, int variation=0);
   void SetupPhotonSF(TString idfile, TString evfile);
   float PhotonSF(LorentzVector pho, int phovar = 0);
   void DoPhoVariation(int tt) { phoVariation = tt;}
   void SetCosThetaStar(float cut, int cat) { cosThetaStarCut = cut; cosThetaStarCutCats = cat;}

};

#endif

#ifdef bbggLTMaker_cxx
bbggLTMaker::bbggLTMaker(TTree *tree) : fChain(0) 
{
   mtotMax = 1200.;
   mtotMin = 230.;
   normalization = 1.;
   photonCR = 0;
   doMX = 1;
   doKinFit = 0;
   outFileName = "LT_output.root";
   btagWP = 0.8;
   btagWP_low = 0.8;
   btagWP_high = 0.8;
   doNoCat = 0;
   doCatMixed = 0;
   doSingleCat = 0;
   bVariation = -999;
   phoVariation = -999;
   cosThetaStarCut = 10;
   cosThetaStarCutCats = 0;
   Init(tree);
}

bbggLTMaker::~bbggLTMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
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
   genWeights = 0;
   leadingPhotonID = 0;
   leadingPhotonISO = 0;
   subleadingPhotonID = 0;
   subleadingPhotonISO = 0;
   CosThetaStar = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
//   fChain->SetMakeClass(1);

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
   fChain->SetBranchAddress("isSignal", &isSignal, &b_isSignal);
   fChain->SetBranchAddress("isPhotonCR", &isPhotonCR, &b_isPhotonCR);
   fChain->SetBranchAddress("leadingJet_flavour", &leadingJet_flavour, &b_leadingJet_flavour);
   fChain->SetBranchAddress("subleadingJet_flavour", &subleadingJet_flavour, &b_subleadingJet_flavour);
   fChain->SetBranchAddress("leadingJet_hadFlavour", &leadingJet_hadFlavour, &b_leadingJet_hadFlavour);
   fChain->SetBranchAddress("subleadingJet_hadFlavour", &subleadingJet_hadFlavour, &b_subleadingJet_hadFlavour);
   fChain->SetBranchAddress("CosThetaStar", &CosThetaStar, &b_CosThetaStar);
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
