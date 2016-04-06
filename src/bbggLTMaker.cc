#define bbggLTMaker_cxx
#include "HiggsAnalysis/bbggLimits/interface/bbggLTMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void bbggLTMaker::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L bbggLTMaker.C
//      Root > bbggLTMaker t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

   o_weight = 0;
   o_bbMass = 0;
   o_ggMass = 0;
   o_bbggMass = 0;
   o_category = -1;
   o_normalization = normalization;

   std::cout << "Output file name: " << outFileName << std::endl;
   std::cout << "Options:\n\t Mtot min: " << mtotMin << "\n\t Mtot max: " <<  mtotMax << "\n\t isPhotonCR: " << isPhotonCR << "\n\t Normalization: " << normalization << std::endl;
 
   outFile = new TFile(outFileName.c_str(), "RECREATE");
   outTree = new TTree("TCVARS", "Limit tree for HH->bbgg analyses");
   outTree->Branch("cut_based_ct", &o_category, "o_category/I"); //0: 2btag, 1: 1btag
   outTree->Branch("evWeight", &o_weight, "o_weight/D");
   outTree->Branch("mjj", &o_bbMass, "o_bbMass/D");
   outTree->Branch("mgg", &o_ggMass, "o_ggMass/D");
   outTree->Branch("mtot", &o_bbggMass, "o_bbggMass/D"); //

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

     if( jentry%1000 == 0 ) std::cout << "[bbggLTMaker::Process] Reading entry #" << jentry << endl;   
   
     o_weight = genTotalWeight*normalization;
     o_bbMass = dijetCandidate->M();
     o_ggMass = diphotonCandidate->M();
     o_bbggMass = diHiggsCandidate->M();
     if(doKinFit)
        o_bbggMass = diHiggsCandidate_KF->M();
     if(doMX)
        o_bbggMass = diHiggsCandidate->M() - dijetCandidate->M() + 125.;
  
     //mtot cut
     bool passedMassCut = 1;
     
     if(o_bbggMass < mtotMin || o_bbggMass > mtotMax)
	passedMassCut = 0;

     if(tilt){
	passedMassCut = 1;
	if( o_bbggMass > mtotMax + o_ggMass - 125.)
	  passedMassCut = 0;
	if( o_bbggMass < mtotMin + o_ggMass - 125.)
	  passedMassCut = 0;
     }

     if(passedMassCut == 0)
	continue;

     if(photonCR == 1 && isPhotonCR == 0)
	continue;
     if(photonCR == 0 && isPhotonCR == 1)
	continue;
   
//   double sumbtag = leadingJet_bDis + subleadingJet_bDis;
//   double upper = 1.83;
//   double lower = 1.11;
//   if ( sumbtag > upper ) o_category = 0;
//   if (sumbtag > lower && sumbtag < upper ) o_category = 1;
//   if (sumbtag < lower ) o_category = -1;
      if ( doNoCat == 0 )
      {
         if ( leadingJet_bDis > btagWP && subleadingJet_bDis > btagWP ) {o_category = 0;}
         else if ( leadingJet_bDis > btagWP && subleadingJet_bDis < btagWP ) { o_category = 1; }
         else if ( leadingJet_bDis < btagWP && subleadingJet_bDis > btagWP ) { o_category = 1; }
         else if ( leadingJet_bDis < btagWP && subleadingJet_bDis < btagWP ) { o_category = -1; }
       } 
       else if ( doNoCat == 1 )
       {
	 if ( leadingJet_bDis > btagWP && subleadingJet_bDis > btagWP ) {o_category = 0;}
	 if ( leadingJet_bDis > btagWP && subleadingJet_bDis < btagWP ) {o_category = 0;}
	 if ( leadingJet_bDis < btagWP && subleadingJet_bDis > btagWP ) {o_category = 0;}
	 if ( leadingJet_bDis < btagWP && subleadingJet_bDis < btagWP ) {o_category = -1;}
       }
      outTree->Fill();
   }

   std::cout << "Finished reading input tree..." << std::endl;

   outTree->Write();
   outFile->Close();
//   delete outTree;
//   delete outFile;
   std::cout << "Output tree/file saved and deleted..." << std::endl;
}
