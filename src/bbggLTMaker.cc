#define bbggLTMaker_cxx
#include "HiggsAnalysis/bbggLimits/interface/bbggLTMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

bool DEBUG = 0;

void bbggLTMaker::Loop()
{
  o_evt = 0;
  o_run = 0;

  for (UInt_t n=0; n<NRWTOT; n++) o_NRWeights[n]=0;

  o_weight = 0;
  o_preweight = 0;
  o_btagweight = 0;
  o_bbMass = 0;
  o_ggMass = 0;
  o_bbggMass = 0;
  o_category = -1;
  o_phoevWeight = 1;
  o_normalization = normalization;
  o_HHTagger = -10;
  o_ttHTagger = -10;
  o_ljet_bdis = -10;
  o_sjet_bdis = -10;
  o_isSignal = -10;
  o_jt1diffweight = -99;
  o_jt2diffweight = -99;
  o_diffweight = -99;

  std::cout << "Output file name: " << outFileName << std::endl;
  std::cout << "Options:\n\t Mtot min: " << mtotMin << "\n\t Mtot max: " <<  mtotMax << "\n\t isPhotonCR: " << photonCR
	    << "\n\t Normalization: " << normalization << std::endl;

  f_event_dump.open(outFileName+".event_dump.txt",ofstream::out);
  
  outFile = new TFile(outFileName.c_str(), "RECREATE");
  outTree = new TTree("TCVARS", "Limit tree for HH->bbgg analyses");
  outTree->Branch("cut_based_ct", &o_category, "o_category/I"); //0: 2btag, 1: 1btag
  outTree->Branch("isSignal", &o_isSignal, "o_isSignal/I"); //0: 2btag, 1: 1btag
  outTree->Branch("evWeight", &o_weight, "o_weight/D");
  outTree->Branch("preWeight", &o_preweight, "o_preweight/D");
  outTree->Branch("btagevWeight", &o_btagweight, "o_btagweight/D");
  outTree->Branch("phoevWeight", &o_phoevWeight, "o_phoevWeight/D");
  outTree->Branch("mjj", &o_bbMass, "o_bbMass/D");
  outTree->Branch("mgg", &o_ggMass, "o_ggMass/D");
  outTree->Branch("mtot", &o_bbggMass, "o_bbggMass/D");
  outTree->Branch("met", &o_met, "o_met/D");
  outTree->Branch("jet1PT", &jet1PT, "jet1PT/D");
  outTree->Branch("jet2PT", &jet2PT, "jet2PT/D");
  outTree->Branch("jet1ETA", &jet1ETA, "jet1ETA/D");
  outTree->Branch("jet2ETA", &jet2ETA, "jet2ETA/D");
  outTree->Branch("HHTagger", &o_HHTagger, "o_HHTagger/D");
  outTree->Branch("ljet_bdis", &o_ljet_bdis, "o_ljet_bdis/D");
  outTree->Branch("sjet_bdis", &o_sjet_bdis, "o_sjet_bdis/D");
  outTree->Branch("isSignal", &o_isSignal, "o_isSignal/I");
  outTree->Branch("jt1diffweight", &o_jt1diffweight, "o_jt1diffweight/D");
  outTree->Branch("jt2diffweight", &o_jt2diffweight, "o_jt2diffweight/D");
  outTree->Branch("diffweight", &o_diffweight, "o_diffweight/D");

  outTree->Branch("ttHTagger", &o_ttHTagger, "o_ttHTagger/D");

  outTree->Branch("evt", &o_evt, "o_evt/l");
  outTree->Branch("run", &o_run, "o_run/i");

  outTree->Branch("gen_mHH", &gen_mHH, "gen_mHH/D");
  outTree->Branch("gen_cosTheta", &gen_cosTheta, "gen_cosTheta/D");

  if (doNonResWeights){

    for (UInt_t n=0; n<NRWTOT; n++)
      outTree->Branch(Form("evWeight_NRW_%d",n), &o_NRWeights[n], Form("o_evWeight_NRW_%d/F",n));
  }
  if (fChain == 0) return;

  //cout<<"We are here"<<endl;
   
  // Get input histogrems for non-resonant re-weighting.
  if (doNonResWeights){

    // First check that input branches exis:
    if (!b_gen_mHH || !b_gen_cosTheta) {
      cout<<" Looks like you're trying to run non-resonant reweighting, but gen_mHH and gen_cosTheta branches don't exist in your input tree. \n"
	  <<" What a shame...  Fix that and then try again!"<<endl;
      exit(1);
    }

    // Now get the Histograms for re-weighting from the root file.
    TString fileNameWei = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/data/weights_v1_1507_points.root");
    NRwFile = new TFile(fileNameWei, "OPEN");

    // This file includes 12 benchmark v3 weights.
    TString fileNameWei2 = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/data/weights_v3_bench12_points.root");
    NRwFile2 = new TFile(fileNameWei2, "OPEN");
 
    if (NRwFile->IsZombie() || NRwFile2->IsZombie()){
      cout<<" Input file for Non-Res weights does not exist!"<<endl;
      exit(1);
    }
    NRwFile->Print();
    NRwFile2->Print();
    
    TList *histList = NRwFile->GetListOfKeys();
    for (UInt_t n=0; n<1507; n++){
      if (histList->Contains(Form("point_%i_weights",n)))
	NR_Wei_Hists[n] = (TH2F*)NRwFile->Get(Form("point_%i_weights",n));
      else
	cout<<"This one does not existe pas: "<<n<<endl;
    }

    TList *histList2 = NRwFile2->GetListOfKeys();
    for (UInt_t n=0; n<12; n++){
      if (histList2->Contains(Form("point_%i_weights",n)))
	NR_Wei_Hists[1507+n] = (TH2F*)NRwFile2->Get(Form("point_%i_weights",n));
      else
	cout<<"This one does not existe pas: "<<1507+n<<endl;
    }
  }

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  //setup btag weight
  if(bVariation > -100){
    cout << "Btagging SF variation: " << bVariation << endl;
    //TString bSF_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/BTag/CSVv2_Moriond17_G_H.csv");
    //TString bSF_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/BTag/CSVv2.csv");
    TString bSF_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/BTag/btagScaleFactors.txt");
    cout << "bSF file: " << bSF_file << endl;
    TString bEffs_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/BTag/btagEffs.root");
    cout << "bEffs file: " << bEffs_file << endl;
    bbggLTMaker::BTagSetup(bSF_file, bEffs_file);
    btag_reader = new BTagCalibrationReader(calib, BTagEntry::OP_RESHAPING, "iterativefit", std::string(myDiffOpt.Data()));
    cout << "btagging differential: " << bEffs_file << " " << myDiffOpt << std::endl;
  }

  if(phoVariation > -100){
    //    TString phoSFID_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/MVAID/mvaIDsf.root");
    TString phoSFID_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/MVAID/egammaEffi.txt_EGM2D.root");
    cout << "phoSFsID file: " << phoSFID_file << endl;
    //    TString phoSFeveto_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/MVAID/csevsf.root");
    TString phoSFeveto_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/MVAID/ScalingFactors_80X_Summer16.root");
    cout << "phoSFsEV file: " << phoSFeveto_file << endl;
    bbggLTMaker::SetupPhotonSF( phoSFID_file, phoSFeveto_file);
  }

  float trig_sf = 1;
  if(trigVariation > -100){
    TString trig_file = TString(std::getenv("CMSSW_BASE")) + TString("/src/HiggsAnalysis/bbggLimits/Weights/TriggerSF/TriggerSFs.root");
    cout << "TriggSF file: " << trig_file << endl;
    bbggLTMaker::SetupTriggerSF(trig_file);
  }

  Long64_t toRun = nentries;
  if (photonCRNormToSig) {
    std::cout << "Doing normalized photon control region!" << std::endl;
    Long64_t nentries_isSignal = fChain->GetEntries("isSignal==0");
    toRun = nentries_isSignal*0.14;
  }

  Long64_t pcrCounter = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if ( photonCRNormToSig ) {
      if (isSignal == 0) {
        pcrCounter++;
      }
      if (pcrCounter > toRun) break;
    }

    o_weight = 0;
    o_preweight = 0;
    o_btagweight = 0;
    o_bbMass = 0;
    o_ggMass = 0;
    o_bbggMass = 0;
    o_category = -1;
    o_phoevWeight = 1;
    o_ljet_bdis = 0;
    o_sjet_bdis = 0;
    o_isSignal = -10;
    o_HHTagger = -10;
    o_ttHTagger = -10;
    o_ljet_bdis = -10;
    o_sjet_bdis = -10;
    o_jt1diffweight = -99;
    o_jt2diffweight = -99;
    o_diffweight = -99;

      
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( jentry%1000 == 0 ) std::cout << "[bbggLTMaker::Process] Reading entry #" << jentry << endl;

    if( GenDiPhotonFilter && nPromptInDiPhoton < 2 ) continue;

    o_evt = event;
    o_run = run;
    
    //o_preweight = 1; //APZ TEMP. Make background MC look like data
    o_preweight = genTotalWeight*normalization;
    o_bbMass = dijetCandidate->M();
    o_ggMass = (*leadingPhoton + *subleadingPhoton).M();
    o_bbggMass = (*dijetCandidate + *leadingPhoton + *subleadingPhoton).M();
    o_ljet_bdis = leadingJet_bDis;
    o_sjet_bdis = subleadingJet_bDis;
    o_isSignal = isSignal;
    o_HHTagger = HHTagger;
    o_ttHTagger = ttHTagger;

    //std::cout<<"ttH tagger = "<<ttHTagger<<std::endl;
    //Double_t mmm = diHiggsCandidate->M();
    //if( jentry%200 == 0 )
    //if (abs(mmm-o_bbggMass)/o_bbggMass > 0.001)
    //std::cout << "Entry #" << jentry << "  run="<<run<<" event="<<event
    //		<<"   diHigg="<<mmm<<"  (diJe+pho+pho)="<<o_bbggMass<<endl;

    o_met = MET->Pt();
    //if( jentry%200 == 0 )
    // std::cout << "Entry #" << jentry << "  run="<<run<<" event="<<event
    //		<<"  Met = "<<o_met<<"  njets="<<njets<<std::endl;

    // Removing ttbar contributions
    //if (o_met > 100 || njets>=6) continue;


    if(doKinFit)
      o_bbggMass = diHiggsCandidate_KF->M();
    if(doMX)
      o_bbggMass = o_bbggMass - o_bbMass - o_ggMass + 250.;

    //mtot cut
    if(o_bbggMass < mtotMin || o_bbggMass > mtotMax) continue;

    // APZ DBG ETH BDT cuts
    if (isETH){
      HHTagger_LM = HHTagger;
      HHTagger_HM = HHTagger;
    }
    //std::cout<<"mbbjj="<<o_bbggMass<<" HHtagger: "<<HHTagger<<"  LM="<<HHTagger_LM<<"  HM="<<HHTagger_HM<<std::endl;
    

    bool isInterestingRegion = 1;
    if(isSignal){
      if(photonCR==0 && photonCRNormToSig==0) isInterestingRegion = 1;
      else if(photonCR==1 || photonCRNormToSig==1) isInterestingRegion = 0;
    } else if (!isSignal){
      if(photonCR==0 && photonCRNormToSig==0) isInterestingRegion = 0;
      else if(photonCR==1 || photonCRNormToSig==1) isInterestingRegion = 1;
    }
    if(!isInterestingRegion) continue;

    o_category = 2;

    //Calculate b-tagging scale factors
    //if(bVariation > -100) btmap = bbggLTMaker::BTagWeight(*leadingJet, leadingJet_hadFlavour, *subleadingJet, subleadingJet_hadFlavour, bVariation);
    double btagweight = -99;
    if(bVariation < -100) btagweight = 1;

    if(bVariation > -100) {
      /* -->>This is no longer used
       double bt1 = 1, bt2 = 1;
       if( leadingJet_bDis < 0.436) bt1 = (1 - btmap[0].first*btmap[0].second)/(1-btmap[0].second);
       if( leadingJet_bDis > 0.436 && leadingJet_bDis < 0.8 ) bt1 = (btmap[0].first*btmap[0].second - btmap[1].first*btmap[1].second)/(btmap[0].second - btmap[1].second);
       if( leadingJet_bDis > 0.8   && leadingJet_bDis < 0.92) bt1 = (btmap[1].first*btmap[1].second - btmap[2].first*btmap[2].second)/(btmap[1].second - btmap[2].second);
       if( leadingJet_bDis > 0.92) bt1 = btmap[2].first;

       if( subleadingJet_bDis < 0.436) bt2 = (1 - btmap[3].first*btmap[3].second)/(1-btmap[3].second);
       if( subleadingJet_bDis > 0.436 && subleadingJet_bDis < 0.8 ) bt2 = (btmap[3].first*btmap[3].second - btmap[4].first*btmap[4].second)/(btmap[3].second - btmap[4].second);
       if( subleadingJet_bDis > 0.8   && subleadingJet_bDis < 0.92) bt2 = (btmap[4].first*btmap[4].second - btmap[5].first*btmap[5].second)/(btmap[4].second - btmap[5].second);
       if( subleadingJet_bDis > 0.92) bt2 = btmap[5].first;

       btagweight = bt1*bt2;
       <<-- */


      // Code implementation following examples at:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration

      o_diffweight = 1.0;

      for (int j=0; j<2; j++){
	double pt,eta, csv;
	int flavor;
	if (j==0){
	  pt  = leadingJet->Pt();
	  eta = leadingJet->Eta();
	  flavor = leadingJet_hadFlavour;
	  csv = leadingJet_bDis;
	}
	else {
	  pt  = subleadingJet->Pt();
	  eta = subleadingJet->Eta();
	  flavor = subleadingJet_hadFlavour;
	  csv = subleadingJet_bDis;
	}
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;
	if( pt > 1000 ) pt = 999.;
	
	bool isBFlav = false;
	bool isCFlav = false;
	bool isLFlav = false;
	if( abs(flavor)==5 )      isBFlav = true;
	else if( abs(flavor)==4 ) isCFlav = true;
	else                      isLFlav = true;

	double my_jet_sf = 1.;

	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
	if( isBFlav )        jf = BTagEntry::FLAV_B;
	else if( isCFlav )   jf = BTagEntry::FLAV_C;
	else                 jf = BTagEntry::FLAV_UDSG;

	if(DEBUG) cout<<"syst:"<<myDiffOpt<<"  Jet flavour: "<<flavor<<" jf="<<jf<<"  Pt="<<pt<<"  eta="<<eta<<"  csv="<<csv<<endl;	
	
	if( myDiffOpt.EqualTo("central") )     my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_jes") && !isCFlav)        my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_jes") && !isCFlav)      my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_lf") && isBFlav )         my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_lf") && isBFlav )       my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_hf") && isLFlav )         my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_hf") && isLFlav )       my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_hfstats1") && isBFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_hfstats1") && isBFlav ) my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_hfstats2") && isBFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_hfstats2") && isBFlav ) my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_lfstats1") && isLFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_lfstats1") && isLFlav ) my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_lfstats2") && isLFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_lfstats2") && isLFlav ) my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_cferr1") && isCFlav )     my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_cferr1") && isCFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("up_cferr2") && isCFlav )     my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	else if( myDiffOpt.EqualTo("down_cferr2") && isCFlav )   my_jet_sf = btag_reader->eval(jf, eta, pt, csv);
	
	else my_jet_sf = 1;
	
	if (my_jet_sf < 0 || my_jet_sf > 10) my_jet_sf = 1.;
	o_diffweight *= my_jet_sf;
      }
    }

    if(doCatNonRes) {
      //       if ( o_category == 2 && leadingJet_bDis > btagWP_tight && subleadingJet_bDis > btagWP_tight) o_category = 0;
      //       if ( o_category == 2 && leadingJet_bDis > btagWP_loose  && leadingJet_bDis < btagWP_tight && subleadingJet_bDis > btagWP_tight) o_category = 1;
      //       if ( o_category == 2 && subleadingJet_bDis > btagWP_loose  && subleadingJet_bDis < btagWP_tight && leadingJet_bDis > btagWP_tight) o_category = 1;
      //       if ( o_category == 2 ) o_category = -1;
      if ( o_category == 2 && leadingJet_bDis > btagWP_medium && subleadingJet_bDis > btagWP_medium) o_category = 0;
      if ( o_category == 2 && leadingJet_bDis > btagWP_medium ) o_category = 1;
      if ( o_category == 2 && subleadingJet_bDis > btagWP_medium ) o_category = 1;
      if ( o_category == 2 && leadingJet_bDis < btagWP_medium && subleadingJet_bDis < btagWP_medium ) o_category = -1;
    }
    else if (doCatLowMass)
      {
       if (o_category == 2 && (leadingJet_bDis > btagWP_tight && subleadingJet_bDis > btagWP_loose)) o_category = 0;
       if (o_category == 2 && (subleadingJet_bDis > btagWP_tight && leadingJet_bDis > btagWP_loose)) o_category = 0;
       if (o_category == 2 && (leadingJet_bDis > btagWP_tight && subleadingJet_bDis < btagWP_loose)) o_category = 1;
       if (o_category == 2 && (subleadingJet_bDis > btagWP_tight && leadingJet_bDis < btagWP_loose)) o_category = 1;
       if (o_category == 2 && (leadingJet_bDis > btagWP_medium && leadingJet_bDis < btagWP_tight && subleadingJet_bDis > btagWP_medium && subleadingJet_bDis < btagWP_tight )) o_category = 1;
       if (o_category == 2 ) o_category = -1;
    }
    else if (doCatHighMass)
    {
       if (o_category == 2 && (leadingJet_bDis > btagWP_tight )) o_category = 0;
       if (o_category == 2 && (subleadingJet_bDis > btagWP_tight )) o_category = 0;
       if (o_category == 2 && (leadingJet_bDis > btagWP_medium && leadingJet_bDis < btagWP_tight && subleadingJet_bDis < btagWP_tight)) o_category = 1;
       if (o_category == 2 && (subleadingJet_bDis > btagWP_medium && subleadingJet_bDis < btagWP_tight && leadingJet_bDis < btagWP_tight)) o_category = 1;
       if ( o_category == 2 ) o_category = -1;
    } 
    else if (doCatMVA)
    {
       if (o_bbggMass > massThreshold ) {
//         if(!isRes && (leadingJet_bDis < HighMassLeadingJetBtagCut || subleadingJet_bDis < HighMassSubLeadingJetBtagCut)) o_category = -1;
         if((leadingJet_bDis < HighMassLeadingJetBtagCut || subleadingJet_bDis < HighMassSubLeadingJetBtagCut)) o_category = -1;
         if(o_category == 2 && HHTagger_HM > mvaCat0_hm) o_category = 0;
         if(o_category == 2 && HHTagger_HM > mvaCat1_hm && HHTagger_HM < mvaCat0_hm) o_category = 1;
         if(o_category == 2 && HHTagger_HM < mvaCat1_hm) o_category = -1;
       } else {
//         if(leadingJet->Pt() < 50) o_category = -1;
//         if(!isRes && (leadingJet_bDis < 0.57 || subleadingJet_bDis < 0.57)) o_category = -1;
//         if(!isRes && (leadingJet_bDis < LowMassLeadingJetBtagCut || subleadingJet_bDis < LowMassSubLeadingJetBtagCut)) o_category = -1;
         if((leadingJet_bDis < LowMassLeadingJetBtagCut || subleadingJet_bDis < LowMassSubLeadingJetBtagCut)) o_category = -1;
//         if(!isRes && (leadingJet_bDis < 0.8484 || subleadingJet_bDis < 0.8484)) o_category = -1;
         if(o_category == 2 && HHTagger_LM > mvaCat0_lm) o_category = 0;
         if(o_category == 2 && HHTagger_LM > mvaCat1_lm && HHTagger_LM < mvaCat0_lm) o_category = 1;
         if(o_category == 2 && HHTagger_LM < mvaCat1_lm) o_category = -1;
       }
    }

    if ( o_category == 2 ) std::cout << "ERROR ERROR ERROR ERROR ERROR ERROR ERROR" << std::endl;

    o_btagweight = btagweight;

    float pho1_sf = 1, pho2_sf = 1;
    if(phoVariation > -100){
      pho1_sf =  bbggLTMaker::PhotonSF(*leadingPhoton, phoVariation);
      pho2_sf =  bbggLTMaker::PhotonSF(*subleadingPhoton, phoVariation);
    }

    if(trigVariation > -100) {
      trig_sf = bbggLTMaker::TriggerSF(*leadingPhoton, leadingPhotonR9full5x5, *subleadingPhoton, subleadingPhotonR9full5x5, trigVariation);
      if(DEBUG) cout << "Trigger sf: " << trig_sf << " lpt " << leadingPhoton->Pt()<< " spt " << subleadingPhoton->Pt() << " leta "<< leadingPhoton->Eta() << " seta "<< subleadingPhoton->Eta() << " lr9 "<< leadingPhotonR9full5x5 << " sr9 "<< subleadingPhotonR9full5x5 << std::endl;
    }

    o_weight = o_preweight*o_btagweight*pho1_sf*pho2_sf*trig_sf;
    //doing differential btagging sfs if using MVA based categorization
    if (doCatMVA) {
      o_weight = o_preweight*o_diffweight*pho1_sf*pho2_sf*trig_sf;
    }

    if( o_preweight == 1) o_weight = 1; //if o_preweight == 1, this is data, no SF
    o_phoevWeight = pho1_sf*pho2_sf;
    jet1PT = leadingJet->pt();
    jet2PT = subleadingJet->pt();
    jet1ETA = leadingJet->eta();
    jet2ETA = subleadingJet->eta();

    
    // ----------------
    // -- Weights for Non-Res samples are added here
    //------------------------------
    
    if (doNonResWeights){
      //std::cout << "Doing Non-Resonant Signal weights " << std::endl;

      for (UInt_t n=0; n<NRWTOT; n++){
	if (n==324 || n==910 || n==985 || n==990){
	  // The points above do not exist in the input file provided by Alexandra
	  o_NRWeights[n] = 1;
	}
	else {
	  if (gen_mHH > 1800) // The re-weighting histograms are limited by this cut
	    o_NRWeights[n]=0;
	  else {
	    const UInt_t binNum = NR_Wei_Hists[n]->FindBin(gen_mHH, fabs(gen_cosTheta));
	    o_NRWeights[n] = NR_Wei_Hists[n]->GetBinContent(binNum);
	    // Just print out for one n:
	    if (DEBUG && n==100 && jentry%1000 == 0)
	      cout<<n<<" **  mHH = "<<gen_mHH<<"   cosT*="<<fabs(gen_cosTheta)<<"  bin="<<binNum<<" wei="<<o_NRWeights[n]<<endl;
	  }
	}
	
	o_NRWeights[n] = o_NRWeights[n]*normalizationNR*genTotalWeight*o_btagweight*pho1_sf*pho2_sf;
	
      }
    }
    // ---------------------------
    // Finished with Non-Res weights
    // ---------------------------
    
    
    //    std::cout << "cosThetaStarCutCats: " << cosThetaStarCutCats << " - - " << cosThetaStarCut << "  - -  " << o_category << std::endl;
    //    if( o_category == 1  && cosThetaStarCutCats > 0 && fabs(CosThetaStar) > fabs(cosThetaStarCut)) continue;
    //    if( o_category == 0  && cosThetaStarCutCats == 2 && fabs(CosThetaStar) > fabs(cosThetaStarCut)) continue;
    if( fabs(CosThetaStar_CS) > fabs(cosThetaStarCutHigh) ) continue;

    // Put Run number and event number in the file
    if (o_category != -1)
      f_event_dump<<run<<" "<<event<<endl;

    outTree->Fill();
  }

  std::cout << "Finished reading input tree..." << std::endl;
  //   effsFile->Close();
  outFile->cd();
  outTree->Write();
  outFile->Close();
  //   delete outTree;
  //   delete outFile;
  std::cout << "Output tree/file saved and deleted..." << std::endl;
}

void bbggLTMaker::SetupPhotonSF(TString idfile, TString evfile)
{
  photonidFile = new TFile(idfile, "READ");
//  photonIDhist = (TH2F*) photonidFile->Get("mva_id_sfs");
  photonIDhist = (TH2F*) photonidFile->Get("EGamma_SF2D;1");
  csevFile = new TFile(evfile, "READ");
//  csevhist = (TH2F*) csevFile->Get("csev_sfs");
  csevhist = (TH2F*) csevFile->Get("Scaling_Factors_CSEV_R9 Inclusive");
}

void bbggLTMaker::SetupTriggerSF(TString sfFile)
{
  triggerFile = new TFile(sfFile, "READ");
  ltriggerHist = (TH3F*) triggerFile->Get("leadingPhotonTSF");
  striggerHist = (TH3F*) triggerFile->Get("subleadingPhotonTSF");
}

float bbggLTMaker::PhotonSF(bbggLTMaker::LorentzVector pho, int phovar)
{
  float sf_id = -99, sf_ev = -99, sf_id_err = -99, sf_ev_err = -99;
  sf_id_err = photonIDhist->GetBinError( photonIDhist->FindBin(pho.eta(), pho.pt()) );
  sf_id = photonIDhist->GetBinContent( photonIDhist->FindBin(pho.eta(), pho.pt()) ) + phovar*sf_id_err;
  sf_ev_err = csevhist->GetBinContent( csevhist->FindBin(pho.eta(), pho.pt()) );
  if(sf_ev_err==1) sf_ev_err=0;
//  sf_ev = csevhist->GetBinContent( csevhist->FindBin(pho.eta(), pho.pt()) ) + phovar*sf_ev_err;
  sf_ev = csevhist->GetBinContent( csevhist->FindBin(fabs(pho.eta()), pho.pt()) ) + phovar*sf_ev_err;
  float totSF = sf_id*sf_ev;
  if (totSF < 1E-1) return 1;
  else return totSF;
}

void bbggLTMaker::BTagSetup(TString btagfile, TString effsfile)
{
  calib = new BTagCalibration("CSVv2", btagfile.Data());
  b_reader_tight = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "central");
  b_reader_tight_up = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "up");
  b_reader_tight_down = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "down");
  b_reader_medium = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "central");
  b_reader_medium_up = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "up");
  b_reader_medium_down = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "down");
  b_reader_loose = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "central");
  b_reader_loose_up = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "up");
  b_reader_loose_down = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "down");

  c_reader_tight = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "central");
  c_reader_tight_up = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "up");
  c_reader_tight_down = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "mujets", "down");
  c_reader_medium = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "central");
  c_reader_medium_up = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "up");
  c_reader_medium_down = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "down");
  c_reader_loose = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "central");
  c_reader_loose_up = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "up");
  c_reader_loose_down = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "mujets", "down");

  l_reader_tight = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "incl", "central");
  l_reader_tight_up = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "incl", "up");
  l_reader_tight_down = new BTagCalibrationReader(calib, BTagEntry::OP_TIGHT, "incl", "down");
  l_reader_medium = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "incl", "central");
  l_reader_medium_up = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "incl", "up");
  l_reader_medium_down = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "incl", "down");
  l_reader_loose = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "incl", "central");
  l_reader_loose_up = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "incl", "up");
  l_reader_loose_down = new BTagCalibrationReader(calib, BTagEntry::OP_LOOSE, "incl", "down");

  effsFile = new TFile(effsfile, "READ");
  b_eff_medium = (TH2F*) effsFile->Get("b_eff_medium");
  b_eff_tight = (TH2F*) effsFile->Get("b_eff_tight");
  b_eff_loose = (TH2F*) effsFile->Get("b_eff_loose");
  c_eff_medium = (TH2F*) effsFile->Get("c_eff_medium");
  c_eff_tight = (TH2F*) effsFile->Get("c_eff_tight");
  c_eff_loose = (TH2F*) effsFile->Get("c_eff_loose");
  l_eff_medium = (TH2F*) effsFile->Get("l_eff_medium");
  l_eff_tight = (TH2F*) effsFile->Get("l_eff_tight");
  l_eff_loose = (TH2F*) effsFile->Get("l_eff_loose");
}

std::vector<std::pair<double,float>> bbggLTMaker::BTagWeight(bbggLTMaker::LorentzVector jet1, int flavour1, bbggLTMaker::LorentzVector jet2, int flavour2, int variation)
{
  double tight1=-99, medium1=-99, loose1=-99, tightup1=-99, tightdown1=-99, mediumup1=-99, mediumdown1=-99, looseup1=-99, loosedown1=-99;
  double tight2=-99, medium2=-99, loose2=-99, tightup2=-99, tightdown2=-99, mediumup2=-99, mediumdown2=-99, looseup2=-99, loosedown2=-99;
  float jet1eff_medium=-99, jet1eff_tight=-99, jet1eff_loose=-99, jet2eff_medium=-99, jet2eff_tight=-99, jet2eff_loose=-99;

  if(flavour1==5){
    float jet1pt = min(max(jet1.Pt(), 30.1), 669.9);
    if(DEBUG) std::cout << "[bbggLTMaker::BTagWeight] Getting info from jet, pt " << jet1.Pt() << " actual pt " << jet1pt << " jet eta " << jet1.eta() << std::endl;
    tight1 = b_reader_tight->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt);
    medium1 = b_reader_medium->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt);
    loose1 = b_reader_loose->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt);

    float extraUnc = 1.;
    if(jet1.Pt() < 30.1 || jet1.Pt() > 669.9) extraUnc = 2.;

    tightup1 = tight1 + extraUnc*( b_reader_tight_up->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - tight1 );
    mediumup1 = medium1 + extraUnc*( b_reader_medium_up->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - medium1);
    looseup1 = loose1 + extraUnc*( b_reader_loose_up->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - loose1);

    tightdown1 = tight1 + extraUnc*(b_reader_tight_down->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - tight1 );
    mediumdown1 = medium1 + extraUnc*(b_reader_medium_down->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - medium1);
    loosedown1 = loose1 + extraUnc*(b_reader_loose_down->eval(BTagEntry::FLAV_B, jet1.Eta(), jet1pt) - loose1);

    if(DEBUG) std::cout << "[bbggLTMaker::BTagWeight] Now getting efficiencies from root file... " << std::endl;
    if(DEBUG) std::cout << "[bbggLTMaker::BTagWeight] Medium eff bin number: " << b_eff_medium->FindBin(fabs(jet1.Eta()), jet1pt) << std::endl;

    jet1eff_medium = b_eff_medium->GetBinContent( b_eff_medium->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_tight = b_eff_tight->GetBinContent( b_eff_tight->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_loose = b_eff_loose->GetBinContent( b_eff_loose->FindBin(fabs(jet1.Eta()), jet1pt) );
  }
  if(flavour1==4){
    float jet1pt = min(max(jet1.Pt(), 30.1), 669.9);

    tight1 = c_reader_tight->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt);
    medium1 = c_reader_medium->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt);
    loose1 = c_reader_loose->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt);

    float extraUnc = 1.;
    if(jet1.Pt() < 30.1 || jet1.Pt() > 669.9) extraUnc = 2.;

    tightup1 = tight1 + extraUnc*( c_reader_tight_up->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - tight1 );
    mediumup1 = medium1 + extraUnc*( c_reader_medium_up->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - medium1);
    looseup1 = loose1 + extraUnc*( c_reader_loose_up->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - loose1);

    tightdown1 = tight1 + extraUnc*(c_reader_tight_down->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - tight1 );
    mediumdown1 = medium1 + extraUnc*(c_reader_medium_down->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - medium1);
    loosedown1 = loose1 + extraUnc*(c_reader_loose_down->eval(BTagEntry::FLAV_C, jet1.Eta(), jet1pt) - loose1);

    jet1eff_medium = c_eff_medium->GetBinContent( c_eff_medium->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_tight = c_eff_tight->GetBinContent( c_eff_tight->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_loose = c_eff_loose->GetBinContent( c_eff_loose->FindBin(fabs(jet1.Eta()), jet1pt) );
  }
  if(flavour1<4){
    float jet1pt = min(jet1.Pt(), 999.9);
    tight1 = l_reader_tight->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt);
    medium1 = l_reader_medium->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt);
    loose1 = l_reader_loose->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt);

    float extraUnc = 1.;
    if(jet1.Pt() > 999.9) extraUnc = 2.;

    tightup1 = tight1 + extraUnc*( l_reader_tight_up->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - tight1 );
    mediumup1 = medium1 + extraUnc*( l_reader_medium_up->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - medium1);
    looseup1 = loose1 + extraUnc*( l_reader_loose_up->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - loose1);

    tightdown1 = tight1 + extraUnc*(l_reader_tight_down->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - tight1 );
    mediumdown1 = medium1 + extraUnc*(l_reader_medium_down->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - medium1);
    loosedown1 = loose1 + extraUnc*(l_reader_loose_down->eval(BTagEntry::FLAV_UDSG, jet1.Eta(), jet1pt) - loose1);

    jet1eff_medium = l_eff_medium->GetBinContent( l_eff_medium->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_tight = l_eff_tight->GetBinContent( l_eff_tight->FindBin(fabs(jet1.Eta()), jet1pt) );
    jet1eff_loose = l_eff_loose->GetBinContent( l_eff_loose->FindBin(fabs(jet1.Eta()), jet1pt) );
  }

  if(flavour2==5){
    float jet2pt = min(max(jet2.Pt(), 30.1), 669.9);
    tight2 = b_reader_tight->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt);
    medium2 = b_reader_medium->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt);
    loose2 = b_reader_loose->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt);

    float extraUnc = 1.;
    if(jet2.Pt() < 30.1 || jet2.Pt() > 669.9) extraUnc = 2.;

    tightup2 = tight2 + extraUnc*(b_reader_tight_up->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - tight2);
    mediumup2 = medium2 + extraUnc*(b_reader_medium_up->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - medium2);
    looseup2 = loose2 + extraUnc*(b_reader_loose_up->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - loose2);

    tightdown2 = tight2 + extraUnc*(b_reader_tight_down->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - tight2);
    mediumdown2 = medium2 + extraUnc*(b_reader_medium_down->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - medium2);
    loosedown2 = loose2 + extraUnc*(b_reader_loose_down->eval(BTagEntry::FLAV_B, jet2.Eta(), jet2pt) - loose2);

    jet2eff_medium = b_eff_medium->GetBinContent( b_eff_medium->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_tight = b_eff_tight->GetBinContent( b_eff_tight->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_loose = b_eff_loose->GetBinContent( b_eff_loose->FindBin(fabs(jet2.Eta()), jet2pt) );
  }
  if(flavour2==4){
    float jet2pt = min(max(jet2.Pt(), 30.1), 669.9);
    tight2 = c_reader_tight->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt);
    medium2 = c_reader_medium->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt);
    loose2 = c_reader_loose->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt);

    float extraUnc = 1.;
    if(jet2.Pt() < 30.1 || jet2.Pt() > 669.9) extraUnc = 2.;

    tightup2 = tight2 + extraUnc*(c_reader_tight_up->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - tight2);
    mediumup2 = medium2 + extraUnc*(c_reader_medium_up->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - medium2);
    looseup2 = loose2 + extraUnc*(c_reader_loose_up->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - loose2);

    tightdown2 = tight2 + extraUnc*(c_reader_tight_down->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - tight2);
    mediumdown2 = medium2 + extraUnc*(c_reader_medium_down->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - medium2);
    loosedown2 = loose2 + extraUnc*(c_reader_loose_down->eval(BTagEntry::FLAV_C, jet2.Eta(), jet2pt) - loose2);

    jet2eff_medium = c_eff_medium->GetBinContent( c_eff_medium->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_tight = c_eff_tight->GetBinContent( c_eff_tight->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_loose = c_eff_loose->GetBinContent( c_eff_loose->FindBin(fabs(jet2.Eta()), jet2pt) );
  }
  if(flavour2<4){
    float jet2pt = min(jet2.Pt(), 999.9);
    tight2 = l_reader_tight->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt);
    medium2 = l_reader_medium->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt);
    loose2 = l_reader_loose->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt);

    float extraUnc = 1.;
    if(jet2.Pt() > 999.9) extraUnc = 2.;

    tightup2 = tight2 + extraUnc*(l_reader_tight_up->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - tight2);
    mediumup2 = medium2 + extraUnc*(l_reader_medium_up->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - medium2);
    looseup2 = loose2 + extraUnc*(l_reader_loose_up->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - loose2);

    tightdown2 = tight2 + extraUnc*(l_reader_tight_down->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - tight2);
    mediumdown2 = medium2 + extraUnc*(l_reader_medium_down->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - medium2);
    loosedown2 = loose2 + extraUnc*(l_reader_loose_down->eval(BTagEntry::FLAV_UDSG, jet2.Eta(), jet2pt) - loose2);

    jet2eff_medium = l_eff_medium->GetBinContent( l_eff_medium->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_tight = l_eff_tight->GetBinContent( l_eff_tight->FindBin(fabs(jet2.Eta()), jet2pt) );
    jet2eff_loose = l_eff_loose->GetBinContent( l_eff_loose->FindBin(fabs(jet2.Eta()), jet2pt) );
  }

  std::vector<std::pair<double, float>> outMap;
  if(variation == 0){
    outMap.push_back( make_pair( loose1, jet1eff_loose ) );
    outMap.push_back( make_pair( medium1, jet1eff_medium ) );
    outMap.push_back( make_pair( tight1, jet1eff_tight ) );
    outMap.push_back( make_pair( loose2, jet2eff_loose ) );
    outMap.push_back( make_pair( medium2, jet2eff_medium ) );
    outMap.push_back( make_pair( tight2, jet2eff_tight ) );
  }
  if(variation == 1){
    outMap.push_back( make_pair( looseup1, jet1eff_loose ) );
    outMap.push_back( make_pair( mediumup1, jet1eff_medium ) );
    outMap.push_back( make_pair( tightup1, jet1eff_tight ) );
    outMap.push_back( make_pair( looseup2, jet2eff_loose ) );
    outMap.push_back( make_pair( mediumup2, jet2eff_medium ) );
    outMap.push_back( make_pair( tightup2, jet2eff_tight ) );
  }
  if(variation == -1){
    outMap.push_back( make_pair( loosedown1, jet1eff_loose ) );
    outMap.push_back( make_pair( mediumdown1, jet1eff_medium ) );
    outMap.push_back( make_pair( tightdown1, jet1eff_tight ) );
    outMap.push_back( make_pair( loosedown2, jet2eff_loose ) );
    outMap.push_back( make_pair( mediumdown2, jet2eff_medium ) );
    outMap.push_back( make_pair( tightdown2, jet2eff_tight ) );
  }
  return outMap;
}

float bbggLTMaker::TriggerSF(LorentzVector lpho, float lr9, LorentzVector spho, float sr9, int var)
{
   int leadingBin = ltriggerHist->FindBin(lr9, fabs(lpho.Eta()), lpho.Pt());
   int subleadingBin = striggerHist->FindBin(sr9, fabs(spho.Eta()), spho.Pt());
   if(DEBUG) cout << "[bbggLTMaker::TriggerSF] leading Bin: " << leadingBin << " subleadingBin " << subleadingBin << std::endl;

   float leadingSF = ltriggerHist->GetBinContent(leadingBin);
   float sleadingSF =striggerHist->GetBinContent(subleadingBin);
   if(DEBUG) cout << "[bbggLTMaker::TriggerSF] leading SF: " << leadingSF << " subleadingSF " << sleadingSF << std::endl;

   float leadingErr = ltriggerHist->GetBinError(leadingBin);
   float subleadingErr = striggerHist->GetBinError(subleadingBin);
   if(DEBUG) cout << "[bbggLTMaker::TriggerSF] leading SFerr: " << leadingErr << " subleadingSFerr " << subleadingErr << std::endl;

   return (leadingSF + var*leadingErr)*(sleadingSF + var*subleadingErr);
}
