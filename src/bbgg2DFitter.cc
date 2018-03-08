#define bbgg2DFitter_cxx
#include "HiggsAnalysis/bbggLimits/interface/bbgg2DFitter.h"
#include "HiggsAnalysis/bbggLimits/interface/bbggFittingTools.h"
//#include "HiggsAnalysis/bbggLimits/interface/Colors.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
//Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

std::ofstream newCout;


std::vector<float> bbgg2DFitter::EffectiveSigma(RooRealVar* mass, RooAbsPdf* binned_pdf, float wmin, float wmax, float step=0.002, float epsilon=1.e-4)
{
  RooAbsReal* binned_cdf = (RooAbsReal*) binned_pdf->createCdf(*mass);

  float point = wmin;
  std::vector<std::pair<float,float>> points;
  while( point <= wmax) {
    mass->setVal( point );
    if ( binned_pdf->getVal() > epsilon ) {
      points.push_back( std::make_pair(point, binned_cdf->getVal() ) );
    }
    point += step;
  }

  float low = wmin;
  float high = wmax;
  float width = wmax - wmin;

  for( unsigned int ip = 0; ip < points.size(); ip++) {
    for( unsigned int jp = ip; jp < points.size(); jp++) {
      float wy = points[jp].second - points[ip].second;
      if ( fabs( wy - 0.683 ) < epsilon ) {
	float wx = points[jp].first - points[ip].first;
	if ( wx < width ) {
	  low = points[ip].first;
	  high = points[jp].first;
	  width = wx;
	}
      }
    }
  }

  if (_verbLvl>1)
    std::cout << "#Sigma effective: xLow: " << low << ", xHigh: " << high << ", width: " << width << std::endl;

  std::vector<float> outVec;
  outVec.push_back(width);
  outVec.push_back(low);
  outVec.push_back(high);
  outVec.push_back(width/2.);

  return outVec;
}


void bbgg2DFitter::PrintWorkspace() {_w->Print("v");}

void bbgg2DFitter::Initialize(RooWorkspace* workspace, Int_t SigMass, float Lumi,std::string folder_name,
			      std::string energy, Bool_t doBlinding, Int_t nCat, bool AddHiggs,
			      float minMggMassFit,float maxMggMassFit,float minMjjMassFit,float maxMjjMassFit,
			      float minSigFitMgg,float maxSigFitMgg,float minSigFitMjj,float maxSigFitMjj,
			      float minHigMggFit,float maxHigMggFit,float minHigMjjFit,float maxHigMjjFit,
			      Int_t doNRW, std::string logFileName, bool doARW)
{
  //std::cout<<"DBG.  We Initialize..."<<std::endl;

  _doblinding = doBlinding;
  _NCAT = nCat;
  _sigMass = SigMass;
  _addHiggs = AddHiggs;
  _w = new RooWorkspace(*workspace);
  _lumi = Lumi;
  _cut = "1";
  _folder_name=folder_name;
  _energy=energy;
  _minMggMassFit=minMggMassFit;
  _maxMggMassFit=maxMggMassFit;
  _minMjjMassFit=minMjjMassFit;
  _maxMjjMassFit=maxMjjMassFit;
  _minSigFitMgg=minSigFitMgg;
  _maxSigFitMgg=maxSigFitMgg;
  _minSigFitMjj=minSigFitMjj;
  _maxSigFitMjj=maxSigFitMjj;
  _minHigMggFit=minHigMggFit;
  _maxHigMggFit=maxHigMggFit;
  _minHigMjjFit=minHigMjjFit;
  _maxHigMjjFit=maxHigMjjFit;
  TGaxis::SetMaxDigits(3);
  _doARW = doARW;
//  doARW = 0;

  _singleHiggsNames = {"ggh_m125_powheg_13TeV","tth_m125_13TeV"};
  _singleHiggsMap = {
    {"ggh_m125_powheg_13TeV",0},
    {"tth_m125_13TeV",1},
    {"vbf_m125_13TeV",2},
    {"wzh_m125_13TeV_zh",3},
    {"bbh_m125_13TeV",4}
  };


  _singleHiggsWSfileNames =
    {
      {"ggh_m125_powheg_13TeV","hgg.hig.mH125_13TeV.ggh"},
      {"tth_m125_13TeV","hgg.hig.mH125_13TeV.tth"},
      {"vbf_m125_13TeV","hgg.hig.mH125_13TeV.vbf"},
      {"wzh_m125_13TeV_zh","hgg.hig.mH125_13TeV.vh"},
      {"bbh_m125_13TeV","hgg.hig.mH125_13TeV.bbh"}
    };


  _nonResWeightIndex = doNRW;

  // Some defaults here are:
  // -2: do Resonant limits
  // -1: Non-resonant limits from Nodes
  // 0-1506 && 1507-1518: Non-resonant limits with re-weighting
  if (_nonResWeightIndex>=0)
    _wName = Form("evWeight_NRW_%d",doNRW);
  else
    _wName = "evWeight";

  if (doARW) {
    _wName = "new_evWeight";
    _nonResWeightIndex = -10;
  }

  //std::cout<<"DBG.  Finished Initialize..."<<std::endl;


  _c1 = new TCanvas("c1","Square Canvas",800,800);
  _c2 = new TCanvas("c2","Rectangular Canvas",800,600);


  _NR_MassRegion=0;
  if (folder_name.find("LowMass")!=std::string::npos)
    _NR_MassRegion=1;
  else if (folder_name.find("HighMass")!=std::string::npos)
    _NR_MassRegion=2;
  else
    _NR_MassRegion=99;


  if (logFileName!=""){
    // If the file for logging specified, redirect all std::out to it:
    newCout.open(logFileName, std::ofstream::out);
    std::cout.rdbuf(newCout.rdbuf());
  }

  std::cout<<"\t Initialized the fitter"<<std::endl;
  std::cout<<"SigMass: "<<SigMass
	   <<"\n NR_MassRegion: "<<_NR_MassRegion
	   <<"\n doNRW: "<< doNRW
	   <<"\n wName:"<<_wName
	   <<"\n "<<std::endl;
}

RooArgSet* bbgg2DFitter::defineVariables(bool swithToSimpleWeight=false)
{
  RooRealVar* mgg  = new RooRealVar("mgg","M(#gamma#gamma)",_minMggMassFit,_maxMggMassFit,"GeV");
  RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  RooRealVar* mjj  = new RooRealVar("mjj","M(jj)",_minMjjMassFit,_maxMjjMassFit,"GeV");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  RooRealVar* evWeight = 0;
  RooRealVar* new_evWeight = 0;

  TString tmp_wName(_wName.c_str());
  // This is to address a specific issue when adding single Higgs samples, while running with --NRW option.
  // In that case, for a signal sample we should take evWeight_NRW_%d, while for single higgs sample, we should use evWeight
  if (swithToSimpleWeight)
    tmp_wName="evWeight";
  // --- //

  if (!_doARW)
    evWeight = new RooRealVar(tmp_wName,"HqT x PUwei",-100000.,100000, "");
  else
  {
    evWeight = new RooRealVar("evWeight","HqT x PUwei",-100000, 100000,"");
    new_evWeight = new RooRealVar("new_evWeight","HqT x PUwei x ARW",-100000,100000,"");
  }


  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  cut_based_ct->defineType("cat4_2",2);
  cut_based_ct->defineType("cat4_3",3);
  //
  RooArgSet* ntplVars = 0;
  if (_doARW)
    ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight, *new_evWeight);
  else if (_nonResWeightIndex>=-1)
  {
    ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight);
    ntplVars->add(*mtot);
  }
  else
    ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight);

  // AP: Why these are here? They are already in the set:
  //ntplVars->add(*mgg);
  //ntplVars->add(*mjj);
  //ntplVars->add(*cut_based_ct);

  return ntplVars;
}

int bbgg2DFitter::AddSigData(float mass, TString signalfile)
{
  if (_verbLvl>1) std::cout << "================= Add Signal========================== " << _wName.c_str() << " " << _doARW << " " << _nonResWeightIndex << std::endl;
  if (_verbLvl>1) std::cout << " File to open:"<<signalfile  << std::endl;
  TFile *sigFile = TFile::Open(signalfile);
  bool opened=sigFile->IsOpen();
  if(opened==false) return -1;
  if (_verbLvl>1) std::cout << " TFile opened:"<<signalfile  << std::endl;

  TTree* sigTree = (TTree*)sigFile->Get("TCVARS");

  //Luminosity
  RooRealVar lumi("lumi","lumi", _lumi);
  _w->import(lumi);
  //Define variables
  RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
  if(sigTree==nullptr)
    {
      if (_verbLvl>1) std::cout<<"TCVARS for AddSigData  not founded in TTree trying with TCVARS"<<std::endl;
      std::exit(1);
      //sigTree = (TTree*) sigFile.Get("TCVARS");
      //if(sigTree==nullptr)
    }
  //Data set

  //Double_t W;
  //ccbar->SetBranchAddress("weight", &wCCBar);
  //ccbar->GetEntry();
  //RooRealVar ccbarweight("NRweight", "NRweight", );

  if (_verbLvl>0) {
    std::cout<<"[DBG]  Prining ntplVars from sig"<<std::endl;
    ntplVars->Print();
  }

  RooDataSet sigScaled("sigScaled","dataset",sigTree,*ntplVars,_cut, _wName.c_str());
//  if(_doARW) sigScaled = RooDataSet("sigScaled","dataset",sigTree,*ntplVars,_cut, "new_evWeight");
//  else sigScaled = RooDataSet("sigScaled","dataset",sigTree,*ntplVars,_cut, _wName.c_str());

  RooDataSet* sigToFit[_NCAT];
  TString cut0 = " && 1>0";

  RooArgList myArgList(*_w->var("mgg"));

  if (_fitStrategy != 1)
    myArgList.add(*_w->var("mjj"));

  if (_nonResWeightIndex>=-1)
    myArgList.add(*_w->var("mtot"));

  myArgList.Print();

  for ( int i=0; i<_NCAT; ++i)
    {

      std::cout << "-- Reducing cat " << i << std::endl;

      sigToFit[i] = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString::Format(" && cut_based_ct==%d ",i)+cut0);

      if (_fitStrategy == 1)
	sigToFit[i] = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString::Format(" && cut_based_ct==%d && mjj < 140 ",i)+cut0);

      this->SetSigExpectedCats(i, sigToFit[i]->sumEntries());

      if (_verbLvl>0) {
	std::cout << "======================================================================" <<std::endl;
	std::cout<<"[DBG]  Cat="<<i<< "\t Sig sumEntries="<<sigToFit[i]->sumEntries()<<std::endl;
	std::cout<<"mGG:  Mean = "<<sigToFit[i]->mean(*_w->var("mgg"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mgg"))<<std::endl;
	if (_fitStrategy != 1)
	  std::cout<<"mJJ:  Mean = "<<sigToFit[i]->mean(*_w->var("mjj"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mjj"))<<std::endl;

	if (_nonResWeightIndex>=-1)
	  std::cout<<"mTot: Mean = "<<sigToFit[i]->mean(*_w->var("mtot"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mtot"))<<std::endl;
      }

      /*This defines each category*/
      std::cout << "-- Importing cat " << i << std::endl;
      _w->import(*sigToFit[i],Rename(TString::Format("Sig_cat%d",i)));
    }
  // Create full signal data set without categorization
  std::cout << "-- Reducing all signal, no cat" << std::endl;
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(myArgList,_cut);
  if (_fitStrategy == 1)
    sigToFitAll = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString(" && mjj < 140 "));

  _w->import(*sigToFitAll,Rename("Sig"));

  // here we print the number of entries on the different categories
  if (_verbLvl>1) {
    std::cout << "======================================================================" <<std::endl;
    std::cout << "========= the number of entries on the different categories ==========" <<std::endl;
    std::cout << "---- one channel: " << sigScaled.sumEntries() <<std::endl;
    for (int c = 0; c < _NCAT; ++c)
      {
	Float_t nExpEvt = sigToFit[c]->sumEntries();
	std::cout<<TString::Format("nEvt exp. cat%d : ",c)<<nExpEvt<<TString::Format(" eff x Acc cat%d : ",c)<< "%"<<std::endl;
    } // close ncat
    std::cout << "======================================================================" <<std::endl;
    sigScaled.Print("v");
    std::cout << "----- DONE With Adding Signal!" << std::endl;
  }
  return 0;
}

std::vector<float> bbgg2DFitter::AddHigData(float mass, TString signalfile, int higgschannel, TString higName)
{
  if (_verbLvl>1) {
    std::cout << "================= Adding Single Higgs ==========================" <<std::endl;
    std::cout<<" \t mass: "<<mass<<" signalfile="<<signalfile<<" higgschannel="<<higgschannel<<" higName="<<higName<<std::endl;
  }

  RooArgSet* ntplVars = defineVariables(1);

  TFile higFile(signalfile);
  TTree* higTree = (TTree*) higFile.Get("TCVARS");
  if(higTree==nullptr)
    {
      if (_verbLvl>1) std::cout<<"TCVARS for AddHigData  not founded in TTree trying with TCVARS"<<std::endl;
      std::exit(1);
      //higTree = (TTree*) higFile.Get("TCVARS");
      //if(higTree==nullptr)std::exit(1);
    }
  RooDataSet higScaled("higScaled1","dataset",higTree, /* all variables of RooArgList*/*ntplVars,_cut,"evWeight");
  //
  RooDataSet* higToFit[_NCAT];
  TString cut0 = "&& 1>0";
  // we take only mtot to fit to the workspace, we include the cuts
  for ( int i=0; i<_NCAT; ++i)
    {
      higToFit[i] = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d ",i)+cut0);
      if(_fitStrategy == 1) higToFit[i] = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg")),_cut+TString::Format(" && cut_based_ct==%d && mjj < 140 ",i)+cut0);
      _w->import(*higToFit[i],Rename(TString::Format("Hig_%s_cat%d",higName.Data(),i)));
    }
  // Create full signal data set without categorization
  RooDataSet* higToFitAll = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut);
  if(_fitStrategy == 1) higToFitAll = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg")),_cut + TString(" && mjj < 140 "));
  _w->import(*higToFitAll,Rename("Hig"));
  // here we print the number of entries on the different categories
  if (_verbLvl>1) {
    std::cout << "========= the number of entries on the different categories (Higgs data) ==========" <<std::endl;
    std::cout << "---- one channel: " << higScaled.sumEntries() <<std::endl;
    for (int c = 0; c < _NCAT; ++c)
      {
	Float_t nExpEvt = higToFit[c]->sumEntries();
        std::cout<<TString::Format("nEvt exp. cat%d : ",c)<<nExpEvt<<TString::Format(" eff x Acc cat%d : ",c)<<"%"<<std::endl;
    }
    std::cout << "======================================================================" <<std::endl;
    higScaled.Print("v");
    std::cout << "===  DONE With Hig Data =="<<std::endl;
  }

  std::vector<float> thisExpHig;
  for (int c = 0; c < _NCAT; ++c)
  {
    Float_t nExpEvt = higToFit[c]->sumEntries();
    thisExpHig.push_back(nExpEvt);
  }
  return thisExpHig;

}

void bbgg2DFitter::AddBkgData(TString datafile)
{
  //Define variables
  RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
  //RooRealVar weightVar("weightVar","",1,0,1000);
  //weightVar.setVal(1.);
  TFile dataFile(datafile);
  TTree* dataTree = (TTree*) dataFile.Get("TCVARS");
  if(dataTree==nullptr)
    {
      if (_verbLvl>1) std::cout<<"TCVARS for AddBkgData  not founded in TTree trying with TCVARS"<<std::endl;
      std::exit(1);
      //dataTree = (TTree*) dataFile.Get("TCVARS");
      //if(dataTree==nullptr)std::exit(1);
    }
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","evWeight");
  RooDataSet* dataToFit[_NCAT];
  RooDataSet* dataToPlot[_NCAT];
  TString cut0 = "&& 1>0";
  TString cut1 = "&& 1>0";

  
  if (_verbLvl>1) std::cout<<"================= Add Bkg ==============================="<<std::endl;
  if (_verbLvl>1) {
    std::cout<<"\t Total events in root file: "<<Data.sumEntries()<<std::endl;
    std::cout<<"\t reduce to category 0: "<<Data.reduce("cut_based_ct==0")->sumEntries()
	     <<"  From the original Tree: "<<dataTree->GetEntries("cut_based_ct==0")<<std::endl;    
    std::cout<<"\t reduce to category 1: "<<Data.reduce("cut_based_ct==1")->sumEntries() 
	     <<"  From the original Tree: "<<dataTree->GetEntries("cut_based_ct==1")<<std::endl;    
    
  }
    
  for( int i=0; i<_NCAT; ++i)
    {

      dataToFit[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d",i));

      this->SetObservedCats(i, dataToFit[i]->sumEntries());

      if (_verbLvl>1) std::cout<<"\t categ="<<i<<"  events="<<dataToFit[i]->sumEntries()<<std::endl;

      if(_doblinding)
	dataToPlot[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d",i)+cut0);
      else
	dataToPlot[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d",i) );

      if (_fitStrategy == 1) {
	dataToFit[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg")),_cut+TString::Format(" && cut_based_ct==%d && mjj < 140",i));

	if(_doblinding)
	  dataToPlot[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d && mjj < 140 ",i)+cut0);
	else
	  dataToPlot[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d && mjj < 140 ",i) );

      }
      _w->import(*dataToFit[i],Rename(TString::Format("Data_cat%d",i)));
      _w->import(*dataToPlot[i],Rename(TString::Format("Dataplot_cat%d",i)));
    }
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut);
  if (_fitStrategy == 1)
    data = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg")),_cut + TString(" && mjj < 140 "));
  _w->import(*data, Rename("Data"));
  if (_verbLvl>1) data->Print("v");
}

void bbgg2DFitter::SigModelFit(float mass)
{
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[_NCAT];
  RooAbsPdf* mggSig[_NCAT];
  RooAbsPdf* mjjSig[_NCAT];
  RooProdPdf* SigPdf[_NCAT];
  RooAbsPdf* SigPdf1[_NCAT];
  // fit range
  //Float_t minSigFitMgg(115),maxSigFitMgg(135); //This should be an option
  //Float_t minSigFitMjj(60),maxSigFitMjj(180); //This should be an option
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  mgg->setRange("SigFitRange",_minSigFitMgg,_maxSigFitMgg);
  mjj->setRange("SigFitRange",_minSigFitMjj,_maxSigFitMjj);
  for (int c = 0; c < _NCAT; ++c)
    {
      // import sig and data from workspace

      sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
      mggSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggSig_cat%d",c));
      mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjSig_cat%d",c));

      if(_fitStrategy == 2) SigPdf[c] = new RooProdPdf(TString::Format("SigPdf_cat%d",c),"",RooArgSet(*mggSig[c], *mjjSig[c]));
      if(_fitStrategy == 1) SigPdf1[c] = (RooAbsPdf*) mggSig[c]->Clone(TString::Format("SigPdf_cat%d",c));

      ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mass);

      //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
      //peak->setVal(MASS);
      if (_verbLvl>1) std::cout << "OK up to now... Mass point: " <<mass<<std::endl;
      if(_fitStrategy == 2) SigPdf[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));
      if(_fitStrategy == 1) SigPdf1[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));
/*
      if (_verbLvl>1) std::cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() <<std::endl;
      double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+(mass-125.0); // shift the peak //This should be an option

      ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak
      if (_verbLvl>1) std::cout << "mPeak = " << mPeak << std::endl;
      if (_verbLvl>1) std::cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() <<std::endl;
*/
      // IMPORTANT: fix all pdf parameters to constant, why?
      RooArgSet* paramsMjj = (RooArgSet*) SigPdf[c]->getParameters(*mjj);
      TIterator* iterMjj = (TIterator*) paramsMjj->createIterator();
      TObject* tempObjMjj = nullptr;
      RooArgSet sigParams;
      while( (tempObjMjj = iterMjj->Next()) ) {
        if ( (TString(tempObjMjj->GetName()).EqualTo("mjj")) || (TString(tempObjMjj->GetName()).EqualTo("mgg"))) continue;
        std::cout << "Signal variables: " << tempObjMjj->GetName() << std::endl;
        sigParams.add(*_w->var(tempObjMjj->GetName()));
      }
      sigParams.Print("v");
/*
      RooArgSet* paramsMgg = (RooArgSet*) SigPdf[c]->getParameters(*mgg);
      TIterator* iterMjj = paramsMjj->createIterator();
      TObject* tempObjMjj=nullptr;

      if (_verbLvl>1) {
        while((tempObjMjj=iterMjj->Next()))
          {
            RooRealVar* var = (RooRealVar*)tempObjMjj;
            std::cout << "Variables after fit = " << tempObjMjj->GetName() << " " << var->getVal() << "+/-" << var->getError() << std::endl;
        }
        std::cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%s_cat%d",higName.Data(),c)))->getVal() <<std::endl;
      }
*/
/*
      RooArgSet sigParams( RooArgSet( *_w->var(TString::Format("mgg_sig_m0_cat%d",c)),
                           *_w->var(TString::Format("mgg_sig_sigma_cat%d",c)) ) );
      if(_fitStrategy==2) {
        sigParams.add( RooArgSet( *_w->var(TString::Format("mjj_sig_m0_cat%d",c)),
                                *_w->var(TString::Format("mjj_sig_sigma_cat%d",c)) ) );
      }
      if (!_useDSCB) {
         sigParams.add(RooArgSet(
			   *_w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_n_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_frac_cat%d",c))));
         if(_fitStrategy == 2) {
   	   sigParams.add(RooArgSet( *_w->var(TString::Format("mjj_sig_alpha_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_n_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_gsigma_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_frac_cat%d",c))));
         }
      }
      if(_useDSCB) {
         sigParams.add(RooArgSet(
			   *_w->var(TString::Format("mgg_sig_alpha1_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_n1_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_alpha2_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_n2_cat%d",c))));
         if(_fitStrategy == 2) {
   	   sigParams.add(RooArgSet(
				*_w->var(TString::Format("mjj_sig_alpha1_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_n1_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_alpha2_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_n2_cat%d",c))));
         }
      }
*/
      _w->defineSet(TString::Format("SigPdfParam_cat%d",c), sigParams);
      _w->set(TString::Format("SigPdfParam_cat%d",c))->Print("v");
      SetConstantParams(_w->set(TString::Format("SigPdfParam_cat%d",c)));
      if (_verbLvl>1) std::cout<<std::endl;
      if(_fitStrategy == 2) _w->import(*SigPdf[c]);
      if(_fitStrategy == 1) _w->import(*SigPdf1[c]);
    }
}

void bbgg2DFitter::HigModelFit(float mass, int higgschannel, TString higName)
{
  // four categories to fit
  RooDataSet* higToFit[_NCAT];
  RooAbsPdf* mggHig[_NCAT];
  RooAbsPdf* mjjHig[_NCAT];
  RooProdPdf* HigPdf[_NCAT];
  // fit range
  //Float_t minHigMggFit(115),maxHigMggFit(135);//This should be an option
  //Float_t minHigMjjFit(60),maxHigMjjFit(180);//This should be an option
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  mgg->setRange("HigFitRange",_minHigMggFit,_maxHigMggFit);
  mjj->setRange("HigFitRange",_minHigMjjFit,_maxHigMjjFit);
  for (int c = 0; c < _NCAT; ++c)
    {
      // import sig and data from workspace
      higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_%s_cat%d",higName.Data(),c));
      mggHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%s_cat%d",higName.Data(),c));
//      mjjHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",higgschannel,c));
//      if(higgschannel == 1 || higgschannel == 3) mjjHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",higgschannel,c));
      if(higName.Contains("ggh") == 1 || higName.Contains("vbf") == 1) {
        mjjHig[c] = new RooBernstein(TString::Format("mjjHig_%s_cat%d",higName.Data(),c),"",*mjj,
			RooArgList( *_w->var( TString::Format("mjj_hig_slope1_%s_cat%d", higName.Data(),c) ),
				    *_w->var( TString::Format("mjj_hig_slope2_%s_cat%d", higName.Data(),c) ),
				    *_w->var( TString::Format("mjj_hig_slope3_%s_cat%d", higName.Data(),c) ) ));
      } else {
        mjjHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%s_cat%d",higName.Data(),c));
      }
      HigPdf[c] = new RooProdPdf(TString::Format("HigPdf_%s_cat%d",higName.Data(),c),"",RooArgSet(*mggHig[c], *mjjHig[c]));
      std::cout << TString::Format("mggHig_%s_cat%d",higName.Data(),c) << std::endl;
      mggHig[c]->Print();
      std::cout << TString::Format("mjjHig_%s_cat%d",higName.Data(),c) << std::endl;
      mjjHig[c]->Print();
      std::cout << TString::Format("HigPdf_%s_cat%d",higName.Data(),c) << std::endl;
      HigPdf[c]->Print();
      std::cout << TString::Format("Hig_%s_cat%d",higName.Data(),c) << std::endl;
      higToFit[c]->Print();
      //((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(MASS);
      if (_verbLvl>1) std::cout << "OK up to now... Mass point: " <<mass<<std::endl;
      HigPdf[c]->fitTo(*higToFit[c],Range("HigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));
      RooArgSet* paramsMjj;
      paramsMjj = (RooArgSet*) mjjHig[c]->getParameters(*mjj);
      TIterator* iterMjj = paramsMjj->createIterator();
      TObject* tempObjMjj=nullptr;

      if (_verbLvl>1) {
	while((tempObjMjj=iterMjj->Next()))
	  {
	    RooRealVar* var = (RooRealVar*)tempObjMjj;
	    std::cout << "Variables after fit = " << tempObjMjj->GetName() << " " << var->getVal() << "+/-" << var->getError() << std::endl;
	}
	std::cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%s_cat%d",higName.Data(),c)))->getVal() <<std::endl;
      }
      //There are very few events in some fits, so adjust the max by a good amount so the MASS-125.0 shift doesn't touch it.
      /*
      ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setMax( ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getMax()+(mass-125.0) );
      double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal()+(mass-125.0); // shift the peak
      ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(mPeak); // shift the peak
      if (_verbLvl>1) std::cout << "mPeak = " << mPeak <<std::endl;
      if (_verbLvl>1) std::cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() <<std::endl;
      */
      // IMPORTANT: fix all pdf parameters to constant

      RooArgSet* paramsMgg = (RooArgSet*) HigPdf[c]->getParameters(*mgg);
      TIterator* iterMgg = (TIterator*) paramsMgg->createIterator();
      TObject* tempObjMgg = nullptr;
      RooArgSet sigParams;
      while( (tempObjMgg = iterMgg->Next()) ) {
        if ( (TString(tempObjMgg->GetName()).EqualTo("mjj")) || (TString(tempObjMgg->GetName()).EqualTo("mgg"))) continue;
        std::cout << "Higgs variables: " << tempObjMgg->GetName() << std::endl;
        sigParams.add(*_w->var(tempObjMgg->GetName()));
      }
      sigParams.Print("v");
/*
      RooArgSet sigParams;
      sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_m0_%s_cat%d", higName.Data(),c))));
      sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_sigma_%s_cat%d", higName.Data(),c))));
      if(higName.Contains("ggh") == 1 || higName.Contains("vbf") == 1) {
        sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_slope1_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_slope2_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_slope3_%s_cat%d", higName.Data(),c))));
      } else {
        sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_m0_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_sigma_%s_cat%d", higName.Data(),c))));
      }
      if(_useDSCB) {
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_alpha1_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_n1_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_alpha2_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_n2_%s_cat%d", higName.Data(),c))));
        if(higName.Contains("ggh") == 0 && higName.Contains("vbf") == 0) {
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_alpha1_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_n1_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_alpha2_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_n2_%s_cat%d", higName.Data(),c))));
        }
      } else {
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_alpha_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_n_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_gsigma_%s_cat%d", higName.Data(),c))));
        sigParams.add(RooArgSet(*_w->var(TString::Format("mgg_hig_frac_%s_cat%d", higName.Data(),c))));
        if(higName.Contains("ggh") == 0 && higName.Contains("vbf") == 0) {
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_alpha_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_n_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_gsigma_%s_cat%d", higName.Data(),c))));
          sigParams.add(RooArgSet(*_w->var(TString::Format("mjj_hig_frac_%s_cat%d", higName.Data(),c))));
        }
      }
*/
/*
			   *_w->var(TString::Format("mgg_hig_sigma_%d_cat%d",higgschannel,c)));
      if(!_useDSCB) {
	sigParams.add( RooArgSet(
			   *_w->var(TString::Format("mgg_hig_alpha_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_n_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_gsigma_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_frac_%d_cat%d",higgschannel,c))));
      } else {}
      if(higgschannel == 1 || higgschannel == 3){
	sigParams.add(RooArgSet(
				*_w->var(TString::Format("mjj_hig_m0_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_sigma_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_alpha_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_n_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_gsigma_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_frac_%d_cat%d",higgschannel,c)) ) );
      }
*/
      _w->defineSet(TString::Format("HigPdfParam_%s_cat%d",higName.Data(),c), sigParams);
      SetConstantParams(_w->set(TString::Format("HigPdfParam_%s_cat%d",higName.Data(),c)));
      if (_verbLvl>1) std::cout<<std::endl;
      _w->import(*HigPdf[c]);
    } // close for ncat
} // close higgs model fit

void bbgg2DFitter::MakePlots(float mass)
{

  _c1->cd();

  std::vector<TString> catdesc;
  int trueSigMass = _sigMass;

  // AP: What a hell is this?! :
  if(_sigMass >= 9000) trueSigMass = _sigMass - 9000;
  //

  if( _NCAT == 2 )catdesc={" High Purity Category"," Med. Purity Category"};
  if( _NCAT == 1 )catdesc={" High Mass Analysis"," High Mass Analysis"};
  //  else catdesc={" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}",
  //												 " #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}"};
  //
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) _w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  std::vector<RooDataSet*>sigToFit(_NCAT,nullptr);
  std::vector<RooAbsPdf*>mggGaussSig(_NCAT,nullptr);
  std::vector<RooAbsPdf*>mggCBSig(_NCAT,nullptr);
  std::vector<RooAbsPdf*>mggSig(_NCAT,nullptr);
  std::vector<RooAbsPdf*> mjjGaussSig(_NCAT,nullptr);
  std::vector<RooAbsPdf*> mjjCBSig(_NCAT,nullptr);
  std::vector<RooAbsPdf*> mjjSig(_NCAT,nullptr);
  //
  std::vector<RooAbsPdf*> mggBkg(_NCAT,nullptr);
  std::vector<RooAbsPdf*> mjjBkg(_NCAT,nullptr);

  std::vector<float> sigma_mjj;
  std::vector<float> sigma_mjj_std;
  std::vector<float> sigma_mgg;
  std::vector<float> sigma_mgg_std;

  std::vector<float> mean_mgg;
  std::vector<float> mean_mjj;


  for (int c = 0; c < _NCAT; ++c)
    {
      // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
      sigToFit[c]    = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
      mggSig[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mggSig_cat%d",c));
      mggBkg[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mggBkg_cat%d",c));
      mjjSig[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mjjSig_cat%d",c));
      mjjBkg[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mjjBkg_cat%d",c));

      if (!_useDSCB){
        mggGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggGaussSig_cat%d",c));
        mggCBSig[c]    = (RooAbsPdf*) _w->pdf(TString::Format("mggCBSig_cat%d",c));
        mjjGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjGaussSig_cat%d",c));
        mjjCBSig[c]    = (RooAbsPdf*) _w->pdf(TString::Format("mjjCBSig_cat%d",c));
      }

      std::vector<float> effSigmaVecMgg = EffectiveSigma( _w->var("mgg"), mggSig[c], _minSigFitMgg, _maxSigFitMgg);
      sigma_mgg.push_back(effSigmaVecMgg[3]);

      double mgg_sigmaSTD = (mggSig[c]->sigma(*_w->var("mgg")))->getVal();
      sigma_mgg_std.push_back(mgg_sigmaSTD);

      double mjj_sigmaSTD = (mggSig[c]->sigma(*_w->var("mjj")))->getVal();
      sigma_mjj_std.push_back(mjj_sigmaSTD);

      std::vector<float> effSigmaVecMjj = EffectiveSigma( _w->var("mjj"), mjjSig[c], _minSigFitMjj, _maxSigFitMjj);
      sigma_mjj.push_back(effSigmaVecMjj[3]);

//      double mgg_mean = (mggSig[c]->mean(*_w->var("mgg")))->getVal();
      double mgg_mean = ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal();//(mggSig[c]->getParameters(*_w->var("mgg")))->find(
      mean_mgg.push_back(mgg_mean);

//      double mjj_mean = (mjjSig[c]->mean(*_w->var("mjj")))->getVal();
      double mjj_mean = ((RooRealVar*) _w->var(TString::Format("mjj_sig_m0_cat%d",c)))->getVal();
      mean_mjj.push_back(mjj_mean);

    } // close categories


  RooRealVar* mgg = _w->var("mgg");
  mgg->setUnit("GeV");
  //RooAbsPdf* mggGaussSigAll = _w->pdf("mggGaussSig");
  //RooAbsPdf* mggCBSigAll = _w->pdf("mggCBSig");
  //RooAbsPdf* mggSigAll = _w->pdf("mggSig");
  RooRealVar* mjj = _w->var("mjj");
  mjj->setUnit("GeV");
  //RooAbsPdf* mjjGaussSigAll = _w->pdf("mjjGaussSig");
  //RooAbsPdf* mjjCBSigAll = _w->pdf("mjjCBSig");
  //RooAbsPdf* mjjSigAll = _w->pdf("mjjSig");
  //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  // SetParamNames(w);
  //  _minSigFitMjj, _maxSigFitMjj
  Float_t minSigPlotMgg(115),maxSigPlotMgg(135);
  Float_t minSigPlotMjj(_minSigFitMjj),maxSigPlotMjj(_maxSigFitMjj);
  mgg->setRange("SigPlotRange",minSigPlotMgg,maxSigPlotMgg);
  mjj->setRange("SigPlotRange",minSigPlotMjj,maxSigPlotMjj);
  Int_t nBinsMass(20); // just need to plot
  RooPlot* plotmggAll = mgg->frame(Range("SigPlotRange"),Bins(nBinsMass));
  signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  //  _c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  //  TLatex *text = new TLatex();
  //  text->SetNDC();
  //  text->SetTextSize(0.04);
  RooPlot* plotmgg[_NCAT];


  RooRealVar* mtot = _w->var("mtot");
  mtot->setUnit("GeV");
  RooPlot* plotmtot[_NCAT];
  Float_t minSigPlotMtot(300),maxSigPlotMtot(1600);
  mtot->setRange("SigPlotRange",minSigPlotMtot,maxSigPlotMtot);

  //TCanvas* ctmp = new TCanvas("c1","Canvas",800,800);

  if (_verbLvl>1) std::cout << "[MakePlots] Doing now sig Mgg  and Mtot plots" << std::endl;
  for (int c = 0; c < _NCAT; ++c)
    {
      if (_nonResWeightIndex>=-1){
	plotmtot[c] = mtot->frame(Range("SigPlotRange"),Bins(30));
	sigToFit[c]->plotOn(plotmtot[c]);
	plotmtot[c]->Draw();
      	_c1->SaveAs(TString::Format("%s/sigMtot_cat%d.png",_folder_name.data(),c),"QUIET");
      }

      plotmgg[c] = mgg->frame(Range("SigPlotRange"),Bins(nBinsMass));
      sigToFit[c]->plotOn(plotmgg[c]);
      mggSig[c] ->plotOn(plotmgg[c]);
      //    double chi2n = plotmgg[c]->chiSquare(0) ;
      //    if (_verbLvl>1) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
      if (!_useDSCB) {
         mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggGaussSig_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
         mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggCBSig_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
      }
      //    mggSig[c] ->paramOn(plotmgg[c]);
      sigToFit[c] ->plotOn(plotmgg[c]);
      // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
      //    TH1F *hist = new TH1F(TString::Format("histMgg_cat%d",c), "hist", 400, minSigPlotMgg, maxSigPlotMgg);
      //plotmgg[c]->SetTitle("CMS preliminary 19.7/fb ");
      plotmgg[c]->SetMinimum(0.0);
      plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
      plotmgg[c]->GetXaxis()->SetTitle("M(#gamma#gamma) [GeV]");
      plotmgg[c]->Draw();

      //    plotmgg[c]->Draw("SAME");
      TLegend *legmc = new TLegend(0.52,0.7,0.92,0.90);
      legmc->AddEntry(plotmgg[c]->getObject(4),"Simulation","LPE");
      legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
      if (!_useDSCB){
         legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian component","L");
         legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
      }
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      //    TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
      //    pt->SetNDC();
      //pt->SetName("title");
      //    pt->SetBorderSize(0);
      //    pt->SetFillColor(0);
      //    pt->SetTextSize(0.25);
      //    pt->AddText("CMS Preliminary Simulation ");
      //    pt->Draw();
      TLatex *tlat0 = new TLatex();
      tlat0->SetNDC();
      tlat0->SetTextAngle(0);
      tlat0->SetTextColor(kBlack);
      tlat0->SetTextFont(63);
      tlat0->SetTextAlign(11);
      tlat0->SetTextSize(27);
      tlat0->DrawLatex(0.50, 0.92, "CMS");
      tlat0->SetTextFont(53);
      tlat0->DrawLatex(0.58, 0.92, "Simulation Preliminary");
      tlat0->SetTextFont(63);
      tlat0->SetTextSize(27);


      TString str_desc;

      if (_nonResWeightIndex==-2){
	// This is resonant case

	// What is this? Don't know:
	if(_sigMass >=9000) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + "");

	else if (_sigMass < 250 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " (Low Mass)");
	else if(_sigMass > 200 && _sigMass < 8000 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " ");
	tlat0->SetTextFont(43);

	str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), trueSigMass);
      }
      else {
	// This is Non-resonant case


	if (_NR_MassRegion==1)
	  tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " (Low Mass)");
	if (_NR_MassRegion==2)
	  tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " (High Mass)");


	if (_nonResWeightIndex==-1){
	  // This is the Nodes
	  if (_sigMass==0) str_desc=TString::Format(" Nonresonant HH, SM");
	  else if (_sigMass==1)	str_desc=TString::Format(" Nonresonant HH, Box Diagram Only");
	  else str_desc=TString::Format(" Nonresonant HH, Node %d", _sigMass);
	}
	else if (_nonResWeightIndex<1519){
	  // This is the points
	  str_desc=TString::Format(" Nonresonant HH, Weight No.%d", _nonResWeightIndex);
	}
	else
	  std::cout<<"Warning this index is BAD: "<<_nonResWeightIndex<<std::endl;
      }


      tlat0->DrawLatex(0.16, 0.82, str_desc);

      str_desc=TString::Format(" #mu = %.2f GeV",mean_mgg[c]);
      tlat0->DrawLatex(0.165, 0.77, str_desc);
      str_desc=TString::Format(" #sigma_{eff} = %.2f GeV",sigma_mgg[c]);
      tlat0->DrawLatex(0.165, 0.72, str_desc);

      _c1->SaveAs(TString::Format("%s/sigmodelMgg_cat%d.pdf",_folder_name.data(),c),"QUIET");
      _c1->SaveAs(TString::Format("%s/sigmodelMgg_cat%d.png",_folder_name.data(),c),"QUIET");
      //_c1->SaveAs(TString::Format("sigmodelMgg_cat%d.C",c));
    } // close categories

  //if (ctmp) ctmp->Close();
  //  _c1 = new TCanvas("cMjj","mgg",800,600);
  //  _c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  if(_fitStrategy==2) {
    TLatex * text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    RooPlot* plotmjj[_NCAT];
    if (_verbLvl>1) std::cout << "[MakePlots] Doing now sig Mjj plot" << std::endl;
    for (int c = 0; c < _NCAT; ++c)
      {
	plotmjj[c] = mjj->frame(Range("SigPlotRange"),Bins(nBinsMass));
	sigToFit[c]->plotOn(plotmjj[c]);
	mjjSig[c] ->plotOn(plotmjj[c]);
	double chi2n = plotmjj[c]->chiSquare(0) ;
	if (_verbLvl>1) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
        if (!_useDSCB) {
   	   mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjGaussSig_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
	   mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjCBSig_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
        }
	//    mjjSig[c] ->paramOn(plotmjj[c]);
	sigToFit[c] ->plotOn(plotmjj[c]);
	plotmjj[c]->SetMinimum(0.0);
	plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
	plotmjj[c]->GetXaxis()->SetTitle("M(jj) [GeV]");
	//    TCanvas* ctmp = new TCanvas(TString::Format("ctmpSigMjj_cat%d",c),"Background Categories",0,0,500,500);
	//TCanvas* ctmp = new TCanvas(TString::Format("ctmpSigMjj_cat%d",c),"Background Categories",800,800);
	plotmjj[c]->Draw();
	plotmjj[c]->Draw("SAME");
	//    TLegend *legmc = new TLegend(0.62,0.75,0.85,0.85);
	//      TLegend *legmc = new TLegend(0.55,0.7,0.95,0.89);
	TLegend *legmc = new TLegend(0.52,0.7,0.92,0.90);
	legmc->AddEntry(plotmgg[c]->getObject(5),"Simulation","LPE");
	legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
        if (!_useDSCB) {
	   legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian component","L");
	   legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
        }
	//    legmc->SetHeader(" ");
	legmc->SetBorderSize(0);
	legmc->SetFillStyle(0);
	legmc->Draw();

	TLatex *tlat0 = new TLatex();
	tlat0->SetNDC();
	tlat0->SetTextAngle(0);
	tlat0->SetTextColor(kBlack);
	tlat0->SetTextFont(63);
	tlat0->SetTextAlign(11);
	tlat0->SetTextSize(27);
	tlat0->DrawLatex(0.50, 0.92, "CMS");
	tlat0->SetTextFont(53);
	tlat0->DrawLatex(0.58, 0.92, "Simulation Preliminary");
	tlat0->SetTextFont(63);
	tlat0->SetTextSize(27);

	if(_sigMass >=9000) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + "");
	if(_sigMass < 250 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " (Low Mass)");
	if(_sigMass > 200 && _sigMass < 8000 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " ");

	tlat0->SetTextFont(43);

	// Ahhhh common... Repeating the same shit here... Who wrote this code?
	// Not gonna change it here. Better re-write from scracth
	TString str_desc;
	if(trueSigMass < 250 )str_desc=TString::Format(" Nonresonant HH, Node %d", trueSigMass);
	else str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), trueSigMass);
	if(trueSigMass == 0 )str_desc=TString::Format(" Nonresonant HH, Box Diagram Only");
	if(trueSigMass == 1 )str_desc=TString::Format(" Nonresonant HH, SM");
	tlat0->DrawLatex(0.16, 0.82, str_desc);

	str_desc=TString::Format(" #mu = %.2f GeV",mean_mjj[c]);
	tlat0->DrawLatex(0.165, 0.77, str_desc);
	str_desc=TString::Format(" #sigma_{eff} = %.2f GeV",sigma_mjj[c]);
	tlat0->DrawLatex(0.165, 0.72, str_desc);

	_c1->SaveAs(TString::Format("%s/sigmodelMjj_cat%d.pdf",_folder_name.data(),c),"QUIET");
	_c1->SaveAs(TString::Format("%s/sigmodelMjj_cat%d.png",_folder_name.data(),c),"QUIET");
	//1->SaveAs(TString::Format("sigmodelMjj_cat%d.C",c));

      } // close categories
  }

  //if (ctmp) ctmp->Close();

} // close makeplots signal

void bbgg2DFitter::MakePlotsHiggs(float mass)
//void bbgg2DFitter::MakePlotsHiggs(float mass,std::vector<std::string>higgstrue,std::map<std::string,int>_singleHiggsMap)
{
  std::vector<TString> catdesc;
  if( _NCAT == 2 )catdesc={" High Purity Category"," Med. Purity Category"};
  if( _NCAT == 1 )catdesc={" High Mass Analysis"," High Mass Analysis"};
  //  else catdesc={" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}",
  //  if ( _NCAT == 2 )catdesc={" High Purity"," Med. Purity"};
  //  else catdesc={" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}",
  //								" #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}"};
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  //RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  std::vector<TString>component={"ggH","ttH","VBF","VH", "bbH"};
  //    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));
  for (unsigned int d = 0; d < _singleHiggsNames.size(); ++d)
    {
      int realint= _singleHiggsMap[_singleHiggsNames[d]];
      std::vector<RooDataSet*> higToFit(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mggGaussSig(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mggCBSig(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mggSig(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mjjGaussSig(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mjjCBSig(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mjjSig(_NCAT,nullptr);
      //
      std::vector<RooAbsPdf*> mggBkg(_NCAT,nullptr);
      std::vector<RooAbsPdf*> mjjBkg(_NCAT,nullptr);

      for (int c = 0; c < _NCAT; ++c)
	{
	  // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
	  higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",realint,c));
	  mggGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggGaussHig_%d_cat%d",realint,c));
	  mggCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggCBHig_%d_cat%d",realint,c));
	  mggSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%d_cat%d",realint,c));
	  mggBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggBkg_%d_cat%d",realint,c));
	  if(d == 1 || d == 3)
	    {
	      mjjGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjGaussHig_%d_cat%d",realint,c));
	      mjjCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjCBHig_%d_cat%d",realint,c));
	      mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",realint,c));
	    }
	  else mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",realint,c));
	  mjjBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjBkg_%d_cat%d",realint,c));
	} // close categories
      RooRealVar* mgg = _w->var("mgg");
      mgg->setUnit("GeV");
      RooRealVar* mjj = _w->var("mjj");
      mjj->setUnit("GeV");
      //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
      //
      //****************************//
      // Plot mgg Fit results
      //****************************//
      // Set P.D.F. parameter names
      // WARNING: Do not use it if Workspaces are created
      // SetParamNames(w);
      Float_t minHigPlotMgg(115);//,maxHigPlotMgg(135);
      Float_t minHigPlotMjj(80);//,maxHigPlotMjj(180);
      Int_t nBinsMass(20); // just need to plot
      //RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
      //higgsAll->plotOn(plotmggAll);
      gStyle->SetOptTitle(0);
      //TCanvas* _c1 = new TCanvas(TString::Format("cMgg_%d",realint),"mgg",0,0,500,500);
      _c1->cd(1);
      //********************************************//
      // Plot Signal Categories
      //****************************//
      TLatex *text = new TLatex();
      text->SetNDC();
      text->SetTextSize(0.04);
      RooPlot* plotmgg[_NCAT];
      for (int c = 0; c < _NCAT; ++c)
	{
	  plotmgg[c] = mgg->frame(Range("HigFitRange"),Bins(nBinsMass));
	  higToFit[c]->plotOn(plotmgg[c]);
	  mggSig[c] ->plotOn(plotmgg[c]);
	  double chi2n = plotmgg[c]->chiSquare(0) ;
	  if (_verbLvl>1) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
	  mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggGaussHig_%d_cat%d",realint,c)),LineStyle(kDashed),LineColor(kGreen));
	  mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggCBHig_%d_cat%d",realint,c)),LineStyle(kDashed),LineColor(kRed));
	  mggSig[c] ->paramOn(plotmgg[c]);
	  higToFit[c] ->plotOn(plotmgg[c]);
	  // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
	  //      TH1F *hist = new TH1F(TString::Format("histMgg_%d_cat%d",d,c), "hist", 400, minHigPlotMgg, maxHigPlotMgg);
	  //plotmgg[c]->SetTitle("CMS Preliminary 19.7/fb ");
	  plotmgg[c]->SetMinimum(0.0);
	  plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
	  plotmgg[c]->GetXaxis()->SetTitle("M(#gamma#gamma) [GeV]");
	  //TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMgg_%d_cat%d",d,c),"Background Categories",0,0,500,500);
	  plotmgg[c]->Draw();
	  plotmgg[c]->Draw("SAME");
	  TLegend *legmc = new TLegend(0.55,0.7,0.9,0.89);
	  legmc->AddEntry(plotmgg[c]->getObject(5),component[d],"LPE");
	  legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
	  legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
	  legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
	  //      legmc->SetHeader(" ");
	  legmc->SetBorderSize(0);
	  legmc->SetFillStyle(0);
	  legmc->Draw();
	  TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
	  //pt->SetName("title");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetTextSize(0.035);
	  pt->AddText("CMS Preliminary Simulation ");
	  //      pt->Draw();
	  TLatex *tlat0 = new TLatex();
	  tlat0->SetNDC();
	  tlat0->SetTextAngle(0);
	  tlat0->SetTextColor(kBlack);
	  tlat0->SetTextFont(63);
	  tlat0->SetTextAlign(11);
	  tlat0->SetTextSize(27);
	  //      tlat0->DrawLatex(0.5, 0.92, "CMS Preliminary Simulation");
	  tlat0->DrawLatex(0.54, 0.92, "CMS");
	  tlat0->SetTextFont(53);
	  tlat0->DrawLatex(0.61, 0.92, "Simulation Preliminary");
	  // float effS = effSigma(hist);
	  TString str_desc;
	  if(_sigMass==0) str_desc=" Nonresonant HH";
	  else str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), _sigMass);
	  TLatex *lat = new TLatex(minHigPlotMgg+0.5,0.85*plotmgg[c]->GetMaximum(),str_desc);
	  lat->Draw();
	  TLatex *lat2 = new TLatex(minHigPlotMgg+0.5,0.70*plotmgg[c]->GetMaximum(),catdesc.at(c));
	  lat2->Draw();
	  ///////
	  TString myChi2buffer=TString::Format("#chi^{2}/ndof = %f",chi2n);
	  TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
	  latex -> SetNDC();
	  latex -> SetTextFont(42);
	  latex -> SetTextSize(0.04);
	  //latex -> Draw("same");
	  _c1->SaveAs(TString::Format("%s/higmodelMgg_%d_cat%d.pdf",_folder_name.data(),realint,c),"QUIET");
	  _c1->SaveAs(TString::Format("%s/higmodelMgg_%d_cat%d.png",_folder_name.data(),realint,c),"QUIET");
	  //1->SaveAs(TString::Format("sigmodelMgg_cat%d.C",c));

	  //if (ctmp) ctmp->Close();

	} // close categories
      //_c1 = new TCanvas(TString::Format("cMjj_%d",realint),"mjj",0,0,500,500);
      //_c1->cd(1);
      //********************************************//
      // Plot Signal Categories
      //****************************//
      text = new TLatex();
      text->SetNDC();
      text->SetTextSize(0.04);
      RooPlot* plotmjj[_NCAT];
      for (int c = 0; c < _NCAT; ++c)
	{
	  plotmjj[c] = mjj->frame(Range("HigFitRange"),Bins(nBinsMass));
	  higToFit[c]->plotOn(plotmjj[c]);
	  mjjSig[c] ->plotOn(plotmjj[c]);
	  double chi2n = plotmjj[c]->chiSquare(0) ;
	  if (_verbLvl>1) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
	  if(d == 1 || d == 3)
	    {
	      mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjGaussHig_%d_cat%d",realint,c)),LineStyle(kDashed),LineColor(kGreen));
	      mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjCBHig_%d_cat%d",realint,c)),LineStyle(kDashed),LineColor(kRed));
	      mjjSig[c] ->paramOn(plotmjj[c]);
	      higToFit[c] ->plotOn(plotmjj[c]);
	      // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
	      //         TH1F *hist = new TH1F(TString::Format("histMjj_%d_cat%d",d,c), "hist", 400, minHigPlotMjj, maxHigPlotMjj);
	      //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
	      plotmjj[c]->SetMinimum(0.0);
	      plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
	      plotmjj[c]->GetXaxis()->SetTitle("M(jj) [GeV]");
	      //TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMjj_%d_cat_%d",realint,c),"Background Categories",0,0,500,500);
	      plotmjj[c]->Draw();
	      plotmjj[c]->Draw("SAME");
	      TLegend *legmc = new TLegend(0.55,0.7,0.95,0.89);
	      legmc->AddEntry(plotmjj[c]->getObject(5),component[realint],"LPE");
	      legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
	      legmc->AddEntry(plotmjj[c]->getObject(2),"Gaussian Outliers","L");
	      legmc->AddEntry(plotmjj[c]->getObject(3),"Crystal Ball component","L");
	      //         legmc->SetHeader(" ");
	      legmc->SetBorderSize(0);
	      legmc->SetFillStyle(0);
	      legmc->Draw();
	      TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
	      //pt->SetName("title");
	      pt->SetBorderSize(0);
	      pt->SetFillColor(0);
	      pt->SetTextSize(0.035);
	      pt->AddText("CMS Preliminary Simulation ");
	      //   pt->Draw();
	      TLatex *tlat0 = new TLatex();
	      tlat0->SetNDC();
	      tlat0->SetTextAngle(0);
	      tlat0->SetTextColor(kBlack);
	      tlat0->SetTextFont(63);
	      tlat0->SetTextAlign(11);
	      tlat0->SetTextSize(27);
	      //      tlat0->DrawLatex(0.5, 0.92, "CMS Preliminary Simulation");
	      tlat0->DrawLatex(0.54, 0.92, "CMS");
	      tlat0->SetTextFont(53);
	      tlat0->DrawLatex(0.61, 0.92, "Simulation Preliminary");
	      // float effS = effSigma(hist);
	      TString str_desc;
	      if(_sigMass==0)str_desc=" Nonresonant HH";
	      else str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), _sigMass);
	      TLatex *lat = new TLatex(minHigPlotMjj+0.5,0.85*plotmjj[c]->GetMaximum(),str_desc);
	      lat->Draw();
	      TLatex *lat2 = new TLatex(minHigPlotMjj+0.5,0.70*plotmjj[c]->GetMaximum(),catdesc.at(c));
	      lat2->Draw();
	      ///////
	      TString myChi2buffer=TString::Format("#chi^{2}/ndof = %f",chi2n);
	      TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      //latex -> Draw("same");
	      _c1->SaveAs(TString::Format("%s/higmodelMjj_%d_cat%d.pdf",_folder_name.data(),realint,c),"QUIET");
	      _c1->SaveAs(TString::Format("%s/higmodelMjj_%d_cat%d.png",_folder_name.data(),realint,c),"QUIET");
	    }
	  else
	    {
	      //mjjSig[c] ->paramOn(plotmjj[c]);
	      higToFit[c] ->plotOn(plotmjj[c]);
	      // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
	      //         TH1F *hist = new TH1F(TString::Format("histMjj_%d_cat%d",d,c), "hist", 400, minHigPlotMjj, maxHigPlotMjj);
	      //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
	      plotmjj[c]->SetMinimum(0.0);
	      plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
	      plotmjj[c]->GetXaxis()->SetTitle("M(jj) [GeV]");
	      //TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMjj_%d_cat_%d",realint,c),"Background Categories",0,0,500,500);
	      plotmjj[c]->Draw();
	      plotmjj[c]->Draw("SAME");
	      TLegend *legmc = new TLegend(0.55,0.7,0.95,0.89);
	      legmc->AddEntry(plotmjj[c]->getObject(5),component[realint],"LPE");
	      legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
	      //         legmc->SetHeader(" ");
	      legmc->SetBorderSize(0);
	      legmc->SetFillStyle(0);
	      legmc->Draw();
	      TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
	      //pt->SetName("title");
	      pt->SetBorderSize(0);
	      pt->SetFillColor(0);
	      pt->SetTextSize(0.035);
	      pt->AddText("CMS Preliminary Simulation ");
	      //         pt->Draw();
	      TLatex *tlat0 = new TLatex();
	      tlat0->SetNDC();
	      tlat0->SetTextAngle(0);
	      tlat0->SetTextColor(kBlack);
	      tlat0->SetTextFont(63);
	      tlat0->SetTextAlign(11);
	      tlat0->SetTextSize(27);
	      //      tlat0->DrawLatex(0.5, 0.92, "CMS Preliminary Simulation");
	      tlat0->DrawLatex(0.54, 0.92, "CMS");
	      tlat0->SetTextFont(53);
	      tlat0->DrawLatex(0.61, 0.92, "Simulation Preliminary");
	      // float effS = effSigma(hist);
	      TString str_desc;
	      if(_sigMass==0)str_desc=" Nonresonant HH";
	      else str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), _sigMass);
	      TLatex *lat = new TLatex(minHigPlotMjj+0.5,0.85*plotmjj[c]->GetMaximum(),str_desc);
	      lat->Draw();
	      TLatex *lat2 = new TLatex(minHigPlotMjj+0.5,0.70*plotmjj[c]->GetMaximum(),catdesc.at(c));
	      lat2->Draw();
	      ///////
	      TString myChi2buffer=TString::Format("#chi^{2}/ndof = %f",chi2n);
	      TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      //latex -> Draw("same");
	      _c1->SaveAs(TString::Format("%s/higmodelMjj_%d_cat%d.pdf",_folder_name.data(),realint,c),"QUIET");
	      _c1->SaveAs(TString::Format("%s/higmodelMjj_%d_cat%d.png",_folder_name.data(),realint,c),"QUIET");
	    }
	  //_c1->SaveAs(TString::Format("sigmodelMjj_cat%d.C",c));
	} // close categories
    } // close to higgs component
} // close makeplots signal

void bbgg2DFitter::MakeSigWS(std::string fileBaseName)
{
  TString wsDir = TString::Format("%s/workspaces/",_folder_name.data());
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save
  // for statistical tests.
  //**********************************************************************//
  std::vector<RooAbsPdf*> SigPdf(_NCAT,nullptr);
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
/*  _w->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  _w->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
  _w->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  _w->factory("CMS_hbb_sig_sigmaScale[1,1,1]");*/
  for (int c = 0; c < _NCAT; ++c)
    {
      int newC = c + _ncat0;

      _w->factory(TString::Format("CMS_hgg_sig_m0_absShift[1,1,1]"));
      _w->factory(TString::Format("CMS_hbb_sig_m0_absShift[1,1,1]"));
      _w->factory(TString::Format("CMS_hgg_sig_sigmaScale[1,1,1]"));
      _w->factory(TString::Format("CMS_hbb_sig_sigmaScale[1,1,1]"));

      //IMPORTHERE
      //
      SigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("SigPdf_cat%d",c));
      RooArgSet* paramsMjj = (RooArgSet*) SigPdf[c]->getParameters(*_w->var("mjj"));
      TIterator* iterMjj = (TIterator*) paramsMjj->createIterator();
      TObject* tempObjMjj = nullptr;
      std::vector<std::pair<TString,TString>> varsToChange;
      while( (tempObjMjj = iterMjj->Next()) ) {
        if ( (TString(tempObjMjj->GetName()).EqualTo("mjj")) || (TString(tempObjMjj->GetName()).EqualTo("mgg"))) continue;
        TString thisVarName(tempObjMjj->GetName());
        TString newVarName = TString(thisVarName).ReplaceAll(TString::Format("cat%d", c), TString::Format("cat%d", newC));
        if ( !newVarName.Contains("m0") && !newVarName.Contains("sigma") ) {
          if ( newVarName.Contains("mgg") ) newVarName.ReplaceAll("mgg_", "CMS_hgg_");
          if ( newVarName.Contains("mjj") ) newVarName.ReplaceAll("mjj_", "CMS_hbb_");
          varsToChange.push_back(std::make_pair(thisVarName, newVarName));
        }
        std::cout << "Importing variable with new name: old - " << thisVarName << " new - " << newVarName << std::endl;
        _w->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));
        wAll->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));
      }
      //Shifts and smearings
/*      _w->factory(TString::Format("prod::CMS_hgg_sig_m0_cat%d(mgg_sig_m0_cat%d, CMS_hgg_sig_m0_absShift)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hbb_sig_m0_cat%d(mjj_sig_m0_cat%d, CMS_hbb_sig_m0_absShift)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hgg_sig_sigma_cat%d(mgg_sig_sigma_cat%d, CMS_hgg_sig_sigmaScale)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hbb_sig_sigma_cat%d(mjj_sig_sigma_cat%d, CMS_hbb_sig_sigmaScale)", newC, newC));*/
      _w->factory(TString::Format("prod::CMS_hgg_sig_m0_cat%d(mgg_sig_m0_cat%d, CMS_hgg_sig_m0_absShift)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hbb_sig_m0_cat%d(mjj_sig_m0_cat%d, CMS_hbb_sig_m0_absShift)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hgg_sig_sigma_cat%d(mgg_sig_sigma_cat%d, CMS_hgg_sig_sigmaScale)", newC, newC));
      _w->factory(TString::Format("prod::CMS_hbb_sig_sigma_cat%d(mjj_sig_sigma_cat%d, CMS_hbb_sig_sigmaScale)", newC, newC));
      if(!_useDSCB) {
/*        _w->factory(TString::Format("prod::CMS_hgg_gsigma_cat%d(mgg_sig_gsigma_cat%d, CMS_hgg_sig_sigmaScale)", newC, newC));
        _w->factory(TString::Format("prod::CMS_hbb_gsigma_cat%d(mjj_sig_gsigma_cat%d, CMS_hbb_sig_sigmaScale)", newC, newC));*/
        _w->factory(TString::Format("prod::CMS_hgg_gsigma_cat%d(mgg_sig_gsigma_cat%d, CMS_hgg_sig_sigmaScale)", newC, newC));
        _w->factory(TString::Format("prod::CMS_hbb_gsigma_cat%d(mjj_sig_gsigma_cat%d, CMS_hbb_sig_sigmaScale)", newC, newC));
      }

      TString EditPDF = TString::Format("EDIT::CMS_sig_cat%d(SigPdf_cat%d,", newC, c);
      for (unsigned int iv = 0; iv < varsToChange.size(); iv++)
        EditPDF += TString::Format("%s=%s,", varsToChange[iv].first.Data(), varsToChange[iv].second.Data());
      //Shifted and smeared vars
      if(!_useDSCB) {
        EditPDF += TString::Format("mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d,", c, newC);
        EditPDF += TString::Format("mjj_sig_gsigma_cat%d=CMS_hbb_sig_gsigma_cat%d)", c, newC);
      }
      EditPDF += TString::Format("mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d,", c, newC);
      EditPDF += TString::Format("mjj_sig_m0_cat%d=CMS_hbb_sig_m0_cat%d,", c, newC);
      EditPDF += TString::Format("mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d,", c, newC);
      EditPDF += TString::Format("mjj_sig_sigma_cat%d=CMS_hbb_sig_sigma_cat%d)", c, newC);
      std::cout << "STRINGTOCHANGE   ---  " << EditPDF << std::endl;
      _w->factory(EditPDF);

      wAll->import(*_w->pdf(TString::Format("CMS_sig_cat%d",newC)));// Rename(TString::Format("SigPdf_cat%d", newC)));
      wAll->import(*_w->data(TString::Format("Sig_cat%d",c)), Rename(TString::Format("Sig_cat%d", newC)));
    }
  wAll->Print("v");
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  if (_verbLvl>1) std::cout << "Write signal workspace in: " << filename << " file" << std::endl;
  return;
} // close make signal WP

void bbgg2DFitter::MakeHigWS(std::string fileHiggsName,int higgschannel, TString higName)
{
  TString wsDir = TString::Format("%s/workspaces/",_folder_name.data());
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  std::vector<RooAbsPdf*> HigPdf(_NCAT,nullptr);
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");

  for (int c = 0; c < _NCAT; ++c)
    {
      int newC = c + _ncat0;
      HigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("HigPdf_%s_cat%d",higName.Data(),c));
      RooArgSet* paramsMjj = (RooArgSet*) HigPdf[c]->getParameters(*_w->var("mjj"));
      TIterator* iterMjj = (TIterator*) paramsMjj->createIterator();
      TObject* tempObjMjj = nullptr;
      std::vector<std::pair<TString,TString>> varsToChange;

      while( (tempObjMjj = iterMjj->Next()) ) {

        if ( (TString(tempObjMjj->GetName()).EqualTo("mjj")) || (TString(tempObjMjj->GetName()).EqualTo("mgg"))) continue;
        TString thisVarName(tempObjMjj->GetName());
        TString newVarName = TString(thisVarName).ReplaceAll(TString::Format("cat%d", c), TString::Format("cat%d", newC));

        if ( !newVarName.Contains("m0") && !newVarName.Contains("sigma") ) {
          if ( newVarName.Contains("mgg") ) newVarName.ReplaceAll("mgg_", "CMS_hgg_");
          if ( newVarName.Contains("mjj") ) newVarName.ReplaceAll("mjj_", "CMS_hbb_");
          varsToChange.push_back(std::make_pair(thisVarName, newVarName));
        }

        std::cout << "Importing variable with new name: old - " << thisVarName << " new - " << newVarName << std::endl;
        _w->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));
        wAll->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));

      }

      //Shifts and smearings
      //CMS_hgg_sig_m0_absShift, CMS_hbb_sig_m0_absShift, CMS_hgg_sig_sigmaScale, and CMS_hbb_sig_sigmaScale have already been defined when doing SigWS
/*      _w->factory(TString::Format("prod::CMS_hgg_hig_m0_%s_cat%d(mgg_hig_m0_%s_cat%d, CMS_hgg_sig_m0_absShift)", higName.Data(), newC, higName.Data(), newC));
      _w->factory(TString::Format("prod::CMS_hgg_hig_sigma_%s_cat%d(mgg_hig_sigma_%s_cat%d, CMS_hgg_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));*/
      _w->factory(TString::Format("prod::CMS_hgg_hig_m0_%s_cat%d(mgg_hig_m0_%s_cat%d, CMS_hgg_sig_m0_absShift)", higName.Data(), newC, higName.Data(), newC));
      _w->factory(TString::Format("prod::CMS_hgg_hig_sigma_%s_cat%d(mgg_hig_sigma_%s_cat%d, CMS_hgg_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));

//      if(!_useDSCB) _w->factory(TString::Format("prod::CMS_hgg_gsigma_%s_cat%d(mgg_hig_gsigma_%s_cat%d, CMS_hgg_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));
      if(!_useDSCB) _w->factory(TString::Format("prod::CMS_hgg_gsigma_%s_cat%d(mgg_hig_gsigma_%s_cat%d, CMS_hgg_sig_sigmaScale)",
                                                               higName.Data(), newC, higName.Data(), newC));
      if (higName.Contains("ggh") == 0 && higName.Contains("vbf") == 0) {
/*        _w->factory(TString::Format("prod::CMS_hbb_hig_m0_%s_cat%d(mjj_hig_m0_%s_cat%d, CMS_hbb_sig_m0_absShift)", higName.Data(), newC, higName.Data(), newC));
        _w->factory(TString::Format("prod::CMS_hbb_hig_sigma_%s_cat%d(mjj_hig_sigma_%s_cat%d, CMS_hbb_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));
        if(!_useDSCB) _w->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%s_cat%d(mjj_hig_gsigma_%s_cat%d, CMS_hbb_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));*/
        _w->factory(TString::Format("prod::CMS_hbb_hig_m0_%s_cat%d(mjj_hig_m0_%s_cat%d, CMS_hbb_sig_m0_absShift)", higName.Data(), newC, higName.Data(), newC));
        _w->factory(TString::Format("prod::CMS_hbb_hig_sigma_%s_cat%d(mjj_hig_sigma_%s_cat%d, CMS_hbb_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));
        if(!_useDSCB) _w->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%s_cat%d(mjj_hig_gsigma_%s_cat%d, CMS_hbb_sig_sigmaScale)", higName.Data(), newC, higName.Data(), newC));
      }

      TString EditPDF = TString::Format("EDIT::CMS_hig_%s_cat%d(HigPdf_%s_cat%d,", higName.Data(), newC, higName.Data(), c);
      for (unsigned int iv = 0; iv < varsToChange.size(); iv++)
        EditPDF += TString::Format("%s=%s,", varsToChange[iv].first.Data(), varsToChange[iv].second.Data());
      //Shifted and smeared vars
      if(higName.Contains("ggh") == 0 && higName.Contains("vbf") == 0) {
        if(!_useDSCB) EditPDF += TString::Format("mjj_hig_gsigma_%s_cat%d=CMS_hbb_hig_gsigma_%s_cat%d,", higName.Data(), c, higName.Data(), newC);
        EditPDF += TString::Format("mjj_hig_m0_%s_cat%d=CMS_hbb_hig_m0_%s_cat%d,", higName.Data(), c, higName.Data(), newC);
        EditPDF += TString::Format("mjj_hig_sigma_%s_cat%d=CMS_hbb_hig_sigma_%s_cat%d,", higName.Data(), c, higName.Data(), newC);
      }
      if(!_useDSCB) EditPDF += TString::Format("mgg_hig_gsigma_%s_cat%d=CMS_hgg_hig_gsigma_%s_cat%d,", higName.Data(), c, higName.Data(), newC);
      EditPDF += TString::Format("mgg_hig_m0_%s_cat%d=CMS_hgg_hig_m0_%s_cat%d,", higName.Data(), c, higName.Data(), newC);
      EditPDF += TString::Format("mgg_hig_sigma_%s_cat%d=CMS_hgg_hig_sigma_%s_cat%d)", higName.Data(), c, higName.Data(), newC);
      std::cout << "STRINGTOCHANGE   ---  " << EditPDF << std::endl;
      _w->factory(EditPDF);

      wAll->import(*_w->pdf(TString::Format("CMS_hig_%s_cat%d",higName.Data(),newC)));
      wAll->import(*_w->data(TString::Format("Hig_%s_cat%d",higName.Data(), c)), Rename(TString::Format("Hig_%s_cat%d", higName.Data(), newC)));

    }
  TString filename(wsDir+fileHiggsName+".inputhig.root");
  wAll->Print("v");
  wAll->writeToFile(filename);
  if (_verbLvl>1) std::cout << "Write signal workspace in: " << filename << " file" << std::endl;

  return;
} // close make higgs WP

void bbgg2DFitter::MakeBkgWS(std::string fileBaseName)
{
  TString wsDir = TString::Format("%s/workspaces/",_folder_name.data());
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  std::vector<RooDataSet*> data(_NCAT,nullptr);
  std::vector<RooAbsPdf*> BkgPdf(_NCAT,nullptr);
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < _NCAT; ++c)
    {
      int newC = c + _ncat0;
      BkgPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("BkgPdf_cat%d",c));
      RooArgSet* paramsMjj = (RooArgSet*) BkgPdf[c]->getParameters(*_w->var("mjj"));
      TIterator* iterMjj = (TIterator*) paramsMjj->createIterator();
      TObject* tempObjMjj = nullptr;
      std::vector<std::pair<TString,TString>> varsToChange;

      while( (tempObjMjj = iterMjj->Next()) ) {

        if ( (TString(tempObjMjj->GetName()).EqualTo("mjj")) || (TString(tempObjMjj->GetName()).EqualTo("mgg"))) continue;
        TString thisVarName(tempObjMjj->GetName());
        TString newVarName = TString(thisVarName).ReplaceAll(TString::Format("cat%d", c), TString::Format("cat%d", newC));
        varsToChange.push_back(std::make_pair(thisVarName, newVarName));
        std::cout << "Importing variable with new name: old - " << thisVarName << " new - " << newVarName << std::endl;
        _w->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));
        wAll->import( *_w->var( tempObjMjj->GetName() ), RenameVariable( thisVarName, newVarName));

      }

      TString EditPDF = TString::Format("EDIT::CMS_Bkg_cat%d(BkgPdf_cat%d", newC,  c);
      for (unsigned int iv = 0; iv < varsToChange.size(); iv++)
        EditPDF += TString::Format(",%s=%s", varsToChange[iv].first.Data(), varsToChange[iv].second.Data());
      EditPDF += ")";
      std::cout << "EDITSTRING: " << EditPDF << std::endl;
      _w->factory(EditPDF);

      wAll->import(*_w->pdf(TString::Format("CMS_Bkg_cat%d", newC)));
      wAll->import(*_w->var(TString::Format("BkgPdf_cat%d_norm", c)), RenameVariable(TString::Format("BkgPdf_cat%d_norm", c) , TString::Format("CMS_Bkg_cat%d_norm",newC)));
      wAll->import(*_w->data(TString::Format("Data_cat%d", c)), Rename(TString::Format("data_obs_cat%d", newC) ));

//      data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
      //RooDataHist* dataBinned = data[c]->binnedClone(); // Uncomment if you want to use wights in the limits
//      BkgPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("BkgPdf_cat%d",c));
//      wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));// Comment if you want to use wights in the limits
      //wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c))); // Uncomment if you want to use wights in the limits
//      wAll->import(*_w->pdf(TString::Format("BkgPdf_cat%d",c)));
//      wAll->import(*_w->var(TString::Format("BkgPdf_cat%d_norm",c)));

    } // close ncat

  TString filename(wsDir+fileBaseName+".root");
  wAll->writeToFile(filename);
  if (_verbLvl>1) std::cout << "Write background workspace in: " << filename << " file" <<std::endl;
  if (_verbLvl>1) std::cout << "observation ";
  for (int c = 0; c < _NCAT; ++c)
    {
      int newC = c + _ncat0;
      if (_verbLvl>1) std::cout << " " << wAll->data(TString::Format("data_obs_cat%d",newC))->sumEntries();
    }
  if (_verbLvl>1) std::cout << std::endl;
  return;
} // close make BKG workspace

void bbgg2DFitter::SetConstantParams(const RooArgSet* params)
{
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next())
    {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      if (rrv) rrv->setConstant(true); if (_verbLvl>1) std::cout << " " << rrv->GetName();
    }
} // close set const parameters

TStyle * bbgg2DFitter::style()
{
  // TODO: Need to get this outside of this class..

  TStyle *hggStyle = new TStyle("hggPaperStyle","Hgg Paper Style");

  //hggStyle->SetCanvasColor(0);
  //hggStyle->SetPadColor(0);
  hggStyle->SetFrameFillColor(0);
  hggStyle->SetStatColor(0);
  hggStyle->SetOptStat(0);
  hggStyle->SetTitleFillColor(0);
  hggStyle->SetCanvasBorderMode(0);
  hggStyle->SetPadBorderMode(0);
  hggStyle->SetFrameBorderMode(0);
  //hggStyle->SetFrameBorderSize(3);
  //hggStyle->SetPadBorderSize(3);
  //hggStyle->SetCanvasBorderSize(3);
  hggStyle->SetPadColor(kWhite);
  hggStyle->SetCanvasColor(kWhite);

  //hggStyle->SetCanvasDefH(600); //Height of canvas
  //hggStyle->SetCanvasDefW(600); //Width of canvas
  hggStyle->SetCanvasDefX(0);   //POsition on screen
  hggStyle->SetCanvasDefY(0);

  hggStyle->SetPadLeftMargin(0.10);//0.16);
  hggStyle->SetPadRightMargin(0.05);//0.02);
  hggStyle->SetPadTopMargin(0.04);//0.02);
  hggStyle->SetPadBottomMargin(0.09);//0.02);

  // For hgg axis titles:
  hggStyle->SetTitleColor(1, "XYZ");
  hggStyle->SetTitleFont(42, "XYZ");
  hggStyle->SetTitleSize(0.04, "XYZ");
  hggStyle->SetTitleXOffset(1.2);//0.9);
  hggStyle->SetTitleYOffset(1.2); // => 1.15 if exponents

  // For hgg axis labels:
  hggStyle->SetLabelColor(1, "XYZ");
  hggStyle->SetLabelFont(42, "XYZ");
  hggStyle->SetLabelOffset(0.007, "XYZ");
  //  hggStyle->SetLabelOffset(0.05, "YX");
  hggStyle->SetLabelSize(0.035, "XYZ");

  // Legends
  hggStyle->SetLegendBorderSize(0);
  hggStyle->SetLegendFillColor(kWhite);
  hggStyle->SetLegendFont(42);
  //  hggStyle->SetLegendTextSize(0.045);

  hggStyle->SetFillColor(10);
  // Nothing for now
  hggStyle->SetTextFont(42);
  hggStyle->SetTextSize(0.03);
  hggStyle->cd();
  return hggStyle;

}

RooFitResult* bbgg2DFitter::BkgModelFit(Bool_t dobands, bool addhiggs)
//RooFitResult* bbgg2DFitter::BkgModelFit(Bool_t dobands, bool addhiggs,std::vector<std::string>_singleHiggsNames,std::map<std::string,int>_singleHiggsMap )
{
  const Int_t ncat = _NCAT;
  std::vector<TString> catdesc;
  if ( _NCAT == 2 )catdesc={" High Purity Category"," Med. Purity Category"};
  if ( _NCAT == 1 )catdesc={" High Mass Analysis", " High Mass Analysis"};
  //  else catdesc={" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}",
  //	        " #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}"," #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}"};
  //
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  std::vector<RooDataSet*> data(ncat,nullptr);
  std::vector<RooDataSet*> dataplot(ncat,nullptr); // the data
  std::vector<RooBernstein*> mggBkg(ncat,nullptr);// the polinomial of 4* order
  std::vector<RooBernstein*> mjjBkg(ncat,nullptr);// the polinomial of 4* order
  std::vector<RooPlot*> plotmggBkg(ncat,nullptr);
  std::vector<RooPlot*> plotmjjBkg(ncat,nullptr);;
  std::vector<RooDataSet*>vecset(ncat,nullptr);
  std::vector<RooAbsPdf*>vecpdf(ncat,nullptr);
  std::vector<std::vector<RooDataSet*>>sigToFitvec(5,vecset);
  std::vector<std::vector<RooAbsPdf*>>mggSigvec(5,vecpdf);
  std::vector<std::vector<RooAbsPdf*>>mjjSigvec(5,vecpdf);
  std::vector<RooAbsPdf*> mggSig(ncat,nullptr);
  std::vector<RooAbsPdf*> mjjSig(ncat,nullptr);
  RooProdPdf* BkgPdf = nullptr;
  //  RooAbsPdf* BkgPdf1 = nullptr;
  //  RooAbsPdf* mjjBkgTmpPow1 = nullptr;
  //  RooAbsPdf* mggBkgTmpPow1 = nullptr;
  //  RooExponential* mjjBkgTmpExp1 = nullptr;
  //  RooExponential* mggBkgTmpExp1 = nullptr;
  RooBernstein* mjjBkgTmpBer1 = nullptr;
  RooBernstein* mggBkgTmpBer1 = nullptr;
  //Float_t minMggMassFit(100),maxMggMassFit(180);
  //Float_t minMjjMassFit(60),maxMjjMassFit(180);
  //  if(_sigMass == 260) _maxMggMassFit = 145;
  //  if(_sigMass == 270) _maxMggMassFit = 155;
  // Fit data with background pdf for data limit
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  //RooRealVar* mtot = _w->var("mtot");
  mgg->setUnit("GeV");
  mjj->setUnit("GeV");
  mgg->setRange("BkgFitRange",_minMggMassFit,_maxMggMassFit);
  mjj->setRange("BkgFitRange",_minMjjMassFit,_maxMjjMassFit);
  RooFitResult* fitresults = new RooFitResult();
  //
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //
  if (_verbLvl>1) std::cout << "[BkgModelFit] Starting cat loop " << std::endl;
  for (int c = 0; c < ncat; ++c) { // to each category
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));

    TH2* data_h2 = 0;
    TH1* data_h11 = 0;
    if(_fitStrategy==2)  data_h2= (TH2*) data[c]->createHistogram("mgg,mjj", 60, 40);
    if(_fitStrategy==1)  data_h11= (TH1*) data[c]->createHistogram("mgg", 60);

    if (_verbLvl>1) {
      std::cout<<"\t categ="<<c<<std::endl;
      if(_doblinding==0 && _fitStrategy==2) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h2->Integral() << std::endl;
      if(_doblinding==0 && _fitStrategy==1) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h11->Integral() << std::endl;
      std::cout<<"\t sumEntries()="<<data[c]->sumEntries()<<std::endl;
    }

    int nEvtsObs = -1;
    if(_fitStrategy == 2) nEvtsObs = data_h2->Integral();
    if(_fitStrategy == 1) nEvtsObs = data_h11->Integral();

    //data_h11->Delete();

    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop 1 - cat" << c << std::endl;

    ////////////////////////////////////
    // these are the parameters for the bkg polinomial
    // one slope by category - float from -10 > 10
    // we first wrap the normalization of mggBkgTmp0, mjjBkgTmp0
    // CMS_hhbbgg_13TeV_mgg_bkg_slope1
    _w->factory(TString::Format("BkgPdf_cat%d_norm[1.0,0.0,100000]",c));
    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop 2 - cat" << c << std::endl;
    /*
      RooFormulaVar *mgg_p0amp = new RooFormulaVar(TString::Format("mgg_p0amp_cat%d",c),"","@0*@0",
      *_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d",c)));
      RooFormulaVar *mgg_p1amp = new RooFormulaVar(TString::Format("mgg_p1amp_cat%d",c),"","@0*@0*@1",
      RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d",c)), *mgg_p0amp ));
      if(nEvtsObs > 10) { RooFormulaVar *mgg_p2amp = new RooFormulaVar(TString::Format("mgg_p2amp_cat%d",c),"","@0*@0*@1",
      RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d",c)), *mgg_p1amp)); }

      RooFormulaVar *mjj_p0amp = new RooFormulaVar(TString::Format("mjj_p0amp_cat%d",c),"","@0*@0",
      *_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat%d",c)));
      RooFormulaVar *mjj_p1amp = new RooFormulaVar(TString::Format("mjj_p1amp_cat%d",c),"","@0*@0*@1",
      RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d",c)), *mjj_p0amp ));
      if(nEvtsObs > 10) { RooFormulaVar *mjj_p2amp = new RooFormulaVar(TString::Format("mjj_p2amp_cat%d",c),"","@0*@0*@1",
      RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d",c)), *mjj_p0amp ));}
    */

    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop 3 - cat" << c << std::endl;

    RooFormulaVar *mgg_p0amp = new RooFormulaVar(TString::Format("mgg_p0amp_cat%d",c),"","@0*@0",
						            *_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d",c)));
    RooFormulaVar *mgg_p1amp = new RooFormulaVar(TString::Format("mgg_p1amp_cat%d",c),"","@0*@0",
						 RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d",c)) ));
    RooFormulaVar *mgg_p2amp = new RooFormulaVar(TString::Format("mgg_p2amp_cat%d",c),"","@0*@0",
						 RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d",c)) ));

    RooFormulaVar *mjj_p0amp = new RooFormulaVar(TString::Format("mjj_p0amp_cat%d",c),"","@0*@0",
						            *_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat%d",c)));
    RooFormulaVar *mjj_p1amp = new RooFormulaVar(TString::Format("mjj_p1amp_cat%d",c),"","@0*@0",
						 RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d",c)) ));
    RooFormulaVar *mjj_p2amp = new RooFormulaVar(TString::Format("mjj_p2amp_cat%d",c),"","@0*@0",
						 RooArgList(*_w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d",c)) ));


    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop 4 - cat" << c << std::endl;

    mggBkgTmpBer1 = new RooBernstein(TString::Format("mggBkgTmpBer1_cat%d",c),"",*mgg,RooArgList(*mgg_p0amp,*mgg_p1amp));
    mjjBkgTmpBer1 = new RooBernstein(TString::Format("mjjBkgTmpBer1_cat%d",c),"",*mjj,RooArgList(*mjj_p0amp,*mjj_p1amp));

    if(nEvtsObs > 15) {
      mggBkgTmpBer1 = new RooBernstein(TString::Format("mggBkgTmpBer1_cat%d",c),"",*mgg,RooArgList(*mgg_p0amp,*mgg_p1amp, *mgg_p2amp));
      mjjBkgTmpBer1 = new RooBernstein(TString::Format("mjjBkgTmpBer1_cat%d",c),"",*mjj,RooArgList(*mjj_p0amp,*mjj_p1amp, *mjj_p2amp));
    }

    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop 5 - cat" << c << std::endl;

    if(_fitStrategy==2) BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpBer1, *mjjBkgTmpBer1));
    //    if(_fitStrategy==1) BkgPdf1 = (RooAbsPdf*) mggBkgTmpBer1->Clone(TString::Format("BkgPdf_cat%d",c));

    if (_verbLvl>1) std::cout << "[BkgModelFit] Cat loop " << c << std::endl;

    RooExtendPdf* BkgPdfExt;

    if(_fitStrategy == 2) {
      BkgPdfExt = new RooExtendPdf(TString::Format("BkgPdfExt_cat%d",c),"", *BkgPdf,*_w->var(TString::Format("BkgPdf_cat%d_norm",c)));
      fitresults = BkgPdfExt->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("BkgFitRange"),SumW2Error(kTRUE), Save(kTRUE),PrintLevel(-1));
      _w->import(*BkgPdfExt);
      //	delete BkgPdfExt;
    }

    if(_fitStrategy == 1) {
      fitresults = mggBkgTmpBer1->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("BkgFitRange"),SumW2Error(kTRUE), Save(kTRUE),PrintLevel(-1));
      RooAbsPdf* BkgPdf1 = (RooAbsPdf*) mggBkgTmpBer1->Clone(TString::Format("BkgPdf_cat%d",c));
      _w->import(*BkgPdf1);
    }

    if (data_h2) data_h2->Delete();
    if (data_h11) data_h11->Delete();
  }

  // Don't make plots, exit here!
  return fitresults;

  for (int c = 0; c < ncat; ++c) { // to each category

    TH2* data_h2 = 0;
    TH1* data_h11 = 0;

    if (_verbLvl>1) std::cout << "[BkgModelFit] Done with fit " << c << std::endl;

    ///////////////////////////////////
    //Calculate 2D chisquare by hand //
    ///////////////////////////////////
    //    TH2* pdf_h2 = BkgPdfExt.createHistogram("mgg vs mjj pdf", mgg, Binning(60), YVar(mjj, Binning(40)));
    //    TH2* data_h2 = (TH2*) data[c]->createHistogram("mgg,mjj", 60, 40);
    //    if (_verbLvl>1) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h2->Integral() << std::endl;
    if(_fitStrategy == 2) {
      TH2* pdf_h2 = (TH2*) BkgPdf->createHistogram("mgg,mjj", 60, 40);
      TH1F* data_h1 = new TH1F("data_h1", "data_h1", 2400, 0, 2400);
      TH1F* pdf_h1 = new TH1F("pdf_h1", "pdf_h1", 2400, 0, 2400);
      unsigned int nbins_data_x = data_h2->GetNbinsX();
      unsigned int nbins_pdf_x  =  pdf_h2->GetNbinsX();
      unsigned int nbins_data_y = data_h2->GetNbinsY();
      unsigned int nbins_pdf_y  =  pdf_h2->GetNbinsY();
      //    float Total2DChiSquare = 0;
      int counterBin = 1;
      if( nbins_data_x != nbins_pdf_x || nbins_data_y != nbins_pdf_y ) {
	if (_verbLvl>1) std::cout << "Number of bins for 2D chi square are different!!!!" << std::endl;
      } else {
	for (unsigned int nbx = 1; nbx < nbins_pdf_x + 1; nbx++) {
	  for ( unsigned int nby = 1; nby < nbins_pdf_y + 1; nby++) {
	    float observed_h2 = data_h2->GetBinContent(nbx, nby);
	    float errObs_h2 = data_h2->GetBinError(nbx, nby);
	    float expected_h2 = pdf_h2->GetBinContent(nbx, nby);
	    float errExp_h2 = pdf_h2->GetBinError(nbx, nby);
	    data_h1->SetBinContent(counterBin, observed_h2);
	    data_h1->SetBinError(counterBin, errObs_h2);
	    pdf_h1->SetBinContent(counterBin, expected_h2);
	    pdf_h1->SetBinError(counterBin, errExp_h2);
	    counterBin++;
	    //if (_verbLvl>1) std::cout << "bin x: " << nbx << " bin y: " << nby << " exp: " << expected_h2 << " obs: " << observed_h2 << std::endl;
	    // float csq = (observed_h2 - expected_h2)*(observed_h2 - expected_h2)/expected_h2;
	    // Total2DChiSquare += csq;
	  }
	}
      }

      float pdfNormalization = pdf_h1->Integral();
      float pdfNormalization2 = pdf_h2->Integral();
      if (!pdf_h1->GetSumw2N()) pdf_h1->Sumw2();
      if (!pdf_h1->GetSumw2N()) pdf_h2->Sumw2();
      pdf_h1->Scale(data_h1->Integral()/pdfNormalization);
      pdf_h2->Scale(data_h2->Integral()/pdfNormalization2);
      if (_verbLvl>1) {
	std::cout << "########NUMBER OF PDFFF EVENTSSSS::: " << pdf_h2->Integral() << std::endl;
	std::cout << "################################################################" << std::endl;
	std::cout << "#############  KolmogorovTests  #############################" << std::endl;
	std::cout << "################ 1DKSTEST:" << data_h1->KolmogorovTest(pdf_h1, "") << " ##################################" << std::endl;
	std::cout << "################ 2DKSTEST:" << data_h2->KolmogorovTest(pdf_h2, "") << " ##################################" << std::endl;
	std::cout << "################################################################" << std::endl;
	std::cout << "################################################################" << std::endl;
      }

      pdf_h1->Delete();
      pdf_h2->Delete();
      data_h1->Delete();

      if (data_h2) data_h2->Delete();
      if (data_h11) data_h11->Delete();

    }


    //BkgPdf.fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("BkgFitRange"),SumW2Error(kTRUE), Save(kTRUE));
    //w->import(BkgPdf);
    //************************************************//
    // Plot mgg background fit results per categories
    //************************************************//
    if (_verbLvl>1) std::cout << "[BkgModelFit] Plotting Mgg - cat" << c << std::endl;
    //TCanvas* ctmp = new TCanvas(TString::Format("ctmpBkgMgg_cat%d",c),"mgg Background Categories",800,600);
    Int_t nBinsMass(80);
    plotmggBkg[c] = mgg->frame(nBinsMass);
    dataplot[c] = (RooDataSet*) _w->data(TString::Format("Dataplot_cat%d",c));
    //    if(_sigMass == 260) plotmggBkg[c]->GetXaxis()->SetRangeUser(100.,145.);
    //    if(_sigMass == 270) plotmggBkg[c]->GetXaxis()->SetRangeUser(100.,155.);
    if(_doblinding) dataplot[c]->plotOn(plotmggBkg[c],Invisible());
    else dataplot[c]->plotOn(plotmggBkg[c]);
    /*    if(_sigMass == 0)
	  {
	  if(c == 0 || c == 2) mggBkgTmpExp1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
	  else if(c == 1 || c == 3) mggBkgTmpPow1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
	  }*/
    if(_sigMass > -2)
      {
	mggBkgTmpBer1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
	//       else if(c == 1) mggBkgTmpPow1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
      }
    if(_doblinding) dataplot[c]->plotOn(plotmggBkg[c], Invisible());
    else dataplot[c]->plotOn(plotmggBkg[c]);
    // now we fit the gaussian on signal
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    //    if(c==0||c==2)plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    //    if(c==1||c==3)plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    plotmggBkg[c]->Draw();
    //plotmggBkg[c]->SetTitle("CMS preliminary 19.7/fb");
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    plotmggBkg[c]->SetMaximum(1.40*plotmggBkg[c]->GetMaximum());
    plotmggBkg[c]->GetXaxis()->SetTitle("M(#gamma#gamma) [GeV]");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
    TPaveText *pt = new TPaveText(0.1,0.92,0.9,0.99, "brNDC");
    // pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextSize(0.035);
    pt->AddText(TString::Format("            CMS Preliminary                     L = %.2f fb^{-1}    #sqrt{s} = %s",_lumi,_energy.c_str()));
    pt->Draw();
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
    if (dobands)
      {
    	RooAbsPdf *cpdf=mggBkgTmpBer1;
	//        cpdf = mggBkgTmpExp1;
	/*
	  if(_sigMass == 0 && (c == 0 || c == 2)) cpdf = mggBkgTmpExp1;
	  if(_sigMass == 0 && (c == 1 || c == 3)) cpdf = mggBkgTmpPow1;
	  if(_sigMass != 0 && c == 0) cpdf = mggBkgTmpBer1;
	  if(_sigMass != 0 && c == 1) cpdf = mggBkgTmpPow1;*/
	//      onesigma = new TGraphAsymmErrors();
	//      twosigma = new TGraphAsymmErrors();
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
	nlim->removeRange();
	RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
	for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i)
	  {
	    double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
	    double upedge = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
	    double center = plotmggBkg[c]->GetXaxis()->GetBinCenter(i);
	    double nombkg = nomcurve->interpolate(center);
	    nlim->setVal(nombkg);
	    mgg->setRange("errRange",lowedge,upedge);
	    RooAbsPdf *epdf = 0;
	    epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	    RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	    RooMinimizer minim(*nll);
	    minim.setStrategy(0);
	    //double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	    double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	    minim.migrad();
	    minim.minos(*nlim);
	    // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	    onesigma->SetPoint(i-1,center,nombkg);
	    onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	    minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
	    // the 0.5 is because qmu is -2*NLL
	    // eventually if cl = 0.95 this is the usual 1.92!
	    minim.migrad();
	    minim.minos(*nlim);
	    twosigma->SetPoint(i-1,center,nombkg);
	    twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	    delete nll;
	    delete epdf;
	  } // close for bin
	mgg->setRange("errRange",_minMggMassFit,_maxMggMassFit);
	twosigma->SetLineColor(kYellow);
	twosigma->SetFillColor(kYellow);
	twosigma->SetMarkerColor(kYellow);
	twosigma->Draw("L3 SAME");
	onesigma->SetLineColor(kGreen);
	onesigma->SetFillColor(kGreen);
	onesigma->SetMarkerColor(kGreen);
	onesigma->Draw("L3 SAME");
	plotmggBkg[c]->Draw("SAME");
      }
    else plotmggBkg[c]->Draw("SAME"); // close dobands
    //plotmggBkg[c]->getObject(1)->Draw("SAME");
    //plotmggBkg[c]->getObject(2)->Draw("P SAME");
    ////////////////////////////////////////////////////////// plot higgs
    /*
    if(addhiggs) {
      for(unsigned int d=0;d!=_singleHiggsNames.size();++d)
    	{
	  static std::vector<int>color{2,3,6,7,4};
	  int realint=_singleHiggsMap[_singleHiggsNames[d]];
	  sigToFitvec[realint][c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",realint,c));
	  double norm = 1.0*sigToFitvec[realint][c]->sumEntries(); //
	  mggSigvec[realint][c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%d_cat%d",realint,c));
	  // we are not constructing signal pdf, this is constructed on sig to fit function...
	  mggSigvec[realint][c]->plotOn(plotmggBkg[c],Normalization(norm,RooAbsPdf::NumEvent),DrawOption("F"),LineColor(color[realint]),FillStyle(1001),FillColor(19));
	  mggSigvec[realint][c]->plotOn(plotmggBkg[c],Normalization(norm,RooAbsPdf::NumEvent),LineColor(color[realint]),LineStyle(1));
	}
    }
    */
    //////////////////////////////////////////////////////////
    plotmggBkg[c]->Draw("SAME");
    if(c==0||c==2)plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0||c==2)plotmggBkg[c]->SetMaximum(5.3); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMaximum(20); // no error bar in bins with zero events
    // plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmggBkg[c]->SetLogy(0);
    TLegend *legmc = new TLegend(0.40,0.72,0.62,0.9);
    TLegend *legmcH = new TLegend(0.66,0.72,0.94,0.9);
    if(_doblinding) legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","");
    else  legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","LPE");
    legmc->AddEntry(plotmggBkg[c]->getObject(1),"Bkg Fit","L");
    if(dobands)legmc->AddEntry(onesigma,"Fit #pm1 #sigma","F");
    if(dobands)legmc->AddEntry(twosigma,"Fit #pm2 #sigma ","F"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(11),"bbH ","LPE"); // not...
    if(_nonResWeightIndex>=-1)
      legmc->SetHeader(" Nonresonant HH");
    else
      legmc->SetHeader(TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), _sigMass));
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmcH->SetBorderSize(0);
    legmcH->SetFillStyle(0);
    legmc->Draw();
    legmcH->Draw();
    TLatex *lat2 = new TLatex(_minMggMassFit+1.5,0.85*plotmggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    //
    _c1->SaveAs(TString::Format("%s/databkgoversigMgg_cat%d.pdf",_folder_name.data(),c),"QUIET");
    _c1->SaveAs(TString::Format("%s/databkgoversigMgg_cat%d.png",_folder_name.data(),c),"QUIET");

    if(c==0||c==2)plotmggBkg[c]->SetMaximum(100); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
    _c1->SetLogy(1);
    _c1->SaveAs(TString::Format("%s/databkgoversigMgg_cat%d_log.pdf",_folder_name.data(),c),"QUIET");
    _c1->SaveAs(TString::Format("%s/databkgoversigMgg_cat%d_log.png",_folder_name.data(),c),"QUIET");
    // ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.C",c));
    _c1->SetLogy(0);

    //if (ctmp) ctmp->Close();

    if(_fitStrategy == 2) {
      if (_verbLvl>1 ) std::cout << "[BkgModelFit] Plotting Mgg - cat" << c << std::endl;
      //************************************************//
      // Plot mjj background fit results per categories
      //************************************************//
      //ctmp = new TCanvas(TString::Format("ctmpBkgMjj_cat%d",c),"mjj Background Categories",0,0,500,500);
      nBinsMass = 60;
      plotmjjBkg[c] = mjj->frame(nBinsMass);
      dataplot[c] = (RooDataSet*) _w->data(TString::Format("Dataplot_cat%d",c));
      if(_doblinding) dataplot[c]->plotOn(plotmjjBkg[c],Invisible());
      else dataplot[c]->plotOn(plotmjjBkg[c]);
      if(_sigMass > -2)
	{
          mjjBkgTmpBer1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
	  //       else if(c == 1) mjjBkgTmpPow1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
	}
      if(_doblinding) dataplot[c]->plotOn(plotmjjBkg[c],Invisible());
      else dataplot[c]->plotOn(plotmjjBkg[c]);
      // now we fit the gaussian on signal
      //plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      if(c==0||c==2)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1||c==3)plotmjjBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
      plotmjjBkg[c]->Draw();
      //plotmjjBkg[c]->SetTitle("CMS preliminary 19.7/fb");
      //plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      plotmjjBkg[c]->SetMaximum(1.40*plotmjjBkg[c]->GetMaximum());
      plotmjjBkg[c]->GetXaxis()->SetTitle("M(jj) [GeV]");
      //double test = sigToFit[c]->sumEntries();
      //cout<<"number of events on dataset "<<test<<endl;
      pt = new TPaveText(0.1,0.94,0.9,0.99, "brNDC");
      // pt->SetName("title");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetTextSize(0.035);
      pt->AddText(TString::Format("            CMS Preliminary                     L = %.2f fb^{-1}    #sqrt{s} = %s",_lumi,_energy.c_str()));
      pt->Draw();
      if (dobands)
	{
	  RooAbsPdf *cpdf = mjjBkgTmpBer1;
          TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
          TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
          RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
          nlim->removeRange();
          RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmjjBkg[c]->getObject(1));
          for (int i=1; i<(plotmjjBkg[c]->GetXaxis()->GetNbins()+1); ++i)
	    {
	      double lowedge = plotmjjBkg[c]->GetXaxis()->GetBinLowEdge(i);
	      double upedge = plotmjjBkg[c]->GetXaxis()->GetBinUpEdge(i);
	      double center = plotmjjBkg[c]->GetXaxis()->GetBinCenter(i);
	      double nombkg = nomcurve->interpolate(center);
	      nlim->setVal(nombkg);
	      mjj->setRange("errRange",lowedge,upedge);
	      RooAbsPdf *epdf = 0;
	      epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	      RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	      RooMinimizer minim(*nll);
	      minim.setStrategy(0);
	      //double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	      double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	      minim.migrad();
	      minim.minos(*nlim);
	      // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	      onesigma->SetPoint(i-1,center,nombkg);
	      onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	      minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
	      // the 0.5 is because qmu is -2*NLL
	      // eventually if cl = 0.95 this is the usual 1.92!
	      minim.migrad();
	      minim.minos(*nlim);
	      twosigma->SetPoint(i-1,center,nombkg);
	      twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	      delete nll;
	      delete epdf;
	    } // close for bin
	  mjj->setRange("errRange",_minMjjMassFit,_maxMjjMassFit);
	  twosigma->SetLineColor(kYellow);
	  twosigma->SetFillColor(kYellow);
	  twosigma->SetMarkerColor(kYellow);
	  twosigma->Draw("L3 SAME");
	  onesigma->SetLineColor(kGreen);
	  onesigma->SetFillColor(kGreen);
	  onesigma->SetMarkerColor(kGreen);
	  onesigma->Draw("L3 SAME");
	  plotmjjBkg[c]->Draw("SAME");
	}
      else plotmjjBkg[c]->Draw("SAME"); // close dobands

      plotmjjBkg[c]->Draw("SAME");
      if(c==0||c==2)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1||c==3)plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      if(c==0||c==2)plotmjjBkg[c]->SetMaximum(5.3); // no error bar in bins with zero events
      if(c==1||c==3)plotmjjBkg[c]->SetMaximum(20); // no error bar in bins with zero events
      // plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      //plotmjjBkg[c]->SetLogy(0);
      legmc = new TLegend(0.40,0.72,0.62,0.9);
      legmcH = new TLegend(0.66,0.72,0.94,0.9);
      if(_doblinding) legmc->AddEntry(plotmjjBkg[c]->getObject(2),"Data ","");
      else legmc->AddEntry(plotmjjBkg[c]->getObject(2),"Data ","LPE");
      legmc->AddEntry(plotmjjBkg[c]->getObject(1),"Fit","L");
      if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); // not...
      if(dobands)legmc->AddEntry(onesigma,"one sigma","F");

      if(addhiggs){
	legmcH->AddEntry(plotmjjBkg[c]->getObject(3),"ggH ","LPE"); // not...
	legmcH->AddEntry(plotmjjBkg[c]->getObject(5),"ttH ","LPE"); // not...
	legmcH->AddEntry(plotmjjBkg[c]->getObject(7),"VBF ","LPE"); // not...
	legmcH->AddEntry(plotmjjBkg[c]->getObject(9),"VH ","LPE"); // not...
	legmcH->AddEntry(plotmjjBkg[c]->getObject(11),"bbH ","LPE"); // not...
      }

      if(_nonResWeightIndex>=-1)legmc->SetHeader(" Nonresonant HH");
      else legmc->SetHeader(TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), _sigMass));
      legmcH->SetHeader(" Higgs");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmcH->SetBorderSize(0);
      legmcH->SetFillStyle(0);
      legmc->Draw();
      legmcH->Draw();
      lat2 = new TLatex(_minMjjMassFit+1.5,0.85*plotmjjBkg[c]->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      //
      _c1->SaveAs(TString::Format("%s/databkgoversigMjj_cat%d.pdf",_folder_name.data(),c),"QUIET");
      _c1->SaveAs(TString::Format("%s/databkgoversigMjj_cat%d.png",_folder_name.data(),c),"QUIET");

      if(c==0||c==2)plotmjjBkg[c]->SetMaximum(100); // no error bar in bins with zero events
      if(c==1||c==3)plotmjjBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
      _c1->SetLogy(1);
      _c1->SaveAs(TString::Format("%s/databkgoversigMjj_cat%d_log.pdf",_folder_name.data(),c),"QUIET");
      _c1->SaveAs(TString::Format("%s/databkgoversigMjj_cat%d_log.png",_folder_name.data(),c),"QUIET");
      // _c1->SaveAs(TString::Format("databkgoversigMjj_cat%d.C",c));
      _c1->SetLogy(0);

      //if (ctmp) ctmp->Close();

    }//close if fit strategy == 2

    //if (ctmp) ctmp->Close();

  } // close to each category
  return fitresults;
} // close berestein 3

void bbgg2DFitter::MakeFitsForBias(std::string biasConfig, std::string outputFile)
{
  std::cout<<"\n ** Doing fits for Bias study** \n\n"
	   <<"biasConfig="<<biasConfig
	   <<"outputFile="<<outputFile<<std::endl;

  //Parameters
  std::string plotsDir=_folder_name;
  std::vector<std::string> vars;
  std::vector<std::string> varNames;
  std::vector<std::string> nbins;
  std::vector<std::string> functions;
  std::vector<std::string> functionsToFit;
  std::vector<std::string> functionsToPlot;
  std::vector<std::string> biasFunctions;
  std::string plotTitle="";
  std::vector<std::string> legends;
  double bkgNorm_up, bkgNorm_down, bkgNorm;

  //Read json
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( biasConfig, pt );
  bkgNorm_up = TString( pt.get_child("bkgNorm_up").data().c_str()).Atof();
  bkgNorm_down = TString( pt.get_child("bkgNorm_down").data().c_str()).Atof();
  bkgNorm = TString( pt.get_child("bkgNorm").data().c_str()).Atof();


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
    return ;
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

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functionLegends" ) ){
    legends.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "biasFunctions" ) ){
    biasFunctions.push_back(rowPair.second.data());
  }

  TFile * tFile = new TFile(outputFile.data(), "RECREATE");

  RooArgSet* TreeVars = new RooArgSet();
  const Int_t ncat = _NCAT;
  for ( int cat = 0; cat < ncat; cat++) {
      RooWorkspace* wBias = new RooWorkspace(TString::Format("wBias_cat%d", cat));

      for( unsigned int i = 0; i < vars.size(); i++){
          // Initialize fit variable
          wBias->factory( vars[i].c_str() );
          std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
          std::cout << "Variable to be fit: " << sVar << std::endl;
          wBias->var( sVar.c_str() )->SetTitle(varNames[i].c_str());
          wBias->var( sVar.c_str() )->setUnit("GeV");
          TreeVars->add( *wBias->var( sVar.c_str() ) );
      }

      //Get data Distribution from main Workspace
      RooDataSet* dataBias = (RooDataSet*) _w->data(TString::Format("Data_cat%d",cat));
      wBias->import(*dataBias);

      //Initialize all needed PDFs for background
      for ( unsigned int f = 0; f < functions.size(); f++){
        TString thisFunction(functions[f]);
        std::cout << "Function added: " << thisFunction << std::endl;
        wBias->factory( thisFunction.Data() );
      }
      //Fit all functions
      std::vector<bbggFittingTools::FitRes> bkgresults = bbggFittingTools::FitFunctions(wBias, functionsToFit, dataBias);
      std::ofstream myFitResults(std::string(TString::Format("%s/bias/biasFitsResults_cat%d.txt", plotsDir.c_str(),cat).Data() ), std::ofstream::out);
      myFitResults << std::string(TString::Format("#Tot events: %d\n",dataBias->numEntries()));
      myFitResults << "#Name \t chi2 \t minNLL\n";
      for(unsigned int i = 0; i < bkgresults.size(); i++){
          myFitResults << std::string(TString::Format("%s\t%f\t%f\n", bkgresults[i].function.c_str(), bkgresults[i].chi2, bkgresults[i].minNLL)).c_str();
      }
      myFitResults.close();

      //plot fitted functions
      for( unsigned int i = 0; i < vars.size(); i++){
          std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
          bbggFittingTools::PlotCurves(plotTitle, wBias, functionsToPlot, legends, bkgresults, dataBias, sVar, nbins[i], plotsDir+"/bias/BIASplot_bkg_"+sVar+ TString(plotTitle).ReplaceAll(" ", "_").Data()+TString::Format("_cat%d",cat).Data(), 0, 1);
      }

      //Do bias study business: create multipdf, etc
      RooCategory pdf_index("pdf_index","Index of Pdf which is active");
      RooArgList mypdfs;
      for( unsigned int bias = 0; bias < biasFunctions.size(); bias++){
          std::cout << "Adding functions to multipdf! " << biasFunctions[bias].c_str() << std::endl;
          const char* modelName = biasFunctions[bias].c_str();
          mypdfs.add(* wBias->pdf( modelName ) );
      }
      RooMultiPdf multipdf("roomultipdf", "All Pdfs", pdf_index, mypdfs);
      RooRealVar norm("roomultipdf_norm", "Number of background events", bkgNorm, bkgNorm_down, bkgNorm_up);
//      RooRealVar norm("roomultipdf_norm", "Number of background events", 0.000001, 5000000);

      //Create bias workspace
      RooWorkspace * bW = new RooWorkspace(TString::Format("BiasWorkspace_cat%d",cat));
      bW->import(pdf_index);
      bW->import(norm);
      bW->import(multipdf);
      tFile->cd();
      bW->Write();
  }
  tFile->Close();

  std::cout<<"\n \t* Finished with Bias study fits *\n"<<endl;

}
