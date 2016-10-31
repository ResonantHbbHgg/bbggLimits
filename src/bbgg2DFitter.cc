#define bbgg2DFitter_cxx
#include "HiggsAnalysis/bbggLimits/interface/bbgg2DFitter.h"
#include "HiggsAnalysis/bbggLimits/interface/Colors.h"

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

  if (_verbLvl>=1 && _verbLvl<4)
    std::cout << "#Sigma effective: xLow: " << low << ", xHigh: " << high << ", width: " << width << std::endl;

  std::vector<float> outVec;
  outVec.push_back(width);
  outVec.push_back(low);
  outVec.push_back(high);

  return outVec;
}


void bbgg2DFitter::PrintWorkspace() {_w->Print("v");}

void bbgg2DFitter::Initialize(RooWorkspace* workspace, Int_t SigMass, float Lumi,std::string folder_name,
			      std::string energy, Bool_t doBlinding, Int_t nCat, bool AddHiggs,
			      float minMggMassFit,float maxMggMassFit,float minMjjMassFit,float maxMjjMassFit,
			      float minSigFitMgg,float maxSigFitMgg,float minSigFitMjj,float maxSigFitMjj,
			      float minHigMggFit,float maxHigMggFit,float minHigMjjFit,float maxHigMjjFit,
			      Int_t doNRW)
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
  // 0-1507: Non-resonant limits with re-weighting
  if (_nonResWeightIndex>=0)
    _wName = Form("evWeight_NRW_%d",doNRW);
  else
    _wName = "evWeight";

  //std::cout<<"DBG.  Finished Initialize..."<<std::endl;


  _c1 = new TCanvas("c1","Square Canvas",800,800);
  _c2 = new TCanvas("c2","Rectangular Canvas",800,600);
  
}

RooArgSet* bbgg2DFitter::defineVariables()
{
  RooRealVar* mgg  = new RooRealVar("mgg","M(#gamma#gamma)",_minMggMassFit,_maxMggMassFit,"GeV");
  RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  RooRealVar* mjj  = new RooRealVar("mjj","M(jj)",_minMjjMassFit,_maxMjjMassFit,"GeV");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;

  RooRealVar* evWeight = new RooRealVar(_wName.c_str(),"HqT x PUwei",0.,100,"");
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  cut_based_ct->defineType("cat4_2",2);
  cut_based_ct->defineType("cat4_3",3);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight);
  if (_nonResWeightIndex>=0)
    ntplVars->add(*mtot);
  // AP: Why these are here? They are already in the set:
  //ntplVars->add(*mgg);
  //ntplVars->add(*mjj);
  //ntplVars->add(*cut_based_ct);

  return ntplVars;
}

int bbgg2DFitter::AddSigData(float mass, TString signalfile)
{
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "================= Add Signal==============================" <<std::endl;
  //Luminosity
  RooRealVar lumi("lumi","lumi", _lumi);
  _w->import(lumi);
  //Define variables
  RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
  //Signal file & tree
  TFile sigFile(signalfile);
  bool opened=sigFile.IsOpen();
  if(opened==false) return -1;
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  if(sigTree==nullptr)
    {
      if (_verbLvl>=1 && _verbLvl<4) std::cout<<"TCVARS for AddSigData  not founded in TTree trying with TCVARS"<<std::endl;
      std::exit(1);
      //sigTree = (TTree*) sigFile.Get("TCVARS");
      //if(sigTree==nullptr)
    }
  //Data set

  //Double_t W;
  //ccbar->SetBranchAddress("weight", &wCCBar);
  //ccbar->GetEntry();
  //RooRealVar ccbarweight("NRweight", "NRweight", );

  if (_verbLvl==4) {
    std::cout<<"[DBG]  Prining ntplVars from sig"<<std::endl;
    ntplVars->Print();
  }

  RooDataSet sigScaled("sigScaled","dataset",sigTree,*ntplVars,_cut, _wName.c_str());


  RooDataSet* sigToFit[_NCAT];
  TString cut0 = " && 1>0";

  RooArgList myArgList(*_w->var("mgg"));

  if (_fitStrategy != 1)
    myArgList.add(*_w->var("mjj"));

  if (_nonResWeightIndex>=0)
    myArgList.add(*_w->var("mtot"));

  for ( int i=0; i<_NCAT; ++i)
    {

      sigToFit[i] = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString::Format(" && cut_based_ct==%d ",i)+cut0);
      if (_fitStrategy == 1)
	sigToFit[i] = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString::Format(" && cut_based_ct==%d && mjj < 140 ",i)+cut0);

      this->SetSigExpectedCats(i, sigToFit[i]->sumEntries());
      if (_verbLvl==4) {
	std::cout << "======================================================================" <<std::endl;
	std::cout<<"[DBG]  Cat="<<i<< "\t Sig sumEntries="<<sigToFit[i]->sumEntries()<<std::endl;
	std::cout<<"mGG:  Mean = "<<sigToFit[i]->mean(*_w->var("mgg"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mgg"))<<std::endl;
	if (_fitStrategy != 1)
	  std::cout<<"mJJ:  Mean = "<<sigToFit[i]->mean(*_w->var("mjj"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mjj"))<<std::endl;

	if (_nonResWeightIndex>=0)
	  std::cout<<"mTot: Mean = "<<sigToFit[i]->mean(*_w->var("mtot"))<<"  sigma = "<<sigToFit[i]->sigma(*_w->var("mtot"))<<std::endl;
      }

      /*This defines each category*/
      _w->import(*sigToFit[i],Rename(TString::Format("Sig_cat%d",i)));
    }
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(myArgList,_cut);
  if (_fitStrategy == 1)
    sigToFitAll = (RooDataSet*) sigScaled.reduce(myArgList,_cut+TString(" && mjj < 140 "));

  _w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  if (_verbLvl>=1 && _verbLvl<4) {
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

void bbgg2DFitter::AddHigData(float mass, TString signalfile, int higgschannel)
{
  RooArgSet* ntplVars = defineVariables();
  TFile higFile(signalfile);
  TTree* higTree = (TTree*) higFile.Get("TCVARS");
  if(higTree==nullptr)
    {
      if (_verbLvl>=1 && _verbLvl<4) std::cout<<"TCVARS for AddHigData  not founded in TTree trying with TCVARS"<<std::endl;
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
      _w->import(*higToFit[i],Rename(TString::Format("Hig_%d_cat%d",higgschannel,i)));
    }
  // Create full signal data set without categorization
  RooDataSet* higToFitAll = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut);
  if(_fitStrategy == 1) higToFitAll = (RooDataSet*) higScaled.reduce(RooArgList(*_w->var("mgg")),_cut + TString(" && mjj < 140 "));
  _w->import(*higToFitAll,Rename("Hig"));
  // here we print the number of entries on the different categories
  if (_verbLvl>=1 && _verbLvl<4) {
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

}

void bbgg2DFitter::AddBkgData(TString datafile)
{
  //Define variables
  RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
  RooRealVar weightVar("weightVar","",1,0,1000);
  weightVar.setVal(1.);
  TFile dataFile(datafile);
  TTree* dataTree = (TTree*) dataFile.Get("TCVARS");
  if(dataTree==nullptr)
    {
      if (_verbLvl>=1 && _verbLvl<4) std::cout<<"TCVARS for AddBkgData  not founded in TTree trying with TCVARS"<<std::endl;
      std::exit(1);
      //dataTree = (TTree*) dataFile.Get("TCVARS");
      //if(dataTree==nullptr)std::exit(1);
    }
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","evWeight");
  RooDataSet* dataToFit[_NCAT];
  RooDataSet* dataToPlot[_NCAT];
  TString cut0 = "&& 1>0";
  TString cut1 = "&& 1>0";
  if (_verbLvl>=1 && _verbLvl<4) std::cout<<"================= Add Bkg ==============================="<<std::endl;

  for( int i=0; i<_NCAT; ++i)
    {

      dataToFit[i] = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut+TString::Format(" && cut_based_ct==%d",i));

      this->SetObservedCats(i, dataToFit[i]->sumEntries());

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
  if (_verbLvl>=1 && _verbLvl<4) data->Print("v");
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
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "OK up to now... Mass point: " <<mass<<std::endl;
      if(_fitStrategy == 2) SigPdf[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));
      if(_fitStrategy == 1) SigPdf1[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));

      if (_verbLvl>=1 && _verbLvl<4) std::cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() <<std::endl;
      double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+(mass-125.0); // shift the peak //This should be an option

      ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "mPeak = " << mPeak << std::endl;
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() <<std::endl;

      // IMPORTANT: fix all pdf parameters to constant, why?
      RooArgSet sigParams( *_w->var(TString::Format("mgg_sig_m0_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_n_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
			   *_w->var(TString::Format("mgg_sig_frac_cat%d",c)));
      if(_fitStrategy == 2) {
	sigParams.add(RooArgSet(
				*_w->var(TString::Format("mjj_sig_m0_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_sigma_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_alpha_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_n_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_gsigma_cat%d",c)),
				*_w->var(TString::Format("mjj_sig_frac_cat%d",c))) );
      }
      _w->defineSet(TString::Format("SigPdfParam_cat%d",c), sigParams);
      SetConstantParams(_w->set(TString::Format("SigPdfParam_cat%d",c)));
      if (_verbLvl>=1 && _verbLvl<4) std::cout<<std::endl;
      if(_fitStrategy == 2) _w->import(*SigPdf[c]);
      if(_fitStrategy == 1) _w->import(*SigPdf1[c]);
    }
}

void bbgg2DFitter::HigModelFit(float mass, int higgschannel)
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
      higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));
      mggHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c));
      if(higgschannel == 1 || higgschannel == 3) mjjHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",higgschannel,c));
      else mjjHig[c] = new RooPolynomial(TString::Format("mjjHig_pol0_%d_cat%d",higgschannel,c),"",*mjj,RooArgList());;
      HigPdf[c] = new RooProdPdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c),"",RooArgSet(*mggHig[c], *mjjHig[c]));
      //((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(MASS);
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "OK up to now... Mass point: " <<mass<<std::endl;
      HigPdf[c]->fitTo(*higToFit[c],Range("HigFitRange"),SumW2Error(kTRUE),PrintLevel(-1));
      RooArgSet* paramsMjj;
      paramsMjj = (RooArgSet*) mjjHig[c]->getParameters(*mjj);
      TIterator* iterMjj = paramsMjj->createIterator();
      TObject* tempObjMjj=nullptr;

      if (_verbLvl>=1 && _verbLvl<4) {
	while((tempObjMjj=iterMjj->Next()))
	  {
	    RooRealVar* var = (RooRealVar*)tempObjMjj;
	    std::cout << "Variables after fit = " << tempObjMjj->GetName() << " " << var->getVal() << "+/-" << var->getError() << std::endl;
	}
	std::cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() <<std::endl;
      }
      //There are very few events in some fits, so adjust the max by a good amount so the MASS-125.0 shift doesn't touch it.
      ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setMax( ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getMax()+(mass-125.0) );
      double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal()+(mass-125.0); // shift the peak
      ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(mPeak); // shift the peak
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "mPeak = " << mPeak <<std::endl;
      if (_verbLvl>=1 && _verbLvl<4) std::cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() <<std::endl;
      // IMPORTANT: fix all pdf parameters to constant
      RooArgSet sigParams( *_w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_sigma_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_alpha_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_n_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_gsigma_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mgg_hig_frac_%d_cat%d",higgschannel,c)) );
      if(higgschannel == 1 || higgschannel == 3){
	sigParams.add(RooArgSet(
				*_w->var(TString::Format("mjj_hig_m0_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_sigma_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_alpha_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_n_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_gsigma_%d_cat%d",higgschannel,c)),
				*_w->var(TString::Format("mjj_hig_frac_%d_cat%d",higgschannel,c)) ) );
      }
      _w->defineSet(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c), sigParams);
      SetConstantParams(_w->set(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c)));
      if (_verbLvl>=1 && _verbLvl<4) std::cout<<std::endl;
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
      mggGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggGaussSig_cat%d",c));
      mggCBSig[c]    = (RooAbsPdf*) _w->pdf(TString::Format("mggCBSig_cat%d",c));
      mggSig[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mggSig_cat%d",c));
      mggBkg[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mggBkg_cat%d",c));
      mjjGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjGaussSig_cat%d",c));
      mjjCBSig[c]    = (RooAbsPdf*) _w->pdf(TString::Format("mjjCBSig_cat%d",c));
      mjjSig[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mjjSig_cat%d",c));
      mjjBkg[c]      = (RooAbsPdf*) _w->pdf(TString::Format("mjjBkg_cat%d",c));

      std::vector<float> effSigmaVecMgg = EffectiveSigma( _w->var("mgg"), mggSig[c], _minSigFitMgg, _maxSigFitMgg);
      sigma_mgg.push_back(effSigmaVecMgg[0]);

      double mgg_sigmaSTD = (mggSig[c]->sigma(*_w->var("mgg")))->getVal();
      sigma_mgg_std.push_back(mgg_sigmaSTD);

      double mjj_sigmaSTD = (mggSig[c]->sigma(*_w->var("mjj")))->getVal();
      sigma_mjj_std.push_back(mjj_sigmaSTD);

      std::vector<float> effSigmaVecMjj = EffectiveSigma( _w->var("mjj"), mjjSig[c], _minSigFitMjj, _maxSigFitMjj);
      sigma_mjj.push_back(effSigmaVecMjj[0]);

      double mgg_mean = (mggSig[c]->mean(*_w->var("mgg")))->getVal();
      mean_mgg.push_back(mgg_mean);

      double mjj_mean = (mjjSig[c]->mean(*_w->var("mjj")))->getVal();
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
  Float_t minSigPlotMgg(115),maxSigPlotMgg(135);
  Float_t minSigPlotMjj(80),maxSigPlotMjj(200);
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

  if (_verbLvl>=1 && _verbLvl<4) std::cout << "[MakePlots] Doing now sig Mgg  and Mtot plots" << std::endl;
  for (int c = 0; c < _NCAT; ++c)
    {
      if (_nonResWeightIndex>=0){
	plotmtot[c] = mtot->frame(Range("SigPlotRange"),Bins(30));
	sigToFit[c]->plotOn(plotmtot[c]);
	plotmtot[c]->Draw();
      	_c1->SaveAs(TString::Format("%s/sigMtot_cat%d.png",_folder_name.data(),c),"QUIET");
      }

      plotmgg[c] = mgg->frame(Range("SigPlotRange"),Bins(nBinsMass));
      sigToFit[c]->plotOn(plotmgg[c]);
      mggSig[c] ->plotOn(plotmgg[c]);
      //    double chi2n = plotmgg[c]->chiSquare(0) ;
      //    if (_verbLvl>=1 && _verbLvl<4) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
      mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggGaussSig_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
      mggSig[c] ->plotOn(plotmgg[c],Components(TString::Format("mggCBSig_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
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
      legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian component","L");
      legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
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

      if(_sigMass >=9000) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + "");
      if(_sigMass < 250 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " (Low Mass)");
      if(_sigMass > 200 && _sigMass < 8000 ) tlat0->DrawLatex(0.16, 0.87, catdesc.at(c) + " ");

      tlat0->SetTextFont(43);
      TString str_desc;
      if(trueSigMass < 250 )str_desc=TString::Format(" Nonresonant HH, Node %d", trueSigMass);
      else str_desc=TString::Format(" %s, M_{X} = %d GeV",_signalType.c_str(), trueSigMass);
      if(trueSigMass == 0 )str_desc=TString::Format(" Nonresonant HH, Box Diagram Only");
      if(trueSigMass == 1 )str_desc=TString::Format(" Nonresonant HH, SM");
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
    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[MakePlots] Doing now sig Mjj plot" << std::endl;
    for (int c = 0; c < _NCAT; ++c)
      {
	plotmjj[c] = mjj->frame(Range("SigPlotRange"),Bins(nBinsMass));
	sigToFit[c]->plotOn(plotmjj[c]);
	mjjSig[c] ->plotOn(plotmjj[c]);
	double chi2n = plotmjj[c]->chiSquare(0) ;
	if (_verbLvl>=1 && _verbLvl<4) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
	mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjGaussSig_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
	mjjSig[c] ->plotOn(plotmjj[c],Components(TString::Format("mjjCBSig_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
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
	legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian component","L");
	legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
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
	  else mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_pol0_%d_cat%d",realint,c));
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
	  if (_verbLvl>=1 && _verbLvl<4) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
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
	  if (_verbLvl>=1 && _verbLvl<4) std::cout << "------------------------- Experimental chi2 = " << chi2n <<std::endl;
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
  for (int c = 0; c < _NCAT; ++c)
    {
      SigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("SigPdf_cat%d",c));
      wAll->import(*_w->pdf(TString::Format("SigPdf_cat%d",c)));
    }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  if ( _NCAT > 2 )
    {
      wAll->factory("prod::CMS_hgg_sig_m0_cat2(mgg_sig_m0_cat2, CMS_hgg_sig_m0_absShift)");
      wAll->factory("prod::CMS_hgg_sig_m0_cat3(mgg_sig_m0_cat3, CMS_hgg_sig_m0_absShift)");
    }
  wAll->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_m0_cat0(mjj_sig_m0_cat0, CMS_hbb_sig_m0_absShift)");
  wAll->factory("prod::CMS_hbb_sig_m0_cat1(mjj_sig_m0_cat1, CMS_hbb_sig_m0_absShift)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hbb_sig_m0_cat2(mjj_sig_m0_cat2, CMS_hbb_sig_m0_absShift)");
    wAll->factory("prod::CMS_hbb_sig_m0_cat3(mjj_sig_m0_cat3, CMS_hbb_sig_m0_absShift)");
  }
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  if ( _NCAT > 2 )
    {
      wAll->factory("prod::CMS_hgg_sig_sigma_cat2(mgg_sig_sigma_cat2, CMS_hgg_sig_sigmaScale)");
      wAll->factory("prod::CMS_hgg_sig_sigma_cat3(mgg_sig_sigma_cat3, CMS_hgg_sig_sigmaScale)");
    }
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  if ( _NCAT > 2 )
    {
      wAll->factory("prod::CMS_hgg_sig_gsigma_cat2(mgg_sig_gsigma_cat2, CMS_hgg_sig_sigmaScale)");
      wAll->factory("prod::CMS_hgg_sig_gsigma_cat3(mgg_sig_gsigma_cat3, CMS_hgg_sig_sigmaScale)");
    }
  wAll->factory("CMS_hbb_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat0(mjj_sig_sigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat1(mjj_sig_sigma_cat1, CMS_hbb_sig_sigmaScale)");
  if ( _NCAT > 2 )
    {
      wAll->factory("prod::CMS_hbb_sig_sigma_cat2(mjj_sig_sigma_cat2, CMS_hbb_sig_sigmaScale)");
      wAll->factory("prod::CMS_hbb_sig_sigma_cat3(mjj_sig_sigma_cat3, CMS_hbb_sig_sigmaScale)");
    }
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat0(mjj_sig_gsigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat1(mjj_sig_gsigma_cat1, CMS_hbb_sig_sigmaScale)");
  if ( _NCAT > 2 )
    {
      wAll->factory("prod::CMS_hbb_sig_gsigma_cat2(mjj_sig_gsigma_cat2, CMS_hbb_sig_sigmaScale)");
      wAll->factory("prod::CMS_hbb_sig_gsigma_cat3(mjj_sig_gsigma_cat3, CMS_hbb_sig_sigmaScale)");
    }
  // (4) do reparametrization of signal
  for (int c = 0; c < _NCAT; ++c) {

    if(_fitStrategy == 2) {
      wAll->factory(TString::Format("EDIT::CMS_sig_cat%d(SigPdf_cat%d,",c,c) +
		    TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d,", c,c) +
		    TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d,", c,c) +
		    TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d,", c,c) +
		    TString::Format(" mjj_sig_m0_cat%d=CMS_hbb_sig_m0_cat%d,", c,c) +
		    TString::Format(" mjj_sig_sigma_cat%d=CMS_hbb_sig_sigma_cat%d,", c,c) +
		    TString::Format(" mjj_sig_gsigma_cat%d=CMS_hbb_sig_gsigma_cat%d)", c,c)
		    );
    } else {
      wAll->factory(TString::Format("EDIT::CMS_sig_cat%d(SigPdf_cat%d,",c,c) +
		    TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d,", c,c) +
		    TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d,", c,c) +
		    TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
		    );
    }
  }
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "Write signal workspace in: " << filename << " file" << std::endl;
  return;
} // close make signal WP

void bbgg2DFitter::MakeHigWS(std::string fileHiggsName,int higgschannel)
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
      HigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c));
      wAll->import(*_w->pdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c)));
    }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat0(mgg_hig_m0_%d_cat0, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat1(mgg_hig_m0_%d_cat1, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  if ( _NCAT > 2 )
    {
      wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat2(mgg_hig_m0_%d_cat2, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat3(mgg_hig_m0_%d_cat3, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
    }
  if(higgschannel == 1 || higgschannel == 3)
    {
      wAll->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
      wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat0(mjj_hig_m0_%d_cat0, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat1(mjj_hig_m0_%d_cat1, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
      if ( _NCAT > 2 )
	{
          wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat2(mjj_hig_m0_%d_cat2, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
          wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat3(mjj_hig_m0_%d_cat3, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
	}
    }
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat0(mgg_hig_sigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat1(mgg_hig_sigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  if ( _NCAT > 2 )
    {
      wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat2(mgg_hig_sigma_%d_cat2, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat3(mgg_hig_sigma_%d_cat3, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
    }
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat0(mgg_hig_gsigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat1(mgg_hig_gsigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  if ( _NCAT > 2 )
    {
      wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat2(mgg_hig_gsigma_%d_cat2, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat3(mgg_hig_gsigma_%d_cat3, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
    }
  if(higgschannel == 1 || higgschannel == 3)
    {
      wAll->factory("CMS_hbb_sig_sigmaScale[1,1,1]");
      wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat0(mjj_hig_sigma_%d_cat0, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat1(mjj_hig_sigma_%d_cat1, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
      if ( _NCAT > 2 )
	{
	  wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat2(mjj_hig_sigma_%d_cat2, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
	  wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat3(mjj_hig_sigma_%d_cat3, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        }
      wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat0(mjj_hig_gsigma_%d_cat0, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
      wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat1(mjj_hig_gsigma_%d_cat1, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
      if ( _NCAT > 2 )
	{
	  wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat2(mjj_hig_gsigma_%d_cat2, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
	  wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat3(mjj_hig_gsigma_%d_cat3, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        }
    }
  // (4) do reparametrization of signal
  if(higgschannel == 1 || higgschannel == 3)
    {
      for (int c = 0; c < _NCAT; ++c) wAll->factory(TString::Format("EDIT::CMS_hig_%d_cat%d(HigPdf_%d_cat%d,",higgschannel,c,higgschannel,c) +
						    TString::Format(" mgg_hig_m0_%d_cat%d=CMS_hgg_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mgg_hig_sigma_%d_cat%d=CMS_hgg_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mgg_hig_gsigma_%d_cat%d=CMS_hgg_hig_gsigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mjj_hig_m0_%d_cat%d=CMS_hbb_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mjj_hig_sigma_%d_cat%d=CMS_hbb_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mjj_hig_gsigma_%d_cat%d=CMS_hbb_hig_gsigma_%d_cat%d)",higgschannel, c,higgschannel,c)
						    );
    }
  else
    {
      for (int c = 0; c < _NCAT; ++c) wAll->factory(TString::Format("EDIT::CMS_hig_%d_cat%d(HigPdf_%d_cat%d,",higgschannel,c,higgschannel,c) +
						    TString::Format(" mgg_hig_m0_%d_cat%d=CMS_hgg_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mgg_hig_sigma_%d_cat%d=CMS_hgg_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
						    TString::Format(" mgg_hig_gsigma_%d_cat%d=CMS_hgg_hig_gsigma_%d_cat%d)",higgschannel, c,higgschannel,c)
						    );
    }
  TString filename(wsDir+fileHiggsName+".inputhig.root");
  wAll->writeToFile(filename);
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "Write signal workspace in: " << filename << " file" << std::endl;
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
      data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
      //RooDataHist* dataBinned = data[c]->binnedClone(); // Uncomment if you want to use wights in the limits
      BkgPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("BkgPdf_cat%d",c));
      wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));// Comment if you want to use wights in the limits
      //wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c))); // Uncomment if you want to use wights in the limits
      wAll->import(*_w->pdf(TString::Format("BkgPdf_cat%d",c)));
      wAll->import(*_w->var(TString::Format("BkgPdf_cat%d_norm",c)));

      //    wAll->factory(TString::Format("CMS_bkg_13TeV_cat%d_norm[%g,0.0,100000.0]",c, _w->var(TString::Format("BkgPdf_cat%d_norm",c))->getVal()));

      //    wAll->factory(TString::Format("CMS_CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("mgg_bkg_13TeV_slope1_cat%d",c))->getVal()));
      //    wAll->factory(TString::Format("CMS_mbb_bkg_13TeV_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat%d",c))->getVal()));

      //    wAll->factory(TString::Format("CMS_CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("mgg_bkg_13TeV_slope2_cat%d",c))->getVal()));
      //    wAll->factory(TString::Format("CMS_CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("mgg_bkg_13TeV_slope3_cat%d",c))->getVal()));

      //    wAll->factory(TString::Format("CMS_mbb_bkg_13TeV_slope2_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d",c))->getVal()));
      //    wAll->factory(TString::Format("CMS_mbb_bkg_13TeV_slope3_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d",c))->getVal()));
      /*
	if(_sigMass == 0 || (_sigMass != 0 && c==1))
	{
    	wAll->factory(TString::Format("CMS_bkg_13TeV_cat%d_norm[%g,0.0,100000.0]",c, _w->var(TString::Format("BkgPdf_cat%d_norm",c))->getVal()));
     	wAll->factory(TString::Format("CMS_hgg_bkg_13TeV_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d",c))->getVal()));
    	wAll->factory(TString::Format("CMS_hbb_bkg_13TeV_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat%d",c))->getVal()));
	}
	else if(_sigMass != 0 && c==0)
	{
    	wAll->factory(TString::Format("CMS_bkg_13TeV_cat%d_norm[%g,0.0,100000.0]",c, _w->var(TString::Format("BkgPdf_cat%d_norm",c))->getVal()));
    	wAll->factory(TString::Format("CMS_hgg_bkg_13TeV_slope2_cat%d[%g,-100,100]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d",c))->getVal()));
    	wAll->factory(TString::Format("CMS_hgg_bkg_13TeV_slope3_cat%d[%g,-100,100]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d",c))->getVal()));
    	wAll->factory(TString::Format("CMS_hbb_bkg_13TeV_slope2_cat%d[%g,-100,100]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d",c))->getVal()));
    	wAll->factory(TString::Format("CMS_hbb_bkg_13TeV_slope3_cat%d[%g,-100,100]",c, _w->var(TString::Format("CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d",c))->getVal()));
	}*/
    } // close ncat
  // (2) do reparametrization of background
  //
  //  for (int c = 0; c < _NCAT; ++c)
  //  {
  //	TString modifications = TString::Format("EDIT::CMS_bkg_13TeV_cat%d(BkgPdf_cat%d)",c,c);// +
  //                        TString::Format("BkgPdf_cat%d_norm=CMS_bkg_13TeV_cat%d_norm, ", c, c) +
  //                        TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d=CMS_mgg_bkg_13TeV_slope1_cat%d, ", c, c) +
  //                        TString::Format("mbb_bkg_13TeV_slope1_cat%d=CMS_mbb_bkg_13TeV_slope1_cat%d, ", c, c) +
  //                        TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d=CMS_mgg_bkg_13TeV_slope2_cat%d, ", c, c) +
  //                        TString::Format("mbb_bkg_13TeV_slope2_cat%d=CMS_mbb_bkg_13TeV_slope2_cat%d, ", c, c) +
  //                        TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d=CMS_mgg_bkg_13TeV_slope3_cat%d, ", c, c) +
  //                        TString::Format("mbb_bkg_13TeV_slope3_cat%d=CMS_mbb_bkg_13TeV_slope3_cat%d )", c, c);
  //        if (_verbLvl>=1 && _verbLvl<4) std::cout << "** MODIFICATIONS: " << modifications << std::endl;

  //        wAll->factory( modifications );

  /*
    wAll->factory(
    TString::Format("EDIT::CMS_bkg_13TeV_cat%d(BkgPdf_cat%d,",c,c) +
    TString::Format("BkgPdf_cat%d_norm=CMS_bkg_13TeV_cat%d_norm, ", c, c) +
    TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d=CMS_mgg_bkg_13TeV_slope1_cat%d, ", c, c) +
    TString::Format("mbb_bkg_13TeV_slope1_cat%d=CMS_mbb_bkg_13TeV_slope1_cat%d, ", c, c) +
    TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d=CMS_mgg_bkg_13TeV_slope2_cat%d, ", c, c) +
    TString::Format("mbb_bkg_13TeV_slope2_cat%d=CMS_mbb_bkg_13TeV_slope2_cat%d, ", c, c) +
    TString::Format("CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d=CMS_mgg_bkg_13TeV_slope3_cat%d, ", c, c) +
    TString::Format("mbb_bkg_13TeV_slope3_cat%d=CMS_mbb_bkg_13TeV_slope3_cat%d, ", c, c) );
  */
  /*

    if(_sigMass == 0 || (_sigMass != 0 && c==1))
    {

    wAll->factory(TString::Format("EDIT::CMS_bkg_13TeV_cat%d(BkgPdf_cat%d,",c,c) +
    TString::Format(" BkgPdf_cat%d_norm=CMS_bkg_13TeV_cat%d_norm,", c,c)+
    TString::Format(" CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat%d=CMS_hgg_bkg_13TeV_slope1_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat%d=CMS_hbb_bkg_13TeV_slope1_cat%d)", c,c)  );

    wAll->factory(TString::Format("EDIT::CMS_bkg_13TeV_cat%d(BkgPdf_cat%d,",c,c) +
    TString::Format(" BkgPdf_cat%d_norm=CMS_bkg_13TeV_cat%d_norm,", c,c)+
    TString::Format(" CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d=CMS_hgg_bkg_13TeV_slope2_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d=CMS_hgg_bkg_13TeV_slope3_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d=CMS_hbb_bkg_13TeV_slope2_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d=CMS_hbb_bkg_13TeV_slope3_cat%d)", c,c)  );
    }
    else if(_sigMass != 0 && c==0)
    {
    wAll->factory(TString::Format("EDIT::CMS_bkg_13TeV_cat%d(BkgPdf_cat%d,",c,c) +
    TString::Format(" BkgPdf_cat%d_norm=CMS_bkg_13TeV_cat%d_norm,", c,c)+
    TString::Format(" CMS_hhbbgg_13TeV_mgg_bkg_slope2_cat%d=CMS_hgg_bkg_13TeV_slope2_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mgg_bkg_slope3_cat%d=CMS_hgg_bkg_13TeV_slope3_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mjj_bkg_slope2_cat%d=CMS_hbb_bkg_13TeV_slope2_cat%d,", c,c) +
    TString::Format(" CMS_hhbbgg_13TeV_mjj_bkg_slope3_cat%d=CMS_hbb_bkg_13TeV_slope3_cat%d)", c,c)  );

    } */
  //  } // close for cat
  TString filename(wsDir+fileBaseName+".root");
  wAll->writeToFile(filename);
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "Write background workspace in: " << filename << " file" <<std::endl;
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "observation ";
  for (int c = 0; c < _NCAT; ++c)
    {
      if (_verbLvl>=1 && _verbLvl<4) std::cout << " " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries();
    }
  if (_verbLvl>=1 && _verbLvl<4) std::cout << std::endl;
  return;
} // close make BKG workspace

void bbgg2DFitter::MakeDataCard(std::string fileBaseName, std::string fileBkgName ,Bool_t useSigTheoryUnc)
{
  if (_verbLvl==4){
    std::cout<<" DBG.  Making Data card"<<std::endl;
    std::cout<<" fileBaseName ="<<fileBaseName<<"\n   fileBkgName="<<fileBkgName<<std::endl;
  }


  TString cardDir = TString::Format("%s/datacards/",_folder_name.data());
  TString wsDir   = TString::Format("%s/workspaces/",_folder_name.data());
  std::vector<RooDataSet*>vec(_NCAT,nullptr);
  std::map<std::string,std::vector<RooDataSet*>>higToFits{{"ggh_m125_powheg_13TeV",vec},{"tth_m125_13TeV",vec},{"vbf_m125_13TeV",vec},{"wzh_m125_13TeV_zh",vec},{"bbh_m125_13TeV",vec}};
  std::vector<RooDataSet*> data(_NCAT,nullptr);
  std::vector<RooDataSet*> sigToFit(_NCAT,nullptr);
  for (int c = 0; c < _NCAT; ++c)
    {
      data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
      sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
      //
      if(_addHiggs){
	for(std::vector<std::string>::iterator it=_singleHiggsNames.begin();it!=_singleHiggsNames.end();++it)
	  {
	    higToFits[*it][c]=(RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",_singleHiggsMap[*it],c));
	    if (_verbLvl>=1 && _verbLvl<4) std::cout<<*it<<" found "<<std::endl;
	  }
      }
    } // close cat
  ////////////////////////////////////////////////////////////////////////////////////
  //RooRealVar* lumi = w->var("lumi");
  if (_verbLvl>=1 && _verbLvl<4) {

    std::cout << "======== Expected Events Number =====================" <<std::endl;
    std::cout << ".........Measured Data for L = " << _lumi << " pb-1 ............................" <<std::endl;
    if(!_doblinding){ std::cout << "#Events data: " << _w->data("Data")->sumEntries() << std::endl;}
    else std::cout << "#Events data: -1 " << std::endl;
    if(!_doblinding)
      for (int c = 0; c < _NCAT; ++c) std::cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() <<std::endl;
    else
      for (int c = 0; c < _NCAT; ++c) std::cout << TString::Format("#Events data cat%d: ",c) << -1 << std::endl;
    std::cout << ".........Expected Signal for L = " << _lumi << " pb-1 ............................" << std::endl;
    if(!_doblinding) std::cout << "#Events Signal: " << _w->data("Data")->sumEntries() << std::endl;
    else  std::cout << "#Events Signal: -1 " << std::endl;
    //std::vector<Float_t>siglikeErr(_NCAT,0.0);

    for (int c = 0; c < _NCAT; ++c) {
      std::cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() <<std::endl;
      //siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
    }
    std::cout << "====================================================" <<std::endl;
  }

  TString filename(cardDir+fileBaseName+".txt");
  std::ofstream outFile(filename);
  // outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << std::endl;
  outFile << "#Run with: combine -d hgg.mH350.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3500 -i 50000 --optimizeSim=1 --tries 30" << std::endl;
  outFile << "# Lumi = " << _lumi << " pb-1" << std::endl;
  outFile << "imax "<<_NCAT << std::endl;
  outFile << "jmax "<<_singleHiggsNames.size()+1<< std::endl; // number of BKG
  outFile << "kmax *" << std::endl;
  outFile << "---------------" << std::endl;
  outFile << "shapes data_obs cat0 " << wsDir+fileBkgName+".root" << " w_all:data_obs_cat0" << std::endl;
  if ( _NCAT > 1 ) outFile << "shapes data_obs cat1 "<< wsDir+fileBkgName+".root" << " w_all:data_obs_cat1" << std::endl;
  if ( _NCAT > 2 )
    {
      outFile << "shapes data_obs cat2 " << wsDir+fileBkgName+".root" << " w_all:data_obs_cat2" << std::endl;
      outFile << "shapes data_obs cat3 "<< wsDir+fileBkgName+".root" << " w_all:data_obs_cat3" << std::endl;
    }
  outFile << "############## shape with reparametrization" << std::endl;
  outFile << "shapes Bkg cat0 " << wsDir+fileBkgName+".root" << " w_all:CMS_bkg_13TeV_cat0" << std::endl;
  if ( _NCAT > 1 ) outFile << "shapes Bkg cat1 "<< wsDir+fileBkgName+".root" << " w_all:CMS_bkg_13TeV_cat1" << std::endl;
  if ( _NCAT > 2 )
    {
      outFile << "shapes Bkg cat2 " << wsDir+fileBkgName+".root" << " w_all:CMS_bkg_13TeV_cat2" << std::endl;
      outFile << "shapes Bkg cat3 "<< wsDir+fileBkgName+".root" << " w_all:CMS_bkg_13TeV_cat3" << std::endl;
    }
  outFile << "# signal" << std::endl;
  outFile << "shapes Sig cat0 " << wsDir+fileBaseName+".inputsig.root" << " w_all:CMS_sig_cat0" << std::endl;
  if ( _NCAT > 1 ) outFile << "shapes Sig cat1 " << wsDir+fileBaseName+".inputsig.root" << " w_all:CMS_sig_cat1" << std::endl;
  if ( _NCAT > 2 )
    {
      outFile << "shapes Sig cat2 " << wsDir+fileBaseName+".inputsig.root" << " w_all:CMS_sig_cat2" << std::endl;
      outFile << "shapes Sig cat3 " << wsDir+fileBaseName+".inputsig.root" << " w_all:CMS_sig_cat3" << std::endl;
    }
  for(std::map<std::string,std::vector<RooDataSet*>>::iterator it=higToFits.begin();it!=higToFits.end();++it)
    {
      std::map<std::string,std::string>name{{"ggh_m125_powheg_13TeV","#ggh"},{"tth_m125_13TeV","#tth"},{"vbf_m125_13TeV","#vbf"},{"wzh_m125_13TeV_zh","#vh"},{"bbh_m125_13TeV","#bbh"}};
      std::map<std::string,std::string>name2{{"ggh_m125_powheg_13TeV","Higggh"},{"tth_m125_13TeV","Higtth"},{"vbf_m125_13TeV","Higvbf"},{"wzh_m125_13TeV_zh","Higvh"},
																       {"bbh_m125_13TeV","Higbbh"}};
      std::vector<std::string>::iterator itt=find(_singleHiggsNames.begin(),_singleHiggsNames.end(),it->first);
      if(itt!=_singleHiggsNames.end())
  	{outFile << name[it->first] << std::endl;
	  for(int hh=0;hh<_NCAT;++hh)
	    {
	      std::string catn="cat"+std::to_string(hh);
	      outFile << "shapes "<<name2[it->first]<<" "<<catn<<" "<< wsDir+_singleHiggsWSfileNames[it->first]<<".inputhig.root"<<" w_all:CMS_hig_"<<_singleHiggsMap[it->first]<<"_"<<catn << std::endl;
	    }
	}
    }
  outFile << "---------------" <<std::endl;
  /////////////////////////////////////
  if(1)
    { //
      outFile << "bin cat0 ";
      if ( _NCAT > 1 ) outFile << "cat1 ";
      if ( _NCAT > 2 ) outFile << "cat2 cat3 ";
      if(!_doblinding)
	{
	  outFile << "\nobservation "<< data[0]->sumEntries() <<" " ;
	  //this->SetObservedCats(0, data[0]->sumEntries());
	  if ( _NCAT > 1 ) {
	    outFile << data[1]->sumEntries() <<" ";
	    //this->SetObservedCats(1, data[1]->sumEntries());
	  }
	  if ( _NCAT > 2 ) {
	    outFile << data[2]->sumEntries() <<" " << data[3]->sumEntries() <<" ";
	    //this->SetObservedCats(2, data[2]->sumEntries());
	  }
	}
      else
	{
	  outFile << "\nobservation -1 ";
	  if ( _NCAT > 1 ) outFile << " -1 ";
	  if ( _NCAT > 2 ) outFile << "-1 -1 ";
	}
      outFile << "\n------------------------------" << std::endl;
      std::map<std::string,std::string>Name{{"ggh_m125_powheg_13TeV","Higggh"},{"tth_m125_13TeV","Higtth"},{"vbf_m125_13TeV","Higvbf"},{"wzh_m125_13TeV_zh","Higvh"},
																      {"bbh_m125_13TeV","Higbbh"}};
      outFile << "bin ";
      for(int f=0;f<_NCAT;++f)
	{
	  std::string catn="cat"+std::to_string(f)+" ";

	  for(unsigned int l=0;l!=_singleHiggsNames.size()+2;++l) outFile<<catn;
	}
      outFile << "\nprocess";
      for(int f=0;f<_NCAT;++f)
	{
	  outFile<<" Sig Bkg ";
	  for(unsigned int l=0;l!=_singleHiggsNames.size();++l) outFile<<Name[_singleHiggsNames[l]]<<" ";
	}
      outFile<< "\nprocess ";
      for(int f=0;f<_NCAT;++f)
	{
	  outFile<< "0 1 ";
	  for(unsigned int l=0;l!=_singleHiggsNames.size();++l) outFile<<std::to_string(l+2)+" ";
	}
      outFile << "\nrate ";
      for(int f=0;f<_NCAT;++f)
	{
	  outFile<<" "<<sigToFit[f]->sumEntries()<<" "<<1<<" ";
	  //this->SetSigExpectedCats(f, sigToFit[f]->sumEntries());
	  for(std::map<std::string,std::vector<RooDataSet*>>::iterator it=higToFits.begin();it!=higToFits.end();++it)
	    {
	      std::vector<std::string>::iterator itt=find(_singleHiggsNames.begin(),_singleHiggsNames.end(),it->first);
	      if(itt!=_singleHiggsNames.end())outFile<<(it->second)[f]->sumEntries()<<" ";
	    }

	}
      outFile << " " << std::endl;
      outFile << "############## Total normalisation" <<std::endl;
      std::vector<std::string>lumi{"1.026","-","1.026","1.026","1.026","1.026","1.026"};
      std::vector<std::string>DiphoTrigger{"1.010","-","1.010","1.010","1.010","1.010","1.010"};
      std::vector<std::string>CMS_hgg_eff_g{"1.010","-","1.010","1.010","1.010","1.010","1.010"};
      std::vector<std::string>pTj_cut_acceptance{"1.010","-","1.010","1.010","1.010","1.010","1.010"};
      std::vector<std::string>btag_eff{"1.050","-","1.0508","1.050","1.050","1.050","1.050"};
      std::vector<std::string>btag_eff2{"0.979","-","0.979","0.979","0.979","0.979","0.979"};
      std::vector<std::string>maajj_cut_acceptance{"1.015","-","1.015","1.015","1.015","1.015","1.015"};
      std::vector<std::string>maajj_cut_acceptance2{"0.995","-","0.995","0.995","0.995","0.995","0.995"};
      std::vector<std::string>PDF{"-","-","0.931/1.075","0.919/1.081","0.972/1.026","0.976/1.024","0.976/1.024"};
      std::vector<std::string>QCD_scale{"-","-","0.922/1.072","0.907/1.038","0.998/1.002","0.980/1.020","0.980/1.020"};
      std::vector<std::string>gg_migration{"-","-","1.25","1.25","1.08","1.08","1.08"};
      std::vector<std::string>gluonSplitting{"-","-","1.40","1.40","1.40","1.40","1.40"};
      std::vector<std::pair<std::string,std::vector<std::vector<std::string>>>>MAP
	{
	  {"lumi_"+_energy,{lumi}},
	    {"DiphoTrigger",{DiphoTrigger}},
	      {"CMS_hgg_eff_g",{CMS_hgg_eff_g}},
		{"pTj_cut_acceptance",{pTj_cut_acceptance}},
		  {"btag_eff",{btag_eff,btag_eff2}},
		    {"maajj_cut_acceptance",{maajj_cut_acceptance,maajj_cut_acceptance2}},
		      {"PDF",{PDF}},
			{"QCD_scale",{QCD_scale}},
			  {"gg_migration",{gg_migration}},
			    {"gluonSplitting",{gluonSplitting}}
	};
      std::map<std::string,std::string>Comments
      {
	{"DiphoTrigger","############## Photon selection normalisation uncertainties "},
	  {"CMS_hgg_eff_g","# Trigger efficiency"},
	    {"pTj_cut_acceptance","# photon selection accep.\n\n############## Jet selection and phase space cuts normalisation uncertainties "},
	      {"btag_eff","# JER and JES "},
		{"maajj_cut_acceptance","# b tag efficiency uncertainty"},
		  {"PDF","# uncertainty on mggjj cut acceptance \n\n ############## Theory uncertainties on SM Higgs production "},
		    };


      for(std::vector<std::pair<std::string,std::vector<std::vector<std::string>>>>::iterator it=MAP.begin();it!=MAP.end();++it)
	{
	  if(Comments.find(it->first)!=Comments.end())outFile<<Comments[it->first]<<std::endl;
	  outFile<<it->first<<" lnN ";
	  bool manycase=false;
	  if((it->second).size()==2)manycase=true;
	  for(int f=0;f<_NCAT;++f)
	    {
	      if(manycase==false)outFile<<(it->second)[0][0]<<"  ";
	      else outFile<<(it->second)[f%2][0]<<"  ";
	      if(manycase==false)outFile<<(it->second)[0][1]<<"  ";
	      else outFile<<(it->second)[f%2][1]<<"  ";
	      for(unsigned int l=0;l!=_singleHiggsNames.size();++l)
		{
		  if(manycase==false)outFile<<(it->second)[0][_singleHiggsMap[_singleHiggsNames[l]]+2]<<"  ";
		  else outFile<<(it->second)[f%2][_singleHiggsMap[_singleHiggsNames[l]]+2]<<"  ";
		}
	    }
	  outFile<<std::endl;
	}
      /*outFile << "lumi_13TeV lnN "<< "1.026 - 1.026 1.026 1.026 1.026 1.026 "<< "1.026 - 1.026 1.026 1.026 1.026 1.026 ";
	if ( _NCAT > 2 )outFile << "1.026 - 1.026 1.026 1.026 1.026 1.026 "<< "1.026 - 1.026 1.026 1.026 1.026 1.026 ";
	outFile << " " << std::endl << std::endl;
	outFile << "############## Photon selection normalisation uncertainties " << std::endl;
	outFile << "DiphoTrigger lnN "<< "1.01 - 1.010 1.010 1.010 1.010 1.010 "<< "1.01 - 1.010 1.010 1.010 1.010 1.010 ";
	if ( _NCAT > 2 )outFile << "1.01 - 1.010 1.010 1.010 1.010 1.010 "<< "1.01 - 1.010 1.010 1.010 1.010 1.010 ";
	outFile << "# Trigger efficiency" << std::endl;
	outFile << "CMS_hgg_eff_g lnN "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
	if ( _NCAT > 2 )outFile << "1.010 - 1.010 1.010 1.010 1.010 1.010 "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
	outFile << "# photon selection accep." << std::endl;
	outFile << " " << std::endl;
	outFile << "############## Jet selection and phase space cuts normalisation uncertainties " <<std::endl;
	outFile << "pTj_cut_acceptance lnN "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
	if ( _NCAT > 2 )outFile << "1.010 - 1.010 1.010 1.010 1.010 1.010 "<< "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
	outFile <<"# JER and JES " << std::endl;
	if ( _NCAT == 2)outFile << "btag_eff lnN "<< "1.050 - 1.0508 1.050 1.050 1.050 1.050 "<< "0.979 - 0.979 0.979 0.979 0.979 0.979 ";
	else if ( _NCAT > 2 )outFile << "btag_eff lnN "<< "1.050 - 1.050 1.050 1.050 1.050 1.050 "<< "0.979 - 0.979 0.979 0.979 0.979 0.979 "<< "1.050 - 1.050 1.050 1.050 1.050 1.050 "<< "0.979 - 0.979 0.979 0.979 0.979 0.979 ";
	outFile <<"# b tag efficiency uncertainty" << std::endl;
	if (_NCAT == 2)outFile << "maajj_cut_acceptance lnN "<< "1.015 - 1.015 1.015 1.015 1.015 1.015 "<< "1.015 - 1.015 1.015 1.015 1.015 1.015 ";
	else if (_NCAT > 2)outFile << "maajj_cut_acceptance lnN "<< "0.995 - 0.995 0.995 0.995 0.995 0.995 "<< "0.995 - 0.995 0.995 0.995 0.995 0.995 "<< "1.015 - 1.015 1.015 1.015 1.015 1.015 "<< "1.015 - 1.015 1.015 1.015 1.015 1.015 ";
	outFile << "# uncertainty on mggjj cut acceptance" <<std::endl;
	outFile << " " <<std::endl <<std::endl;
	outFile << "############## Theory uncertainties on SM Higgs production " << std::endl;
	outFile << "PDF lnN "<< " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 "<< " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 ";
	if ( _NCAT > 2 )outFile << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 "<< " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 ";
	outFile << "\nQCD_scale lnN "<< " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 "<< " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 ";
	if ( _NCAT > 2 )outFile << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 "<< " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 ";
	outFile << "\ngg_migration lnN "<< " - - 1.25 1.25 1.08 1.08 1.08 "<< " - - 1.25 1.25 1.08 1.08 1.08 ";
	if ( _NCAT > 2 )outFile << " - - 1.25 1.25 1.08 1.08 1.08 "<< " - - 1.25 1.25 1.08 1.08 1.08 ";
	outFile << "# UEPS" << std::endl;
	outFile << "gluonSplitting lnN "<< " - - 1.40 1.40 1.40 1.40 1.40 "<< " - - 1.40 1.40 1.40 1.40 1.40 ";
	if ( _NCAT > 2 )outFile << " - - 1.40 1.40 1.40 1.40 1.40 "<< " - - 1.40 1.40 1.40 1.40 1.40 ";
	outFile << " " << std::endl<<endl;*/
      if(useSigTheoryUnc)
	{
	  outFile << "############## Theory uncertainty on SM diHiggs production " << std::endl;
	  outFile << "SM_diHiggs_Theory lnN "<< " 0.857/1.136 - - - - - - ";
	  if ( _NCAT > 1 ) outFile << " 0.857/1.136 - - - - - - ";
	  if ( _NCAT > 2 ) outFile << " 0.857/1.136 - - - - - - "<< " 0.857/1.136 - - - - - - ";
	  outFile << " # from 9.96 + 1.35 - 1.42 fb " << std::endl << std::endl;
	}
      outFile << "############## Signal parametric shape uncertainties " << std::endl;
      if ( _sigMass > 0)outFile << "CMS_hgg_sig_m0_absShift param 1 0.0045 # displacement of the dipho mean error = sqrt(0.4^ 2 + 0.2^ 2) " << std::endl;
      else outFile << "CMS_hgg_sig_m0_absShift param 1 0.0054 # displacement of the dipho mean error = sqrt(0.5^ 2 + 0.2^ 2) " << std::endl;
      outFile << "CMS_hgg_sig_sigmaScale param 1 0.05 # optimistic estimate of resolution uncertainty " << std::endl;
      //
      outFile << "CMS_hbb_sig_m0_absShift param 1 0.026 # displacement of the dijet mean error " << std::endl;
      outFile << "CMS_hbb_sig_sigmaScale param 1 0.10 # optimistic estimate of resolution uncertainty " << std::endl;
      //
      outFile << "############## for mggxmjj fit - slopes" << std::endl;
      outFile << "CMS_bkg_13TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << std::endl;
      if ( _NCAT > 1 ) outFile << "CMS_bkg_13TeV_cat1_norm flatParam # Normalization uncertainty on background slope" << std::endl;
      if ( _NCAT > 2 )
	{
	  outFile << "CMS_bkg_13TeV_cat2_norm flatParam # Normalization uncertainty on background slope" << std::endl;
	  outFile << "CMS_bkg_13TeV_cat3_norm flatParam # Normalization uncertainty on background slope" << std::endl;
	  outFile << "CMS_hgg_bkg_13TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hgg_bkg_13TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hgg_bkg_13TeV_slope1_cat2 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hgg_bkg_13TeV_slope1_cat3 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope1_cat2 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope1_cat3 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	}
      else
	{
	  outFile << "CMS_hgg_bkg_13TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hgg_bkg_13TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  if ( _NCAT > 1 ) outFile << "CMS_hgg_bkg_13TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  outFile << "CMS_hbb_bkg_13TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	  if ( _NCAT > 1 ) outFile << "CMS_hbb_bkg_13TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << std::endl;
	}
    } // if ncat == 2 or 4
  /////////////////////////////////////
  outFile.close();
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "Write data card in: " << filename << " file" << std::endl;
} // close write full datacard

void bbgg2DFitter::SetConstantParams(const RooArgSet* params)
{
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next())
    {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      if (rrv) rrv->setConstant(true); if (_verbLvl>=1 && _verbLvl<4) std::cout << " " << rrv->GetName();
    }
} // close set const parameters

TStyle * bbgg2DFitter::style()
{

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

  hggStyle->SetCanvasDefH(600); //Height of canvas
  hggStyle->SetCanvasDefW(600); //Width of canvas
  hggStyle->SetCanvasDefX(0);   //POsition on screen
  hggStyle->SetCanvasDefY(0);

  hggStyle->SetPadLeftMargin(0.13);//0.16);
  hggStyle->SetPadRightMargin(0.1);//0.02);
  hggStyle->SetPadTopMargin(0.085);//0.02);
  hggStyle->SetPadBottomMargin(0.12);//0.02);

  // For hgg axis titles:
  hggStyle->SetTitleColor(1, "XYZ");
  hggStyle->SetTitleFont(42, "XYZ");
  hggStyle->SetTitleSize(0.04, "XYZ");
  hggStyle->SetTitleXOffset(1.2);//0.9);
  hggStyle->SetTitleYOffset(1.75); // => 1.15 if exponents

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
  if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Starting cat loop " << std::endl;
  for (int c = 0; c < ncat; ++c) { // to each category
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));

    TH2* data_h2 = 0;
    TH1* data_h11 = 0;
    if(_fitStrategy==2)  data_h2= (TH2*) data[c]->createHistogram("mgg,mjj", 60, 40);
    if(_fitStrategy==1)  data_h11= (TH1*) data[c]->createHistogram("mgg", 60);

    if (_verbLvl>=1 && _verbLvl<4) {
      if(_doblinding==0 && _fitStrategy==2) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h2->Integral() << std::endl;
      if(_doblinding==0 && _fitStrategy==1) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h11->Integral() << std::endl;
    }
      int nEvtsObs = -1;
    if(_fitStrategy == 2) nEvtsObs = data_h2->Integral();
    if(_fitStrategy == 1) nEvtsObs = data_h11->Integral();

    //data_h11->Delete();

    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop 1 - cat" << c << std::endl;

    ////////////////////////////////////
    // these are the parameters for the bkg polinomial
    // one slope by category - float from -10 > 10
    // we first wrap the normalization of mggBkgTmp0, mjjBkgTmp0
    // CMS_hhbbgg_13TeV_mgg_bkg_slope1
    _w->factory(TString::Format("BkgPdf_cat%d_norm[1.0,0.0,100000]",c));
    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop 2 - cat" << c << std::endl;
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

    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop 3 - cat" << c << std::endl;

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


    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop 4 - cat" << c << std::endl;

    mggBkgTmpBer1 = new RooBernstein(TString::Format("mggBkgTmpBer1_cat%d",c),"",*mgg,RooArgList(*mgg_p0amp,*mgg_p1amp));
    mjjBkgTmpBer1 = new RooBernstein(TString::Format("mjjBkgTmpBer1_cat%d",c),"",*mjj,RooArgList(*mjj_p0amp,*mjj_p1amp));

    if(nEvtsObs > 10) {
      mggBkgTmpBer1 = new RooBernstein(TString::Format("mggBkgTmpBer1_cat%d",c),"",*mgg,RooArgList(*mgg_p0amp,*mgg_p1amp, *mgg_p2amp));
      mjjBkgTmpBer1 = new RooBernstein(TString::Format("mjjBkgTmpBer1_cat%d",c),"",*mjj,RooArgList(*mjj_p0amp,*mjj_p1amp, *mjj_p2amp));
    }

    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop 5 - cat" << c << std::endl;

    if(_fitStrategy==2) BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpBer1, *mjjBkgTmpBer1));
    //    if(_fitStrategy==1) BkgPdf1 = (RooAbsPdf*) mggBkgTmpBer1->Clone(TString::Format("BkgPdf_cat%d",c));

    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Cat loop " << c << std::endl;

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

    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Done with fit " << c << std::endl;

    ///////////////////////////////////
    //Calculate 2D chisquare by hand //
    ///////////////////////////////////
    //    TH2* pdf_h2 = BkgPdfExt.createHistogram("mgg vs mjj pdf", mgg, Binning(60), YVar(mjj, Binning(40)));
    //    TH2* data_h2 = (TH2*) data[c]->createHistogram("mgg,mjj", 60, 40);
    //    if (_verbLvl>=1 && _verbLvl<4) std::cout << "########NUMBER OF OBSERVED EVENTSSSS::: " << data_h2->Integral() << std::endl;
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
	if (_verbLvl>=1 && _verbLvl<4) std::cout << "Number of bins for 2D chi square are different!!!!" << std::endl;
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
	    //if (_verbLvl>=1 && _verbLvl<4) std::cout << "bin x: " << nbx << " bin y: " << nby << " exp: " << expected_h2 << " obs: " << observed_h2 << std::endl;
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
      if (_verbLvl>=1 && _verbLvl<4) {
	std::cout << "########NUMBER OF PDFFF EVENTSSSS::: " << pdf_h2->Integral() << std::endl;
	std::cout << "################################################################" << std::endl;
	std::cout << "################################################################" << std::endl;
	std::cout << "################ 1DKSTEST:" << data_h1->KolmogorovTest(pdf_h1, "D") << " ##################################" << std::endl;
	std::cout << "################ 2DKSTEST:" << data_h2->KolmogorovTest(pdf_h2, "D") << " ##################################" << std::endl;
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
    if (_verbLvl>=1 && _verbLvl<4) std::cout << "[BkgModelFit] Plotting Mgg - cat" << c << std::endl;
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
    if(_sigMass > -1)
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
    if(_sigMass==0)
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
      if (_verbLvl>=1 && _verbLvl<4 ) std::cout << "[BkgModelFit] Plotting Mgg - cat" << c << std::endl;
      //************************************************//
      // Plot mjj background fit results per categories
      //************************************************//
      //ctmp = new TCanvas(TString::Format("ctmpBkgMjj_cat%d",c),"mjj Background Categories",0,0,500,500);
      nBinsMass = 60;
      plotmjjBkg[c] = mjj->frame(nBinsMass);
      dataplot[c] = (RooDataSet*) _w->data(TString::Format("Dataplot_cat%d",c));
      if(_doblinding) dataplot[c]->plotOn(plotmjjBkg[c],Invisible());
      else dataplot[c]->plotOn(plotmjjBkg[c]);
      if(_sigMass > -1)
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
      //plotmjjBkg[c]->getObject(1)->Draw("SAME");
      //plotmjjBkg[c]->getObject(2)->Draw("P SAME");
      ////////////////////////////////////////////////////////// plot higgs
      if(addhiggs){
	for(unsigned int d=0;d!=_singleHiggsNames.size();++d)
	  {
	    static std::vector<int>color{2,3,6,7,4};
	    int realint=_singleHiggsMap[_singleHiggsNames[d]];
	    sigToFitvec[realint][c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",realint,c));
	    double norm = 1.0*sigToFitvec[realint][c]->sumEntries(); //
	    //norm0 = 0.0000001;
	    mjjSigvec[realint][c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",realint,c));
	    // we are not constructing signal pdf, this is constructed on sig to fit function...
	    mjjSigvec[realint][c] ->plotOn(plotmjjBkg[c],Normalization(norm,RooAbsPdf::NumEvent),DrawOption("F"),LineColor(color[realint]),FillStyle(1001),FillColor(19));
	    mjjSigvec[realint][c]->plotOn(plotmjjBkg[c],Normalization(norm,RooAbsPdf::NumEvent),LineColor(color[realint]),LineStyle(1));
	    //
	  }
      }
      //////////////////////////////////////////////////////////
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

      if(_sigMass==0)legmc->SetHeader(" Nonresonant HH");
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
