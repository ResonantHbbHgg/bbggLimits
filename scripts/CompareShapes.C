void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000); 
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.15);
  defaultStyle->SetPadRightMargin(0.04);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.10,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  // defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

  defaultStyle->SetTitleColor(1, "XYZ");
  defaultStyle->SetTitleFont(42, "XYZ");
  defaultStyle->SetTitleSize(0.06, "XYZ");
  defaultStyle->SetTitleXOffset(0.9);
  defaultStyle->SetTitleYOffset(1.05);
  
  // For the axis labels:
  defaultStyle->SetLabelColor(1, "XYZ");
  defaultStyle->SetLabelFont(42, "XYZ");
  defaultStyle->SetLabelOffset(0.007, "XYZ");
  defaultStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(310, "XYZ");
    defaultStyle->SetPadTickX(1);
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
  return;
}



std::vector<float> EffectiveSigma(RooRealVar* mass, RooAbsPdf* binned_pdf, float wmin=115, float wmax=135, float step=0.002, float epsilon=1.e-4)
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

 
  std::cout << "#Sigma effective: xLow: " << low << ", xHigh: " << high << ", width/2: " << width/2 << std::endl;
    
  std::vector<float> outVec;
  outVec.push_back(width);
  outVec.push_back(low);
  outVec.push_back(high);
  outVec.push_back(width/2.);

  return outVec;
}





void CompareShapes(){

  using namespace RooFit;
  using namespace RooStats ;


 style();

 TCanvas* ctmp = new TCanvas("ctmp","Shapes Compare",0,0,600,600);



   RooPlot *plot;
   RooAbsPdf *pdf_2015, *pdf_2016;
   RooRealVar *mass; 
   string datacard;

   TLegend *legmc = new TLegend(0.60,0.55,0.99,0.80);
   legmc->SetHeader("Comparison for HM HP cat");

   TFile* f_hgg_2015 = TFile::Open("hh_ggbb_SM_13TeV.inputsig.root");  
   TFile* f_hgg_2016 = TFile::Open(Form("/afs/cern.ch/work/m/mgouzevi/private/LIMITS/CLEAN_EPS_GGBB_UPDATE/CMSSW_7_4_7/src/HiggsAnalysis/bbggLimits/outDir_m5/HighMass_ARW_/workspaces/hhbbgg.mH125_13TeV.inputsig.root"));

   
   RooWorkspace *w_hgg_2015 = (RooWorkspace*)f_hgg_2015->Get("w_all");
   RooWorkspace *w_hgg_2016 = (RooWorkspace*)f_hgg_2016->Get("w_all");

   // Make a plot (data is a toy dataset)
   mass = w_hgg_2015->var("mgg");
   
   mass->setRange("fitrange",115,135);
   plot = mass->frame(Range("fitrange")); //  data->plotOn(plot);
    

   
   pdf_2015 = w_hgg_2015->pdf("mggSig_cat0_CMS_sig_cat0");
   
   std::vector<float> effsig_2015 = EffectiveSigma(mass, pdf_2015);

   pdf_2015->plotOn(plot,RooFit::LineColor(kRed), RooFit::LineStyle(1), Range("fitrange"),NormRange("fitrange"), Normalization(1,RooAbsReal::NumEvent));
   
   legmc->AddEntry(plot->getObject(0),"HIG-16-032 - 2015","L");
   legmc->AddEntry(plot->getObject(0),Form("#sigma_{eff} = %.2f", effsig_2015[3]),"L");

   mass = w_hgg_2016->var("mgg");
   pdf_2016 = w_hgg_2016->pdf("mggSig_cat0_CMS_sig_cat2");

   std::vector<float> effsig_2016 = EffectiveSigma(mass, pdf_2016);


   pdf_2016->plotOn(plot,RooFit::LineColor(kBlue), RooFit::LineStyle(1), Range("fitrange"),NormRange("fitrange"), Normalization(1,RooAbsReal::NumEvent));

   legmc->AddEntry(plot->getObject(1),"HIG-17-008 - 2016","L");
   legmc->AddEntry(plot->getObject(1),Form("#sigma_{eff} = %.2f", effsig_2016[3]),"L");
   

   plot->SetTitle("; M_{#gamma#gamma} [GeV]; Fraction/GeV");
   plot->Draw();


   legmc->SetBorderSize(0);
   legmc->SetFillStyle(0);
   legmc->SetTextFont(62);
   legmc->SetTextSize(0.03);
    
    
   legmc->Draw();
   ctmp->SaveAs("2015-2016-mgg-Comparison.png");

}
