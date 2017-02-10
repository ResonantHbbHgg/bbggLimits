using namespace RooFit ;

void Correlations(){

  // tool to look how well background model or signal model describe signal/background in 2D
  
  int icat = 1; // category you want to explore

  bool bSig = 1; // Signal or background

  string label("R300");

  TFile *f_hgg = new TFile();
  if (bSig) f_hgg = TFile::Open("hhbbgg.mH125_13TeV.inputsig.root");
  else f_hgg = TFile::Open("hhbbgg.inputbkg_13TeV.root");


  RooWorkspace *w_hgg = (RooWorkspace*)f_hgg->Get("w_all;1");

  RooRealVar *mgg =  w_hgg->var("mgg");
  RooRealVar *mjj =  w_hgg->var("mjj");
  if(bSig)   {mgg->setMin(120);    mgg->setMax(130);}
  mjj->setMax(170);     mjj->setMin(80);
  

  RooDataSet *signal;
  if (bSig) signal = (RooDataSet*) w_hgg->data(Form("Sig_cat%d",icat));
  else signal = (RooDataSet*) w_hgg->data(Form("data_obs_cat%d",icat));
   

  cout << "num entries = " << signal->numEntries() << endl;

  int numEntries =  signal->numEntries(); // Used to label the final file
   
  // Build histogram


  TH2F* hh_signal = new TH2F();
  if (bSig) hh_signal = signal->createHistogram(*mgg, *mjj, 20, 10);
  else hh_signal = signal->createHistogram(*mgg, *mjj, 10, 10);

   hh_signal->SetStats(0);

   double NormSig = hh_signal->Integral();

   cout << "NormSig = " << NormSig << endl;
   //numEntries =  50;//signal->numEntries(); // Possibility to change signal normalisation

   hh_signal->Scale(numEntries/NormSig);
   NormSig = hh_signal->Integral(); // Variable used to normalise the signal


   
   TCanvas* c_signal = new TCanvas("c_signal");

   hh_signal->Draw("COLZ");

   RooAbsPdf *pdf;
   if (bSig) pdf = w_hgg->pdf(Form("SigPdf_cat%d",icat));
   else pdf = w_hgg->pdf(Form("BkgPdf_cat%d",icat));

   RooArgSet* set = new RooArgSet(*mgg, *mjj);


    TH2F* hpdf_signal = new TH2F();
    hh_signal->Copy(*hpdf_signal);
    hpdf_signal->Scale(0);

    for (int i = 1; i < hpdf_signal->GetNbinsX()+1; i++)
      for (int j = 1; j < hpdf_signal->GetNbinsY()+1; j++){
	double mggmin = hpdf_signal->GetXaxis()->GetBinLowEdge(i);
	double mggmax = hpdf_signal->GetXaxis()->GetBinUpEdge(i);

	double mjjmin = hpdf_signal->GetYaxis()->GetBinLowEdge(j);
	double mjjmax = hpdf_signal->GetYaxis()->GetBinUpEdge(j);

	mgg->setRange("signal",mggmin,mggmax) ;
	mjj->setRange("signal",mjjmin,mjjmax) ;

	RooAbsReal* igx_sig = pdf->createIntegral(*set,*set, "signal") ;
	double contrib = igx_sig->getVal()*NormSig;

	hpdf_signal->SetBinContent(i, j, contrib);

      }

    TCanvas* f_signal = new TCanvas("f_signal");
      
    hpdf_signal->Draw("COLZ");

    cout << "signal integral = " << hpdf_signal->Integral() << endl;


    // calculate the difference

    TH2F* hdiff = new TH2F();
    hh_signal->Copy(*hdiff);
    hdiff->Add(hpdf_signal, -1);
    const double alpha = 1 - 0.6827;

    TH2F* herror = new TH2F();
    hh_signal->Copy(*herror);

    double chi2 = 0;

    for (int i = 1; i < herror->GetNbinsX()+1; i++)
      for (int j = 1; j < herror->GetNbinsY()+1; j++){
	double N = herror->GetBinContent(i, j);
	double L = (N==0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
	double U = ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
	double ddiff = hdiff->GetBinContent(i, j);

	double err = 0;

	if (ddiff > 0) err = N-L;
	else err = U-N;
	hdiff->SetBinContent(i, j, ddiff/err);
	double pull =  ddiff/err;

	cout << "pull " << pull << " err = " << err << " i = " << i << " j = " << j << endl; 

	chi2 += pull*pull;

      }


    cout << "chi2 = " << chi2 << endl;
    //    hdiff->Divide(herror);
   
    TCanvas* f_diff = new TCanvas("f_diff");
      
    hdiff->Draw("COLZ");

    if (bSig){
      c_signal->SaveAs(Form("DataSignal_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      c_signal->SaveAs(Form("DataSignal_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));
      
      f_signal->SaveAs(Form("PDFSignal_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      f_signal->SaveAs(Form("PDFSignal_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));
      

      f_diff->SaveAs(Form("PullSignal_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      f_diff->SaveAs(Form("PullSignal_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));

    } else {
      c_signal->SaveAs(Form("DataBackground_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      c_signal->SaveAs(Form("DataBackground_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));
      
      f_signal->SaveAs(Form("PDFBackground_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      f_signal->SaveAs(Form("PDFBackground_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));
      

      f_diff->SaveAs(Form("PullBackground_%s_cat%d_norm%d.pdf", label.c_str(), icat, numEntries));
      f_diff->SaveAs(Form("PullBackground_%s_cat%d_norm%d.png", label.c_str(), icat, numEntries));
    }


    /*
    c_signal->SaveAs(Form("DataSignal_cat%d_norm%d.pdf", icat, numEntries));
    c_signal->SaveAs(Form("DataSignal_cat%d_norm%d.png", icat, numEntries));

    f_signal->SaveAs(Form("PDFSignal_cat%d_norm%d.pdf", icat, numEntries));
    f_signal->SaveAs(Form("PDFSignal_cat%d_norm%d.png", icat, numEntries));


    f_diff->SaveAs(Form("PullSignal_cat%d_norm%d.pdf", icat, numEntries));
    f_diff->SaveAs(Form("PullSignal_cat%d_norm%d.png", icat, numEntries));
    */

   /*
   string file;
   if (category == 0) file = string("CMS_bkg_fit_CMS_jj_4btag_cat0");
   if (category == 1) file = string("CMS_bkg_fit_CMS_jj_3btag_HPHP_cat1");

   // Get three of the functions inside, exponential, linear polynomial, power law
   RooAbsPdf *pdf_exp = w_hgg->pdf(file.c_str());

   // Fit the functions to the data to set the "prefit" state (note this can and should be redone with combine when doing 
   // bias studies as one typically throws toys from the "best-fit"
   //RooDataSet *data = w_hgg->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
   //pdf_exp->fitTo(*data);  // index 0
   //   pdf_pow->fitTo(*data); // index 1 
   //   pdf_pol->fitTo(*data);   // index 2

   // Make a plot (data is a toy dataset)
   RooPlot *plot = mass->frame(); //  data->plotOn(plot);
   pdf_exp->plotOn(plot,RooFit::LineColor(kGreen));
   //   pdf_pol->plotOn(plot,RooFit::LineColor(kBlue));
   //  pdf_pow->plotOn(plot,RooFit::LineColor(kRed));
   plot->SetTitle("PDF fits to toy data");
   plot->Draw();

   cout << "Norm = " << pdf_exp->getNorm() << endl;

   // Make a RooCategory object. This will control which of the pdfs is "active"
   RooCategory cat("pdf_index","Index of Pdf which is active");

   // Make a RooMultiPdf object. The order of the pdfs will be the order of their index, ie for below 
   // 0 == exponential
   // 1 == linear function
   // 2 == powerlaw
   RooArgList mypdfs;
   mypdfs.add(*pdf_exp);
   //  mypdfs.add(*pdf_pol);
   //  mypdfs.add(*pdf_pow);
   
   RooMultiPdf multipdf("roomultipdf","All Pdfs",cat,mypdfs);
   
   // As usual make an extended term for the background with _norm for freely floating yield
   RooRealVar norm("roomultipdf_norm","Number of background events",0,10000);
   
   // Save to a new workspace
   TFile *fout = new TFile(Form("background_pdfs_cat%d.root",category),"RECREATE");
   RooWorkspace wout("backgrounds","backgrounds");
   wout.import(cat);
   wout.import(norm);
   wout.import(multipdf);
   wout.Print();
   wout.Write();
   */
}
