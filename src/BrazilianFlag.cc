#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TNtuple.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "HiggsAnalysis/bbggLimits/interface/Colors.h"
double br=0;
float x01=0;
float y01=0;
float z01=0;
using namespace std;
bool PrintRadion13TeV=true;
TStyle* CreateStyle(std::string name1,std::string name2)
{
	TStyle *defaultStyle = new TStyle(name1.c_str(),name2.c_str());
	//  defaultStyle->SetOptStat(0000);
	//  defaultStyle->SetOptFit(000); 
	//  defaultStyle->SetPalette(1);
  	/////// pad ////////////
  	defaultStyle->SetPadBorderMode(1);
  	defaultStyle->SetPadBorderSize(1);
  	defaultStyle->SetPadColor(0);
  	defaultStyle->SetPadTopMargin(0.08);
  	defaultStyle->SetPadBottomMargin(0.13);
  	defaultStyle->SetPadLeftMargin(0.15);
  	defaultStyle->SetPadRightMargin(0.02);
 	/////// canvas /////////
  	defaultStyle->SetCanvasBorderMode(0);
  	defaultStyle->SetCanvasColor(0);
	//  defaultStyle->SetCanvasDefH(600);
	//  defaultStyle->SetCanvasDefW(600);
  	/////// frame //////////
  	defaultStyle->SetFrameBorderMode(0);
  	defaultStyle->SetFrameBorderSize(1);
  	defaultStyle->SetFrameFillColor(0); 
  	defaultStyle->SetFrameLineColor(1);
  	/////// label //////////
	//  defaultStyle->SetLabelOffset(0.005,"XY");
	//  defaultStyle->SetLabelSize(0.05,"XY");
  	defaultStyle->SetLabelFont(42,"XY");
 	/////// title //////////
	//  defaultStyle->SetTitleOffset(1.1,"X");
	//  defaultStyle->SetTitleSize(0.01,"X");
	//  defaultStyle->SetTitleOffset(1.25,"Y");
	//  defaultStyle->SetTitleSize(0.05,"Y");
  	defaultStyle->SetTitleFont(42, "XYZ");
  	/////// various ////////
  	defaultStyle->SetNdivisions(303,"Y");
  	//defaultStyle->SetTitleFillStyle(10, "Z");

	//  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

	//    defaultStyle->SetTitleColor(1, "XYZ");
	//    defaultStyle->SetTitleFont(42, "XYZ");
	//    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    	// defaultStyle->SetTitleYSize(Float_t size = 0.02);
    	//defaultStyle->SetTitleXOffset(0.9);
    	//defaultStyle->SetTitleYOffset(1.05);
    	// defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    	// For the axis labels:
    	defaultStyle->SetLabelColor(1, "XYZ");
    	defaultStyle->SetLabelFont(42, "XYZ");
   	// defaultStyle->SetLabelOffset(0.007, "XYZ");
    	defaultStyle->SetLabelSize(0.04, "XYZ");

    	// For the axis:
	//    defaultStyle->SetAxisColor(1, "XYZ");
    	defaultStyle->SetStripDecimals(kTRUE);
    	defaultStyle->SetTickLength(0.03, "XYZ");
    	defaultStyle->SetNdivisions(510, "XYZ");
	//    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
	//    defaultStyle->SetPadTickY(1);
    	defaultStyle->cd();
        return defaultStyle;
}

void BrazilianFlag(std::string path_dir,bool HH,bool base,bool low,bool obs,bool twobtag,std::string energy, float lumi) 
{
	std::cout<<green<<"CREATING BRAZILIAN FLAG"<<normal<<std::endl;
	std::vector<double>radMASS;
  	std::vector<std::string>dirs;
  	boost::filesystem::path dir(path_dir);
  	if(boost::filesystem::is_directory(path_dir)==false)
  	{
  		std::cout<<red<<"Directory "<<path_dir<<" doesn't exist"<<normal<<std::endl;
  		std::exit(1);
  	}
	for (boost::filesystem::directory_iterator itr(path_dir); itr!=boost::filesystem::directory_iterator(); ++itr)
	{
		if (boost::filesystem::is_directory(itr->status()))
		{
			std::string folder_name=itr->path().filename().string();
			
      			const char * pattern = "\\d+";
      			boost::regex re(pattern);
      			boost::sregex_iterator it(folder_name.begin(), folder_name.end(), re);
      			boost::sregex_iterator end;
      			std::vector<std::string>number;
			for( ; it != end; ++it)
    			{
         			number.push_back(it->str()); 
    			}
                        if(number.size()!=0)
			{
				radMASS.push_back(std::stod(number[0]));
				dirs.push_back(path_dir+"/"+folder_name+"/combine/higgsCombineTest.Asymptotic.mH125_m"+number[0]);
			}
		}
		
	}
	std::string name1="defaultStyle";
	std::string name2="Default Style";
	TStyle *defaultStyle=CreateStyle(name1,name2);
	/*bool HH=false;
	bool base = true;
	bool low=false;
	bool obs=false;//radlim${radion[$j]}_CSV
	bool twobtag=false;
	int combo = 0;
 	if (combo == 0) HH=0, low=0, twobtag=0;
 	if (combo == 1) HH=0, low=1, twobtag=1;
 	if (combo == 2) HH=0, low=1, twobtag=0;
 	if (combo == 3) HH=1, low=0, twobtag=0;
 	if (combo == 4) HH=1, low=1, twobtag=1;
 	if (combo == 5) HH=1, low=1, twobtag=0;*/
	//  TLegend *leg = new TLegend(0.65,0.5,0.99,0.94);
  	TLegend *leg = new TLegend(0.23,0.65,0.45,0.92);
  	if (low) leg = new TLegend(0.23,0.69,0.45,0.92);
  	leg->SetFillColor(kWhite);
  	leg->SetFillStyle(0);
  	leg->SetBorderSize(0);
  	leg->SetTextSize(0.037);

  	TLegend *leg1 = new TLegend(0.55,0.65,0.99,0.86);
  	leg1->SetFillColor(kWhite);
  	leg1->SetFillStyle(0);
  	leg1->SetBorderSize(0);
  	leg1->SetTextSize(0.036);
	
	//std::vector<double>radMASS{{260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,400,450,500,550,600,650,700,800,900,1000,1100}};
        //std::vector<double>radMASS2{{260,270,300,350,400}};
        //std::vector<double>radMASS3{{260,270,300,350,400}};
	//radMASS={260,270,300,350,400};
        std::vector<double>radCX(int(radMASS.size()),1.0);
        //std::vector<double>radCX2(int(radMASS2.size()),1.0);
        //std::vector<double>radCX3(int(radMASS3.size()),1.0);
    if(!HH) br=1;//1./(0.577*0.00228); 
    if(HH) br=1./(2*0.577*0.00228); // HH
    //////////////////////////////////////////
    // draw the radion line // MR      rad_CX(fb)       grav_CX(fb)
    //ntuple = new TNtuple("ntuple","NTUPLE","x:y:z");
    TNtuple *ntuple = new TNtuple("ntuple","NTUPLE","x:y"); 
    //float x0,y0,z0;
    std::string line="";
    std::string line_bfhh="";
    //char line[127], line_bfhh[127];
    std::ifstream fp("CX_radion_for_paper.data");
    std::ifstream fp_bfhh("CX_radion_for_paper_BFHH.data");
    while(getline(fp,line)&&getline(fp_bfhh,line_bfhh))
    {
      sscanf(&line[0],"%f %f",&x01,&y01);
      sscanf(&line_bfhh[0],"%f %f",&x01,&z01);
      printf("x0=%f, y0=%f, z0=%f\n",x01,y01,z01);
      ntuple->Fill(x01,y01*z01);
    }
    std::cout<<ntuple->GetNvar()<<" "<<ntuple->GetVar1()<<ntuple->GetVar2()<<std::endl;
    // graviton
    std::cout<<"Bulk graviton"<<std::endl;
    std::ifstream fpp("CX_grav_asfeaseability.data");
    TNtuple *ntupleg = new TNtuple("ntupleg","NTUPLE","x:y");
    while (getline(fpp,line))
    {
      sscanf(&line[0],"%f %f",&x01,&y01);
      printf("x0=%f, y0=%f\n",x01,y01);
      ntupleg->Fill(x01,y01); 
      //cout<<x0<<endl;
      //ntuple->Fill(x0,y0); 
    }
    std::cout<<ntupleg->GetNvar()<<" "<<ntupleg->GetVar1()<<ntupleg->GetVar2()<<std::endl;
    // RS graviton
    std::cout<<"RS graviton"<<std::endl;
    std::ifstream fppp("CX_RSgrav_asfeaseability_low.data");
    TNtuple *ntuplegrs = new TNtuple("ntuplegrs","NTUPLE","x:y"); 
    while (getline(fppp,line)) 
    {
      sscanf(&line[0],"%f %f",&x01,&y01);
      printf("x0=%f, y0=%f\n",x01,y01);
      ntuplegrs->Fill(x01,y01); 
      //cout<<x0<<endl;
      //ntuple->Fill(x0,y0); 
    }
    std::cout<<ntuplegrs->GetNvar()<<" "<<ntupleg->GetVar1()<<ntupleg->GetVar2()<<std::endl;
    // draw the BR flag -- Main results
    //std::string path="/home/maxime/Documents/POUGNE_20111106/CERN/BOOSTED/RADION/ANPAS2014/papers/HIG-13-032/trunk/figures/fit_functions_paper/radlim_v44_fitTo2D_resSearch_withRegKinFit/radlim";
    std::string path="";
    //std::string pathh="";
    //std::string pathhh="";
    /*if(argc==1)
    {
	std::cout<<"Please give the folder in argument"<<std::endl;
        std::exit(-2);
    }
    int sz = sizeof(argv[1]);
    path.assign(argv[1],sz);
    path+="/radlim";*/
   // path="./radlim_fitTo2D_resSearch_withRegKinFit_interpolation/radlim";
    //pathh="./radlim_fitTo2D_resSearch_withRegKinFit/radlim";
    //pathhh="./radlim_fitTo2D_resSearch_withRegKinFit_justshape/radlim";
    //std::string path2="/higgsCombineTest.Asymptotic.mH125.mR";
    //std::string path4=".root";
   // std::vector<std::string>Mass{"260","270","280","290","300","310","320","330","340","350","360","370","380","390","400","410","420","430","440","450"};//,"400","450","500","550","600","650","700","800","900","1000","1100"};
    //Mass={"260","270","300","350","400"};
    //std::vector<std::string>Mass2{"260","270","300","350","400"};
    //std::vector<std::string>Mass3{"260","270","300","350","400"};
    std::vector<TFile*>TFileVec;
    std::vector<TTree*>TTrees;
    //std::vector<TFile*>TFileVec2;
    //std::vector<TTree*>TTrees2;
      //std::vector<TFile*>TFileVec3;
    //std::vector<TTree*>TTrees3;
    TBranch *b_limit=nullptr; 
    //TBranch *b_limit2=nullptr;
    //TBranch *b_limit3=nullptr;
    Double_t limit=0;
    //Double_t limit2=0;
    //Double_t limit3=0;
    std::vector<std::vector<float>>rad(6,std::vector<float>(int(radMASS.size()),0));
    //std::vector<std::vector<float>>rad2(6,std::vector<float>(int(Mass2.size())));
    //std::vector<std::vector<float>>rad3(6,std::vector<float>(int(Mass3.size())));
    TFile *file=nullptr;
    for(unsigned int i =0;i!=radMASS.size();++i)
    {
        std::string name=dirs[i];//path+radMASS[i]+path2+radMASS[i];
        //if(i<20) 
        //{
          //if (twobtag==false) name+="_onecatnohiggs.root";
          /*else*/ name+="_higgs.root";
        //}
        //std::cout<<i<<std::endl;
        //if (twobtag) name+="_onecat";
        //name+=path4;
        std::cout<<name<<std::endl;
        file = new TFile(name.c_str());
        if(!(file->IsOpen()))
				{
					std::cout<<"File : "<<name<<" not found"<<std::endl;
					std::exit(-1);
				}
        TFileVec.push_back(file);
        TTrees.push_back((TTree*)TFileVec[i]->Get("limit;1"));
        TTrees[i]->SetMakeClass(1);
        TTrees[i]->SetBranchAddress("limit", &limit, &b_limit);
        for (int k = 0; k<6; k++)
        {
          TTrees[i]->GetTree()->GetEntry(k);
          rad[k][i]=limit*radCX[i]*br;
          std::cout <<"MX = "<< radMASS[i]<<" centrality "<< i<<" limit = " << limit << std::endl; 
        }
    }
    /*for(unsigned int i =0;i!=Mass2.size();++i)
    {
        std::string name=pathh+Mass2[i]+path2+Mass2[i];
        std::cout<<name<<std::endl;
        if(i<20) 
        {
          if (twobtag) name+="_onecatnohiggs";
          else name+="_higgs";
        }
        std::cout<<i<<std::endl;
        if (twobtag) name+="_onecat";
        name+=path4;
        std::cout<<name<<std::endl;
        file = new TFile(name.c_str());
        if(!(file->IsOpen()))std::exit(-1);
        TFileVec2.push_back(file);
        TTrees2.push_back((TTree*)TFileVec2[i]->Get("limit;1"));
        TTrees2[i]->SetMakeClass(1);
        TTrees2[i]->SetBranchAddress("limit", &limit2, &b_limit2);
        for (int k = 0; k<6; k++)
        {
          TTrees2[i]->GetTree()->GetEntry(k);
          rad2[k][i]=limit2*radCX2[i]*br;
          std::cout <<"MX = "<< radMASS2[i]<<" centrality "<< i<<" limit = " << limit2 << std::endl; 
        }
    }*/
    /*for(unsigned int i =0;i!=Mass3.size();++i)
    {
        std::string name=pathhh+Mass3[i]+path2+Mass3[i];
        std::cout<<name<<std::endl;
        if(i<20) 
        {
          if (twobtag) name+="_onecatnohiggs";
          else name+="_higgs";
        }
        std::cout<<i<<std::endl;
        if (twobtag) name+="_onecat";
        name+=path4;
        std::cout<<name<<std::endl;
        file = new TFile(name.c_str());
        if(!(file->IsOpen()))std::exit(-1);
        TFileVec3.push_back(file);
        TTrees3.push_back((TTree*)TFileVec3[i]->Get("limit;1"));
        TTrees3[i]->SetMakeClass(1);
        TTrees3[i]->SetBranchAddress("limit", &limit3, &b_limit3);
        for (int k = 0; k<6; k++)
        {
          TTrees3[i]->GetTree()->GetEntry(k);
          rad3[k][i]=limit3*radCX3[i]*br;
          std::cout <<"MX = "<< radMASS3[i]<<" centrality "<< i<<" limit = " << limit3 << std::endl; 
        }
    }*/
    delete file;
    TPaveText *pt = new TPaveText(0.15,0.93,0.95,0.97, "brNDC");
    //   pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    //   pt->SetShadowColor(kWhite);
    if (!HH) pt->AddText(TString::Format("CMS                                    L = %.2f fb^{-1}    #sqrt{s} = %s",lumi,energy.c_str()));
    else pt->AddText(TString::Format("CMS (Unpublished)                                L = %.2f fb^{-1}    #sqrt{s} = %s",lumi,energy.c_str()));
    pt->SetTextSize(0.04);
    TPaveText *Purity = new TPaveText(0.78,0.53,0.88,0.58, "brNDC");
    //   pt->SetName("title");
    Purity->SetBorderSize(0);
    Purity->SetFillColor(0);
    //   pt->SetShadowColor(kWhite);
    if(twobtag) Purity->AddText("High Purity");
    else Purity->AddText("High+Medium Purity");
    Purity->SetTextSize(0.04);
    Purity->SetTextColor(kBlue);

    TPaveText *BF = new TPaveText(0.74,0.60,0.88,0.63, "brNDC");
    //   pt->SetName("title");
    BF->SetBorderSize(0);
    BF->SetFillColor(0);
    //   pt->SetShadowColor(kWhite);
    BF->AddText("X #rightarrow HH #rightarrow #gamma#gammab#bar{b}");
    BF->SetTextSize(0.04);
    BF->SetTextColor(kBlue);

    //  TPaveText *BF = new TPaveText(0.7,0.93,0.8,0.99, "brNDC");
    //   pt->SetName("title");
    //pt->SetBorderSize(0);
    //pt->SetFillColor(0);
    //   pt->SetShadowColor(kWhite);
    //pt->AddText("CMS Preliminary           L = 19.7 fb^{-1}          #sqrt{s} = 8 TeV");
    //pt->SetTextSize(0.04);

    // we do a plot r*MR
    TMultiGraph *mg = new TMultiGraph();
    if(!HH) mg->SetMinimum(0.01);
    if(!HH && low) mg->SetMinimum(0);
    if(HH) mg->SetMinimum(10); // HH
    if(HH && low) mg->SetMinimum(0); // HH 1000000
    if(HH && low) mg->SetMaximum(5200); // HH 1000000
    if(!HH && low) mg->SetMaximum(39.99); // 10000
    if(HH && !low) mg->SetMaximum(200000); // HH 
    if(!HH && !low) mg->SetMaximum(600); // 
    if(!HH && !low) mg->SetMaximum(40);

    TLine* T = new TLine(400, 0.01, 400, 1);
    T->SetLineStyle(2);
    //  if(!HH && twobtag) {mg->SetTitle("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb - 2 btag only }}");
    //    TText *text = pt->AddText("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb - 2 btag only}}");
    //  }
    //  if(!HH && !twobtag) {mg->SetTitle("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb }}");
    //    TText *text = pt->AddText("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb }}");
    //  }

    ////    mg->SetTitle("CMS Preliminary     L = 19.7 fb^{-1}    #sqrt{s} = 8 TeV; M_{X} (GeV); #sigma(pp -> X)*BR(X -> HH -> #gamma #gamma bb) (fb)");
    if (!HH) mg->SetTitle("; m_{X} (GeV); #sigma(pp #rightarrow X) #times BR(X #rightarrow HH #rightarrow #gamma#gammab#bar{b}) (fb)");
    if (HH) mg->SetTitle("; m_{X} (GeV); #sigma(pp #rightarrow X) #times BR(X #rightarrow HH) (fb)");


    //  if(HH) {mg->SetTitle("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)(fb) - SM Higgs BR}}{#scale[0.8]{CMS preliminary 19.7/fb  }}");
    //    TText *text = pt->AddText("#splitline{#scale[1.0]{#sigma(pp -> X)*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb  }}");
    //}
    //if(HH) mg->SetTitle("Radion > HH > #gamma #gamma bb~");
    //  float radMASS[nmass]={300,500,700,1000};
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    int nmax=0;
    //int nmax2=0;
    //int nmax3=0;
    if(low)
    {
	if(radMASS.size()<5) nmax= radMASS.size();//nmass;//nmass
	else nmax=5;
    }
    if(!low) 
    {
			nmax= radMASS.size();//nmass
       // nmax2= Mass2.size();//nmass
        // nmax3= Mass3.size();
    }
    Double_t yobs[nmax], y2up[nmax], y1up[nmax], y1down[nmax], y2down[nmax], ymean[nmax];
    for(int i=0;i!=nmax;++i)
    {
	yobs[i]=0;
	y2up[i]=0; 
	y1up[i]=0;
	y1down[i]=0;
	y2down[i]=0;
	ymean[i]=0;
    }
    //Double_t yobs2[nmax2], y2up2[nmax2], y1up2[nmax2], y1down2[nmax2], y2down2[nmax2], ymean2[nmax2];
    //Double_t yobs3[nmax3], y2up3[nmax3], y1up3[nmax3], y1down3[nmax3], y2down3[nmax3], ymean3[nmax3];
    //
    
    c1->cd();
    pt->Draw();

    if(low) c1->SetLogy(0);
    if(!low) c1->SetLogy(1);
    c1->SetLogy(0);
    //c1->SetGrid();

    //  if(HH) ntuple->Draw("y*(9/1)*0.25 : x");
    //  if(!HH)  ntuple->Draw("y*0.25*0.577*0.00228*2*9 : x");

    if(HH) ntuple->Draw("y*9 : x");
    if(!HH)  ntuple->Draw("y*0.577*0.00228*2*9 : x");
    TGraphErrors *radion10 = new TGraphErrors(ntuple->GetSelectedRows(), ntuple->GetV2(), ntuple->GetV1());
    radion10->SetLineColor(kRed); 
    radion10->SetLineWidth(3);
    if(HH)  ntuple->Draw("y : x");
    if(!HH)  ntuple->Draw("y*0.577*0.00228*2 : x");
    TGraphErrors *radion = new TGraphErrors(ntuple->GetSelectedRows(), ntuple->GetV2(), ntuple->GetV1());//point,radionx,radiony);
    radion->SetLineWidth(3);
    TGraphErrors *rsgrav =nullptr;
    TGraphErrors *bulk3=nullptr;
    if(base)
    {
      // RS graviton lambda =1

      if(HH)  ntupleg->Draw("y : x");
      if(HH)  ntuplegrs->Draw("y : x");
      if(!HH)  ntuplegrs->Draw("y*0.577*0.00228*2 : x");
      rsgrav = new TGraphErrors(ntuplegrs->GetSelectedRows(), ntuplegrs->GetV2(), ntuplegrs->GetV1());
      rsgrav->SetLineColor(kRed); 
      rsgrav->SetLineWidth(3);
      rsgrav->SetLineStyle(2);
      // bulk graviton
      if(HH)  ntupleg->Draw("y : x");
      if(!HH)  ntupleg->Draw("y*0.577*0.00228*2 : x");
      bulk3 = new TGraphErrors(ntupleg->GetSelectedRows(), ntupleg->GetV2(), ntupleg->GetV1());//
      bulk3->SetLineColor(kBlack); 
      bulk3->SetLineWidth(3);
      bulk3->SetLineStyle(2);
    } //if base plt graviton
    c1->Clear();
    for (Int_t i=0;i<nmax;i++)
    {
      std::cout << "" << radMASS[i] << " &  & " << Form("%.2f",rad[2][i]) << " &  &  \\\\" <<std::endl;
      //yobs[i] = rad[5][i];
      y2up[i] = rad[0][i];
      y1up[i] = rad[1][i];
      ymean[i] = rad[2][i]; 
      y1down[i] = rad[3][i];
      y2down[i] =  rad[4][i];  
      std::cout<<yellow<<y2down[i]<<normal<<std::endl;
      std::cout<<green<<y1down[i]<<normal<<std::endl;
      std::cout<< ymean[i]<<std::endl; 
      std::cout<<green<<y1up[i]<<normal<<std::endl;
      std::cout<<yellow<<y2up[i]<<normal<<std::endl;
      
    }



    /*for (Int_t i=0;i<nmax2;i++)
    {
      std::cout << "" << radMASS2[i] << " &  & " << Form("%.2f",rad2[2][i]) << " &  &  \\\\" <<std::endl;
      //yobs[i] = rad[5][i];
      y2up2[i] = rad2[0][i];
      y1up2[i] = rad2[1][i];
      ymean2[i] = rad2[2][i]; 
      y1down2[i] = rad2[3][i];
      y2down2[i] =  rad2[4][i];  
      std::cout<<yellow<<y2down2[i]<<normal<<std::endl;
      std::cout<<green<<y1down2[i]<<normal<<std::endl;
      std::cout<< ymean2[i]<<std::endl; 
      std::cout<<green<<y1up2[i]<<normal<<std::endl;
      std::cout<<yellow<<y2up2[i]<<normal<<std::endl;
      
    }*/
   /* for (Int_t i=0;i<nmax3;i++)
    {
      std::cout << "" << radMASS3[i] << " &  & " << Form("%.2f",rad3[2][i]) << " &  &  \\\\" <<std::endl;
      //yobs[i] = rad[5][i];
      y2up3[i] = rad3[0][i];
      y1up3[i] = rad3[1][i];
      ymean3[i] = rad3[2][i]; 
      y1down3[i] = rad3[3][i];
      y2down3[i] =  rad3[4][i];  
      std::cout<<yellow<<y2down3[i]<<normal<<std::endl;
      std::cout<<green<<y1down3[i]<<normal<<std::endl;
      std::cout<< ymean3[i]<<std::endl; 
      std::cout<<green<<y1up3[i]<<normal<<std::endl;
      std::cout<<yellow<<y2up3[i]<<normal<<std::endl;
      
    }*/
    TGraphErrors *grobs = new TGraphErrors(1);
    grobs->SetMarkerStyle(kFullDotLarge); 
    grobs->SetMarkerColor(kBlue); 
    grobs->SetLineColor(kBlue);
    grobs->SetLineWidth(2);
    grobs->SetLineStyle(2);

    TGraphErrors *grmean = new TGraphErrors(1);
    grmean->SetLineColor(1);
    grmean->SetLineWidth(3);
    grmean->SetLineStyle(3);
    grmean->SetMarkerSize(0);
    for(int j=0;j<nmax;j++)
    {
      grobs->SetPoint(j, radMASS[j], yobs[j]);
      grmean->SetPoint(j, radMASS[j], ymean[j]);
    }
 
  //  TGraphErrors *grobs2 = new TGraphErrors(1);
  //  grobs2->SetMarkerStyle(kStar); 
  //  grobs2->SetMarkerColor(kBlue); 
  //  grobs2->SetLineColor(kBlue);
  //  grobs2->SetLineWidth(2);
  //  grobs2->SetLineStyle(2);
   // 
  //  TGraphErrors *grmean2 = new TGraphErrors(1);
  //  grmean2->SetLineColor(kBlue);
   // grmean2->SetMarkerStyle(kStar);
  //  grmean2->SetLineWidth(3);
  //  grmean2->SetLineStyle(3);
  //  grmean2->SetMarkerSize(0.01);

  //  TGraphErrors *grgreen2 = new TGraphErrors(1);
   // grgreen2->SetLineColor(kGreen);
  //  grgreen2->SetMarkerStyle(kStar);
   // grgreen2->SetLineWidth(3);
  // grgreen2->SetLineStyle(3);
   // grgreen2->SetMarkerSize(1);
  //  grgreen2->SetMarkerColor(kOrange+10);

   // TGraphErrors *gryellow2 = new TGraphErrors(1);
   // gryellow2->SetLineColor(kRed);
  //  gryellow2->SetMarkerStyle(kStar);
  //  gryellow2->SetLineWidth(3);
  //  gryellow2->SetLineStyle(3);
  //  gryellow2->SetMarkerSize(1);
  //  gryellow2->SetMarkerColor(kGreen);



  //  TGraphErrors *grgreen3 = new TGraphErrors(1);
  //  grgreen3->SetLineColor(kGreen);
  //  grgreen3->SetMarkerStyle(kStar);
   // grgreen3->SetLineWidth(3);
   // grgreen3->SetLineStyle(3);
   // grgreen3->SetMarkerSize(1);
   // grgreen3->SetMarkerColor(kOrange+10);

    //TGraphErrors *gryellow3 = new TGraphErrors(1);
   // gryellow3->SetLineColor(kRed);
   // gryellow3->SetMarkerStyle(kStar);
    //gryellow3->SetLineWidth(3);
  //  gryellow3->SetLineStyle(3);
   // gryellow3->SetMarkerSize(1);
   // gryellow3->SetMarkerColor(kGreen);




    
















   //TGraphErrors *grobs3 = new TGraphErrors(1);
  //  grobs3->SetMarkerStyle(kFullSquare); 
  //  grobs3->SetMarkerColor(kBlue); 
 //   grobs3->SetLineColor(kBlue);
  //  grobs3->SetLineWidth(2);
  //  grobs3->SetLineStyle(2);
    
   // TGraphErrors *grmean3 = new TGraphErrors(1);
  //  grmean3->SetLineColor(kMagenta);
  //  grmean3->SetMarkerStyle(22);
   // grmean3->SetMarkerColor(kMagenta);
   // grmean3->SetLineWidth(3);
   // grmean3->SetLineStyle(3);
   // grmean3->SetMarkerSize(1);

   // TGraphErrors *grgreen4 = new TGraphErrors(1);
   // grgreen4->SetLineColor(kGreen);
   // grgreen4->SetMarkerStyle(kFullSquare);
   // grgreen4->SetLineWidth(3);
   // grgreen4->SetLineStyle(3);
   // grgreen4->SetMarkerSize(1);
   // grgreen4->SetMarkerColor(kOrange+2);

   /// TGraphErrors *gryellow4 = new TGraphErrors(1);
   //gryellow4->SetLineColor(kRed);
    //gryellow4->SetMarkerStyle(kFullSquare);
    //gryellow4->SetLineWidth(3);
   // gryellow4->SetLineStyle(3);
    //gryellow4->SetMarkerSize(1);
    //gryellow4->SetMarkerColor(kGreen+2);



   // TGraphErrors *grgreen5 = new TGraphErrors(1);
   // grgreen5->SetLineColor(kGreen);
   // grgreen5->SetMarkerStyle(kFullSquare);
    //grgreen5->SetLineWidth(3);
    //grgreen5->SetLineStyle(3);
   // grgreen5->SetMarkerSize(1);
    //grgreen5->SetMarkerColor(kOrange+2);

   // TGraphErrors *gryellow5 = new TGraphErrors(1);
  //  gryellow5->SetLineColor(kRed);
  //  gryellow5->SetMarkerStyle(kFullSquare);
   /// gryellow5->SetLineWidth(3);
   // gryellow5->SetLineStyle(3);
   // gryellow5->SetMarkerSize(1);
   // gryellow5->SetMarkerColor(kGreen+2);










   /*  for(int j=0;j<nmax2;j++)
    {
      grobs2->SetPoint(j, radMASS2[j], yobs2[j]);
      grmean2->SetPoint(j, radMASS2[j], ymean2[j]);
      grgreen2->SetPoint(j,radMASS2[j],y2up2[j]);
      grgreen3->SetPoint(j,radMASS2[j],y2down2[j]);
      gryellow2->SetPoint(j,radMASS2[j],y1up2[j]);
      gryellow3->SetPoint(j,radMASS2[j],y1down2[j]);
    }*/
   /*  for(int j=0;j<nmax3;j++)
    {
      grobs3->SetPoint(j, radMASS3[j], yobs3[j]);
      grmean3->SetPoint(j, radMASS3[j], ymean3[j]);
      grgreen4->SetPoint(j,radMASS3[j],y2up3[j]);
      grgreen5->SetPoint(j,radMASS3[j],y2down3[j]);
      gryellow4->SetPoint(j,radMASS3[j],y1up3[j]);
      gryellow5->SetPoint(j,radMASS3[j],y1down3[j]);
    }*/
   mg->Add(grmean,"L*");//->Draw("same,AC*");
   if(obs) mg->Add(grobs,"L,P");//->Draw("AC*");



  mg->Draw("AP");
  mg->GetXaxis()->SetRangeUser(100,1400);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.042);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.3);
  mg->GetYaxis()->CenterTitle(true);
  mg->GetXaxis()->SetTitleOffset(1.01);
  mg->GetXaxis()->CenterTitle(true);
  // histo to shade
   int n=nmax;
   TGraph *grgreen = new TGraph(2*n);
   TGraph *gryellow = new TGraph(2*n);
   for (int i=0;i<n;i++) 
   {
    grgreen->SetPoint(i,radMASS[i],y2up[i]);
    grgreen->SetPoint(n+i,radMASS[n-i-1],y2down[n-i-1]);
    //
    gryellow->SetPoint(i,radMASS[i],y1up[i]);
    gryellow->SetPoint(n+i,radMASS[n-i-1],y1down[n-i-1]);

    std::cout<<" observed "<<radMASS[i]<<" "<<y2down[i]<<std::endl; 
  }

   /*TGraph *grgreen2 = new TGraph(2*nmax2);
   TGraph *gryellow2 = new TGraph(2*nmax2);
   for (int i=0;i<nmax2;i++) 
   {
    grgreen2->SetPoint(i,radMASS2[i],y2up2[i]);
    grgreen2->SetPoint(n+i,radMASS2[i],y2down2[i]);
    //
    gryellow2->SetPoint(i,radMASS2[i],y1up2[i]);
    gryellow2->SetPoint(i,radMASS2[i],y1down2[i]);

    std::cout<<" observed "<<radMASS2[i]<<" "<<y2down2[i]<<std::endl; 
  }*/


  grgreen->SetMarkerSize(0.00001);
  grgreen->SetMarkerColor(kOrange);
  grgreen->SetLineColor(kOrange);
  grgreen->SetFillColor(kOrange);
  grgreen->Draw("f"); 
  gryellow->SetFillColor(kGreen+2);
  gryellow->SetLineColor(kGreen+2);
  gryellow->SetMarkerColor(kGreen+2);
  gryellow->SetMarkerSize(0.00001);
  gryellow->Draw("f"); 
  grmean->Draw("L,same");


/*
  grgreen2->SetMarkerSize(0.00001);
  grgreen2->SetMarkerColor(kRed);
  grgreen2->SetMarkerStyle(kStar);
  grgreen2->SetLineColor(kRed);
  grgreen2->SetFillColor(kRed);*/
  //grgreen2->Draw("*,same"); 
  //grgreen3->Draw("*,same");
  //grgreen4->Draw("*,same"); 
  //grgreen5->Draw("*,same");
  /*gryellow2->SetFillColor(kBlue+2);
  gryellow2->SetMarkerStyle(kStar);
  gryellow2->SetLineColor(kBlue+2);
  gryellow2->SetMarkerColor(kBlue+2);
  gryellow2->SetMarkerSize(0.00001);*/
 
  //gryellow2->Draw("*,same");
  //gryellow3->Draw("*,same");
  //gryellow4->Draw("*,same");
  //gryellow5->Draw("*,same");
  //grmean2->SetMarkerStyle(kStar);
 // grmean2->SetMarkerStyle(kBlue);
 // grmean2->Draw("*,same");
 // grmean3->Draw("P,same");
  //grobs3->Draw("*,same");
  //grmean2->SetMarkerSize(1);
  if(obs) grobs->Draw("L,P,E,same");
  
  /////////////////////////////////////////////////////////////
  // cross checks
  TGraphErrors *gcheck1=nullptr;
  TGraphErrors *gcheck2=nullptr;
  TGraphErrors *gcheck4=nullptr;
  TGraphErrors *gcheck7=nullptr;
  if(!HH)
  {
    int check1mass[3] = {400,450,500};
    Double_t check1limit[3]={1.66,1.25,0.96};
    gcheck1 = new TGraphErrors(1);
    gcheck1->SetMarkerStyle(kFullDotLarge); 
    gcheck1->SetLineColor(kRed+2);
    gcheck1->SetMarkerColor(kRed+2);
    gcheck1->SetLineWidth(3);

    for(int j=0;j<3;j++) gcheck1->SetPoint(j, check1mass[j], check1limit[j]*br);
    if(low && !base)gcheck1->Draw("L,P");
    // sidebands
    int check2mass[5] = {270,300,350,400,450};
    Double_t check2limit[5]={2.49, 2.63, 1.96, 1.37, 1.09};
    gcheck2 = new TGraphErrors(1);
    gcheck2->SetMarkerStyle(kFullDotLarge); 
    gcheck2->SetLineColor(kRed-7);
    gcheck2->SetMarkerColor(kRed-7);
    gcheck2->SetLineWidth(3);
    for(int j=0;j<5;j++) gcheck2->SetPoint(j, check2mass[j], check2limit[j]*br);
    if(low && !base)gcheck2->Draw("L,P");
    //
    int check4mass[11] = {400,450,500,550,600,650,700,800,900,1000,1100};
    Double_t check4limit[11]={2.65,1.88,1.14,1.13,0.93,0.82,0.76,0.68,0.68,0.74,0.93};
    gcheck4 = new TGraphErrors(1);
    gcheck4->SetMarkerStyle(kFullDotLarge); 
    gcheck4->SetLineColor(kViolet);
    gcheck4->SetMarkerColor(kViolet);
    gcheck4->SetLineWidth(3);
    for(int j=0;j<11;j++) gcheck4->SetPoint(j, check4mass[j], check4limit[j]*br);
    if(!low && !base) gcheck4->Draw("L,P");
    //
    //
    int check7mass[11] = {400,450,500,550,600,650,700,800,900,1000,1100}; // 4body wo/ kinfit
    Double_t check7limit[11]={2.90,2.46,2.29,1.63,1.40,1.16,1.11,0.75,0.68,0.67,0.74 };
    gcheck7 = new TGraphErrors(1);
    gcheck7->SetMarkerStyle(kFullDotLarge); 
    gcheck7->SetLineColor(kRed+3);
    gcheck7->SetMarkerColor(kRed+3);
    gcheck7->SetLineWidth(3);
    for(int j=0;j<11;j++) gcheck7->SetPoint(j, check7mass[j], check7limit[j]*br);
    if(!low && !base) gcheck7->Draw("L,P");
  } // close cross checks
  
  //////////////////////////////////////////////////////////////////////
  // 2btag wo higgs
  int check3mass[6] = {270,300,350,400,450,500};
  Double_t check3limit[6]={2.85,3.14,2.71,2.01,1.53,1.22};
  TGraphErrors *gcheck3 = new TGraphErrors(1);
  gcheck3->SetMarkerStyle(kFullDotLarge); 
  gcheck3->SetLineColor(kCyan);
  gcheck3->SetMarkerColor(kCyan);
  gcheck3->SetLineWidth(3);

  for(int j=0;j<3;j++) gcheck3->SetPoint(j, check3mass[j], check3limit[j]*br);
  
  // 2btag w higgs
  int check5mass[6] = {270,300,350,400,450,500};
  Double_t check5limit[6]={3.1406,3.6094,2.7266,2.3359,1.7734,1.3945};
  TGraphErrors *gcheck5 = new TGraphErrors(1);
  gcheck5->SetMarkerStyle(kFullDotLarge); 
  gcheck5->SetLineColor(kMagenta);
  gcheck5->SetMarkerColor(kMagenta);
  gcheck5->SetLineWidth(3);
  for(int j=0;j<6;j++) gcheck5->SetPoint(j, check5mass[j], check5limit[j]*br);
  if(!low && !base) gcheck5->Draw("L,P");
  // wo higgs
  int check6mass[6] = {270,300,350,400,450,500};
  Double_t check6limit[6]={3.1406,3.6094,2.7266,2.3359,1.7734,1.3945};
  TGraphErrors *gcheck6 = new TGraphErrors(1);
  gcheck6->SetMarkerStyle(kFullDotLarge); 
  gcheck6->SetLineColor(kBlue);
  gcheck6->SetMarkerColor(kBlue);
  gcheck6->SetLineWidth(3);
  for(int j=0;j<3;j++) gcheck6->SetPoint(j, check6mass[j], check6limit[j]*br);
  ///////////////////////////////////////////
  radion->Draw("L,same");
  if(!low)  radion10->Draw("L,same");
  if(base)  rsgrav->Draw("L,same");
  if(base && !low)  bulk3->Draw("L,same");

  if(base)  leg->SetHeader("WED: kl = 35, k/Mpl = 0.2, elementary top, no r/H mixing");

  if(obs) leg1->AddEntry(grobs, "Observed 95% upper limit", "L,P");  
  if(!base && low) {leg1->AddEntry(grmean, "Expected (m#gamma #gamma fit w/ Higgs)", "L"); }
  else leg1->AddEntry(grmean, "Expected 95% upper limit", "L");

  leg1->AddEntry(gryellow, "Expected limit #pm 1 #sigma", "f");
  leg1->AddEntry(grgreen, "Expected limit #pm 2 #sigma", "f");
  leg->AddEntry(radion, "radion (#Lambda_{R} = 3 TeV)", "L");
//  if(!low)  leg->AddEntry(radion10, "radion (#Lambda_{R} = 1 TeV)", "L");
//  if(base)  leg->AddEntry(rsgrav, "RS1 KK-graviton  ", "L");
//  if(base && !low)  leg->AddEntry(bulk3, "Bulk KK-graviton ", "L");


  if(!base)
  {
    if(low) leg->AddEntry(gcheck1, "m#gamma #gamma  fit (w/ higgs)", "LP");
    if(!low)leg->AddEntry(gcheck4, " 4 body fit - 2 btag only (w/ kinfit)", "LP");
    if(!low )leg->AddEntry(gcheck7, " 4 body fit - wo/ kinfit ", "LP");
    if(!low)leg->AddEntry(gcheck5, "m#gamma #gamma  fit - 2 btag only (wo/ Higgs)", "LP");
    if(low)leg->AddEntry(gcheck2, "#gamma #gamma Sidebands (stat only)", "LP");
  }
  leg->Draw();
  leg1->Draw();
  pt->Draw();
  if (!low) T->Draw();
  if (twobtag) Purity->Draw();
  if (HH) BF->Draw();
  
  if( !base)
  {
    TLine l(400,0.01,400,4.6);
    l.Draw();
  }
  c1->Update();
  std::string folderr=path_dir+"/BrazilianFlag";
  std::string name="/WP4_cutbased";
	if(HH)name+="_HH";
  if(low)name+="_low";
  if(base)name+="_base";
  if(twobtag)name+="_onecat";
  boost::filesystem::create_directory(folderr);
	folderr+=name;
 

 TGraph* ThisIsRadion;
  if(PrintRadion13TeV){
        Double_t RadToHH[] = {2238.55, 2376.22, 2104.53, 1804.82, 1540.86, 1154.92, 896.161, 733.914, 495.567, 362.581, 277.513, 219.251, 177.292, 146.142, 122.775, 104.76 , 90.6293, 79.4959};
        Double_t RadTobbgg[18];
        for ( int mm = 0; mm < 18; mm++){
                RadTobbgg[mm] = RadToHH[mm]*2.6e-03;
        }

        Double_t RadMasses[] = {260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900};

        ThisIsRadion = new TGraph(18, RadMasses, RadTobbgg);
        c1->cd();
        ThisIsRadion->Draw("Csame");
        ThisIsRadion->SetLineWidth(3);
        ThisIsRadion->SetLineColor(1);
        c1->Update();
  } 

  c1->SaveAs((folderr+".png").c_str());
        c1->SaveAs((folderr+".pdf").c_str());
        c1->SaveAs((folderr+".root").c_str());

  mg->SetMaximum(200);
  mg->SetMinimum(0.01);
  c1->SetLogy();
  c1->Update();
  c1->SaveAs((folderr+"LOG.png").c_str());
  c1->SaveAs((folderr+"LOG.pdf").c_str());

 //if(HH && low ) c1->SaveAs("/WP4_cutbased_HH_low.png"); // HH
  //if(HH && low ) c1->SaveAs("/WP4_cutbased_HH_low.pdf"); // HH
//  if(HH && low ) c1->SaveAs("/WP4_cutbased_HH_low.root"); // HH
 // if(!HH && low && !base) c1->SaveAs("/WP4_cutbased_low.png");
 // if(!HH && low && !base) c1->SaveAs("/WP4_cutbased_low.pdf");
 // if(!HH && low && !base) c1->SaveAs("/WP4_cutbased_low.root");
 //if(HH && !low ) c1->SaveAs("/WP4_cutbased_HH.png"); // HH
 // if(HH && !low ) c1->SaveAs("/WP4_cutbased_HH.pdf"); // HH
//if(HH && !low ) c1->SaveAs("/WP4_cutbased_HH.root"); // HH
 //if(!HH && !low && !base) c1->SaveAs("/WP4_cutbased.png");
  //if(!HH && !low && !base) c1->SaveAs("/WP4_cutbased.pdf");
  //if(!HH && !low && !base) c1->SaveAs("/WP4_cutbased.root");
  //
  //if(!HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_base_onecat.png"); // HH
  //if(!HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_base_onecat.pdf"); // HH
 // if(!HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_base_onecat.root"); // HH
  
 // if(HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_HH_onecat.png"); // HH
 // if(HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_HH_onecat.pdf"); // HH
 // if(HH && low && base && twobtag) c1->SaveAs("/WP4_cutbased_low_HH_onecat.root"); // HH
 // if(!HH && low && base && !twobtag) c1->SaveAs("/WP4_cutbased_low_base.png"); // HH
 // if(!HH && low && base && !twobtag) c1->SaveAs("/WP4_cutbased_low_base.pdf"); // HH
 // if(!HH && low && base && !twobtag) c1->SaveAs("/WP4_cutbased_low_base.root"); // HH
 // if(!HH && !low && base) c1->SaveAs("/WP4_cutbased_base.png");
 // if(!HH && !low && base) c1->SaveAs("/WP4_cutbased_base.pdf");
 // if(!HH && !low && base) c1->SaveAs("/WP4_cutbased_base.root");
  //return c1;

  delete defaultStyle;
  /*delete leg;
  delete leg1;
  delete ntuple;
  delete ntupleg;
  delete ntuplegrs;
  delete pt;
  delete Purity;
  delete BF;
  delete mg;
  delete T;
  delete c1;
  delete radion10;
  delete radion;
  delete rsgrav;
  delete bulk3;
  delete grobs;
  delete grgreen;
  delete gryellow;
  delete gcheck1;
  delete gcheck2;
  delete gcheck4;
  delete gcheck7;
  delete gcheck3;
  delete gcheck5;*/
  //for(unsigned int i=0;i!=TTrees.size();++i) delete TTrees[i];
  //for(unsigned int i=0;i!=TFileVec.size();++i)delete TFileVec[i];
}
