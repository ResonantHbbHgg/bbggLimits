#define bbggFittingTools_cxx
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TString.h"
#include "TChain.h"
#include "TLegend.h"
#include "RooMinuit.h"
#include <TStyle.h>
#include "TH1F.h"
#include "RooChi2Var.h"
#include "RooFitResult.h"
#include "RooCategory.h"
//#include "RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "TFile.h"
#include <vector>
#include <iostream>
#include <string>
#include "TObjArray.h"
#include "HiggsAnalysis/bbggLimits/interface/bbggFittingTools.h"

using namespace RooFit;

int colors[] = {kBlue, 632, 417, 616, 432, 800, 820, 840, 860};
int styles[] = {2, 3, 0, 5, 6, 7, 8, 9};

void bbggFittingTools::PlotCurves(std::string plotTitle, RooWorkspace* w, std::vector<std::string> functionsToFit, std::vector<std::string> legends, std::vector<bbggFittingTools::FitRes> fitResults, RooDataSet* data, std::string sVar, std::string nbins, std::string plotName, int error, int isBkg, float totNorm)
{
    RooPlot* frame = w->var(sVar.c_str())->frame(Bins(TString( nbins ).Atoi() ));
    data->plotOn(frame, DataError(RooAbsData::SumW2));
    int funcCounter = 0;

//    TStyle * gStyle = new TStyle("a", "a");
//    gStyle->SetLegendTextSize(0.03);
    TLegend* leg = new TLegend(0.5, 0.75, 0.89, 0.89);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetBorderSize(0);
    for ( unsigned int ff = 0; ff < functionsToFit.size(); ff++) {
 	std::cout << "PRINTING FUNCTION: " << functionsToFit[ff] << std::endl;       
        TObjArray* funcNames = (TObjArray*) TString(functionsToFit[ff]).Tokenize(":");
        const char* modelName = ((TObjString*) funcNames->At(0) )->String().Data();
        
        if( TString(modelName).Contains(TString(sVar)) == 0) continue;
        
        TObjArray* sComponents = 0;
        if(funcNames->GetEntries() > 1)
            sComponents  = (TObjArray*) ((TObjString*) funcNames->At(1) )->String().Tokenize(",");

        std::vector<const char*> components;
        if(sComponents != 0){
            for (int comp = 0; comp < sComponents->GetEntries(); comp++) {
                components.push_back( ((TObjString*) sComponents->At(comp) )->String().Data() );
            }
        }
        std::cout << "Fitting model: " << modelName << std::endl;
        if(components.size() > 0) std::cout << "\t including " << components.size() << " components" << std::endl;
        
        RooFitResult* fitResult = 0;
        for ( unsigned int fr = 0; fr < fitResults.size(); fr++){
            if ( fitResults[fr].function ==  modelName ){
                fitResult = fitResults[fr].result;
                break;
            }
        }
        //If there's only one, then it's signal
        if( fitResults.size() == 1)
		fitResult = fitResults[0].result;

        if(!fitResult && error){
            std::cout << "Problem finding fit result!" << std::endl;
            return;
        }
        
        w->pdf( modelName )->plotOn(frame, LineColor( colors[funcCounter] ), LineStyle(styles[funcCounter]), Name(modelName), MoveToBack(), Precision(0.00001));
        
        if(error == 1){
            w->pdf( modelName )->plotOn(frame, FillColor( colors[funcCounter]-9 ), VisualizeError(*fitResult, 2, kFALSE));
            w->pdf( modelName )->plotOn(frame, FillColor( colors[funcCounter]-168 ), VisualizeError(*fitResult, 1, kFALSE));
            w->pdf( modelName )->plotOn(frame, LineColor( colors[funcCounter] ), Name(modelName));
            data->plotOn(frame);
        }
        
        for (unsigned int comp = 0; comp < components.size(); comp++){
            std::cout << "Plotting model component: " << components[comp] << std::endl;
            w->pdf( modelName )->plotOn(frame, Components( components[comp] ), LineColor( colors[funcCounter] ), LineStyle( styles[comp] ));
        }
	funcCounter++;

        leg->SetHeader(plotTitle.c_str());
	leg->AddEntry( frame->findObject( modelName), legends[ff].c_str(), "l");
        
    }
    
    TCanvas* c = new TCanvas("a", "a", 1200, 1000);
    frame->SetTitle("");
    frame->Draw();
    if(TString(sVar).Contains("mgg") && !isBkg){
        frame->GetXaxis()->SetRangeUser(115, 135);
    }
    leg->Draw("same");
    c->SaveAs(TString(plotName) + ".pdf");
    c->SaveAs(TString(plotName) + ".png");
}

//std::vector<bbggFittingTools::FitRes> bbggFittingTools::FitFunctions(RooWorkspace* w, std::vector<std::string> functionsToFit, RooDataSet* data, std::string sVar, TH1F* dataHist)
std::vector<bbggFittingTools::FitRes> bbggFittingTools::FitFunctions(RooWorkspace* w, std::vector<std::string> functionsToFit, RooDataSet* data, int isSignal)
{
    std::vector<bbggFittingTools::FitRes> fitresults;
    
    for ( unsigned int ff = 0; ff < functionsToFit.size(); ff++) {
        TObjArray* funcNames = (TObjArray*) TString(functionsToFit[ff]).Tokenize(":");
        const char* modelName = ((TObjString*) funcNames->At(0) )->String().Data();

        TObjArray* sComponents = 0;
        if(funcNames->GetEntries() > 1)
            sComponents  = (TObjArray*) ((TObjString*) funcNames->At(1) )->String().Tokenize(",");

        std::vector<const char*> components;
        if(sComponents != 0){
            for (int comp = 0; comp < sComponents->GetEntries(); comp++) {
                components.push_back( ((TObjString*) sComponents->At(comp) )->String().Data() );
            }
        }
        std::cout << "Fitting model: " << modelName << std::endl;
        if(components.size() > 0) std::cout << "\t including " << components.size() << " components" << std::endl;

        if(w->pdf(modelName) == 0 ){
            std::cout << "Model " << modelName << " not found in workspace! Are you sure you have added it to the list of funtions?" << std::endl;
            continue;
        }
	
	RooFitResult * fitResult = 0;
	if(!isSignal){
        	RooAbsReal* nll = w->pdf( modelName )->createNLL(*data, NumCPU(8));
        	RooMinuit minuit(*nll);
        	minuit.migrad();
        	minuit.hesse();
//                minuit.minos();
        	fitResult = minuit.save();
	}
	if(isSignal){
		std::cout << "Fitting signal!" << std::endl;
//		fitResult = w->pdf( modelName )->fitTo(*data, SumW2Error(kTRUE), Save(kTRUE), Strategy(2), Minos(kTRUE));
		fitResult = w->pdf( modelName )->fitTo(*data, SumW2Error(kTRUE), Save(kTRUE) );
	}

        int nParams = w->pdf( modelName )->getParameters(data)->getSize();

        RooChi2Var Chi2 ("chi2", "chi2", *(w->pdf( modelName )), *(data->binnedClone()));
        double chi2_val = Chi2.getVal(); 
        double finalchi2_val = chi2_val/((float) nParams);
        double minNLL = fitResult->minNll();
        
        bbggFittingTools::FitRes thisResult;

        thisResult.function = modelName;
//        thisResult.kolmogorov = kolmo;
        thisResult.chi2 = finalchi2_val;
        thisResult.minNLL = minNLL;
        thisResult.result = fitResult;
        fitresults.push_back(thisResult);
        
//        delete nll;
//        delete histo;
    }
    return fitresults;
}
