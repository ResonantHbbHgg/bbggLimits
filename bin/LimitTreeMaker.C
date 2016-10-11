#include "HiggsAnalysis/bbggLimits/interface/bbggLTMaker.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int Process(string file, string outFile, string mtotMin, string mtotMax, string cosThetaStar, string CTScats, string scale, 
	    string photonCR, string doKinFit, string doMX, string doTilt, string tiltWindow, 
	    string doNoCat, string btagWP, string doCatMixed, string btagHigh, string btagLow, string singleCat,
	    string doVariation, string doPhoVariation, string doNRWeights ){
    TFile* iFile = new TFile(TString(file), "READ");
    TTree* iTree = (TTree*) iFile->Get("bbggSelectionTree");
    cout << "[LimitTreeMaker:Process] Processing tree with " << iTree->GetEntries() << " entries." << endl;
    bbggLTMaker t(iTree);
    t.SetMax( TString(mtotMax).Atof() );
    t.SetMin( TString(mtotMin).Atof() );
    t.SetNormalization( TString(scale).Atof() );
    t.IsPhotonCR( TString(photonCR).Atoi() );
    t.IsMX( TString(doMX).Atoi() );
    t.IsKinFit( TString(doKinFit).Atoi() );
    t.SetOutFileName( outFile );
    t.SetTilt( TString(doTilt).Atoi());
    t.DoNoCat( TString(doNoCat).Atoi());
    t.SetBTagWP( TString(btagWP).Atof());
    t.DoCatMixed( TString(doCatMixed).Atoi());
    t.SetBTagWP_High( TString(btagHigh).Atof());
    t.SetBTagWP_Low( TString(btagLow).Atof());
    t.DoSingleCat( TString(singleCat).Atoi());
    t.DoBVariation( TString(doVariation).Atoi());
    t.DoPhoVariation( TString(doPhoVariation).Atoi());
    t.SetCosThetaStar( TString(cosThetaStar).Atof(), TString(CTScats).Atoi());
    t.DoNRWeights( TString(doNRWeights).Atoi());
    
    t.Loop();

    return 1;
}

int main(int argc, char *argv[]) {
    string outputFile = "";
    string inputFile = "";
    string outputLocation = "";
    string mtotMax = "10000";
    string mtotMin = "0";
    string scale = "1";
    string photonCR = "0";
    string doKinFit = "0";
    string doMX = "0";
    string doTilt = "0";
    string doNoCat = "0";
    string tiltWindow = "0";
    string btagWP = "0.8";
    string doCatMixed = "0";
    string btagHigh = "0.8";
    string btagLow = "0.435";
    string singleCat = "0";
    string inputRootFile = "";
    string doVariation = "-999";
    string doPhoVariation = "-999";
    string cosThetaStar = "100";
    string cosThetaStarCats = "0";

    string doNRWeights = "0";
    
    for( int i = 1; i < argc; i++){
	if ( std::string(argv[i]) == "-i"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		inputFile = string(std::string(argv[i+1]));
		i++;
	}
    	else if ( std::string(argv[i]) == "-inputFile") {
        	if ( (i+1) == argc) {
            		std::cout << "Invaliv number of arguments!" << std::endl;
            		break;
           	}
        	inputRootFile = string(std::string(argv[i+1]));
        	i++;
    	}
	else if ( std::string(argv[i]) == "-tilt"){
		doTilt = "1";
	}
	else if ( std::string(argv[i]) == "-tiltWindow"){
		if( (i+1) == argc){
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		tiltWindow = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-o"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		outputLocation = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-max"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		mtotMax = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-min"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		mtotMin = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-btagWP"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		btagWP = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-scale"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		scale = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-btagHigh"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		btagHigh = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-btagLow"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		btagLow = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-doBVariation"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		doVariation = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-doPhoVariation"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		doPhoVariation = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-cosThetaStar"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		cosThetaStar = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-cosThetaStarCats"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		cosThetaStarCats = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-photonCR"){
		photonCR = "1";
	}
    	else if ( std::string(argv[i]) =="-KF"){
        	doKinFit = "1";
    	}
    	else if ( std::string(argv[i]) =="-MX"){
        	doMX = "1";
    	}
	else if ( std::string(argv[i]) =="-doNoCat"){
		doNoCat = "1";
	}
	else if ( std::string(argv[i]) =="-doCatMixed"){
		doCatMixed = "1";
	}
	else if ( std::string(argv[i]) =="-singleCat"){
		singleCat = "1";
	}
	else if ( std::string(argv[i]) == "-NRW"){
	  doNRWeights = "1";
	}
	else {
		cout << "UNKNOWN OPTION! " << std::string(argv[i]) << endl;
		cout << "Usage: LimitTreeMaker -i <input list of files> ( or -inputFile <single root file> ) -o <output location>" << endl;// \n
		cout << "options: " << endl;//\n
		cout << "\t -min <min mtot> -max <max mtot> " << endl;//\n
		cout << "\t -scale <scale factor> " << endl;//\n
		cout << "\t -photonCR (do photon control region) " << endl;//\n
		cout << "\t -KF (use Mtot_KF to cut on mass window) " << endl;//\n
		cout << "\t -MX (use MX to cut on mass window) (choose either -MX or -KF!) " << endl;//\n
		cout << "\t -tilt (select tilted mass window) " << endl;//\n
		cout << "\t -doNoCat (dont cut on categories) " << endl;//\n
		cout << "\t -btagWP <WP> (set btagging working point for categories." << endl;
		cout << "\t -doCatMixed (do categories with mixed btagging - cat0: 2>low, cat1: 1<low+1>high)" << endl;
		cout << "\t -btagHigh (for mixed cat, highest value)" << endl;
		cout << "\t -btagLow (for mixed cat, lowest value)" << endl;
		cout << "\t -singleCat (only one category)" << endl;
		cout << "\t -doBVariation (Apply b-tagging SF factors: 1 or -1)" << endl;
		cout << "\t -NRW (add non-resonant weights)" << endl;
		return -1;
	}	
    }

    if( (inputFile == "" && inputRootFile == "") || outputLocation == "") {
		cout << "YOU DIDNT SPECIFY INPUT AND OR OUTPUT" << endl;
		cout << "Usage: LimitTreeMaker -i <input list of files> ( or -inputFile <single root file> ) -o <output location>" << endl;// \n
		cout << "options: " << endl;//\n
		cout << "\t -min <min mtot> -max <max mtot> " << endl;//\n
		cout << "\t -scale <scale factor> " << endl;//\n
		cout << "\t -photonCR (do photon control region) " << endl;//\n
		cout << "\t -KF (use Mtot_KF to cut on mass window) " << endl;//\n
		cout << "\t -MX (use MX to cut on mass window) (choose either -MX or -KF!) " << endl;//\n
		cout << "\t -tilt (select tilted mass window) " << endl;//\n
		cout << "\t -doNoCat (dont cut on categories) " << endl;//\n
		cout << "\t -btagWP <WP> (set btagging working point for categories." << endl;
		cout << "\t -doCatMixed (do categories with mixed btagging - cat0: 2>low, cat1: 1<low+1>high)" << endl;
		cout << "\t -btagHigh (for mixed cat, highest value)" << endl;
		cout << "\t -btagLow (for mixed cat, lowest value)" << endl;
		cout << "\t -singleCat (only one category)" << endl;
	return 0;
    }

    if(inputRootFile != "") {
        cout << "Processing file: " << inputRootFile << endl;
        if (inputRootFile.find(".root")==std::string::npos) {
            std::cout << inputRootFile << " is not a valid root file!" << std::endl;
            return 0;
        }
        size_t sLoc = inputRootFile.find_last_of("/");
        string outF = outputLocation;
        outF.append("/LT_");
        outF.append(inputRootFile.substr(sLoc+1));
        cout << "Output file: " << outF << endl;
        Process(inputRootFile, outF, mtotMin, mtotMax, cosThetaStar, cosThetaStarCats, scale, photonCR, doKinFit, doMX, doTilt, tiltWindow, doNoCat, btagWP, doCatMixed, btagHigh, btagLow, singleCat, doVariation, doPhoVariation, doNRWeights);
        return 0;
    }
    
    ifstream infile(inputFile);
    string line;
    while(getline(infile,line)){
        string rootFileName = "";
        cout << "Processing file: " << line << endl;
        string token;
        stringstream ss(line);
        while(getline(ss, token, '/')){
            if (token.find(".root")!=std::string::npos) rootFileName = token;
        }
        if(rootFileName.find(".root")==std::string::npos) {
            cout << "[LimitTreeMaker] The following line in file is not a root file: " << line << endl;
            continue;
        }
        size_t sLoc = line.find_last_of("/");
        string outF = outputLocation;
        outF.append("/LT_");
        outF.append(line.substr(sLoc+1));
        cout << "Output file: " << outF << endl;
        Process(line, outF, mtotMin, mtotMax, cosThetaStar, cosThetaStarCats, scale, photonCR, doKinFit, doMX, doTilt, tiltWindow, doNoCat, btagWP, doCatMixed, btagHigh, btagLow, singleCat, doVariation, doPhoVariation, doNRWeights);
    }
    
    return 0;
}

