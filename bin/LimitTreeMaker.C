#include "HiggsAnalysis/bbggLimits/interface/bbggLTMaker.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int Process(string file, string outFile, string mtotMin, string mtotMax,string scale, string photonCR, string doKinFit, string doMX, string doTilt, string tiltWindow, string doNoCat){
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
    t.SetBTagWP( 0.8 );
//    t.SetTilt( TString(doTilt).Atoi(), TString(tiltWindow).Atof());
    t.SetTilt( TString(doTilt).Atoi());
    t.DoNoCat( TString(doNoCat).Atoi());
    t.Loop();

//    delete iFile;
//    delete iTree;
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

    for( int i = 1; i < argc; i++){
	if ( std::string(argv[i]) == "-i"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		inputFile = string(std::string(argv[i+1]));
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
	else if ( std::string(argv[i]) =="-scale"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		scale = string(std::string(argv[i+1]));
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
	else {
		cout << "Usage: LimitTreeMaker -i <input list of files> -o <output location> [optional: -min <min mtot> -max <max mtot> -scale <scale factor> -photonCR (do photon control region) -KF (use Mtot_KF to cut on mass window) -MX (use MX to cut on mass window) (choose either -MX or -KF!)" << endl;
		return -1;
	}	
    }

    if(inputFile == "" || outputLocation == "") {
	cout << "Usage: LimitTreeMaker -i <input list of files> -o <output location> [optional: -min <min mtot> -max <max mtot> -scale <scale factor> -photonCR (do photon control region) -KF (use Mtot_KF to cut on mass window) -MX (use MX to cut on mass window) (choose either -MX or -KF!)" << endl;
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
        Process(line, outF, mtotMin, mtotMax, scale, photonCR, doKinFit, doMX, doTilt, tiltWindow, doNoCat);
    }
    
    return 0;
}

