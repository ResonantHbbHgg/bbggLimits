#include<iostream>
#include<vector>
#include"HiggsAnalysis/bbggLimits/interface/Colors.h"
#include"HiggsAnalysis/bbggLimits/interface/BrazilianFlag.h"
#include<cstring>
bool ConvertToBool(const char* a)
{
	if(std::strcmp(a, "true") == 0||std::strcmp(a, "1") == 0||std::strcmp(a, "true") == 0) return true;
	else return false;
	return false;
}
int main(int argc, const char* argv[])
{
	bool HH=true;
	bool base=true;
	bool low=false;
	bool obs=false;
	bool twotag=false;
        std::string energy="";
	float lumi=0.0;
	std::string path_dir="";
	if(argc < 9)
	{
		std::cout <<red<< "Please provide :"<<normal << std::endl;
		std::cout <<red<< "* The path of the folder"<<normal << std::endl;
		std::cout <<red<< "* HH (true,false)"<<normal << std::endl;
		std::cout <<red<< "* base (true,false)"<<normal << std::endl;
		std::cout <<red<< "* low (true,false)"<<normal << std::endl;
		std::cout <<red<< "* obs (true,false)"<<normal << std::endl;
		std::cout <<red<< "* twotag (true,false)"<<normal << std::endl;
		std::cout <<red<< "* energy"<<normal << std::endl;
		std::cout <<red<< "* lumi"<<normal << std::endl;
		std::exit(1);
	}
	else
	{
		path_dir=argv[1];
		HH=ConvertToBool(argv[2]);
		base=ConvertToBool(argv[3]);
		low=ConvertToBool(argv[4]);
		obs=ConvertToBool(argv[5]);
		twotag=ConvertToBool(argv[6]);
                energy=argv[7];
		lumi=std::atof(argv[8]);
		BrazilianFlag(path_dir,HH,base,low,obs,twotag,energy,lumi); 
	}
	return 0;
}
