#include<iostream>
#include<vector>
#include"HiggsAnalysis/bbggLimits/interface/Colors.h"
#include"HiggsAnalysis/bbggLimits/interface/BrazilianFlag.h"

int main(int argc, const char* argv[])
{
	bool HH=false;
	bool base=true;
	bool low=false;
	bool obs=false;
	bool twotag=false;
	std::string path_dir="";
	if(argc < 7)
	{
		std::cout <<red<< "Please provide :"<<normal << std::endl;
		std::cout <<red<< "* The path of the folder"<<normal << std::endl;
		std::cout <<red<< "* HH (true,false)"<<normal << std::endl;
		std::cout <<red<< "* base (true,false)"<<normal << std::endl;
		std::cout <<red<< "* low (true,false)"<<normal << std::endl;
		std::cout <<red<< "* obs (true,false)"<<normal << std::endl;
		std::cout <<red<< "* twotag (true,false)"<<normal << std::endl;
	}
	else
	{
		path_dir=argv[1];
		HH=argv[2];
		base=argv[3];
		low=argv[4];
		obs=argv[5];
		twotag=argv[6];
		BrazilianFlag(path_dir,HH,base,low,obs,twotag); 
	}
	return 0;
}
