#include<iostream>
#include<vector>
#include"HiggsAnalysis/bbggLimits/interface/Colors.h"
#include "HiggsAnalysis/bbggLimits/interface/RunCombine.h"
int main(int argc, const char* argv[])
{	
	bool doblinding=true;
	std::string path_dir="";
	if(argc < 3) 
	{
		std::cout <<red<< "Please provide the folder name as first argument and the blinding flag"<<normal << std::endl;
	}
	else
	{
		path_dir=argv[1];
		doblinding=argv[2];
		RunCombine(path_dir,doblinding);
	}
	return 0;
}
