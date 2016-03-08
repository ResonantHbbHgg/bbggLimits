#include<iostream>
#include<vector>
#include"HiggsAnalysis/bbggLimits/interface/Colors.h"
#include"HiggsAnalysis/bbggLimits/interface/BrazilianFlag.h"
#include<cstring>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/optional/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
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
	if(argc==3)
	{
		if ( !boost::filesystem::exists( argv[1] ) || !boost::filesystem::exists( argv[1] ) )
		{
			std::cout<<red<<"File "<<argv[1]<<" or/and "<<argv[2]<<" doesn't exist "<<normal<<std::endl;
			std::exit(2);
		}
		boost::property_tree::ptree pt;
  		boost::property_tree::read_json( argv[1], pt );
  		std::cout <<green<< "Reading input configuration file..."<<normal << std::endl;
  		BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "" ) )
  		{
    			if (rowPair.first == "other") 
    			{
      				HH = rowPair.second.get<bool>("HH");
				base = rowPair.second.get<bool>("base");
				low = rowPair.second.get<bool>("low");
				obs = rowPair.second.get<bool>("obs");
				twotag = rowPair.second.get<bool>("twotag");
				lumi = rowPair.second.get<float>("integratedLumi");
      				energy = rowPair.second.get<std::string>("energy");
      			}
		}
	}
	else
	{
		std::cout <<red<< "Please provide the json file and the directory you want to use"<<normal << std::endl;
		std::exit(2);
	}
	std::cout<<"Running BrazilianFlag with path_dir:"<<path_dir<<" HH:"<<HH<<" base:"<<base<<" low:"<<low<<" obs:"<<obs<<" twotag:"<<twotag<<" energy:"<<energy<<" lumi:"<<lumi<<normal<<std::endl;
	path_dir=argv[2];
	BrazilianFlag(path_dir,HH,base,low,obs,twotag,energy,lumi); 
	return 0;
}
