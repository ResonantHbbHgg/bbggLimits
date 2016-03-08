#include<iostream>
#include<vector>
#include"HiggsAnalysis/bbggLimits/interface/Colors.h"
#include "HiggsAnalysis/bbggLimits/interface/RunCombine.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/optional/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
bool ConvertToBool(const char* a)
{
	if(std::strcmp(a, "true") == 0||std::strcmp(a, "1") == 0||std::strcmp(a, "True") == 0) return true;
	else return false;
	return false;
}

int main(int argc, const char* argv[])
{	
	bool doblinding=true;
	std::string path_dir="";
	if(argc < 2) 
	{
		std::cout <<red<< "Please provide the json file and the directory you want to use"<<normal << std::endl;
		std::exit(2);
	}
	else if(argc==3)
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
      				doblinding = rowPair.second.get<bool>("doBlinding");
      			}
		}
	}
        path_dir=argv[2];
	RunCombine(path_dir,doblinding);
	return 0;
}
