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
int main(int argc, const char* argv[])
{	
	bool doblinding=true;
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
      				doblinding = rowPair.second.get<bool>("doBlinding");
      			}
		}
	}
	else
	{
		std::cout <<red<< "Please provide the json file and the directory you want to use"<<normal << std::endl;
		std::exit(2);
	}
        path_dir=argv[2];
	RunCombine(path_dir,doblinding);
	return 0;
}
