#include<string>
#include<iostream>
#include <boost/filesystem.hpp>
#include "HiggsAnalysis/bbggLimits/interface/Colors.h"
#include<boost/regex.hpp>


void RunCombine(std::string path_dir,bool doblinding)
{
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
            if (folder_name == "SeparateCategories") continue;
            const char * pattern = "\\d+";
            boost::regex re(pattern);
            boost::sregex_iterator it(folder_name.begin(), folder_name.end(), re);
            boost::sregex_iterator end;
            std::vector<std::string>number;
            for( ; it != end; ++it)
            {
                number.push_back(it->str()); 
            }
			std::string blinded="";
			if(doblinding==1)blinded="--run blind ";
//			std::string txtname="hgg.mH125_8TeV";
			std::string txtname="hhbbgg_13TeV_DataCard";
		 	std::string location=path_dir+"/"+folder_name+"/datacards/";
		 	std::string logs=path_dir+"/"+folder_name+"/logs";
			std::string cresult=path_dir+"/"+folder_name+"/combine"; 
			std::string mvv="mv higgsCombineTest.Asymptotic.mH120.root "+cresult+"/higgsCombineTest.Asymptotic.mH125_m"+number[0];        				
			boost::filesystem::create_directory(logs);
			boost::filesystem::create_directory(cresult);
		  std::string command1="combine -M Asymptotic "+blinded+location+txtname+".txt --X-rtd TMCSO_AdaptivePseudoAsimov=50 >>"+logs+"/higgsCombineTest.Asymptotic.mH125.0._m"+number[0]+"_higgs.txt";
			//std::string command2="combine -M Asymptotic "+blinded+location+txtname+"onecatnohiggs.txt >>"+logs+"/higgsCombineTest.Asymptotic.mH125.0._m"+number[0]+"_onecatnohiggs.txt";
			//std::string command3="combine -M Asymptotic "+blinded+location+txtname+".txt -S 0 >>"+logs+"/higgsCombineTest.Asymptotic.mH125.0._m"+number[0]+"_nosyst_higgs.txt";
			std::cout<<green<<"RUNNING COMBINE FOR "<<folder_name<<" : "<<normal<<std::endl;
      			std::cout<<command1<<std::endl;
			std::system(command1.c_str());
			std::string mvvv=mvv+"_higgs.root";
			std::system(mvvv.c_str());
      			//std::cout<<command2<<std::endl;                                    			                        		        			
			//std::system(command2.c_str());
			//mvvv=mvv+"_onecatnohiggs.root";
			//std::system(mvvv.c_str());
			//std::cout<<command3<<std::endl;						              			                        		        						
			//std::system(command3.c_str());
			//std::system("mv higgsCombineTest.Asymptotic.mH120.root"+cresult+"/higgsCombineTest.Asymptotic.mH125_nosyst_higgs.root")
		}
	}
}
