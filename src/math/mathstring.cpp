/**
 * \file mathstring.cpp
 * \brief implements mathematical functions, mainly associated with string retrieval.
 * E.G. if you  want to read an array from a file...
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathstring.cpp 1526 2012-07-03 20:44:54Z gronoff $
 *
 *
 */



#include "mathstring.hpp"
using namespace std;
using namespace boost;

namespace mathessentials
{
	double modulo(double vA, double vB)
	{
		double tmp=fmod(vA,vB);
		if(tmp<0)
			tmp+=vB;
		return tmp;
	}



	string trim(string vText,string vChars)
	{
		if("" == vText)
			return vText;
		regex my_regex;
		string my_reg="(^";
		my_reg+=vChars+"|"+vChars+"$)";

		try
		{
			my_regex=my_reg;
		}
		catch(boost::regex_error& err_reg)
		{
			cout<<"L'expression rationnelle /"<<my_reg<<"/ n'est pas valide"<<endl;
			exit(1);
			return vText;
		}

		vText=regex_replace(vText,my_regex,"",boost::match_default | boost::format_all);


		return vText;

	}



	std::string StrReplace(std::string vStr,std::string vAig,std::string vRep)
	{
		string resu=vStr;
		size_t i=resu.find(vAig,0);
		while(i!=string::npos)
		{
			resu.replace(i,vAig.size(),vRep);
			i=resu.find(vAig,i+1);

		}
		return resu;



	}


	void NoNegative(ublas::vector<double>& vMat)
	{
		for(size_t i=0;i<vMat.size();++i)
		{
			if(vMat[i]<0)
				vMat[i]=0;
		}
	}
	void MinValue(ublas::vector<double>& vMat,double vMin)
	{
		for(size_t i=0;i<vMat.size();++i)
		{
			if(vMat[i]<vMin)
				vMat[i]=vMin;
		}
	}
};
