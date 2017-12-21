/** 
 * \file main.cpp
 * \brief
 * Copyright G Gronoff Feb 2014
 *
 */


#include "main.hpp"
namespace bub = boost::numeric::ublas;
using namespace std;
int main(int argc,char** argv)
{
	cout<<"Welcome to the Photodissociation Branching Ratio computation version "<<VERSION<<endl;


	string logname="Log.txt";

	if(argc>2)
	{

		string suffix(argv[2]);
		logname+=suffix;
	}

	Log::Init(logname);


	//cout<<"hello world"<<endl;
	//cout<<"nombre de parametres: "<<argc<<" les parametres:"<<argv[0]<<endl;
	if(argc<2)
	{
		//	Log::SetPriority(Log::INFO);
		Log::mS<<"A control file is needed: please use the correct xml file"<<endl;
		return 1;
	}
	Log::mS<<"We will use your file "<<argv[1]<<endl;



	try{

		//	Log::SetPriority(Log::INFO,"main()");
		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We set up the atmosphere"<<endl;
		Log::mS<<"____________________________________________"<<endl;
		Atmo atmosphere(argv[1]);


		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We proceed the computation and writing "<<endl;
		Log::mS<<"____________________________________________"<<endl;

		if(argc>2)
		{
			string suffix(argv[2]);
			atmosphere.ComputePhotodissociation(suffix);
		}else
		{
			atmosphere.ComputePhotodissociation();
		}
	}
	catch(Error &err)
	{
		err.Affiche();
		return 1;
	}
	catch(MathError &err)
	{
		err.Affiche();
		return 1;
	}
	catch(...)
	{
		//Log::SetPriority(Log::ERROR,"main()");
		Log::mE<<"OUUUUPS, unknown error"<<endl;
	}

	Log::mS<<"Program exited correctly!!!!!"<<endl;
	cout<<"bye"<<endl;
	return 0;
}


