/** 
 * \file main.cpp
 * \brief
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: main.cpp 1595 2012-11-03 01:30:26Z gronoff $
 *
 */


#include "main.hpp"
namespace bub = boost::numeric::ublas;
using namespace std;





/**
  Reads the files and start the computation
	\todo 	
	      - 3)  FLUID MODEL (easy! No, just kidding)

	\todo logarithmic uncertainties ??? 
  */

int main(int argc,char** argv)
{
	cout<<"Welcome to Aeroplanets version "<<VERSION<<endl;


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


	//	Log::SetPriority(Log::INFO,"main()");
		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We launch the computation"<<endl;
		Log::mS<<"____________________________________________"<<endl;
		atmosphere.Compute();

	//	Log::SetPriority(Log::INFO,"main()");
		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We proceed the outputs "<<endl;
		Log::mS<<"____________________________________________"<<endl;

		if(argc>2)
		{
			string suffix(argv[2]);
			Log::mL<<"You have a second argument : "<<suffix<<endl;
			Log::mL<<suffix<<" is used now as a prefix for the output"<<endl;
			atmosphere.ProceedEmissions(suffix);
			atmosphere.ProceedOutputs(suffix);

		}else
		{
			atmosphere.ProceedEmissions();
			atmosphere.ProceedOutputs();
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


