/** 
 * \file paramfit.cpp
 * \brief For fitting the atmosphere with pre defined functions and limb observations
 * Copyright G Gronoff July 2010
 * Last Modification $Id$
 *
 */
#include "paramfit.hpp"
using namespace std;

namespace bub = boost::numeric::ublas;

int main(int argc,char** argv)
{
	cout<<"Welcome to AeroFitParam version "<<VERSION<<endl;


	string logname="Log.txt";

	if(argc>2)
	{
		string suffix(argv[2]);
		logname+=suffix;
	}
	Log::Init(logname);
	if(argc<2)
	{
		Log::mS<<"A control file is needed: please use the correct xml file"<<endl;
		return 1;
	}
	Log::mS<<"We will use your file "<<argv[1]<<endl;


	try
	{
		XmlParameters parameter(argv[1]);

		Log::mI<<"Check aerofile"<<endl;
		parameter.ExistsOrDie("/aerofit/aerofile","We need the aerofile to be able to compute your comparison data ");

		string aerofile=parameter.GetSubFileName("/aerofit/aerofile/");
		Atmo* atmosphere = new Atmo(aerofile);


		// The Monte Carlo System is activated after the call to the subfile
		Log::mI<<"Check if MC is activated"<<endl;
		if(parameter.Exists("/aerofit/SetMonteCarloActive"))
		{
			Log::mI<<"Monte Carlo analysis activated!!!"<<endl;
			parameter.SetMonteCarloActive();
		}




		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We launch the computation"<<endl;
		Log::mS<<"____________________________________________"<<endl;
		//atmosphere.Compute(); // here, we have to change for a most adapted model

		atmosphere->InitMultipleComp();// Initialise the objects for multiple computations


		Adjustements adjust(atmosphere,&parameter);
		Log::mS<<"End of the initializations..."<<endl;
		/*
		// test: launch a multiple computation to check if it is correctly working
		atmosphere->MCompute();
		atmosphere->MCompute();
		atmosphere->MCompute();
		atmosphere->MCompute();

		std::deque<double>chapparams;
		chapparams.push_back(120.);
		chapparams.push_back(10000.);
		chapparams.push_back(20.);

		atmosphere->ResetElectronDens(1,chapparams);
		//atmosphere->ResetSpecieDens("CO2",1,chapparams);

		atmosphere->MCompute();
		*/
		adjust.StartAdjustements();




		Log::mS<<"____________________________________________"<<endl;
		Log::mS<<"We proceed the outputs "<<endl;
		Log::mS<<"____________________________________________"<<endl;


		if(argc>2)
		{
			string suffix(argv[2]);
			Log::mI<<"You have a second argument : "<<suffix<<endl;
			Log::mI<<suffix<<" is used now as a prefix for the output"<<endl;
			atmosphere->ProceedEmissions(suffix);
			atmosphere->ProceedOutputs(suffix);
			adjust.PrintReport(suffix);

		}else
		{
			atmosphere->ProceedEmissions();
			atmosphere->ProceedOutputs();
			adjust.PrintReport();
		}

		atmosphere->FinishMultipleComp();// Destroy the objects of multiple computations
		delete atmosphere;



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
		Log::mE<<"OUUUUPS, unknown error"<<endl;
	}

	Log::mS<<"Program exited correctly!!!!!"<<endl;
	cout<<"bye"<<endl;
	return 0;
}




