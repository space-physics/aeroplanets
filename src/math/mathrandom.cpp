/**
 * \ingroup Math 
 * \file mathrandom.cpp
 * \brief Implements MathRandom
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathrandom.cpp 1205 2011-01-20 03:33:04Z gronoff $
 *
 */

#include "mathrandom.hpp"

using namespace std;

boost::mt19937 MathRandom::msRng;
unsigned long MathRandom::msSeed=0;
bool MathRandom::mIsInitialized=false;

void MathRandom::Initialize()
{
	MathRandom::mIsInitialized=true;

	ifstream in("/dev/random",ios::in|ios::binary);

	//Log::SetPriority(Log::CONFIG);
	if(in)
	{
		in.read((char*) &MathRandom::msSeed,sizeof MathRandom::msSeed);
		in.close();
	}else
	{
		MathRandom::msSeed=static_cast<long unsigned int>(clock());
		Log::mL<<"Your Random seed is initialized by the clock"<<endl;
		Log::mL<<"It means that is is a very bad idea to launch it twice at the same time. (Bye bye your super multicore, you had to use A REAL OS!"<<endl;
	}
	Log::mD<<"Your Random SEED : "<<MathRandom::msSeed<<endl;

	MathRandom::msRng.seed( static_cast<boost::mt19937::result_type>(MathRandom::msSeed));
}



double MathRandom::GetNormal(double vMean,double vSigma)
{// We suppose that the seed and msRng are correctly defined.

	if(!MathRandom::mIsInitialized)
	{
		MathRandom::Initialize();
	}

	boost::normal_distribution<> ndistro(vMean,vSigma);
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > nrnd(MathRandom::msRng,ndistro);
	return nrnd();
}

double MathRandom::GetUniformReal(double vMin,double vMax)
{// We suppose that the seed and msRng are correctly defined.

	if(!MathRandom::mIsInitialized)
	{
		MathRandom::Initialize();
	}
	boost::uniform_real<> ndistro(vMin,vMax);
	boost::variate_generator<boost::mt19937&,boost::uniform_real<> > nrnd(MathRandom::msRng,ndistro);
	return nrnd();
}

