/**
 * \file protoncrs.cpp
 * \brief Implements the proton cross section class
 * Copyright G Gronoff November 2011
 * Last Modification : $Id$
 */

#include "protoncrs.hpp"
using namespace std;

/*/////////////////////////////////
 * The initializations
 */////////////////////////////////


ProtCrossSection::ProtCrossSection(ublas::vector<double> vGreV, std::string vParticleName): CrossSection(vGreV), mParticleName(vParticleName)
{
	mbIsElossElasComputed=false;
}


/*/////////////////////////////////
 * Load Cross Sections
 */////////////////////////////////


bool ProtCrossSection::LoadCrs(std::string vFilename, std::string vName, bool bIsMonteCarlo, bool bLogInterp)
{

	mpParam=XmlParameters::AttachFileParameter(vFilename);
	mParamLoaded=true;
	mSpecies=vName;
	bool resu=LoadCrs(mpParam, vName, bIsMonteCarlo, bLogInterp);
	return resu;
}

bool ProtCrossSection::LoadCrs(XmlParameters *pParams, std::string vName, bool bIsMonteCarlo, bool bLogInterp)
{
	Log::mL<<"Load electron impact cross section for"<<vName<<endl;
	string name=trim(vName);
	if(!	CrossSection::LoadCrs(pParams,name,bIsMonteCarlo,bLogInterp))
	{
		return false;
	}
	if(pParams->Exists("/crs/"+name+"/ZeroCrs"))
	{
		mIsDefinedCrs=false;
		return true;
	}
	// We check for the inelastic cross section
	if(!pParams->Exists("/crs/"+name+"/ElasticCrs"))
	{
		Error err("ElecCrossSection::LoadCrs","ElasticCrs","Your elastic cross section has not been defined");
		throw err;
	}
	mbIsDefinedElastic=true;

	if(pParams->Exists("/crs/"+name+"/ElasticCrs/use_hydrogen_function"))
	{
		HydrogenElasticFunction();
	}else if(pParams->Exists("/crs/"+name+"/ElasticCrs/use_proton_function"))
	{
		ProtonElasticFunction();
	}else if(pParams->Exists("/crs/"+name+"/ElasticCrs/zero"))
	{

		mElasticCrscm2.resize(mGrideV.size());
		mElasticCrscm2.clear();
	}
	else
	{
		TiXmlNode* elasticnode=pParams->GetNode("/crs/"+name+"/ElasticCrs");
		if(bLogInterp)
		{// nb no threshold for elastic cross section
			mElasticCrscm2=ExtractCrsLog(pParams,elasticnode,0.,mMaxDefinedEnergyElasticCrseV);

		}else
		{
			mElasticCrscm2=ExtractCrs(pParams,elasticnode,0.,mMaxDefinedEnergyElasticCrseV);
		}
	}
	// We check if we have excitation cross sections

	if(mbTotalSumIneel)
	{
		mIsDefinedTotalCrs=true;
		mTotalCrscm2.resize(mElasticCrscm2.size());
		mTotalCrscm2.clear();
		mTotalCrscm2=mElasticCrscm2;

	}
	mTotIonizationCrscm2.resize(mElasticCrscm2.size());
	mTotIonizationCrscm2.clear();
	mTotExcitationCrscm2.resize(mElasticCrscm2.size());
	mTotExcitationCrscm2.clear();
	mTotExchangeCrscm2.resize(mElasticCrscm2.size());
	mTotExchangeCrscm2.clear();
	assert(mExcitationThresholdeV.size()==mExcitationCrsPosition.size());//mExcitationCrscm2.size());

	for(unsigned j=0;j<mExcitationCrsPosition.size();++j)
	{
		unsigned k=mExcitationCrsPosition[j];
		//		Log::mL<<"Super important pour le check de k "<<k<<endl;
		mExcitationCrscm2.push_back(&mCrscm2[k]);
		if(mbTotalSumIneel)
		{
			mTotalCrscm2+=mCrscm2[k];
		}
		mTotExcitationCrscm2+=mCrscm2[k];
	}
	if(mExcitationCrscm2.size()==0)
	{// Nothing was defined : error
		Error err("ProtCrossSection::LoadCrs","No excitation cross section","The code was unable to find excitation cross section (even null). It is considered as an error");
		throw err;
	}

	Log::mD<<"C est partit pour une bonne ionization"<<endl;

	assert(mIonizationThresholdeV.size()==mIonizationCrsPosition.size());

	for(unsigned j=0;j<mIonizationCrsPosition.size();++j)
	{
		unsigned k=mIonizationCrsPosition[j];
		mIonizationCrscm2.push_back(&(mCrscm2[k]));
		if(mbTotalSumIneel)
		{
			mTotalCrscm2+=mCrscm2[k];
		}
		mTotIonizationCrscm2+=mCrscm2[k];
	}


	if(mIonizationCrscm2.size()==0)
	{// Nothing was defined : error
		Error err("ProtCrossSection::LoadCrs","No ionization cross section","The code was unable to find ionization cross section (even null). It is considered as an error");
		throw err;
	}

	assert(mExchangeThresholdeV.size()==mExchangeCrsPosition.size());//mExcitationCrscm2.size());

	for(unsigned j=0;j<mExchangeCrsPosition.size();++j)
	{
		unsigned k=mExchangeCrsPosition[j];
		//		Log::mL<<"Super important pour le check de k "<<k<<endl;
		mExchangeCrscm2.push_back(&mCrscm2[k]);
		if(mbTotalSumIneel)
		{
			mTotalCrscm2+=mCrscm2[k];
		}
		mTotExchangeCrscm2+=mCrscm2[k];
	}
	if(mExchangeCrscm2.size()==0)
	{// Nothing was defined : error
		Error err("ProtCrossSection::LoadCrs","No exchange cross section","The code was unable to find exchange cross section (even null). It is considered as an error");
		throw err;
	}


	////////////////////////////
	//
	// mWexcitation = minimum mExcitationThresholdeV

	mWexcitation = * (min_element(mExcitationThresholdeV.begin(), mExcitationThresholdeV.end()));

	//
	// mW10 = minimum mExchangeThresholdeV


	mW10 = * (min_element(mExchangeThresholdeV.begin(), mExchangeThresholdeV.end()));

	//
	//
	// vion = minimum mIonizationThresholdeV
	double vion =  * (min_element(mIonizationThresholdeV.begin(), mIonizationThresholdeV.end()));
	mVion = vion;
	mW10 = vion - mW10;
	double alpha = 0.91;
	// mElossIoni = vion  + pow(vion * mGrideV[i] * GAMMA_PROTON), 0.5/alpha // similar for both proton and H
	// mElossHioni = mW10 + 1 + mGrideV[i] * GAMMA_PROTON // in this context (H) W10 does not have any sense; in the p context, mElossHioni does not have any sense

	mElossIoni.resize(mGrideV.size());
	mElossIoni.clear();
	mElossHioni.resize(mGrideV.size());
	mElossHioni.clear();

	for(size_t i = 0; i < mGrideV.size(); ++i)
	{
		mElossIoni[i] = vion  + sqrt(vion * mGrideV[i] * ELECTRON_PROTON_MASS_RATIO)/ alpha; // similar for both proton and H
		mElossHioni[i] = mW10 + 1 + mGrideV[i] * ELECTRON_PROTON_MASS_RATIO; // in this context (H) W10 does not have any sense; in the p context, mElossHioni does not have any sense
	}

	// mElossElas -> depends on the species mass and the number of angles... Will be computed in the initialization of the proton model.
	//
	return true;
}


void ProtCrossSection::PrintCrs(std::string vFilename,std::string vRedisFilename)
{
	CrossSection::PrintCrs(vFilename);


	if(vRedisFilename != "")
	{
		// We do the same processus for the redistribution


		Log::mL<<"Redistribution file name: "<<vRedisFilename<<endl;
		if(FileExists(vRedisFilename))
		{
			Log::mD<<"Low level warning : We overwrite"<<vRedisFilename<<endl;
		}

		ofstream orf(vRedisFilename.c_str());
		orf.precision(9);
		orf.setf(ios::scientific);
		orf<<Log::msMessageLog<<endl;
		orf.close();
	}
}



void ProtCrossSection::HydrogenElasticFunction()
{
	double a[5];
	a[0]=4.183;
	a[1]=-0.39348;
	a[2]=0.65156E-2;
	a[3]=0.30656E-2;
	a[4]=-0.59441E-3;
	mElasticCrscm2.resize(mGrideV.size());
	mElasticCrscm2.clear();

	for(size_t i = 0; i < mGrideV.size(); ++i)
	{
		double kd = 0;
		if(mGrideV[i] < 4E5)
		{
			for(size_t j = 0; j < 5; ++j)
			{
				kd += a[j] * pow(log(mGrideV[i]*1E-3), static_cast<double>(j));
			}
			mElasticCrscm2[i] = 1E-16*0.529*0.529*exp(kd);
		}else
		{
			mElasticCrscm2[i] = 1E-16*0.529*0.529*exp(6.147 - 0.7 * log(mGrideV[i]*1E-3));
		}
	}

}

void ProtCrossSection::ProtonElasticFunction()
{
	double a[5];
	a[0]=5.3746;
	a[1]=-0.20842;
	a[2]=0.40561E-2;
	a[3]=-0.33036E-2;
	a[4]=-0.67680E-3;
	mElasticCrscm2.resize(mGrideV.size());
	mElasticCrscm2.clear();

	for(size_t i = 0; i < mGrideV.size(); ++i)
	{
		double kd = 0;
		if(mGrideV[i] < 4E5)
		{
			for(size_t j = 0; j < 5; ++j)
			{
				kd += a[j] * pow(log(mGrideV[i]*1E-3), static_cast<double>(j));
			}
			mElasticCrscm2[i] = 1E-16*0.529*0.529*exp(kd);
		}else
		{
			mElasticCrscm2[i] = 1E-16*0.529*0.529*exp(8.658 -  log(mGrideV[i]*1E-3));
		}
	}

}



void ProtCrossSection::ComputeElasticLoss(const double& vSpeciesMass,const MathFunction::GaussianAngle& vGa, bool vbIsRedis, bool vbForce)
{
	// We do not compute again if it is done, and not forced
	if(mbIsElossElasComputed and ! vbForce)
		return;

	ublas::matrix<double> tmp(vGa.mNbAngles, vGa.mNbAngles);
	tmp.clear();
	mElossElas.resize(mGrideV.size());


	if(vbIsRedis)
	{ // if the redistribution is taken into account: we compute it

		for(size_t i =0 ; i< static_cast<unsigned>(vGa.mNbAngles); ++i)
		{
			for(size_t j =0 ; j< static_cast<unsigned>(vGa.mNbAngles); ++j)
			{// The 1 in the mass equation is the proton/hydrogen mass
				
				double val = 1 - pow(1 / ( 1 + vSpeciesMass) , 2) * 
					pow( cos(vGa.mDThetaR(i,j)) +
							sqrt( pow(vSpeciesMass, 2) // vSpeciesMass/MassProton
								- pow(sin(vGa.mDThetaR(i,j)) , 2)	     
							    )
							,2);
				//if(val  > 1)
					tmp(i,j) = val;
			}
		}
		for(size_t i=0; i< mGrideV.size(); ++i)
		{
			mElossElas[i] = tmp * mGrideV[i];
		}
	}else
	{

		for(size_t i=0; i< mGrideV.size(); ++i)
		{
			mElossElas[i] = tmp;
		}
	}
	mbIsElossElasComputed = true;

}
