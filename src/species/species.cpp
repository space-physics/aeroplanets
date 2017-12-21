 /** 
 * \file species.cpp
 * \brief Implements the specie class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: species.cpp 1922 2014-02-26 22:38:52Z gronoff $
 *
 */


#include "species.hpp"
using namespace std;
typedef std::deque< ublas::vector<double>* > DequeUblas;
typedef std::deque< std::deque< ublas::vector<double> > > DequeDeque;



Specie::Specie(std::string vName,
		ublas::vector<double> *pPhgrideV,
		ublas::vector<double> *pElgrideV,
		ublas::vector<double> *pElEngddeV,
		ublas::vector<double> *pPrgrideV,
		std::string vParamfilename,
		bool bIsMonteCarlo)
{
	mName=trim(vName); // Clear the beginning and ending spaces, just to verify...
//	mpParams=new XmlParameters(vParamfilename);
	mpParams=XmlParameters::AttachFileParameter(vParamfilename);
	if(bIsMonteCarlo)
	{
		mpParams->SetMonteCarloActive();
	}
	mbPorterLoaded=false;
	mIsmpParamsHere=true;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;

	mProdPicPosition=0;
	mIsPhotoCrsCopied=false;
	mIsElecCrsCopied=false;
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	mpProtonEgrideV=pPrgrideV;
	mpPhotonEgrideV=pPhgrideV;
	mpElectronEgrideV=pElgrideV;
	mpElectronEngddeV=pElEngddeV;

	mProtonRedistributionFunction.resize(3);
	mProtonRedistributionFunction.clear();
	mProtonSigma.resize(3);
	mProtonSigma.clear();
	LoadSpecie(bIsMonteCarlo);
}

Specie::Specie(std::string vName,XmlParameters* spParams)
{
	mName=trim(vName); // Clear the beginning and ending spaces, just to verify...
//	cout<<"We load the specie "<<mName<<endl;
	mpParams=spParams;
	mbPorterLoaded=false;
	mIsmpParamsHere=false;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;

	mProtonRedistributionFunction.resize(3);
	mProtonRedistributionFunction.clear();
	mProtonSigma.resize(3);
	mProtonSigma.clear();
	
	mProdPicPosition=0;
	mIsPhotoCrsCopied=false;
	mIsElecCrsCopied=false;
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	LoadResuSpecie();
}

Specie::Specie(std::string vName,std::string vParamfilename)
{
	mName=trim(vName); // Clear the beginning and ending spaces, just to verify...
//	cout<<"We load the specie "<<mName<<endl;
	mpParams=XmlParameters::AttachFileParameter(vParamfilename);
	mbPorterLoaded=false;
	mIsmpParamsHere=false;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;

	mProtonRedistributionFunction.resize(3);
	mProtonRedistributionFunction.clear();
	mProtonSigma.resize(3);
	mProtonSigma.clear();

	mProdPicPosition=0;
	mIsPhotoCrsCopied=false;
	mIsElecCrsCopied=false;
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	LoadResuSpecie();
}



Specie::Specie(std::string vName,
		ublas::vector<double> *pPhgrideV,
		ublas::vector<double> *pElgrideV,
		ublas::vector<double> *pElEngddeV,
		ublas::vector<double> *pPrgrideV,
		XmlParameters* spParams,
		bool bIsMonteCarlo)
{

	mName=trim(vName); // Clear the beginning and ending spaces, just to verify...
	mpParams=spParams;
	
	mProtonRedistributionFunction.resize(3);
	mProtonRedistributionFunction.clear();
	mProtonSigma.resize(3);
	mProtonSigma.clear();
	
	if(bIsMonteCarlo)
	{
		mpParams->SetMonteCarloActive();
	}
	mProdPicPosition=0;
	mIsmpParamsHere=false;
	mbPorterLoaded=false;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;
	mIsPhotoCrsCopied=false;
	mIsElecCrsCopied=false;
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	mpProtonEgrideV=pPrgrideV;
	mpPhotonEgrideV=pPhgrideV;
	mpElectronEngddeV=pElEngddeV;
	mpElectronEgrideV=pElgrideV;
	LoadSpecie(bIsMonteCarlo);
}





Specie::~Specie()
{
	if(mIsPhotoCrsLoaded)
	{
		delete mpPhotoCrs;
	}

	if(mIsElecCrsLoaded)
	{
		delete mpElecCrs;
	}
	if(mIsProtCrsLoaded)
	{
		delete mpProtCrs;
	}
	if(mIsHCrsLoaded)
	{
		delete mpHCrs;
	}
	if(mbPorterLoaded)
	{
		//Log::SetPriority(Log::DEBUGG,"~Specie");
		Log::mD<<"We delete the Porter arrays"<<endl;
		delete mpEPortereV;
		delete mpGPorter;
		delete mpBPorter;
		delete mpAPorter;
		Log::mD<<"We finished deleting the Porter arrays"<<endl;
	}
	if(mIsmpParamsHere)
	{// When attached, this value is true!
		//delete mpParams;
		mpParams->DetachFileParameter();
	}


}



void Specie::PrintDensities(std::string vFilename)
{
	Log::mL<<vFilename<<endl;
}

void Specie::PrintCrossSections(std::string vFileprefix)
{ 
	string extname=vFileprefix+mName+"_";
	// Photon -> PrintCrs
	if(mIsPhotoCrsLoaded)
	{
		mpPhotoCrs->PrintCrs(extname+"Photon.dat");

	}
	// Electron -> PrintCrs(filename,redisfilename)
	if(mIsElecCrsLoaded)
	{
		mpElecCrs->PrintCrs(extname+"Electron.dat",extname+"ElectronRedist.dat");
	}
	// proton -> printcrs (now)
	if(mIsProtCrsLoaded)
	{
		mpProtCrs->PrintCrs(extname+"Proton.dat");
	}
	if(mIsHCrsLoaded)
	{
		mpHCrs->PrintCrs(extname+"Hydrogen.dat");
	}
}


void Specie::PrintProductions(std::string vFilename,const ublas::vector<double>& vAltGridKm)
{
	//cout<<vFilename<<endl;

	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Species production "<<mName<<endl;

	of<<"# Productions in cm-3"<<endl;
	of<<"# Altitude in km"<<endl;

	for(unsigned i=0;i<mTotSpecies.size();++i)
	{
		of<<"# "<<mTotSpecies[i].StandardName()<<endl;
	}

	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		of<<vAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mTotSpeciesProductioncm_3s_1.size();++j)
		{
			//of<<mTotSpeciesProductioncm_3s_1(j,i)<<"\t";
			of<<mTotSpeciesProductioncm_3s_1[j][i]<<"\t";
		}
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
	// Names std::vector<SpecieID> TotSpecies
	// Productions vector<vector<double>> TotSpeciesProduction
	// Altitudes 
}


void Specie::MergingProd()
{
	for(unsigned i=0;i<mCreatedSpecies.size();++i)
	{
		for(unsigned sp=0;sp<mCreatedSpecies[i].size();++sp)
		{
			unsigned position=0;
			if(PosInVector(mTotSpecies,mCreatedSpecies[i][sp],position))
			{
				if(mTotSpeciesProductioncm_3s_1[position].size()==0)
				{
					//mTotSpeciesProductioncm_3s_1[position]= MultDeQue( mSpeciesProductioncm_3s_1[i], mMultiplicativeFactor[i][sp] );
					mTotSpeciesProductioncm_3s_1[position]= mSpeciesProductioncm_3s_1[i]  *  mMultiplicativeFactor[i][sp] ;
				}else
				{
					mTotSpeciesProductioncm_3s_1[position]=( ( mSpeciesProductioncm_3s_1[i] * mMultiplicativeFactor[i][sp] ) + mTotSpeciesProductioncm_3s_1[position]);

				}
			}else
			{
				mTotSpecies.push_back( mCreatedSpecies[i][sp]);
				mTotSpeciesProductioncm_3s_1.push_back((mSpeciesProductioncm_3s_1[i] * mMultiplicativeFactor[i][sp] ));

			}
		}
	}
	// we compute the total electron
	mElecProductioncm_3s_1=mPhotoElecProductioncm_3s_1;
	mElecProductioncm_3s_1=AddOrFullVec(mElecProductioncm_3s_1,mElecElecProductioncm_3s_1);
	mElecProductioncm_3s_1=AddOrFullVec(mElecProductioncm_3s_1,mProtElecProductioncm_3s_1);
	mElecProductioncm_3s_1=AddOrFullVec(mElecProductioncm_3s_1,mHElecProductioncm_3s_1);
	mElecProductioncm_3s_1=AddOrFullVec(mElecProductioncm_3s_1,mCosmoElecProductioncm_3s_1);


}



void Specie::LoadResuSpecie()
{
	// Ok, nothing here now, but could be necessary in the future!

	/* Not really logical. Be careful, there is a bug here
	if(!spParams->exists("/"+name))
	{
		cout<<"Impossible to find your species"<<endl;
		cout<<"Species : "<<name<<endl;
		exit(1);
	}
	
	if(!spParams->exists("/"+name+"/mass"))
	{
		cout<<"Impossible to find your species mass"<<endl;
		cout<<"Species : "<<name<<endl;
		exit(1);
	}

	spParams->get_value("/"+name+"/mass",mass);
	*/
// To simplify the work, there is no predefined result states in the result species
// It could be added later, but complexifies a lot the merging functions
// as some vectors can be full or not.
/*	
#ifdef DEBUG
	if(!spParams->exists("/"+name+"/states"))
	{
		cout<<"Impossible to find your species states"<<endl;
		cout<<"Species : "<<name<<endl;
		exit(1);
	}
#endif
	states=spParams->get_param_array("/"+name+"/states");
*/

}

void Specie::LoadSpecie(bool bIsMonteCarlo)
{
	string name=StrReplace(mName,"+","_PLUS");
	if(!mpParams->Exists("/aero_species/"+name))
	{
		Log::mE<<"Impossible to find your species"<<endl;
		Log::mE<<"Species : "<<name<<endl;
		Error err("LoadSpecies","Debug: search species","Impossible to find your species");
		throw err;
	}
	
	if(!mpParams->Exists("/aero_species/"+name+"/mass"))
	{
		Log::mE<<"Impossible to find your species mass"<<endl;
		Log::mE<<"Species : "<<name<<endl;
		Error err("LoadSpecies","Debug: search species mass","Impossible to find your species mass");
		throw err;
	}

	mpParams->GetValue("/aero_species/"+name+"/mass",mMass);
	
#ifdef DEBUG
	if(!mpParams->Exists("/aero_species/"+name+"/states"))
	{
		Log::mE<<"Impossible to find your species states"<<endl;
		Log::mE<<"Species : "<<name<<endl;
		Error err("LoadSpecie","Debug: search specie state","Impossible to find your specie states");
		throw err;
	
	}
#endif
	mStates=mpParams->GetParamDeQue("/aero_species/"+name+"/states");


	if(!mpParams->Exists("/aero_species/"+name+"/ioni_photon"))
	{
		Log::mI<<"Impossible to find your species photoionization cross section"<<endl;
		Log::mI<<"Species : "<<name<<endl;
	//	exit(1);
	}else
	{
#ifdef DEBUG
		if((*mpPhotonEgrideV).size()==0)
		{
			Log::mW<<"Warning photon energy grid not defined"<<endl;
		}
#endif
		mIsPhotoCrsLoaded=true;
		mpPhotoCrs=new CrossSection(*mpPhotonEgrideV);
		string nfilename=mpParams->GetSubFileName("/aero_species/"+name+"/ioni_photon");
		mpPhotoCrs->LoadCrs(nfilename,name,bIsMonteCarlo);
		
	}


	if(!mpParams->Exists("/aero_species/"+name+"/ioni_electron"))
	{
		Log::mI<<"Impossible to find your species electron impact ionization cross section"<<endl;
		Log::mI<<"Species : "<<name<<endl;
	//	exit(1);
	}else
	{
#ifdef DEBUG
		if((*mpElectronEgrideV).size()==0)
		{
			Log::mW<<"Warning electron energy grid not defined"<<endl;
		}
#endif

		mIsElecCrsLoaded=true;
		mpElecCrs=new ElecCrossSection(*mpElectronEgrideV,*mpElectronEngddeV);
		string nfilename=mpParams->GetSubFileName("/aero_species/"+name+"/ioni_electron");

		bool is_opal = false;
		double opal_ebareV=0;
		double reesalpha = 1/31.5;
		double reesbeta = 339.;
		double reesgamma = 1/2.49;

		if(mpParams->Exists("/aero_species/"+name+"/Opal"))
		{
			mpParams->GetValue("/aero_species/"+name+"/Opal",opal_ebareV);
			assert(opal_ebareV!=0);
			is_opal=true;
		}else{
			mpParams->GetValueOrDefault("/aero_species/"+name+"/ReesAlpha",reesalpha);
			mpParams->GetValueOrDefault("/aero_species/"+name+"/ReesBeta",reesbeta);
			mpParams->GetValueOrDefault("/aero_species/"+name+"/ReesGamma",reesgamma);
		}



		mpElecCrs->LoadCrs(nfilename,name,bIsMonteCarlo,true,reesalpha,reesbeta,reesgamma,
				is_opal,opal_ebareV);
		
	}

	if(!mpParams->Exists("/aero_species/"+name+"/ioni_proton"))
	{
		Log::mI<<"Impossible to find your species proton impact cross section"<<endl;
		Log::mI<<"Species : "<<name<<endl;
	//	exit(1);
	}else
	{
#ifdef DEBUG
		if((*mpProtonEgrideV).size()==0)
		{
			Log::mI<<"Warning proton energy grid not defined"<<endl;
		}
#endif

		mIsProtCrsLoaded=true;
		mpProtCrs=new ProtCrossSection(*mpProtonEgrideV, "Proton");
		string nfilename=mpParams->GetSubFileName("/aero_species/"+name+"/ioni_proton");
		mpProtCrs->LoadCrs(nfilename,name,bIsMonteCarlo);
		if(mpParams->Exists("/aero_species/"+name+"/proton_ioni_redistribution"))
		{
			mpParams->GetValue("/aero_species/"+name+"/proton_ioni_redistribution", mProtonRedistributionFunction[0]);
			if(1 == mProtonRedistributionFunction[0])
			{
				mpParams->ExistsOrDie("/aero_species/"+name+"/proton_ioni_redistribution_sigma","You need a sigma for the redistribution of the protons");
				mpParams->GetValue("/aero_species/"+name+"/proton_ioni_redistribution_sigma", mProtonSigma[0]);
			}
		}
		if(mpParams->Exists("/aero_species/"+name+"/proton_elas_redistribution"))
		{
			mpParams->GetValue("/aero_species/"+name+"/proton_elas_redistribution", mProtonRedistributionFunction[1]);
			if(1 == mProtonRedistributionFunction[1])
			{
				mpParams->ExistsOrDie("/aero_species/"+name+"/proton_elas_redistribution_sigma","You need a sigma for the redistribution of the protons");
				mpParams->GetValue("/aero_species/"+name+"/proton_elas_redistribution_sigma", mProtonSigma[1]);
			}
		}
		if(mpParams->Exists("/aero_species/"+name+"/proton_exch_redistribution"))
		{
			mpParams->GetValue("/aero_species/"+name+"/proton_exch_redistribution", mProtonRedistributionFunction[2]);
			if(1 == mProtonRedistributionFunction[2])
			{
				mpParams->ExistsOrDie("/aero_species/"+name+"/proton_exch_redistribution_sigma","You need a sigma for the redistribution of the protons");
				mpParams->GetValue("/aero_species/"+name+"/proton_exch_redistribution_sigma", mProtonSigma[2]);
			}
		}
	}

	if(!mpParams->Exists("/aero_species/"+name+"/ioni_hydrogen"))
	{
		Log::mI<<"Impossible to find your species hydrogen impact cross section"<<endl;
		Log::mI<<"Species : "<<name<<endl;
	//	exit(1);
	}else
	{
#ifdef DEBUG
		if((*mpProtonEgrideV).size()==0)
		{
			Log::mI<<"Warning proton/hydrogen energy grid not defined"<<endl;
		}
#endif

		mIsHCrsLoaded=true;
		mpHCrs=new ProtCrossSection(*mpProtonEgrideV, "Hydrogen");
		string nfilename=mpParams->GetSubFileName("/aero_species/"+name+"/ioni_hydrogen");
		mpHCrs->LoadCrs(nfilename,name,bIsMonteCarlo);
		
	}
}


void Specie::InitPhotoIonization(unsigned nb_alt)
{
	mProcessNames=mpPhotoCrs->mProcessNames;
	mCreatedSpecies=mpPhotoCrs->mCreatedSpecies;
	ublas::vector<double> tmp(nb_alt);
	tmp.clear();
	mSpeciesProductioncm_3s_1.resize(mCreatedSpecies.size(),tmp);
	mSpeciesProductionProbabilitys_1.resize(mCreatedSpecies.size(), 0);
	mPhotoElecProductionProbabilitys_1 = 0;
	mPhotoDissociationProbabilitys_1 = 0;
	mMultiplicativeFactor=mpPhotoCrs->mCreatedSpeciesMultiplicator;
	mPhotoElecProductioncm_3s_1.resize(nb_alt);
	mPhotoElecProductioncm_3s_1.clear();
	InitDummyStates(nb_alt);

}

void Specie::InitElectronImpact(unsigned nb_alt)
{
	mProcessNames=mpElecCrs->mProcessNames;
	mCreatedSpecies=mpElecCrs->mCreatedSpecies;
	ublas::vector<double> tmp(nb_alt);
	tmp.clear();
	mSpeciesProductioncm_3s_1.resize(mCreatedSpecies.size(),tmp);
	mMultiplicativeFactor=mpElecCrs->mCreatedSpeciesMultiplicator;
	mElecElecProductioncm_3s_1.resize(nb_alt);
	mElecElecProductioncm_3s_1.clear();
	InitDummyStates(nb_alt);
}

void Specie::InitProtonImpact(unsigned nb_alt, boost::shared_ptr<MathFunction::GaussianAngle> vGa, bool vIsRedis, double vBeamSpreading)
{// We add the proton and hydrogen to that initialization, because they are computed together.
	assert(mIsProtCrsLoaded && mIsHCrsLoaded);
	
	if(not mIsProtCrsLoaded && mIsHCrsLoaded)
	{
		Error err("Proton impact ionization", "Proton::init","There was an error in the cross sections for the proton/hydrogen ionization");
		throw err;
	}
	mProcessNames=mpProtCrs->mProcessNames;
	ConcatenateDeQue(mProcessNames,mpHCrs->mProcessNames);
	mCreatedSpecies=mpProtCrs->mCreatedSpecies;
	ConcatenateDeQue(mCreatedSpecies, mpHCrs->mCreatedSpecies);
	ublas::vector<double> tmp(nb_alt);
	tmp.clear();
	mSpeciesProductioncm_3s_1.resize(mCreatedSpecies.size(),tmp);

	mMultiplicativeFactor=mpProtCrs->mCreatedSpeciesMultiplicator;
	ConcatenateDeQue(mMultiplicativeFactor, mpHCrs->mCreatedSpeciesMultiplicator);
	
	mProtElecProductioncm_3s_1.resize(nb_alt);
	mProtElecProductioncm_3s_1.clear();
	mHElecProductioncm_3s_1.resize(nb_alt);
	mHElecProductioncm_3s_1.clear();
	InitDummyStates(nb_alt);

	InitProtonPhase(vGa, vIsRedis, vBeamSpreading);
}

void Specie::InitProtonPhase(boost::shared_ptr<MathFunction::GaussianAngle> vGa, bool vIsRedis, double vBeamSpreading)
{
	int nba = vGa->mNbAngles;
	ublas::matrix<double> tmp(nba,nba);
	tmp.clear();
	mProtonPhase.resize(3);
	for(unsigned k = 0; k < 3; ++k)
	{
	//	if(mIsProtCrsLoaded)
	//	{
	//		Log::mL<<"Position "<<k<<" and "<< mProtonRedistributionFunction[k]<<endl;
	//	}else
		if(!mIsProtCrsLoaded)
		{
			Log::mL<<"Problem: you proton CRS is not loaded!!!"<<endl;
		}
		mProtonPhase[k] = tmp;
		if( vIsRedis && 1 == mProtonRedistributionFunction[k])
		{ // Maxwell function with sigma
//			Log::mL<<"Maxwell function, for position"<<k<<endl;
			if(mProtonSigma[k] > 1E-30)
			{
				for(int i = 0; i < nba ; ++ i)
				{
					double sum = 0;
					for(int j = 0; j < nba ; ++ j)
					{
						mProtonPhase[k](i,j) = 1. / sqrt(2*PI) / mProtonSigma[k] * exp(vGa->mDThetaR(i,j) * vGa->mDThetaR(i,j) / (2 * mProtonSigma[k] * mProtonSigma[k]));
						sum += mProtonPhase[k](i,j) * vGa->mWeight[j];
					}
					for(int j = 0; j < nba ; ++ j)
					{
						mProtonPhase[k](i,j) /= sum;
					}
				}
				continue;
			} // Else, we initialize as phaseforward
		}
		if( vIsRedis && 2 == mProtonRedistributionFunction[k])
		{ // Comb function (I learnt a word)
//			Log::mL<<"Comb function, for position"<<k<<endl;
			for(int i = 1; i < nba-1 ; ++ i)
			{
				int j = i - 1 ;
				mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 4);
				j = i + 1;
				mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 4);
				j = i;
				mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 2);
			}
			int i = 0;
			int j = i;
			mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 2);
			j = i + 1;
			mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 2);
			i = nba - 1;
			j = i - 1;
			mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 2);
			j = i;
			mProtonPhase[k](i,j) = 1 / (vGa->mWeight[j] * 2);
			continue;
		}
		if( vIsRedis && 3 == mProtonRedistributionFunction[k])
		{ // Screened Rutherford
//			double eps = 1E-3; // this is the beam spreading parameter
			double eps = vBeamSpreading;
//			Log::mL<<"Screened Rutherford, for position"<<k<<endl;
			for(int i = 0; i < nba ; ++ i)
			{
				double sum = 0.;
				for(int j = 0; j < nba ; ++ j)
				{
					double ruthcos = vGa->mXmu[i] * vGa->mXmu[j] + sqrt(1 - vGa->mXmu[j] * vGa->mXmu[j]) * sqrt(1 - vGa->mXmu[i] * vGa->mXmu[i]) * cos( vGa->mDThetaR(i, j))  ;
					double tmpr = 1 + 2 * eps - ruthcos;
					mProtonPhase[k](j,i) = 4. * eps * (1. + eps) / (tmpr * tmpr);
					sum += mProtonPhase[k](j,i) * vGa->mWeight[j];
				}
				for(int j = 0; j < nba ; ++ j)
				{
					mProtonPhase[k](j,i) /= sum;
				}
			}
			continue;
		}
		if(not vIsRedis)
		{
			Log::mL<<"The redistribution is not selected"<<endl;
		}
//		Log::mL<<"Standard, for position"<<k<<endl;
		for(int i = 0; i < nba ; ++ i)
		{ 
			if(vGa->mWeight[i] > 1E-30)
			{
				mProtonPhase[k](i,i) = 1. / vGa->mWeight[i];
			}
		}

	}
}

void Specie::InitCosmoImpact(unsigned vNbAlt)
{
	
	
	string name=StrReplace(mName,"+","_PLUS");
	mpParams->ExistsOrDie("/aero_species/"+name+"/cosmic/Specie","You should define the cosmic part of your species to work!");

	std::vector<TiXmlNode*> messp=mpParams->GetNodes("/aero_species/"+name+"/cosmic/Specie");
	mCreatedSpecies.resize(1);// We start with a non null vector
	mCreatedSpecies[0].resize(0);// But with one species.
	mCosmoNumberOfElectrons.resize(0);// We start at 0
	std::deque<double> fractions;
	std::vector<TiXmlNode*>::iterator jt;
	for(jt=messp.begin();jt!=messp.end();++jt)
	{
		string spname=mpParams->GetKey(*jt,"/","name");
		string spstate=mpParams->GetKey(*jt,"/","state");
		//			cout<<"spstate : "<<spstate<<endl;
		double fraction=0;
		mpParams->GetNKey(*jt,"/","fraction",fraction);
		fractions.push_back(fraction);
		double elnumber=0;
		mpParams->GetNKey(*jt,"/","electrons",elnumber);
		mCosmoNumberOfElectrons.push_back(elnumber*fraction);
		
		SpecieId tmpsp(spname,spstate);
		mCreatedSpecies[0].push_back(tmpsp);
	}
	mProcessNames.push_back(mName+" + Cosmic Ray -> Ions ");
	ublas::vector<double> tmp(vNbAlt);
	tmp.clear();
	mSpeciesProductioncm_3s_1.resize(mCreatedSpecies.size(),tmp);
	mMultiplicativeFactor.resize(1);// 1 process!
	mMultiplicativeFactor[0].resize(mCreatedSpecies[0].size());
	for(unsigned i=0;i<mMultiplicativeFactor[0].size();++i)
	{
		mMultiplicativeFactor.at(0)[i]=fractions.at(i);/////;
	}
	mSpeciesProductioncm_3s_1.resize(1);// 1 process!
	mSpeciesProductioncm_3s_1.at(0)=tmp;
	mSpeciesProductioncm_3s_1.at(0).clear();
	mCosmoElecProductioncm_3s_1.resize(vNbAlt);
	mCosmoElecProductioncm_3s_1.clear();
	InitDummyStates(vNbAlt);
}



void Specie::InitDummyStates(unsigned nb_alt)
{
	unsigned nbstates=mStates.size();
	mStateDensitycm_3.resize(nbstates);
	mStateProductioncm_3s_1.resize(nbstates);
	mStateLosscm_3s_1.resize(nbstates);
	mStateEscapecm_2s_1.resize(nbstates);
	ublas::zero_vector<double> zero(nb_alt);
	for(unsigned i=0;i<nbstates;++i)
	{
		mStateDensitycm_3.at(i)=zero;
		mStateDensitycm_3.at(i).clear();
		mStateProductioncm_3s_1.at(i)=zero;
		mStateProductioncm_3s_1.at(i).clear();
		mStateLosscm_3s_1.at(i)=zero;
		mStateLosscm_3s_1.at(i).clear();
		mStateEscapecm_2s_1.at(i)=0;
	}
}


/********************************************************/
/*                                                      */
/* void Specie::print_productions(std::string filename) */
/* {                                                    */
/*         cout<<"Hello"<<endl;                         */
/*                                                      */
/*                                                      */
/* }                                                    */
/********************************************************/


namespace SpecieUtils
{

	void SpeciesToResu(std::deque<Specie*>& rInitspec,std::deque<Specie*>& rResuspec)
	{
		//	cout<<"Salut les gars, c'est partit pour le travail"<<endl;
		//	cout<<"Taille de initspec "<<rInitspec.size()<<endl;
		deque<Specie*>::iterator it;
		for(it=rInitspec.begin();it!=rInitspec.end();++it)
		{
			//		cout<<"Merging specie "<<(*it)->mName<<endl;
			(*it)->FinishProduction();


			for(unsigned i=0;i<(*it)->mTotSpecies.size();++i)
			{
				//			cout<<"position "<<i<<" dans mise en place des differentes choses"<<endl;
				int pos=PosOfSp((*it)->mTotSpecies.at(i).mName,rResuspec);
				if(pos==-1)
				{// Create a new species and push it into resuspec
					Specie*tmp=new Specie((*it)->mTotSpecies.at(i).mName,(*it)->mpParams);
					tmp->mTotProductioncm_3s_1.clear();
					tmp->PutStateProduction((*it)->mTotSpecies.at(i).mState,(*it)->mTotSpeciesProductioncm_3s_1.at(i));
					rResuspec.push_back(tmp);

				}else
				{
					rResuspec[pos]->PutStateProduction((*it)->mTotSpecies.at(i).mState,(*it)->mTotSpeciesProductioncm_3s_1.at(i));
				}
			}
			//		//cout<<"Fin merging specie"<<endl;
		}
		//	cout<<"Fin tous merging, statestototprod"<<endl;
		for(it=rResuspec.begin();it!=rResuspec.end();++it)
		{
			(*it)->ComputePicProd();
		}
		//	cout<<rResuspec.size()<<endl;


	}

	int PosOfSp(std::string vName,const std::deque<Specie*>& vSp)
	{
		vName=trim(vName);
		for(unsigned i=0;i<vSp.size();++i)
		{
	//		cout<<i<<endl;
	//		cout<<vSp[i]->mName<<endl;
			if(vSp[i]->mName==vName)
			{
				return static_cast<int>(i);
			}
		}
		return -1;
	}

	bool GetSpecieDensity(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double> & rDensity)
	{
		int specieposition=SpecieUtils::PosOfSp(vName,vSp);
		if(specieposition<0)
		{
			return false;
		}
		if(""==trim(vState))
		{
			rDensity=vSp[specieposition]->mTotDensitycm_3;
			return true;
		}

		return vSp[specieposition]->GetStateDensity(vState,rDensity);
	}
	bool PutSpecieDensity(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double> vDensity)
	{
		int specieposition=SpecieUtils::PosOfSp(vName,vSp);
		if(specieposition<0)
		{
			return false;
		}
		vSp[specieposition]->PutStateDensity(vState,vDensity);
		return  true;
	}


	bool GetSpecieProduction(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double> & rDensity)
	{
		int specieposition=SpecieUtils::PosOfSp(vName,vSp);
		if(specieposition<0)
		{
			return false;
		}
		return vSp[specieposition]->GetStateProduction(vState,rDensity);
		return true;
	}
};

/*
Specie Specie::operator=(const Specie&s1)
{

	Log::SetPriority(Log::DEBUGG);
	Log::mL<<"copy constructor"<<endl;
	Specie resu(s1.mName,s1.mpParams);
	// The new species does not own the different vectors
	// so, it does not delete them when it is closed


*/
Specie::Specie(const Specie&s1)
{
	//	Log::SetPriority(Log::DEBUGG);
	//	Log::mL<<"copy constructor for "<<s1.mName<<endl;
	mName=s1.mName;
	mpParams=s1.mpParams;
	mpParams=s1.mpParams;
	mbPorterLoaded=false;
	mIsmpParamsHere=false;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;

	mIsPhotoCrsCopied=false;
	mIsElecCrsCopied=false;
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	if(s1.CheckPhotCrs())
	{
		mIsPhotoCrsCopied=true;
	}

	if(s1.CheckElecCrs())
	{
		mIsElecCrsCopied=true;
	}
	if(s1.CheckProtCrs())
	{
		mIsProtCrsCopied=true;
		mIsHCrsCopied=true;
	}
	LoadResuSpecie();
	mCosmoNumberOfElectrons=s1.mCosmoNumberOfElectrons;

	mProtonSigma = s1.mProtonSigma;
	mProtonRedistributionFunction = s1.mProtonRedistributionFunction;

	mpPhotonEgrideV=s1.mpPhotonEgrideV;
	mpElectronEgrideV=s1.mpElectronEgrideV;
	mpProtonEgrideV=s1.mpProtonEgrideV;
	mMass=s1.mMass;
	mStates=s1.mStates;
	mTotDensitycm_3=s1.mTotDensitycm_3;
	mColDenscm_2=s1.mColDenscm_2;
	mScaleHcm=s1.mScaleHcm;
	mTotProductioncm_3s_1=s1.mTotProductioncm_3s_1;
	mTotLosscm_3s_1=s1.mTotLosscm_3s_1;
	mTotEscapecm_2s_1=s1.mTotEscapecm_2s_1;
	mStateDensitycm_3=s1.mStateDensitycm_3;
	mStateProductioncm_3s_1=s1.mStateProductioncm_3s_1;
	mStateLosscm_3s_1=s1.mStateLosscm_3s_1;
	mStateEscapecm_2s_1=s1.mStateEscapecm_2s_1;

	mpPhotoCrs=s1.mpPhotoCrs;
	mpElecCrs=s1.mpElecCrs;
	mpProtCrs=s1.mpProtCrs;
	mpHCrs=s1.mpHCrs;
	mPhotoElecProductioncm_3s_1=s1.mPhotoElecProductioncm_3s_1;
	mElecElecProductioncm_3s_1=s1.mElecElecProductioncm_3s_1;
	mProtElecProductioncm_3s_1=s1.mProtElecProductioncm_3s_1;
	mHElecProductioncm_3s_1=s1.mHElecProductioncm_3s_1;
	mCosmoElecProductioncm_3s_1=s1.mCosmoElecProductioncm_3s_1;
	mElecProductioncm_3s_1=s1.mElecProductioncm_3s_1;


	mProcessNames=s1.mProcessNames;
	mCreatedSpecies=s1.mCreatedSpecies;
	mMultiplicativeFactor=s1.mMultiplicativeFactor;
	mSpeciesProductioncm_3s_1=s1.mSpeciesProductioncm_3s_1;

	mTotSpecies=s1.mTotSpecies;
	mTotSpeciesProductioncm_3s_1=s1.mTotSpeciesProductioncm_3s_1;
	mProdPicPosition=s1.mProdPicPosition;
}

Specie::Specie(std::deque< Specie* > vSpecies,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath)
{
	if(vSpecies.size()==0)
	{
		Error err("Copy of species for a bent atmosphere","No species in the list","There is no species in the input list... You must compute an atmosphere before doing any work");
		throw err;
	}
	size_t nbsp=vSpecies.size();

	Specie* s1=vSpecies[0];
	mName=s1->mName;
	mpParams=s1->mpParams;
	mbPorterLoaded=false;
	mIsmpParamsHere=false;
	mIsPhotoCrsLoaded=false;
	mIsElecCrsLoaded=false;
	mIsProtCrsLoaded=false;
	mIsHCrsLoaded=false;
	LoadResuSpecie();
	mpPhotonEgrideV=s1->mpPhotonEgrideV;
	mpElectronEgrideV=s1->mpElectronEgrideV;
	mpElectronEngddeV=s1->mpElectronEngddeV;
	mpProtonEgrideV=s1->mpProtonEgrideV;
	mpEPortereV=s1->mpEPortereV;
	mpGPorter=s1->mpGPorter;
	mpBPorter=s1->mpBPorter;
	mpAPorter=s1->mpAPorter;
	mCosmoNumberOfElectrons=s1->mCosmoNumberOfElectrons;
	mMass=s1->mMass;
	mStates=s1->mStates;

	mIsPhotoCrsCopied=false;
	if(s1->CheckPhotCrs())
	{
		mIsPhotoCrsCopied=true;
	}

	mIsElecCrsCopied=false;
	if(s1->CheckElecCrs())
	{
		mIsElecCrsCopied=true;
	}
	mIsProtCrsCopied=false;
	mIsHCrsCopied=false;
	if(s1->CheckProtCrs())
	{
		mIsProtCrsCopied=true;
		mIsHCrsCopied=true;

	}
	mProtonSigma = s1->mProtonSigma;
	mProtonRedistributionFunction = s1->mProtonRedistributionFunction;
	mpPhotoCrs=s1->mpPhotoCrs;
	mpElecCrs=s1->mpElecCrs;
	mpProtCrs=s1->mpProtCrs;
	mpHCrs=s1->mpHCrs;
	mProcessNames=s1->mProcessNames;
	mCreatedSpecies=s1->mCreatedSpecies;
	mMultiplicativeFactor=s1->mMultiplicativeFactor;
	mTotSpecies=s1->mTotSpecies;

	DequeUblas tmptotdensity;
	Log::mD<<"Bend density"<<endl;
	if(vSpecies[0]->mTotDensitycm_3.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmptotdensity.push_back(&(vSpecies[i]->mTotDensitycm_3));
		}
		mTotDensitycm_3=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmptotdensity,vPath.mAltitudeKm,vPath.mSZADegree);
		MinValue(mTotDensitycm_3,1E-42);
	}
	DequeUblas tmpscaleh;

	Log::mD<<"Bend scale H"<<endl;
	if(vSpecies[0]->mScaleHcm.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpscaleh.push_back(&(vSpecies[i]->mScaleHcm));
		}
		mScaleHcm=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmpscaleh,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	Log::mD<<"Col dens"<<endl;
	if(mScaleHcm.size()>0)
	{
		assert(mScaleHcm.size()==mTotDensitycm_3.size());
		mColDenscm_2.resize(mTotDensitycm_3.size());
		mColDenscm_2[0]=mTotDensitycm_3[0]*mScaleHcm[0];
		double coldensinit=mColDenscm_2[0];

		for(unsigned i=1;i<mTotDensitycm_3.size();++i)
		{// We make an integration along the path, so mLengthKm is integrated
			mColDenscm_2[i]=(MathFunction::TrapzInt(vPath.mLengthKm,mTotDensitycm_3,i)*1E5+coldensinit);
		}
	}

	DequeUblas tmptp;

	Log::mD<<"Tot production"<<endl;
	if(vSpecies[0]->mTotProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mTotProductioncm_3s_1),1E-42);
			tmptp.push_back(&(vSpecies[i]->mTotProductioncm_3s_1));
		}
		mTotProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmptp,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	DequeUblas tmplp;

	Log::mD<<"Tot loss"<<endl;
	if(vSpecies[0]->mTotLosscm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{

			MinValue((vSpecies[i]->mTotLosscm_3s_1),1E-42);
			tmplp.push_back(&(vSpecies[i]->mTotLosscm_3s_1));
		}
		mTotLosscm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmplp,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	DequeDeque tmpstdens;

	Log::mD<<"Bend state density"<<endl;
	if(vSpecies[0]->mStateDensitycm_3.size()>0)
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpstdens.push_back((vSpecies[i]->mStateDensitycm_3));
		}
		mStateDensitycm_3=MathFunction::SetIntLogPath(vAltitudesKm,vSZADegree,tmpstdens,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	DequeDeque tmpstprod;

	Log::mD<<"Bend state production "<<endl;
	if(vSpecies[0]->mStateProductioncm_3s_1.size()>0)
	{
		size_t tailleo=(vSpecies[0]->mStateProductioncm_3s_1).size();
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpstprod.push_back((vSpecies[i]->mStateProductioncm_3s_1));
			size_t tailleo2=(vSpecies[i]->mStateProductioncm_3s_1).size();
			assert(tailleo==tailleo2);
		}
		try
		{
			mStateProductioncm_3s_1=MathFunction::SetIntLogPath(vAltitudesKm,vSZADegree,tmpstprod,vPath.mAltitudeKm,vPath.mSZADegree);
		}
		catch(MathError& err)
		{
			Log::mE<<"Math error"<<endl;
			err.Affiche();
		}
		catch(Error& err)
		{
			Log::mE<<"Standard Error"<<endl;
			err.Affiche();
		}
		catch(...)
		{
			Log::mE<<"Unknown error"<<endl;
		}
	}

	DequeDeque tmpstloss;
	Log::mD<<"Bend stateloss "<<endl;
	if(vSpecies[0]->mStateLosscm_3s_1.size()>0)
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpstloss.push_back((vSpecies[i]->mStateLosscm_3s_1));
		}
		mStateLosscm_3s_1=MathFunction::SetIntLogPath(vAltitudesKm,vSZADegree,tmpstloss,vPath.mAltitudeKm,vPath.mSZADegree);
	}



	DequeUblas tmppephot;
	Log::mD<<"Bend photelecproduction "<<endl;
	if(vSpecies[0]->mPhotoElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mPhotoElecProductioncm_3s_1),1E-42);
			tmppephot.push_back(&(vSpecies[i]->mPhotoElecProductioncm_3s_1));
		}
		mPhotoElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppephot,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	DequeUblas tmppee;
	Log::mD<<"Bend elec elec production "<<endl;
	if(vSpecies[0]->mElecElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mElecElecProductioncm_3s_1),1E-42);
			tmppee.push_back(&(vSpecies[i]->mElecElecProductioncm_3s_1));
		}
		mElecElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppee,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	DequeUblas tmppeprot;

	Log::mD<<"Bend protelec production "<<endl;
	if(vSpecies[0]->mProtElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mProtElecProductioncm_3s_1),1E-42);
			tmppeprot.push_back(&(vSpecies[i]->mProtElecProductioncm_3s_1));
		}
		mProtElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppeprot,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	DequeUblas tmppeh;

	Log::mD<<"Bend protelec production "<<endl;
	if(vSpecies[0]->mHElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mHElecProductioncm_3s_1),1E-42);
			tmppeh.push_back(&(vSpecies[i]->mHElecProductioncm_3s_1));
		}
		mHElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppeh,vPath.mAltitudeKm,vPath.mSZADegree);
	}



	DequeUblas tmppecosmo;

	Log::mD<<"Bend cosmo elec production "<<endl;
	if(vSpecies[0]->mCosmoElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mCosmoElecProductioncm_3s_1),1E-42);
			tmppecosmo.push_back(&(vSpecies[i]->mCosmoElecProductioncm_3s_1));
		}
		mCosmoElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppecosmo,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	DequeUblas tmppe;

	Log::mD<<"Bend elec production "<<endl;
	if(vSpecies[0]->mElecProductioncm_3s_1.size()==vAltitudesKm.size())
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			MinValue((vSpecies[i]->mElecProductioncm_3s_1),1E-42);
			tmppe.push_back(&(vSpecies[i]->mElecProductioncm_3s_1));
		}
		mElecProductioncm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmppe,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	Log::mD<<"bend escape"<<endl;
	// HERE WE WORK ON THE ESCAPE
	// we consider the escape at the upper level...
	ublas::vector<double> tmpescape(nbsp);
	ublas::vector<double> tmpsza(nbsp);
	for(unsigned i=0;i<nbsp;++i)
	{
		tmpescape(i)=vSpecies[i]->mTotEscapecm_2s_1;
		tmpsza(i)=vSZADegree[i];
	}
	ublas::vector<double> outsza(1);
	outsza[0]=vPath.mSZADegree[0];
	ublas::vector<double> resu=MathFunction::IntLin(tmpsza,tmpescape,outsza);
	mTotEscapecm_2s_1=resu[0];

	for(unsigned j=0;j<s1->mStateEscapecm_2s_1.size();++j)
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpescape(i)=vSpecies[i]->mStateEscapecm_2s_1[j];
		}
		ublas::vector<double> resu=MathFunction::IntLin(tmpsza,tmpescape,outsza);
		mStateEscapecm_2s_1.push_back(resu[0]);

	}

	Log::mD<<"End of bend escape"<<endl;

	// End of escape considerations
	DequeDeque tmpspprod;
	if(vSpecies[0]->mSpeciesProductioncm_3s_1.size()>0)
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmpspprod.push_back((vSpecies[i]->mSpeciesProductioncm_3s_1));
		}
		mSpeciesProductioncm_3s_1=MathFunction::SetIntLogPath(vAltitudesKm,vSZADegree,tmpspprod,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	DequeDeque tmptotspprod;
	if(vSpecies[0]->mTotSpeciesProductioncm_3s_1.size()>0)
	{
		for(unsigned i=0;i<nbsp;++i)
		{
			tmptotspprod.push_back((vSpecies[i]->mTotSpeciesProductioncm_3s_1));
		}
		mTotSpeciesProductioncm_3s_1=MathFunction::SetIntLogPath(vAltitudesKm,vSZADegree,tmptotspprod,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	Log::mD<<"Compute pic prod"<<endl;
	ComputePicProd();
	Log::mD<<"End of compute pic prod"<<endl;

}

Specie& Specie::operator+=(const Specie& s2)
{
	assert(mName==s2.mName);
	//assert(mMass==s2.mMass);
	//assert(mpPhotonEgrideV==s2.mpPhotonEgrideV);
	//assert(mpProtonEgrideV==s2.mpProtonEgrideV);
	// This is intented to add the result species: we only add production, losses, and escape for the 1st point
	mTotProductioncm_3s_1 = AddOrFullVec(mTotProductioncm_3s_1, s2.mTotProductioncm_3s_1);
	mTotLosscm_3s_1 = AddOrFullVec(mTotLosscm_3s_1, s2.mTotLosscm_3s_1);
	mTotEscapecm_2s_1 = mTotEscapecm_2s_1 + s2.mTotEscapecm_2s_1;

	if(mStates.size()!=mStateProductioncm_3s_1.size())
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Probleme sur les tailles au niveau de l'initialisation"<<endl;
		Log::mE<<"Mstates size : "<<mStates.size();

	}

	if(mCosmoNumberOfElectrons.size()==0)
	{
		mCosmoNumberOfElectrons=s2.mCosmoNumberOfElectrons;
	}


	// First difficult part : the states
	for(unsigned i=0;i<s2.mStates.size();++i)
	{
		unsigned pos=0;
		if(!PosInVector(mStates,s2.mStates[i],pos))
		{// We can add, because it is not the s1 which is read by the loop
			mStates.push_back(s2.mStates[i]);
			if(s2.mStateDensitycm_3.size()!=0)
			{
				mStateDensitycm_3.push_back(s2.mStateDensitycm_3[i]);
			}
			mStateProductioncm_3s_1.push_back(s2.mStateProductioncm_3s_1[i]);
			if(s2.mStateLosscm_3s_1.size()!=0)
			{
				mStateLosscm_3s_1.push_back(s2.mStateLosscm_3s_1[i]);
			}
			if(s2.mStateEscapecm_2s_1.size()!=0)
			{
				mStateEscapecm_2s_1.push_back(s2.mStateEscapecm_2s_1[i]);
			}
		} else
		{// we can merge density, production and loss
			if(mStateDensitycm_3.size()!=0)
			{
				mStateDensitycm_3[pos]= AddOrFullVec(mStateDensitycm_3[pos],s2.mStateDensitycm_3[i]);
			}
			mStateProductioncm_3s_1[pos]= AddOrFullVec(mStateProductioncm_3s_1[pos],s2.mStateProductioncm_3s_1[i]);
			if(mStateLosscm_3s_1.size()!=0)
			{
				mStateLosscm_3s_1[pos]= AddOrFullVec(mStateLosscm_3s_1[pos],s2.mStateLosscm_3s_1[i]);
			}
			if(mStateEscapecm_2s_1.size()!=0)
			{
				mStateEscapecm_2s_1[pos] += s2.mStateEscapecm_2s_1[i];
			}
		}

	}

	if(mStates.size()!=mStateProductioncm_3s_1.size())
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Probleme sur les tailles apres l'ajout"<<endl;
		Log::mE<<"Mstates size : "<<mStates.size();

	}
	assert(mStates.size()==mStateProductioncm_3s_1.size());

	mPhotoElecProductioncm_3s_1=AddOrFullVec(mPhotoElecProductioncm_3s_1,s2.mPhotoElecProductioncm_3s_1);
	mElecElecProductioncm_3s_1=AddOrFullVec(mElecElecProductioncm_3s_1,s2.mElecElecProductioncm_3s_1);
	mProtElecProductioncm_3s_1=AddOrFullVec(mProtElecProductioncm_3s_1,s2.mProtElecProductioncm_3s_1);
	mHElecProductioncm_3s_1=AddOrFullVec(mHElecProductioncm_3s_1,s2.mHElecProductioncm_3s_1);
	mCosmoElecProductioncm_3s_1=AddOrFullVec(mCosmoElecProductioncm_3s_1,s2.mCosmoElecProductioncm_3s_1);
	mElecProductioncm_3s_1=AddOrFullVec(mElecProductioncm_3s_1,s2.mElecProductioncm_3s_1);
	
	// Second difficult part : the created species
	for(unsigned i=0;i<s2.mProcessNames.size();++i)
	{
		unsigned pos=0;
		//PosInVector(s2.mProcessNames[i],mProcessNames,pos);
		if(PosInVector(mProcessNames,s2.mProcessNames[i],pos))
		{// we can merge density, production and loss
			// Theoretically this case should never append
			// In fact, if we add two species, that's because
			// the processes are different
			// so I raise a warning

			Log::mW<<"WARNING"<<endl;
			Log::mW<<"WARNING"<<endl;
			Log::mW<<"WARNING adding two species with same ionization process_names!!!!!!!!!!!"<<endl;
			Log::mW<<"WARNING if you do not understand this message, it must be considered as critical!"<<endl;
			Log::mW<<"WARNING"<<endl;
			Log::mW<<"WARNING"<<endl;
			mMultiplicativeFactor[pos]=AddDeQue(mMultiplicativeFactor[pos],s2.mMultiplicativeFactor[i]);
			mSpeciesProductioncm_3s_1[pos]=AddOrFullVec(mSpeciesProductioncm_3s_1[pos],s2.mSpeciesProductioncm_3s_1[i]);
		}else
		{// We can add, because it is not the s1 which is read by the loop
			mProcessNames.push_back(s2.mProcessNames[i]);
			mCreatedSpecies.push_back(s2.mCreatedSpecies[i]);
			mMultiplicativeFactor.push_back(s2.mMultiplicativeFactor[i]);
			mSpeciesProductioncm_3s_1.push_back(s2.mSpeciesProductioncm_3s_1[i]);
		}
	}


	for(unsigned i=0;i<s2.mTotSpecies.size();++i)
	{ // Nice SpecieId has a == parameter!
		unsigned pos=0;
		
		if(PosInVector(mTotSpecies,s2.mTotSpecies[i],pos))
		{
			mTotSpeciesProductioncm_3s_1[pos]=AddOrFullVec(mTotSpeciesProductioncm_3s_1[pos],s2.mTotSpeciesProductioncm_3s_1[i]);
		}else
		{
			mTotSpecies.push_back(s2.mTotSpecies[i]);
			mTotSpeciesProductioncm_3s_1.push_back(s2.mTotSpeciesProductioncm_3s_1[i]);
		}
	}
	ComputePicProd();
	return *this;
}

	/*
	Log::SetPriority(Log::DEBUGG);
	Log::mL<<"Addition of two species"<<endl;
	assert(s1.mName==s2.mName);
	assert(s1.mMass==s2.mMass);
	assert(s1.mpPhotonEgrideV==s2.mpPhotonEgrideV);
	assert(s1.mpProtonEgrideV==s2.mpProtonEgrideV);
	// This is intented to add the result species: we only add production, losses, and escape for the 1st point
	s1.mTotProductioncm_3s_1=AddVec(s1.mTotProductioncm_3s_1,s2.mTotProductioncm_3s_1);
	s1.mTotLosscm_3s_1=AddVec(s1.mTotLosscm_3s_1,s2.mTotLosscm_3s_1);
	s1.mTotEscapecm_2s_1+=s2.mTotEscapecm_2s_1;

	// First difficult part : the states
	for(unsigned i=0;i<s2.mStates.size();++i)
	{
		int pos=PosInVector(s2.mStates[i],s1.mStates);
		if(pos!=-1)
		{// we can merge density, production and loss
			s1.mStateDensitycm_3[pos]=AddVec(s1.mStateDensitycm_3[i],s2.mStateDensitycm_3[i]);
			s1.mStateProductioncm_3s_1[pos]=AddVec(s1.mStateProductioncm_3s_1[i],s2.mStateProductioncm_3s_1[i]);
			s1.mStateLosscm_3s_1[pos]=AddVec(s1.mStateLosscm_3s_1[i],s2.mStateLosscm_3s_1[i]);
			s1.mStateEscapecm_2s_1[pos]+=s2.mStateEscapecm_2s_1[i];
		}else
		{// We can add, because it is not the s1 which is read by the loop
			s1.mStates.push_back(s2.mStates[i]);
			s1.mStateDensitycm_3.push_back(s2.mStateDensitycm_3[i]);
			s1.mStateProductioncm_3s_1.push_back(s2.mStateProductioncm_3s_1[i]);
			s1.mStateLosscm_3s_1.push_back(s2.mStateLosscm_3s_1[i]);
			s1.mStateEscapecm_2s_1.push_back(s2.mStateEscapecm_2s_1[i]);
		}
	}

	s1.mPhotoElecProductioncm_3s_1=AddVec(s1.mPhotoElecProductioncm_3s_1,s2.mPhotoElecProductioncm_3s_1);
	s1.mElecElecProductioncm_3s_1=AddVec(s1.mElecElecProductioncm_3s_1,s2.mElecElecProductioncm_3s_1);
	s1.mProtElecProductioncm_3s_1=AddVec(s1.mProtElecProductioncm_3s_1,s2.mProtElecProductioncm_3s_1);
	s1.mElecProductioncm_3s_1=AddVec(s1.mElecProductioncm_3s_1,s2.mElecProductioncm_3s_1);
	
	
	// Second difficult part : the created species
	for(unsigned i=0;i<s2.mProcessNames.size();++i)
	{
		int pos=PosInVector(s2.mProcessNames[i],s1.mProcessNames);
		if(pos!=-1)
		{// we can merge density, production and loss
			// Theoretically this case should never append
			// In fact, if we add two species, that's because
			// the processes are different
			// so I raise a warning

			cout<<"WARNING"<<endl;
			cout<<"WARNING"<<endl;
			cout<<"WARNING adding two species with same ionization process_names!!!!!!!!!!!"<<endl;
			cout<<"WARNING if you do not understand this message, it must be considered as critical!"<<endl;
			cout<<"WARNING"<<endl;
			cout<<"WARNING"<<endl;
			s1.mMultiplicativeFactor[pos]=AddVec(s1.mMultiplicativeFactor[i],s2.mMultiplicativeFactor[i]);
			s1.mSpeciesProductioncm_3s_1[pos]=AddVec(s1.mSpeciesProductioncm_3s_1[i],s2.mSpeciesProductioncm_3s_1[i]);
		}else
		{// We can add, because it is not the s1 which is read by the loop
			s1.mProcessNames.push_back(s2.mProcessNames[i]);
			s1.mCreatedSpecies.push_back(s2.mCreatedSpecies[i]);
			s1.mMultiplicativeFactor.push_back(s2.mMultiplicativeFactor[i]);
			s1.mSpeciesProductioncm_3s_1.push_back(s2.mSpeciesProductioncm_3s_1[i]);
		}
	}


	for(unsigned i=0;i<s2.mTotSpecies.size();++i)
	{ // Nice SpecieId has a == parameter!
		int pos=PosInVector(s2.mTotSpecies[i],s1.mTotSpecies);
		if(pos!=-1)
		{
			s1.mTotSpeciesProductioncm_3s_1[pos]=AddVec(s1.mTotSpeciesProductioncm_3s_1[pos],s2.mTotSpeciesProductioncm_3s_1[i]);
		}else
		{
			s1.mTotSpecies.push_back(s2.mTotSpecies[i]);
			s1.mTotSpeciesProductioncm_3s_1.push_back(s2.mTotSpeciesProductioncm_3s_1[i]);
		}
	}
	return s1;


}
*/


void SpecieUtils::PrintProduction(std::deque<Specie*> vSpecies,ublas::vector<double> vAltGridKm,std::string vFilename,bool bPrintStateProduction,bool vIsLength)
{
	//cout<<"we print the production : "<<vFilename<<endl;

	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Species production"<<endl;

	of<<"# Productions in cm-3"<<endl;
	if(vIsLength)
	{
		of<<"# Length in km (bent case)"<<endl;
	}else
	{
		of<<"# Altitude in km"<<endl;
	}

	deque<Specie*>::iterator sp;


	// Firstly, we print the different names
	for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
	{
		of<<"# "<<(*sp)->mName<<" (peak "<<ntostr(vAltGridKm[(*sp)->mProdPicPosition])<<" Km)"<<endl; 
		if(bPrintStateProduction)
		{
			for(unsigned i=0;i<(*sp)->mStates.size();++i)
			{
				of<<"# "<<(*sp)->mName<<"("<<(*sp)->mStates[i]<<")"<<endl;
			}
		}
	}

	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		of<<vAltGridKm[i]<<"\t";
		for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
		{
			if( ((*sp)->mTotProductioncm_3s_1).size() == vAltGridKm.size())
			{
				of<<(*sp)->mTotProductioncm_3s_1[i]<<"\t";
			}else
			{
				of<<0.<<"\t";
			}
			//of<<(*sp)->mTotProductioncm_3s_1[i]<<"\t";
			if(bPrintStateProduction)
			{
				//assert((*sp)->mStates.size()==(*sp)->mStateProductioncm_3s_1.size());
				for(unsigned j=0;j<(*sp)->mStates.size();++j)
				{
					//	of<<"# "<<(*sp)->name<<"("<<(*sp)->states[i]<<")"<<endl;
					if( ((*sp)->mStateProductioncm_3s_1[j]).size() == vAltGridKm.size())
					{
						of<<(*sp)->mStateProductioncm_3s_1[j][i]<<"\t";
					}else
					{
						of<<0.<<"\t";
					}
				}
			}
		}	
		of<<endl;
	}
	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}

void SpecieUtils::PrintElectronProduction(std::deque<Specie*> vSpecies,ublas::vector<double> vAltGridKm,std::string vFilename,bool vIsLength)
{
	//cout<<"we print the production : "<<vFilename<<endl;

	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Electron production"<<endl;

	of<<"# Productions in cm-3"<<endl;
	if(vIsLength)
	{
		of<<"# Length in km (bent case)"<<endl;
	}else
	{
		of<<"# Altitude in km"<<endl;
	}

	deque<Specie*>::iterator sp;

	ublas::vector<double> tote,phote,prote,cosme,elece,hydrogene;


	for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
	{
		tote=AddOrFullVec((*sp)->mElecProductioncm_3s_1,tote);
		phote=AddOrFullVec((*sp)->mPhotoElecProductioncm_3s_1,phote);
		prote=AddOrFullVec((*sp)->mProtElecProductioncm_3s_1,prote);
		hydrogene=AddOrFullVec((*sp)->mHElecProductioncm_3s_1,hydrogene);
		cosme=AddOrFullVec((*sp)->mCosmoElecProductioncm_3s_1,cosme);
		elece=AddOrFullVec((*sp)->mElecElecProductioncm_3s_1,elece);
	}
	bool isphot=false;
	bool isprot=false;
	bool ishydrogene=false;
	bool iscosm=false;
	bool iselec=false;

	assert(tote.size()==vAltGridKm.size());
	of<<"# Tot e- (peak "<<ntostr( vAltGridKm[PosOfMax(tote)])<<")"<<endl;
	if(phote.size()>0)
	{
		assert(phote.size()==tote.size());
		of<<"# Phot e- (peak "<<ntostr( vAltGridKm[PosOfMax(phote)])<<")"<<endl;
		isphot=true;
	}
	if(prote.size()>0)
	{
		assert(prote.size()==tote.size());
		of<<"# Prot e- (peak "<<ntostr( vAltGridKm[PosOfMax(prote)])<<")"<<endl;
		isprot=true;
	}	
	if(hydrogene.size()>0)
	{
		assert(hydrogene.size()==tote.size());
		of<<"# H e- (peak "<<ntostr( vAltGridKm[PosOfMax(hydrogene)])<<")"<<endl;
		ishydrogene=true;
	}
	if(cosme.size()>0)
	{
		assert(cosme.size()==tote.size());
		of<<"# Cosmo e- (peak "<<ntostr( vAltGridKm[PosOfMax(cosme)])<<")"<<endl;
		iscosm=true;
	}
	if(elece.size()>0)
	{
		assert(elece.size()==tote.size());
		of<<"# Elec e- (peak "<<ntostr( vAltGridKm[PosOfMax(elece)])<<")"<<endl;
		iselec=true;
	}
/*
	for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
	{
		of<<"# "<<(*sp)->mName<<" (peak "<<ntostr(vAltGridKm[(*sp)->mProdPicPosition])<<" Km)"<<endl; 
		if(bPrintStateProduction)
		{
			for(unsigned i=0;i<(*sp)->mStates.size();++i)
			{
				of<<"# "<<(*sp)->mName<<"("<<(*sp)->mStates[i]<<")"<<endl;
			}
		}
	}
*/
	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		of<<vAltGridKm[i]<<"\t"<<tote[i]<<"\t";
		if(isphot)
			of<<phote[i]<<"\t";
		if(isprot)
			of<<prote[i]<<"\t";
		if(ishydrogene)
			of<<hydrogene[i]<<"\t";
		if(iscosm)
			of<<cosme[i]<<"\t";
		if(iselec)
			of<<elece[i]<<"\t";
		// Secondly, we print the different names
		/*		for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
				{
				of<<(*sp)->mTotProductioncm_3s_1[i]<<"\t";
				if(bPrintStateProduction)
				{
				assert((*sp)->mStates.size()==(*sp)->mStateProductioncm_3s_1.size());
				for(unsigned j=0;j<(*sp)->mStates.size();++j)
				{
		//	of<<"# "<<(*sp)->name<<"("<<(*sp)->states[i]<<")"<<endl;
		of<<(*sp)->mStateProductioncm_3s_1[j][i]<<"\t";
		}
		}
		}	
		*/
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}



std::deque<Specie*> SpecieUtils::MergeResult(std::deque<Specie*> s1,std::deque<Specie*> s2)
{// We want the possibility to use effectiveley the new species created
// Therefore, we copy all the species

	deque<Specie*> resu;
//	Log::mI<<"Start merge "<<endl;

	std::deque<Specie*>::iterator sp;
	for(sp=s1.begin();sp!=s1.end();++sp)
	{
		Specie* tmp=new Specie(*(*sp));
		resu.push_back(tmp);
	}
	for(sp=s2.begin();sp!=s2.end();++sp)
	{
		int pos=PosOfSp((*sp)->mName,resu);
//		Log::mI<<"Species "<<(*sp)->mName<<endl;
		if(pos==-1)
		{
//		Log::mI<<"Add "<<(*sp)->mName<<endl;
			Specie* tmp=new Specie(*(*sp));
		//	*tmp=*(*sp); // Normally, there is a copy here
			resu.push_back(tmp);
		}else
		{
//		Log::mI<<"Sum "<<(*sp)->mName<<endl;
			*(resu[pos])+=*(*sp); // We add!

		}
//		Log::mI<<"End "<<(*sp)->mName<<endl;

	}

	for(sp=resu.begin();sp!=resu.end();++sp)
	{
//		Log::mI<<"Compute pic "<<(*sp)->mName<<endl;
		(*sp)->ComputePicProd();
	}
//	Log::mI<<"End merge "<<endl;
	return resu;
}

std::deque<Specie*> SpecieUtils::BendAtmosphere(std::deque< std::deque<Specie*> > vAtmoSpecies,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath)
{
	std::deque<Specie*> outspecies;
	if(vAtmoSpecies.size()<1)
	{
		return outspecies;
	}
	size_t nb_d_sp=vAtmoSpecies[0].size();

	for(unsigned j=0;j<nb_d_sp;++j)
	{
		std::deque<Specie*> tmpspecies;
		for(unsigned i=0;i<vAtmoSpecies.size();++i)
		{
			// We copy only the pointer to the species.
			tmpspecies.push_back(vAtmoSpecies[i][j]);
		}
		assert(vSZADegree.size()==tmpspecies.size());

		Specie* mysp=new Specie(tmpspecies,vAltitudesKm,vSZADegree,vPath);
		outspecies.push_back(mysp);

	}
	return outspecies;
}


std::deque<Specie*>  SpecieUtils::CopyAtmo(std::deque<Specie*> vAtmo)
{
	std::deque<Specie*> resultat;
	for(unsigned i=0;i<vAtmo.size();++i)
	{
		Specie* tmp=new Specie(*(vAtmo[i]));
		resultat.push_back(tmp);
	}
	return resultat;
}

void SpecieUtils::SelectedPrintProduction(std::deque<SpecieId>& rStates, std::deque<Specie*> rSpecies, ublas::vector<double>  vAltKmGrid, std::string rFilename,bool vIsLength)
{
//	cout<<"we print the production : "<<rFilename<<endl;

	if(FileExists(rFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<rFilename<<endl;
	}
	ofstream of(rFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Species production "<<endl;

	of<<"# Productions in cm-3"<<endl;
	if(vIsLength)
	{
		of<<"# Length in km (bent case)"<<endl;
	}else
	{
		of<<"# Altitude in km"<<endl;
	}
	deque<Specie*>::iterator sp;


	deque< ublas::vector<double>* > the_productions;
	ublas::vector<double> zero_production(vAltKmGrid.size());
	zero_production.clear();

	// Firstly, we print the different names.
	// While doing that, we find the production vectors
	// and we store their pointeur in the_production
	// this vector will be used to print
		
	for(unsigned i=0;i<rStates.size();++i)
	{
		of<<"# "<<rStates[i].StandardName()<<"    ";//endl;;
	//	of<<"# "<<(*sp)->mName<<" (peak "<<ntostr(vAltGridKm[(*sp)->mProdPicPosition])<<" Km)"<<endl; 
		int pos=PosOfSp(rStates[i].mName,rSpecies);
		if(pos==-1)
		{
			the_productions.push_back(&zero_production);
		
		}else
		{
			if(rStates[i].mState=="")
			{
				the_productions.push_back(&(rSpecies[pos]->mTotProductioncm_3s_1));
				of<<" (peak "<<ntostr(vAltKmGrid[rSpecies[pos]->mProdPicPosition])<<" Km)"; 
			}else
			{
				unsigned sta=0;
				if(!PosInVector(rSpecies[pos]->mStates,rStates[i].mState,sta))
				{
					the_productions.push_back(&zero_production);
				}else
				{
					the_productions.push_back(&(rSpecies[pos]->mStateProductioncm_3s_1[sta]));
				}
			}
		}
		of<<endl;
	}
	for(unsigned i=0;i<vAltKmGrid.size();++i)
	{
		of<<vAltKmGrid[i]<<"\t";
		for(unsigned j=0;j<the_productions.size();++j)
		{
			of<<(*(the_productions[j]))[i]<<"\t";
		}
		of<<endl;


	/*	for(sp=species.begin();sp!=species.end();++sp)
		{
			of<<(*sp)->TotProduction[i]<<"\t";
			if(print_state_production)
			{
				for(unsigned j=0;j<(*sp)->states.size();++j)
				{
					//	of<<"# "<<(*sp)->name<<"("<<(*sp)->states[i]<<")"<<endl;
					of<<(*sp)->StateProduction[j][i]<<"\t";
				}
			}
		}	
		of<<endl;
	*/
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();

}



void Specie::PutStateProduction(std::string vName,ublas::vector<double> vProductioncm_3s_1)
{
	unsigned pos=0;



	bool is_in_total=true;// True if the state contributes to the total production of the species; it happens if the end is different to -NOTOT
	if(vName.size()>6)
	{
		size_t siz=vName.size();
		string endstr=vName.substr(siz-6);
		if(endstr=="-NOTOT")
		{
			vName=vName.substr(0,siz-6);
			is_in_total=false;
		}
	}

	
	if(!PosInVector(mStates,vName,pos))
	{
		mStates.push_back(vName);
		mStateProductioncm_3s_1.push_back(vProductioncm_3s_1);
		if(is_in_total)
		{
			mTotProductioncm_3s_1=AddOrFullVec(mTotProductioncm_3s_1,vProductioncm_3s_1);
		}
	}else
	{
		if(mStateProductioncm_3s_1[pos].size()==0)
		{
			mStateProductioncm_3s_1[pos]=vProductioncm_3s_1;
		}else
		{
			mStateProductioncm_3s_1[pos]=(vProductioncm_3s_1+mStateProductioncm_3s_1[pos]);
		}
		if(is_in_total)
		{
			mTotProductioncm_3s_1=AddOrFullVec(mTotProductioncm_3s_1,vProductioncm_3s_1);
		}
	}


}
// Deprecated, the state is now in the main system 
void Specie::StatesToTotProd()
{
	Error err("deprecated function","","");
	throw err;
	// We can find strange things here. 
	// I prefer to erase it even if I know it is empty
	//mTotProductioncm_3s_1.erase(mTotProductioncm_3s_1.begin(),mTotProductioncm_3s_1.end());
	mTotProductioncm_3s_1.clear();
	if(mStateProductioncm_3s_1.size()==0) // Quasi-impossible case
		return;

	mTotProductioncm_3s_1=mStateProductioncm_3s_1[0];

	for(unsigned i=1;i<mStateProductioncm_3s_1.size();++i)
	{
		mTotProductioncm_3s_1=(mTotProductioncm_3s_1+mStateProductioncm_3s_1[i]);

	}
}




// Returns the porter energy parameter
ublas::vector<double>* Specie::ReturnEPorter()
{
	return mpEPortereV;
}
// Returns the porter gamm parameter
ublas::vector<double>* Specie::ReturnGPorter()
{
	return mpGPorter;
}
// Returns the porter beta parameter
ublas::vector<double>* Specie::ReturnBPorter()
{
	return mpBPorter;
}
// Returns the porter alpha parameter
ublas::vector<double>* Specie::ReturnAPorter()
{
	return mpAPorter;
}

// Load porter
void Specie::LoadPorter()
{
	mbPorterLoaded=true;
	mpEPortereV=new ublas::vector<double>();
	mpGPorter=new ublas::vector<double>();
	mpBPorter=new ublas::vector<double>();
	mpAPorter=new ublas::vector<double>();
	string name=StrReplace(mName,"+","_PLUS");
	if(mpParams->Exists("/aero_species/"+name+"/EPorter")&&
	   mpParams->Exists("/aero_species/"+name+"/GPorter")&&
	   mpParams->Exists("/aero_species/"+name+"/BPorter")&&
	   mpParams->Exists("/aero_species/"+name+"/APorter"))
	{

		mpParams->Get1DArray("/aero_species/"+name+"/EPorter",*mpEPortereV);
		mpParams->Get1DArray("/aero_species/"+name+"/GPorter",*mpGPorter);
		mpParams->Get1DArray("/aero_species/"+name+"/BPorter",*mpBPorter);
		mpParams->Get1DArray("/aero_species/"+name+"/APorter",*mpAPorter);
	}

}


void Specie::ClearProd()
{
	mTotProductioncm_3s_1.clear();
	mTotLosscm_3s_1.clear();
	mStateDensitycm_3.clear();
	mStateProductioncm_3s_1.clear();
	mStateLosscm_3s_1.clear();
	mProcessNames.clear();
	mCreatedSpecies.clear();
	mMultiplicativeFactor.clear();
	mSpeciesProductioncm_3s_1.clear();
	mTotSpecies.clear();
	mTotSpeciesProductioncm_3s_1.clear();
}

void Specie::ComputePicProd()
{
	if(mTotProductioncm_3s_1.size()==0)
	{
		//Log::SetPriority(Log::WARNING);
		Log::mW<<"The tot production has been found empty during the computation of the maximum position"<<endl;
		mProdPicPosition=0;
		return;
	}
	double maxval=mTotProductioncm_3s_1[0];
	for(size_t i=1;i<mTotProductioncm_3s_1.size();++i)
	{
		if(mTotProductioncm_3s_1[i]>maxval)
		{	
			maxval=mTotProductioncm_3s_1[i];
			mProdPicPosition=i;
		}
	}


}


bool Specie::GetStateDensity(std::string vName,ublas::vector<double>& rDensitycm_3)
{
	unsigned pos=0;
	rDensitycm_3.clear();
	if(vName=="")
	{
		rDensitycm_3=mTotDensitycm_3;
		return true;
	}
	if(PosInVector(mStates,vName,pos))
	{
		if(mStateDensitycm_3.size()>pos)
		{
			rDensitycm_3=mStateDensitycm_3.at(pos);
			return true;
		}else
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"The position of the density does not matche the position of your state!!! "<<mName<<" ("<<vName<<")"<<endl;
			//Log::SetPriority(Log::DEBUGG);
		}
	}
	return false;
}


void Specie::PutStateDensity(std::string vName,ublas::vector<double> vDensitycm_3)
{
	Log::mD<<" Put state density for "<<mName<<" with the state |"<<vName<<"|"<<endl;
	unsigned pos=0;
	if(vName=="")
	{
		Log::mE<<"Error: impossible to overload the main density"<<endl;
		Error err("Impossible to overload the main density with PutStateDensity","Specie::PutStateDensity","");
		throw err;
	}
	if(!PosInVector(mStates,vName,pos))
	{
		pos=mStates.size();
		mStates.push_back(vName);

	}
	Log::mD<<"Insertion of your state density, position : "<<pos<<" with the size of mstate : "<<mStateDensitycm_3.size()<<endl;
	if(mStateDensitycm_3.size()>pos)
	{
		mStateDensitycm_3.at(pos)=vDensitycm_3;
		Log::mD<<"bye bye putstatedensity"<<endl;
		return;
	}else
	{
		Log::mD<<"Resize"<<endl;
		mStateDensitycm_3.resize(pos+1);
		mStateDensitycm_3.at(pos)=vDensitycm_3;
	}
	Log::mD<<"bye bye putstatedensity"<<endl;
}





bool Specie::GetStateProduction(std::string vName,ublas::vector<double>& rProductioncm_3s_1)
{
	unsigned pos=0;
	rProductioncm_3s_1.clear();
	if(PosInVector(mStates,vName,pos))
	{
		if(mStateProductioncm_3s_1.size()>pos)
		{
			rProductioncm_3s_1=mStateProductioncm_3s_1.at(pos);
			return true;
		}
	}
	return false;
}



std::string Specie::ReturnPhotoProbabilityInfo()
{
	string resu = "\t<" + mName + ">\n";

	resu += "\t\t<Electron> " + ntostr(mPhotoElecProductionProbabilitys_1) + "</Electron> \n";
	resu += "\t\t<Total> " + ntostr(mPhotoDissociationProbabilitys_1) + "</Total> \n";
	for(unsigned i = 0; i < mProcessNames.size(); ++i)
	{
		resu += "\t\t<Process" + ntostr(i) + "> \n \t\t<!-- " + mProcessNames[i] + "-->\n";
		resu += ntostr(mSpeciesProductionProbabilitys_1[i]);
		resu += "\t\t</Process" + ntostr(i) + "> \n";
	}
	resu += "\t</" + mName + ">\n";
	return resu;
}
