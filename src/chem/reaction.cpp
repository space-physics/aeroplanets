/**
 * \file reaction.cpp
 * \brief Implements the chemical reaction class ChemReact
 * This class allows to work with the chemical reaction, and the possibility to add the uncertainties, or to modify their value!
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: reaction.cpp 1111 2010-08-12 19:43:33Z gronoff $
 */


#include "reaction.hpp"
using namespace std;

/*
ChemReact::ChemReact(XmlParameters* vpParam): mpParam(vpParam)
{
	//mbIsEmit=false;
//	this->SetId();
	ReadParameters();
}
*/

ChemReact::~ChemReact()
{
}

void ChemReact::ReadParameters()
{
	// Search for the node matching id
	mpParam->ExistsOrDie("/aero_main/chem","The chemistry section /aero_main/chem must be defined in your xml file");
	vector<TiXmlNode*> fluxes = mpParam->GetNodes("/aero_main/chem/chem_reac");

	unsigned i=0;
	bool found = false;

	if(mpParam->Exists("/aero_main/chem/MC_for_all"))
	{
		mbIsMC=true;
	}

	while(!found && i< fluxes.size())
	{
		unsigned val_id=0;
		mpParam->GetNKey(fluxes[i],"/","id",val_id);
		if(val_id==mId)
		{
			if(mpParam->NbParams(fluxes[i],"/","value")>0)
			{// We overload the main value
				mpParam->GetNKey(fluxes[i],"/","value",mMainValue);
				Log::mD<<"Value replacing for "<<mId<<" : "<<mMainValue<<endl;
			}

			if(mpParam->NbParams(fluxes[i],"/","setMC")>0)
			{
				if(mpParam->GetKey(fluxes[i],"/","setMC")=="active")
				{
					mbIsMC=true;
				}
			}

			if(mpParam->NbParams(fluxes[i],"/","uncertainty")>0)
			{
				mUncertainty=mpParam->GetKey(fluxes[i],"/","uncertainty");
			}
		}
		++i;
	}
	mOrigMainValue=mMainValue;

	if(mbIsMC&&mpParam->GetMonteCarlo())
	{
		// Modification of the main value
		double val=0;
		strton(mUncertainty,val);
		if(mUncertainty.find("%")!=string::npos)
		{
			double sigma=val*mMainValue/100.;
			mMainValue=MathRandom::GetNormal(mMainValue,sigma);


		}else
		{
			mMainValue=MathRandom::GetNormal(mMainValue,val);
		}
	}
}


std::string ChemReact::GetInfo()
{
	string info="Reaction "+ntostr(mId)+"\t:\t";
	for(size_t i=0;i<mReactant.size();++i)
	{
		info+=mReactant[i].StandardName();
		if(i!=mReactant.size()-1)
		{
			info+=" + ";
		}
	}

	for(size_t i=0;i<mCatalys.size();++i)
	{
		if(i==0)
		{
			info+=" + ";
		}
		info+=mCatalys[i].StandardName();
		if(i!=mCatalys.size()-1)
		{
			info+=" + ";
		}
	}
	info+="\t\t -> \t ";

	if(mCatalys.size()+mProducts.size() > 0)
	{
		for(size_t i=0;i<mProducts.size();++i)
		{
			info+=ntostr(mEfficiency[mProducts[i]]);
			info+=" x "+mProducts[i].StandardName();
			if(i+1!=mProducts.size())
			{
				info+=" + ";
			}
		}	
		for(size_t i=0;i<mCatalys.size();++i)
		{
			if(i==0)
			{
				info+=" + ";
			}
			info+=mCatalys[i].StandardName();
			if(i!=mCatalys.size()-1)
			{
				info+=" + ";
			}
		}

		if(mbIsEmit)
		{
			info+=" + hv("+ntostr(mEmitFreqnm)+" nm)";
		}


	}
	info+="\t MainVal : "+ntostr(mOrigMainValue)+" ; MainVal used (monte carlo) : "+ntostr(mMainValue)+" ; uncertainty : "+ntostr(mUncertainty)+"\t"+mSupplInfo;
	return info;
}




double ChemReact::GetEfficiency(std::string vName,std::string vState)
{
	double resu=0;
	try
	{
		resu=mEfficiency[SpecieId(vName,vState)];

	}
	catch(...)
	{
		return 0;
	}
	return resu;
}


void CalcDensPhEq::Init( std::deque<Specie*> vNeutralDenscm_3,
		std::deque<Specie*> vTotProdcm_3s_1,
		std::deque<Specie*> vPhotProdcm_3s_1,
		std::deque<Specie*> vElecProdcm_3s_1,
		std::deque<Specie*> vProtProdcm_3s_1,
		std::deque<Specie*> vCosmoProdcm_3s_1,
		ublas::vector<double> vTempNeutreK,
		ublas::vector<double> vTempElecK,
		ublas::vector<double> vTempIonK,
		ublas::vector<double>* vpElecDenscm_3,
		std::deque< boost::shared_ptr<ChemReact> > vChemList
		)
{

	mSize=vTempNeutreK.size();
	mNeutralDenscm_3=vNeutralDenscm_3;
	mTotProdcm_3s_1=vTotProdcm_3s_1;
	mPhotProdcm_3s_1=vPhotProdcm_3s_1;
	mElecProdcm_3s_1=vElecProdcm_3s_1;
	mProtProdcm_3s_1=vProtProdcm_3s_1;
	mCosmoProdcm_3s_1=vCosmoProdcm_3s_1;
	mTempNeutreK=vTempNeutreK;
	mTempElecK=vTempElecK;
	mTempIonK=vTempIonK;
	mpElecDenscm_3=vpElecDenscm_3;
	mChemList=vChemList; // Ok, perfect with shared_ptr


}



ublas::vector<double> CalcDensPhEq::GetDens(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieDensity(vSpecieName,vSpecieState,mNeutralDenscm_3,resu);
	if(resu.size()==0)
	{
		Log::mD<<"The density of the species "<<vSpecieName<<"("<<vSpecieState<<")"<<" has not been found or has not been computed "<<endl;
		resu.resize(mSize);
		resu.clear();
	}
	assert(resu.size()==mSize);
	return resu;
}

ublas::vector<double> CalcDensPhEq::GetTotProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mTotProdcm_3s_1,resu);
	if(resu.size()==0)
	{
		resu.resize(mSize);
		resu.clear();
	}
	assert(resu.size()==mSize);
	return resu;
}


