/**
 * \file emit.cpp
 * \brief Implement the emit class: it allows to compute and write the spectrum.
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emit.cpp 1491 2012-05-07 21:35:29Z gronoff $
 */

#include "emit.hpp"
using namespace std;




void Emit::EmitSpectrumIntegrate()
{
#ifdef DEBUG
	Log::mD<<"altgridsize : "<<mAltGridKm.size()<<endl;
	Log::mD<<"colemitsize : "<<mColEmitphcm_3s_1.size()<<endl;
#endif
	if(mAltGridKm.size()!=mColEmitphcm_3s_1.size())
	{
		mColEmitphcm_3s_1.resize(mAltGridKm.size());
		mColEmitphcm_3s_1.clear();
	}
	mIntegratedColEmitR=MathFunction::TrapzInt(mAltGridKm,mColEmitphcm_3s_1)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2)


	std::map<double, ublas::vector<double> >::iterator it;
	for(it=mEmitSpectrumcm_3s_1.begin();it!=mEmitSpectrumcm_3s_1.end();++it)
	{
		if(mAltGridKm.size()==it->second.size())
		{
			mEmitSpectrumIntegratedR[it->first]=MathFunction::TrapzInt(mAltGridKm,it->second)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2)
		}else{
			mEmitSpectrumIntegratedR[it->first]=0;
		}
		mLambdanm.push_back(it->first);
	}

	if(trim(mColumnFileName)!="")
	{
		WriteColumn();
	}
}


void Emit::EmitSpectrumToLimb(boost::shared_ptr<GeoObs> vPath)
{
	if(vPath->IsPath())
	{
		mbPathUsed=true;
		mProfileIntegratedR=vPath->PathIntegrate(mAltGridKm,mColEmitphcm_3s_1)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2)

		if(mbUseAbsorption && (mColumnAbsorption_cm[-1]).size()==mAltGridKm.size())
		{
			if(mAltGridKm.size()==mColEmitphcm_3s_1.size())
				mProfileSpR[-1]=vPath->PathIntegrateAbsorbed(mAltGridKm,mColEmitphcm_3s_1,mColumnAbsorption_cm[-1]*1E5)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2); the 1E5 in absorption is the transformation _cm to km
		}

		std::map<double, ublas::vector<double> >::iterator it;
		for(it=mEmitSpectrumcm_3s_1.begin();it!=mEmitSpectrumcm_3s_1.end();++it)
		{
			if(mbUseAbsorption)
			{
				Log::mL<<"Search if the absorption is defined for "<<it->first<<endl;
			}


			if(!mbUseAbsorption   || (mColumnAbsorption_cm[it->first]).size()!=mAltGridKm.size())
			{
				if(mAltGridKm.size()==it->second.size())
					mProfileSpR[it->first]=vPath->PathIntegrate(mAltGridKm,it->second)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2)
			}else
			{
				Log::mL<<"Absorption taken into account for "<<it->first<<endl;
				if(mAltGridKm.size()==it->second.size())
					mProfileSpR[it->first]=vPath->PathIntegrateAbsorbed(mAltGridKm,it->second,mColumnAbsorption_cm[it->first]*1E5)/10.; // 0.1 = 1E5 (km to cm) / 1E6 (1R = 1E6 ph cm-2); the 1E5 in absorption is the transformation _cm to km

			}
		}
	}

	if(trim(mProfileFileName)!="")
	{
		WriteProfile(vPath);
	}
}



void Emit::WriteColumn()
{
	
	if(FileExists(mColumnFileName))
	{
		Log::mD<<"Low level warning : We overwrite"<<mColumnFileName<<endl;
	}
	ofstream of(mColumnFileName.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Volume profile of "<<mName<<" emissions"<<endl;
	of<<"# Emission in ph cm-3s-1 "<<endl;
	of<<"# Altitude in Km"<<endl;
	// Put the name of the total (+ integrate)
	of<<"# Total : "<<mIntegratedColEmitR<<" R"<<endl;
	// Put the name of the wavelength (+ integrate)
	
	std::map<double,double>::iterator it;
        for(it=mEmitSpectrumIntegratedR.begin();it!=mEmitSpectrumIntegratedR.end();++it)
	{
		of<<"# "<<it->first<<" : "<<it->second<<" R"<<endl;
	}
	// print the profile cm-3s-1
	
	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t"<<mColEmitphcm_3s_1[i]<<"\t";
		for(std::map<double, ublas::vector<double> >::iterator jt=mEmitSpectrumcm_3s_1.begin();jt!=mEmitSpectrumcm_3s_1.end();++jt)
		{
			if(mAltGridKm.size()==jt->second.size())
			{
				of<<(jt->second)[i]<<"\t";
			}else{
				of<<0<<"\t";
			}
		}

		of<<endl;
	}
	

	of<<Log::msMessageLog<<endl;
	of.close();

}

void Emit::WriteProfile(boost::shared_ptr<GeoObs> vPath)
{
	if(!vPath->IsPath())
		return;

	if(FileExists(mProfileFileName))
	{
		Log::mD<<"Low level warning : We overwrite"<<mProfileFileName<<endl;
	}
	ofstream of(mProfileFileName.c_str());
	of<<"# Satellite simulation of the "<<mName<<" emission"<<endl;
	of<<"# Flux in R"<<endl;
	std::deque<double> parameter;
	bool print_parameter=false;
	if(vPath->IsTanDefined())
	{
		of<<"# Tangent altitude in Km"<<endl;
		parameter=vPath->GetTanList();
		print_parameter=true;
	}
	if(vPath->IsDecDefined())
	{
		of<<"# Observation declination in Degree"<<endl;
		parameter=vPath->GetDecList();
		print_parameter=true;
	}

	// Put the name of the total
	of<<"# Total Flux"<<endl;
	// Put the name of the wavelength

	for(std::map<double, ublas::vector<double> >::iterator it=mProfileSpR.begin();it!=mProfileSpR.end();++it)
	{
		if(it->first < 0)
		{
			of<<"# Total (absorbed)"<<endl;
		}else
		{
			of<<"# "<<it->first<<endl;
		}
	}
	// print the integrated profile R
	of.precision(9);
	of.setf(ios::scientific);

	for(unsigned i=0;i< mProfileIntegratedR.size();++i)
	{
		if(print_parameter)
			of<<parameter[i]<<"\t";
		of<<mProfileIntegratedR[i]<<"\t";

		for(std::map<double, ublas::vector<double> >::iterator it=mProfileSpR.begin();it!=mProfileSpR.end();++it)
		{
			of<<(it->second)[i]<<"\t";
		}
		of<<endl;

	}


	of<<Log::msMessageLog<<endl;
	of.close();
}



ublas::vector<double> Emit::GetVER(const int vPos)
{
	if(vPos<0)
	{
		return mColEmitphcm_3s_1;
	}//else
	if(! mLambdanm.size()>vPos)
	{
		Error err("GetVER","Position not found","The emission number "+ntostr(vPos)+" was not found for your species "+mEmittingSpecie.StandardName()+" and at the position "+ntostr(vPos)+" (Check the position before!)");
		throw err;
	}
	return mEmitSpectrumcm_3s_1[mLambdanm[vPos]];
}

ublas::vector<double> Emit::GetLimb(const int vPos)
{


	if(vPos<-1)
	{
		return mProfileIntegratedR;
	}//else
	if(vPos<0)
	{
		if(mProfileSpR[-1].size() > 0)
		{
			return mProfileSpR[-1];
		}
		return mProfileIntegratedR;
	}
	//else
	if(! mLambdanm.size()>vPos)
	{
		Error err("GetVER","Position not found","The emission number "+ntostr(vPos)+" was not found for your species "+mEmittingSpecie.StandardName());
		throw err;
	}
	return mProfileSpR[mLambdanm[vPos]];
}


void Emit::ComputeAbsorption(std::map< std::string, std::map<double,double> > vAbsorptions, boost::shared_ptr<Chem> vChem)
{
	std::map< std::string, std::map<double,double> >::iterator absit;
	std::map<double,double>::iterator lambit;
	
	for(absit = vAbsorptions.begin(); absit != vAbsorptions.end(); ++absit)
	{
		string species = absit->first;
		ublas::vector<double> densite = vChem->GetDens(species,"");
		for(lambit = (absit->second).begin(); lambit != (absit->second).end(); ++lambit)
		{
			if(mColumnAbsorption_cm[lambit->first].size() != mAltGridKm.size())
			{
				mColumnAbsorption_cm[lambit->first] = densite * lambit->second;
			}else
			{
				mColumnAbsorption_cm[lambit->first] += densite * lambit->second;
			}
		}
	}


	//mColumnAbsorption_cm[297.229]=vChem->GetDens("CO2","")*1E-25;
}

