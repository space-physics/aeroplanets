/**
 * \file geoobs.cpp
 * \brief  Implements the  experiment class : the set of observations. Can be seen as a collection of line of sight
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geoobs.cpp 1111 2010-08-12 19:43:33Z gronoff $
 */

#include "geoobs.hpp"
using namespace std;
GeoObs::GeoObs(XmlParameters* vpParams,double vRpKm,double vOutAtmoKm):mpParam(vpParams),mRpKm(vRpKm),mOutAtmoKm(vOutAtmoKm)
{
	mbIsSatFixed=false;
	mbIsTanDefined=false;
	mbIsDecDefined=false;
	mNbPoints=0;
	mbIsGeoUsed=false;
	if(mpParam->Exists("/aero_main/geometry/use_geometry"))
	{
		mbIsGeoUsed=true;
		ReadParam();
		CreatePath();
	}
}



void GeoObs::ReadParam()
{

	mpParam->ExistsOrDie("/aero_main/geometry/satellite","You must define an observation  position if you want to work with the geometry");
	mpParam->ExistsOrDie("/aero_main/geometry/satellite/altitude","You must define the altitude of your satellite");
	mpParam->ExistsOrDie("/aero_main/geometry/satellite/latitude","You must define the latitude of your satellite");
	mpParam->ExistsOrDie("/aero_main/geometry/satellite/longitude","You must define the longitude of your satellite");

	mSatAltitudes=ReadNode(mpParam->GetNode("/aero_main/geometry/satellite/altitude"));
	mSatLatDegree=ReadNode(mpParam->GetNode("/aero_main/geometry/satellite/latitude"));
	mSatLonDegree=ReadNode(mpParam->GetNode("/aero_main/geometry/satellite/longitude"));

	if(mpParam->Exists("/aero_main/geometry/atmosphere_limit"))
	{
		mpParam->GetValue("/aero_main/geometry/atmosphere_limit",mOutAtmoKm);
	}

	if(mSatAltitudes.size()==1 && mSatLatDegree.size()==1 && mSatLonDegree.size()==1)
	{
		mbIsSatFixed=true;
	}else
	{
		unsigned size=nmax(mSatAltitudes.size(),mSatLatDegree.size());
		size=nmax(size,static_cast<unsigned>(mSatLonDegree.size()));
		if(mSatAltitudes.size()==1)
		{
			for(size_t i=1;i<size;++i)
			{
				mSatAltitudes.push_back(mSatAltitudes[0]);
			}
		}
		if(mSatLatDegree.size()==1)
		{
			for(size_t i=1;i<size;++i)
			{
				mSatLatDegree.push_back(mSatLatDegree[0]);
			}
		}
		if(mSatLonDegree.size()==1)
		{
			for(size_t i=1;i<size;++i)
			{
				mSatLonDegree.push_back(mSatLonDegree[0]);
			}
		}
		if(mSatAltitudes.size()!=size)
		{
			Error err("Satellite altitude","Size mismatch"," your satellite altitude size is not  in correspondance with the other sizes");
			throw err;
		}	
		if(mSatLatDegree.size()!=size)
		{
			Error err("Satellite latitude","Size mismatch"," your satellite latitude size is not  in correspondance with the other sizes");
			throw err;
		}	
		if(mSatLonDegree.size()!=size)
		{
			Error err("Satellite longitude","Size mismatch"," your satellite longitude size is not  in correspondance with the other sizes");
			throw err;
		}	
	}


	mpParam->ExistsOrDie("/aero_main/geometry/los","You must define an observation  direction if you want to work with the geometry");
	mpParam->ExistsOrDie("/aero_main/geometry/los/azimut","You must define an observation azimut if you want to work with the geometry");

	if(mpParam->Exists("/aero_main/geometry/los/tangent_alt"))
	{
		mbIsTanDefined=true;
	}
	if(mpParam->Exists("/aero_main/geometry/los/declination"))
	{
		mbIsDecDefined=true;
	}
	if(mbIsTanDefined==mbIsDecDefined)
	{
		Error err("Read los","too much defined","the los/tangent_alt and los/declination cannot be both defined. But one has to be defined. Either the tangent altitude of the line of sight or the declination of the line of Sight");
		throw err;
	}


	mAzimut=ReadNode(mpParam->GetNode("/aero_main/geometry/los/azimut"));
	if(mbIsTanDefined)
	{
		mTanAltKm=ReadNode(mpParam->GetNode("/aero_main/geometry/los/tangent_alt"));
	}
	if(mbIsDecDefined)
	{
		mDec=ReadNode(mpParam->GetNode("/aero_main/geometry/los/declination"));
	}
	unsigned size=nmax(mDec.size(),mTanAltKm.size());
	size=nmax(size,static_cast<unsigned>(mAzimut.size()));
	size=nmax(size,static_cast<unsigned>(mSatLonDegree.size()));
	mNbPoints=size;


	if(mAzimut.size()==1)
	{
		for(size_t i=1;i<size;++i)
		{
			mAzimut.push_back(mAzimut[0]);
		}
	}
	if(mbIsTanDefined&&mTanAltKm.size()==1)
	{
		for(size_t i=1;i<size;++i)
		{
			mTanAltKm.push_back(mTanAltKm[0]);
		}
	}
	if(mbIsDecDefined&&mDec.size()==1)
	{
		for(size_t i=1;i<size;++i)
		{
			mDec.push_back(mDec[0]);
		}
	}
	if(mAzimut.size()!=size)
	{
		Error err("Azimut","Size mismatch"," your azimut size is not  in correspondance with the other sizes");
		throw err;
	}	
	if(mbIsTanDefined&&mTanAltKm.size()!=size)
	{
		Error err("Tangent altitude","Size mismatch"," your tangent altitude size is not  in correspondance with the other sizes");
		throw err;
	}	
	if(mbIsDecDefined&&mDec.size()!=size)
	{
		Error err("Declination","Size mismatch"," your declination size is not  in correspondance with the other sizes");
		throw err;
	}	

}


std::deque<double> GeoObs::ReadNode(TiXmlNode* vNode)
{
	std::deque<double> resu;

	if(mpParam->Exists(vNode,"//const"))
	{
		double value=0;
		mpParam->GetValue(vNode,"//const",value);
		resu.push_back(value);
	}else if(mpParam->Exists(vNode,"//list"))
	{
		ublas::vector<double> tmp;
		mpParam->Get1DArray(vNode,"//list",tmp);
		for(size_t i=0;i<tmp.size();++i)
		{
			resu.push_back(tmp[i]);
		}
	}else if(mpParam->Exists(vNode,"//range"))
	{
		return ReadRange(mpParam->GetNode(vNode,"//range"));
	}else
	{
		Error err("Value undefined","Your value in the node is undefined","Please check your satellite and los definitions");
	}
	return resu;
}


std::deque<double> GeoObs::ReadRange(TiXmlNode* vNode)
{

	std::deque<double> resu;
	
	mpParam->ExistsOrDie(vNode,"//start","You must define the starting value");
	mpParam->ExistsOrDie(vNode,"//end","You must define the ending value");
	mpParam->ExistsOrDie(vNode,"//number","You must define the number of points");

	double start=0,end=0;
	unsigned number=0;
	mpParam->GetValue(vNode,"//start",start);
	mpParam->GetValue(vNode,"//end",end);
	mpParam->GetValue(vNode,"//number",number);

	assert(number>1);
	for(unsigned i=0;i<number;++i)
	{
		double tmpval=(end-start)/(static_cast<double>(number)-1)*static_cast<double>(i)+start;

		resu.push_back(tmpval);
	}


	return resu;
}

void GeoObs::CreatePath()
{

#ifdef DEBUG
	Log::mD<<mAzimut.size()<<"\t"<<mTanAltKm.size()<<"\t"<<mDec.size()<<std::endl;
#endif 
	assert(mbIsTanDefined == (mAzimut.size()==mTanAltKm.size()));
	assert(mbIsDecDefined == (mAzimut.size()==mDec.size()));
	assert(mbIsTanDefined!=mbIsDecDefined);
	assert(mNbPoints>0);
	assert(mSatAltitudes.size()==mSatLatDegree.size());
	assert(mSatAltitudes.size()==mSatLonDegree.size());
	assert(mbIsSatFixed!=(mSatAltitudes.size()==mAzimut.size()));
	for(size_t i=0;i<mNbPoints;++i)
	{
		GeoPoint tmpsatpos(mRpKm);
		if(!mbIsSatFixed)
		{
			tmpsatpos.SetGeoPosition(mSatAltitudes[i],mSatLatDegree[i],mSatLonDegree[i]);
		}else
		{
			tmpsatpos.SetGeoPosition(mSatAltitudes[0],mSatLatDegree[0],mSatLonDegree[0]);
		}
		GeoPath tmpgeopath(tmpsatpos);
		if(mbIsDecDefined)
		{
			tmpgeopath.InitAzDecPath(mAzimut[i],mDec[i],mOutAtmoKm,mNbPoints);
		}
		if(mbIsTanDefined)
		{
			tmpgeopath.InitAzTanPath(mAzimut[i],mTanAltKm[i],mOutAtmoKm,mNbPoints);
		}
		mPaths.push_back(tmpgeopath);
	}



}



ublas::vector<double> GeoObs::PathIntegrate(ublas::vector<double> vAltGridKm,ublas::vector<double> vValue)
{
	assert(vAltGridKm.size()==vValue.size());
	assert(mPaths.size()==mNbPoints);
	ublas::vector<double> resu(mPaths.size());
	resu.clear();
	for(size_t i=0;i<mPaths.size();++i)
	{
		resu[i]=mPaths[i].PathIntegrate(vAltGridKm,vValue);
	}
	return resu;
	
}
ublas::vector<double> GeoObs::PathIntegrateAbsorbed(ublas::vector<double> vAltGridKm,ublas::vector<double> vValue,ublas::vector<double> vAbsorb_Km)
{
	assert(vAltGridKm.size()==vValue.size());
	assert(vAltGridKm.size()==vAbsorb_Km.size());
	assert(mPaths.size()==mNbPoints);
	ublas::vector<double> resu(mPaths.size());
	resu.clear();
	for(size_t i=0;i<mPaths.size();++i)
	{
		resu[i]=mPaths[i].PathIntegrateAbsorbed(vAltGridKm,vValue,vAbsorb_Km);
	}
	return resu;
	
}
