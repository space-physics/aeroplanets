/**
 * \file path.cpp
 * \brief Implements the path class
 * Copyright G Gronoff January 2010
 * Last Modification : $Id: path.cpp 832 2010-01-09 17:07:07Z gronoff $
 */
#include "path.hpp"

Path::Path(ublas::vector<double> vAltKm,ublas::vector<double> vLenKm,ublas::vector<double> vSZADeg):mAltitudeKm(vAltKm),mLengthKm(vLenKm),mSZADegree(vSZADeg)
{
}


void Path::InitAltLen(ublas::vector<double> vAltitudeKm, ublas::vector<double> vPathKm,double vSZADegree)
{

	if(vAltitudeKm.size()!=vPathKm.size())
	{
		Error err("Path::InitAltLen","Error in the size of the vectors","Your altitude vector and your length vector have not the same size");
		throw err;
	}
	mAltitudeKm=vAltitudeKm;
	mLengthKm=vPathKm;
	mSZADegree.resize(vPathKm.size());
	for(unsigned i=0;i<vPathKm.size();++i)
	{
		mSZADegree[i]=vSZADegree;
	}

}

void Path::InitAltLenSZA(ublas::vector<double> vAltitudeKm, ublas::vector<double> vPathKm,ublas::vector<double> vSZADegree)
{

	if((vAltitudeKm.size()!=vPathKm.size()) || (vAltitudeKm.size()!=vSZADegree.size()))
	{
		Error err("Path::InitAltLenSZA","Error in the size of the vectors","Your altitude, length, and SZA vectors have not the same size");
		throw err;
	}
	mAltitudeKm=vAltitudeKm;
	mLengthKm=vPathKm;
	mSZADegree=vSZADegree;
}



void Path::CartesianToPath(ublas::vector<double> vXRp,ublas::vector<double> vYRp,ublas::vector<double> vZRp,double vRpKm)
{
//	assert(vXRp.size()==vYRp.size());
//	assert(vXRp.size()==vZRp.size());

	if((vXRp.size()!=vYRp.size()) || (vXRp.size()!=vZRp.size()))
	{
		Error err("Path::CartesianToPath","Error in the size of the vectors","Your X,Y, and Z vectors have not the same size");
		throw err;
	}
	size_t siz=vXRp.size();
	
	mAltitudeKm.resize(siz);
	mLengthKm.resize(siz);
	mSZADegree.resize(siz);


	mLengthKm[0]=0;
	mSZADegree[0]=atan2(sqrt(vYRp[0]*vYRp[0]+vZRp[0]*vZRp[0]),vXRp[0])*180/PI;
	mAltitudeKm[0]=(sqrt(vYRp[0]*vYRp[0]+vZRp[0]*vZRp[0]+vXRp[0]*vXRp[0])-1)*vRpKm;

	for(size_t i=1;i<siz;++i)
	{
		mLengthKm[i]=sqrt((vYRp[i]-vYRp[i-1])*(vYRp[i]-vYRp[i-1])
				+(vZRp[i]-vZRp[i-1])*(vZRp[i]-vZRp[i-1])
				+(vXRp[i]-vXRp[i-1])*(vXRp[i]-vXRp[i-1])
				)*vRpKm;
		mSZADegree[i]=atan2(sqrt(vYRp[i]*vYRp[i]+vZRp[i]*vZRp[i]),vXRp[i])*180/PI;
		mAltitudeKm[i]=(sqrt(vYRp[i]*vYRp[i]+vZRp[i]*vZRp[i]+vXRp[i]*vXRp[i])-1)*vRpKm;
	}

	for(size_t i=0;i<siz;++i)
	{
		mLengthKm[i]-=mLengthKm[siz-1];
	}
}

