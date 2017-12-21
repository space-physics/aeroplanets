/**
 * \file geopath.cpp
 * \brief Implements the path of a line of sight
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geopath.cpp 1595 2012-11-03 01:30:26Z gronoff $
 */


#include "geopath.hpp"



void GeoPath::PointsPathToArrays()
{
	assert(mNbPoints==mPointsPath.size());
	mAltitudeKm.resize(mNbPoints);
	mAltitudeKm.clear();
	mLengthKm.resize(mNbPoints);
	mLengthKm.clear();
	mSZADegree.resize(mNbPoints);
	mSZADegree.clear();

	mLengthKm[0]=0;
	mAltitudeKm[0]=mPointsPath[0].GetAlt();
	mSZADegree[0]=mPointsPath[0].GetSZADegree();

	for(size_t i=1;i<mNbPoints;++i)
	{
		GeoPoint pas=mPointsPath[i]-mPointsPath[i-1];
		mLengthKm[i]=mLengthKm[i-1]+pas.GetDist();
		mAltitudeKm[i]=mPointsPath[i].GetAlt();
		mSZADegree[i]=mPointsPath[i].GetSZADegree();
	}


}

void GeoPath::StraightExtremas(unsigned vNbPoints)
{
	mNbPoints=vNbPoints;
	assert(vNbPoints>1);
//	GeoPoint tmp=;

	mAltitudeKm.resize(vNbPoints);
	mAltitudeKm.clear();
	mLengthKm.resize(vNbPoints);
	mLengthKm.clear();
	mSZADegree.resize(vNbPoints);
	mSZADegree.clear();

	GeoPoint pas=(mEndPoint-mStartPoint)*(1/static_cast<double>(vNbPoints-1));
	mPointsPath.push_back(mStartPoint);
	mLengthKm[0]=0;
	mAltitudeKm[0]=mStartPoint.GetAlt();
	mSZADegree[0]=mStartPoint.GetSZADegree();

	double paslong=pas.GetDist();
	for(size_t i=1;i<vNbPoints;++i)
	{
		mPointsPath.push_back(mPointsPath[i-1]+pas);
		mLengthKm[i]=mLengthKm[i-1]+paslong;
		mAltitudeKm[i]=mPointsPath[i].GetAlt();
		mSZADegree[i]=mPointsPath[i].GetSZADegree();
	}
}

GeoPath::GeoPath(GeoPoint vStartPoint,GeoPoint vEndPoint,unsigned vNbPoints): mStartPoint(vStartPoint),mEndPoint(vEndPoint),mNbPoints(vNbPoints)
{
	mbIsPathDefined=true;
	StraightExtremas(vNbPoints);
}

GeoPath::GeoPath(ublas::vector<double> vAltGridKm, double vSZAdeg, double vRplanetKm) : mStartPoint(GeoPoint(vRplanetKm)), mEndPoint(GeoPoint(vRplanetKm))
{
	assert(vAltGridKm.size() > 1);
	assert(vAltGridKm[0] > vAltGridKm[vAltGridKm.size() - 1]);
	
//	GeoPoint end(vRplanetKm);
//	mEndPoint = GeoPoint(vRplanetKm);
	mEndPoint.SetGeoPosition(vAltGridKm[vAltGridKm.size() - 1], 0, 0);
	mStartPoint = mEndPoint.ReturnAllSortieAtmo(0, 90 - vSZAdeg, vAltGridKm[0]); // We can have problems here with SZA > 90... We will have to check that!
	// First, the length should not be inverted (if SZA > 90 we have 90 - SZA)
	// Second, if the altitude in the path goes below 0km, we should go dark
	// Third, we should have some warning if the path goes below the minimum altitude
	mbIsPathDefined=true;
	StraightExtremas(vAltGridKm.size() * 2);
	ResetGrid(vAltGridKm);
}

GeoPath::GeoPath(GeoPoint vStartPoint): mStartPoint(vStartPoint),mEndPoint(vStartPoint)
{
	mbIsPathDefined=false;
}

void GeoPath::InitAzTanPath(double vAzimuthDegree,double vTanAltKm,double vOutAltKm,unsigned vNbPoints)
{
	
	mEndPoint=mStartPoint.ReturnSortieAtmoHo(vAzimuthDegree,vTanAltKm,vOutAltKm);
	mbIsPathDefined=true;
	StraightExtremas(vNbPoints);

}

void GeoPath::InitAzDecPath(double vAzimuthDegree,double vDecDegree,double vOutAltKm,unsigned vNbPoints)
{
	
	mEndPoint=mStartPoint.ReturnAllSortieAtmo(vAzimuthDegree,vDecDegree,vOutAltKm);
	mbIsPathDefined=true;
	StraightExtremas(vNbPoints);
}




double GeoPath::PathIntegrate(ublas::vector<double> vAltKm, ublas::vector<double> vVals)
{
	if(!mbIsPathDefined || 0==mLengthKm.size())
	{
		Error err("GeoPath::PathIntegrate","Path not defined","You try to make an integration while your path is not correctly defined");
		throw err;
	}

	ublas::vector<double> values=MathFunction::IntLin(vAltKm,vVals,mAltitudeKm);
	/*
#define NOTESTING
#ifndef NOTESTING
	Log::mL<<" You length values : "<<mLengthKm<<std::endl;
	// juste un test a la noix
	ublas::vector<double> test(mLengthKm.size());
	ublas::vector<double> ones(mLengthKm.size());
	for(unsigned i=0;i<mLengthKm.size();++i)
	{
		test[i]=0.;
		ones[i]=1.;
	}

	for(unsigned i=0;i<mLengthKm.size();++i)
	{
		//unsigned nvals=mLengthKm.size()-i-1;
		//test[nvals]=MathFunction::TrapzInt(mLengthKm,ones,nvals);
		test[i]=exp(-MathFunction::TrapzInt(mLengthKm,ones,i));
	}
	Log::mL<<" Your test values : "<<test<<std::endl;
	Log::mL<<" Your values : "<<values<<std::endl;
	Log::mL<<" Your funny value : "<<test*values<<std::endl;

#endif
*/
	return MathFunction::TrapzInt(mLengthKm,values);
}

double GeoPath::PathIntegrateAbsorbed(ublas::vector<double> vAltKm, ublas::vector<double> vVals,ublas::vector<double> vAbsorption_Km)
{
	assert(vAltKm.size()==vVals.size());
	assert(vAltKm.size()==vAbsorption_Km.size());
	
	if(!mbIsPathDefined || 0==mLengthKm.size())
	{
		Error err("GeoPath::PathIntegrate","Path not defined","You try to make an integration while your path is not correctly defined");
		throw err;
	}

	ublas::vector<double> values=MathFunction::IntLin(vAltKm,vVals,mAltitudeKm);
	ublas::vector<double> absorb_values=MathFunction::IntLin(vAltKm,vAbsorption_Km,mAltitudeKm);
	ublas::vector<double> absorb(mLengthKm.size());
	
	// absorb_values contains sigma*N
	for(unsigned i=0;i<mLengthKm.size();++i)
	{ // we compute exp(- \int sigma N dl) where l is the path position
		absorb[i]=exp(-MathFunction::TrapzInt(mLengthKm,absorb_values,i));
	}
	// This absorption is multiplied by the values
	values=values*absorb; 
	// The integration of the values*absorption gives the final value
	return MathFunction::TrapzInt(mLengthKm,values);
}


void GeoPath::TestNewPath(ublas::vector<double> vAltListKm)
{
	GeoPoint pas=(mEndPoint-mStartPoint);
//	double lengthpas = pas.GetDist();
	Log::mL<< " Altitude list length"<<vAltListKm.size()<<std::endl;
	Log::mL<< "Initial altitude list length"<< mAltitudeKm.size()<<std::endl;
	Log::mL<< mLengthKm<<std::endl;
	Log::mL<< "Mlength "<< mLengthKm[mLengthKm.size() - 1]<<std::endl;
	ublas::vector<double> newLength = MathFunction::IntLin(mAltitudeKm, mLengthKm, vAltListKm) /  mLengthKm[mLengthKm.size() - 1]; //(*(mLengthKm.end()));
	Log::mL<< " Newlength list length"<<newLength.size()<<std::endl;
	Log::mL<<"==============================="<<std::endl;
	Log::mL<<"==============================="<<std::endl;
	Log::mL<<" We work on the altitude list"<<std::endl;
	Log::mL<<"==============================="<<std::endl;
	for(unsigned i = 0; i < vAltListKm.size(); ++i)
	{
		GeoPoint newpoint = mStartPoint + pas * newLength[i];
		Log::mL<< vAltListKm[i] << "\t" << newpoint.GetAlt() << "\t"<< vAltListKm[i] - newpoint.GetAlt() << "\t" << (newpoint - mStartPoint).GetDist() << std::endl; 
	}
	Log::mL<<"==============================="<<std::endl;

//		mLengthKm[i]=mLengthKm[i-1]+paslong;
//		mAltitudeKm[i]=mPointsPath[i].GetAlt();
}


void GeoPath::ResetGrid(ublas::vector<double> vAltListKm)
{
	mNbPoints=vAltListKm.size();
	assert(mNbPoints>1);
	//	GeoPoint tmp=;
	for(unsigned i = 0; i < mAltitudeKm.size(); ++i)
	{
		Log::mL<< mAltitudeKm[i] << "\t" << mLengthKm[i] <<  std::endl;
	}
	
	ublas::vector<double> newLength = MathFunction::IntLin(mAltitudeKm, mLengthKm, vAltListKm) /  mLengthKm[mLengthKm.size() - 1];
	mAltitudeKm.resize(mNbPoints);
	mAltitudeKm.clear();
	mLengthKm.resize(mNbPoints);
	mLengthKm.clear();
	mSZADegree.resize(mNbPoints);
	mSZADegree.clear();

	GeoPoint pas=(mEndPoint-mStartPoint);
	mPointsPath.clear();
	for(unsigned i = 0; i < vAltListKm.size(); ++i)
	{
		GeoPoint newpoint = mStartPoint + pas * newLength[i];
		Log::mL<< vAltListKm[i] << "\t" << newpoint.GetAlt() << "\t"<< vAltListKm[i] - newpoint.GetAlt() << "\t" << (newpoint - mStartPoint).GetDist() << std::endl;
		mPointsPath.push_back(newpoint);
		mLengthKm[i]=(newpoint - mStartPoint).GetDist() ;
		mAltitudeKm[i]=newpoint.GetAlt();
		mSZADegree[i]=newpoint.GetSZADegree();
	}
	assert(mSZADegree.size() == mAltitudeKm.size());
	assert(mSZADegree.size() == mPointsPath.size());
	assert(mSZADegree.size() == mLengthKm.size());
}
