/**
 * \file geopoint.cpp
 * \brief Implements the geographical point 
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geopoint.cpp 1595 2012-11-03 01:30:26Z gronoff $
 */

#include "geopoint.hpp"

/*       _\|/_
         (o o)
 +----oOO-{_}-OOo---------------------------------------+
 |                                                      |
 |Automatic conversion between the different coordinates|
 |                                                      |
 +-----------------------------------------------------*/
void GeoPoint::GeoToCar()
{
	mXKm=mRKm*cos(mLatrad)*cos(mLonrad);
	mYKm=mRKm*cos(mLatrad)*sin(mLonrad);
	mZKm=mRKm*sin(mLatrad);
}
void GeoPoint::CarToGeo()
{
	mRKm=sqrt(mXKm*mXKm+mYKm*mYKm+mZKm*mZKm);
	mLatrad=asin(mZKm/mRKm);
	mLonrad=atan2(mYKm/mRKm,mXKm/mRKm);
	mAltKm=mRKm-mRpKm;
}
void GeoPoint::PositionToSZA()
{
	mSZArad=acos(cos(mLatrad)*cos(mLonrad));
}
/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |Initialization  |
 |                |
 +---------------*/
GeoPoint::GeoPoint(double vRpKm): mRpKm(vRpKm)
{
}

GeoPoint GeoPoint::Copy()
{
	GeoPoint resu(mRpKm);
	resu.SetCartPosition(mXKm,mYKm,mZKm);
	return resu;
}
void GeoPoint::SetGeoPosition(double vAltKm,double vLatDegree,double vLonDegree)
{
	mAltKm=vAltKm;
	mLatrad=vLatDegree*PI/180.;
	mLonrad=vLonDegree*PI/180.;
	mRKm=mRpKm+mAltKm;
	GeoToCar();
	PositionToSZA();
}
void GeoPoint::SetCartPosition(double vXKm, double vYKm, double vZKm)
{
	mXKm=vXKm;
	mYKm=vYKm;
	mZKm=vZKm;
	CarToGeo();
	PositionToSZA();
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo------------------------------------------------+
 |                                                               |
 |Mainly internal function, let public because it could be useful|
 |                                                               |
 +--------------------------------------------------------------*/
double GeoPoint::EcartAngleDegree(GeoPoint vC)
{
	std::deque<double> mp=GetCartesianReduced();
	std::deque<double> cp=vC.GetCartesianReduced();

	double d=sqrt((mp[0]-cp[0])*(mp[0]-cp[0])+(mp[1]-cp[1])*(mp[1]-cp[1])+ (mp[2]-cp[2])*(mp[2]-cp[2]));
	return 2.*asin(d/2.)*180./PI;
}
double GeoPoint::ReturnTangentPointAltitude(double vLatDegree, double vLongDegree)
{
	GeoPoint tmppoint(mRpKm);
	tmppoint.SetGeoPosition(0,vLatDegree,vLongDegree);
	double angle=EcartAngleDegree(tmppoint);
	double newr=cos(angle*PI/180.)*mRKm;
	return newr-mRpKm;
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-------------------+
 |                                  |
 |Returns the different LOS geoPoint|
 |                                  |
 +---------------------------------*/

GeoPoint GeoPoint::ReturnTrajPoint(double vAzimuthDegree)
{
	double az=vAzimuthDegree*PI/180.;
	double newlat=asin(cos(mLatrad)*cos(az));
	double delong=0;
	if(0==newlat)
	{
		delong=vAzimuthDegree;
	}else
	{
		delong=sin(az)/cos(newlat);
		if(nabs( nabs(delong)-1 )<1E-5)
		{
			if(delong>0)
			{// in degree
				delong=90;
			}else
			{
				delong=-90;
			}
		}else
		{
			delong=asin(delong)*180/PI;
		}
	}
	newlat*=180./PI;
	GeoPoint resu(mRpKm);
	resu.SetGeoPosition(mAltKm,newlat,delong);
	return resu;
}

GeoPoint GeoPoint::ReturnSubPoint(double vAzimuthDegree,double vHoKm)
{
	GeoPoint resu(mRpKm);
	double hdr = (vHoKm + mRpKm) / mRKm;
	if(hdr>1.)
	{
		Error err("GeoPoint::ReturnSubPoint","hdr>1.: trying to find a tangent point while looking upside","The altitude requested for the tangent point is above the altitude of the observer, which is impossible");
		throw err;
	}
	double newlat=asin( sin(mLatrad)*hdr+cos(mLatrad)*sqrt(1-hdr*hdr)*cos(vAzimuthDegree*PI/180.));
	double delong=mLonrad;
	if(newlat!=0.)
	{// quasi-always true: the subsolar point is unlikely to be exactly at the pole
		delong+=asin(sin(vAzimuthDegree*PI/180.)/cos(newlat)*sqrt(1-hdr*hdr));
	}
	resu.SetGeoPosition(vHoKm,newlat*180./PI,delong*180./PI);
	return resu;
}

GeoPoint GeoPoint::ReturnSubViz(double vAzimuthDegree, double vDecDegree)
{
	if(vDecDegree>0)
	{

		Error err("GeoPoint::ReturnSubViz","vDecDegree>0.: trying to find a tangent point while looking upside","The declination for the line of sight should be negative (looking downside) or = 0");
		throw err;
	}

	double Ho=cos(vDecDegree * PI / 180.) * mRKm - mRpKm;
	if(nabs(vDecDegree) < 1E-5)
		Ho=mAltKm;
	// GeoPoint resu(mRpKm);
	//return resu.ReturnSubPoint(vAzimuthDegree,Ho);
	return ReturnSubPoint(vAzimuthDegree, Ho);
}



GeoPoint GeoPoint::ReturnSortieAtmo(double vAzimuthDegree, double vDecDegree, double vAltOutKm)
{
	GeoPoint tmp(mRpKm);
	double dout = 0;
	double rout = mRpKm + vAltOutKm;
	if(nabs(vDecDegree) < 1E-5)
	{ // Case: we look at 90 Degree
		tmp = ReturnSubViz(vAzimuthDegree, -45);
		double newaltitude = tmp.GetDist() * 2 - tmp.GetRp();
		tmp.SetGeoPosition(newaltitude, tmp.GetLatDegree(), tmp.GetLonDegree());
		dout = sqrt(rout * rout - mRKm * mRKm);
	}else
	{
		tmp = ReturnSubViz(vAzimuthDegree, vDecDegree);
		double rho = tmp.GetDist();
		dout = sqrt(mRKm * mRKm - rho * rho) + sqrt(rout * rout - rho * rho);
	}
	GeoPoint vectordirection = tmp - *this;
	double dis = vectordirection.GetDist();
	Log::mL<<"DIS "<< dis <<std::endl; 
	Log::mL<<"DOUT "<< dout <<std::endl; 
	double multiplicator = dout / dis;
	//GeoPoint resu=Copy();
	return *this + vectordirection * multiplicator;
}

GeoPoint GeoPoint::ReturnAllSortieAtmo(double vAzimuthDegree, double vDecDegree, double vAltOutKm)
{
	if(vDecDegree<=0)
	{// Standard case
		return ReturnSortieAtmo(vAzimuthDegree,vDecDegree,vAltOutKm);
	}// We are looking above: we must invert the computation

	GeoPoint tmp = ReturnSubViz(vAzimuthDegree + 180., -vDecDegree); // Tangent point in the other direction
	GeoPoint vectordirection= *this - tmp;// The direction of the vector is inverted with respect to the tangent point
	double rout = mRpKm + vAltOutKm;
	double rho = tmp.GetDist();
	double dout = sqrt(rout * rout - rho * rho)-sqrt(mRKm * mRKm - rho * rho);// This is the last difference with the dec<0 situation, the point to the tangent point distance is substracted to the tangent point to outside distance
	double dis = vectordirection.GetDist();
	double multiplicator = dout / dis;
	return *this + vectordirection * multiplicator;
}



GeoPoint GeoPoint::ReturnSortieAtmoHo(double vAzimuthDegree, double vHoKm, double vAltOutKm)
{
	
	double rho=mRpKm+vHoKm;
	double rout=mRpKm+vAltOutKm;
	GeoPoint tmp= ReturnSubPoint(vAzimuthDegree,vHoKm);
	double dout=sqrt(mRKm*mRKm-rho*rho)+sqrt(rout*rout-rho*rho);
	GeoPoint vectordirection=tmp-*this;
	double dis=vectordirection.GetDist();
	double multiplicator=dout/dis;
	return *this+vectordirection*multiplicator;
}
