/**
 * \file planet.cpp
 * \brief Implements the planet abstract class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: planet.cpp 1451 2012-03-20 14:31:42Z gronoff $
 */

#include "planet.hpp"

using namespace std;


Planete::Planete(XmlParameters* pParam)
{
	cout<<"Hello, I am a new planet"<<endl;
	mName="";
	mUA=-1;
	mGms_2=0;
	mLatDegree=0;
	mLoDegree=0;
	mHrLoc=0;
	mSZADegree=0;
	mIsTemperatureDefined=false;

	mpParameter=pParam;
	
	mIsLoadCoordError=false;


}

Planete::Planete(XmlParameters* pParam,double vDistSun)
{
	mUA=vDistSun;
	mGms_2=0;
	mLatDegree=0;
	mLoDegree=0;
	mHrLoc=0;
	mSZADegree=0;
	mpParameter=pParam;
	mName="";
	mIsLoadCoordError=false;
	mIsTemperatureDefined=false;
}
Planete::~Planete()
{
	Log::mL<<"Bye Bye planete"<<endl;
}



void Planete::ShowDistance()
{
	if(mUA<0)
	{
		Log::mL<<"La distance de votre planete au soleil n'a pas ete definie"<<endl;
	}else
	{
		Log::mL<<"Distance au soleil : "<<mUA<<endl;
	}
}

void Planete::ShowName()
{
	if(mName=="")
	{
		Log::mL<<"Le nom de votre planete  n'a pas ete defini"<<endl;
	}else
	{
		Log::mL<<"Votre planete : "<<mName<<endl;
	}
}

void Planete::ShowData()
{
	ShowName();
	ShowDistance();
	Log::mL<<"Gravity : "<<mGms_2<<endl;
	Log::mL<<"Radius : "<<mRKm<<endl;
	Log::mL<<"Latitude : "<<mLatDegree<<endl;
	Log::mL<<"Longitude : "<<mLoDegree<<endl;
	Log::mL<<"SZA :"<<mSZADegree<<endl;
	Log::mL<<"Local hour : "<<mHrLoc<<endl;
}




void Planete::SetCoord(double vLatitudeDegree,double vLongitudeDegree)
{
	mLatDegree=vLatitudeDegree;
	mLoDegree=vLongitudeDegree;


	double lat=mLatDegree;
	double lo=mLoDegree;
	double rlat=lat*PI/180;
	double rlo=lo*PI/180;


	mSZADegree=acos(cos(rlat)*cos(rlo))*180./PI;
	if(lat>90||lat<-90||lo>90||lo<-90)
	{
		/// We are in the nightside
		mSZADegree=360.-mSZADegree;
	}
}


bool Planete::LoadCoords()
{
	bool mIsLat=mpParameter->Exists("/aero_main/planet/lat");
	bool mIsLo=mpParameter->Exists("/aero_main/planet/long");
	bool mIsSZA=mpParameter->Exists("/aero_main/planet/SZA");
	bool mIsHrLoc=mpParameter->Exists("/aero_main/planet/local_hour");

	// This one, copied, is very important!
	// mSZADegree=modulo(mSZADegree,180);

	if(mIsLat&&mIsLo&&!mIsSZA&&!mIsHrLoc)
	{
		mpParameter->GetValue("/aero_main/planet/lat",mLatDegree);
		mpParameter->GetValue("/aero_main/planet/long",mLoDegree);

		if(mLatDegree>90||mLatDegree<-90)
		{
			cout<<"Error : latitude should be in the 90:-90 range"<<endl;
			exit(1);
		}
		// Longitude is 0 at subsolar point, but hour is 0 at midnight
		mHrLoc=(modulo((mLoDegree+180),360))/15;
		mSZADegree=acos(cos(mLatDegree*PI/180.)*cos(mLoDegree*PI/180))*180/PI; // SZA is in the range 0 180
		mSZADegree=modulo(mSZADegree,180);
		return true;
	}

	if(mIsLat&&mIsHrLoc&&!mIsLo&&!mIsSZA)
	{
		mpParameter->GetValue("/aero_main/planet/lat",mLatDegree);
		mpParameter->GetValue("/aero_main/planet/local_hour",mHrLoc);

		mLoDegree=(modulo((mHrLoc+12),24))*15;
		mSZADegree=acos(cos(mLatDegree*PI/180.)*cos(mLoDegree*PI/180))*180/PI; // SZA is in the range 0 180
		mSZADegree=modulo(mSZADegree,180);
		return true;
	}

	if(mIsLat&&mIsSZA&&!mIsLo&&!mIsHrLoc)
	{
		mpParameter->GetValue("/aero_main/planet/lat",mLatDegree);
		mpParameter->GetValue("/aero_main/planet/SZA",mSZADegree);
		if(mLatDegree>89.9 || mLatDegree<-89.9)
		{
			mLoDegree=0;
			mHrLoc=0;
		}

		double tmpcos=cos(mSZADegree*PI/180)/cos(mLatDegree*PI/180);
		if(tmpcos>1.||tmpcos<-1.)
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"Error : your SZA is incompatible with your latitude"<<endl;
			mIsLoadCoordError=true;
			return false;
		}
		//mLoDegree=acos(cos(mSZADegree*PI/180)/cos(mLatDegree*PI/180))*180/PI;
		mLoDegree=acos(tmpcos)*180/PI;
		mHrLoc=(modulo((mLoDegree+180),360))/15;
		mSZADegree=modulo(mSZADegree,180);
		return true;
	}

	if(mIsLo&&mIsSZA&&!mIsLat&&!mIsHrLoc)
	{
		mpParameter->GetValue("/aero_main/planet/long",mLoDegree);
		mpParameter->GetValue("/aero_main/planet/SZA",mSZADegree);
		mSZADegree=modulo(mSZADegree,180);
		mHrLoc=(modulo((mLoDegree+180),360))/15;
		double tmpcos=cos(mSZADegree*PI/180)/cos(mLoDegree*PI/180);
		if(tmpcos>1.||tmpcos<-1.)
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"Error : your SZA is incompatible with your latitude"<<endl;
			mIsLoadCoordError=true;
			return false;
		}
		mLatDegree=acos(tmpcos)*180/PI;
		mSZADegree=modulo(mSZADegree,180);
		return true;
	}

	if(mIsHrLoc&&mIsSZA&&!mIsLat&&!mIsLo)
	{
		mpParameter->GetValue("/aero_main/planet/local_hour",mHrLoc);
		mpParameter->GetValue("/aero_main/planet/SZA",mSZADegree);
		mSZADegree=modulo(mSZADegree,180);
		mLoDegree=(modulo((mHrLoc+12),24))*15;
		double tmpcos=cos(mSZADegree*PI/180)/cos(mLoDegree*PI/180);
		if(tmpcos>1.||tmpcos<-1.)
		{
			cout<<"Error : your SZA is incompatible with your latitude"<<endl;
			mIsLoadCoordError=true;
			return false;
		}
		mLatDegree=acos(tmpcos)*180/PI;
		mSZADegree=modulo(mSZADegree,180);
		return true;
	}

	mIsLoadCoordError=true;
	return false;



}



int Planete::ReaddB_B(const ublas::vector<double>& vAltGridKm, ublas::vector<double>& vrdB_B)
{

	vrdB_B.resize(vAltGridKm.size());
	vrdB_B.clear();
	if(! mpParameter->Exists("/aero_main/planet/dB_B"))
	{
		return 0;
	}
	int type = 0;

	mpParameter->GetNKey("/aero_main/planet/dB_B", "type", type);

	if(type != 0)
		return type;

	mpParameter->ExistsOrDie("/aero_main/planet/dB_B/alt","You have to define the altitude for dB_B");
	mpParameter->ExistsOrDie("/aero_main/planet/dB_B/dB_B","You have to define the value for dB_B");

	ublas::vector<double> lalt;
	mpParameter->Get1DArray("/aero_main/planet/dB_B/alt",lalt);
	ublas::vector<double> lb;
	mpParameter->Get1DArray("/aero_main/planet/dB_B/dB_B",lb);
	vrdB_B = MathFunction::IntLin(lalt, lb, vAltGridKm);

	return 0;
}




