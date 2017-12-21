/**
 * \file earth.cpp
 * \brief Implements the earth planet
 * Copyright G Gronoff Nov 2009
 * Last Modifcation $Id: earth.cpp 1451 2012-03-20 14:31:42Z gronoff $
 */



#include "earth.hpp"
using namespace std;
Earth::Earth(XmlParameters* pParam) : Planete(pParam)
{
	mName="Earth";
	mUA=1.;
	mGms_2=9.81;
	mRKm=6370;

	mIsTemperatureDefined=false;
	mbIsModelLoaded=false;
	LoadCoords();
	



}

Earth::Earth(XmlParameters* pParam,double vUA) : Planete(pParam,vUA)
{
	mName="Earth";
	mGms_2=9.81;
	mRKm=6370;
	mIsTemperatureDefined=false;

	LoadCoords();
	mbIsModelLoaded=false;
}
Earth::~Earth()
{
	if(mbIsModelLoaded)
	{
		delete mEarth;
	}
}
std::map< std::string, ublas::vector<double> > Earth::AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType) 
{
	map<string, ublas::vector<double> > resultat;

	if(vType==1)
	{
		if(!mbIsModelLoaded)
		{
			LoadEarth();
		}

		
		deque< ublas::vector<double> > allatmo=mEarth->ReturnNeutralAtmo(vAltitudeGridKm);
		deque< ublas::vector<double> > alliono=mEarth->ReturnIono(vAltitudeGridKm);
		mTemperatureModelGridK=mEarth->ReturnNTemp(vAltitudeGridKm);

		mIsTemperatureDefined=true;

		deque<string>::iterator it;
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			if(*it=="N2")
			{
				resultat[*it]=allatmo[0];
			}else if(*it=="O2")
			{
				resultat[*it]=allatmo[1];
			}else if(*it=="O")
			{
				resultat[*it]=allatmo[2];
			}else if(*it=="H")
			{
				resultat[*it]=allatmo[3];
			}else if(*it=="He")
			{
				resultat[*it]=allatmo[4];
			}else if(*it=="N")
			{
				resultat[*it]=allatmo[5];
			}else if(*it=="A")
			{
				resultat[*it]=allatmo[5];
			}else if(*it=="NO")
			{
				resultat[*it]=allatmo[5];
			}else if(*it=="O+")
			{
				resultat[*it]=alliono[1];
			}else if(*it=="H+")
			{
				resultat[*it]=alliono[2];
			}else if(*it=="He+")
			{
				resultat[*it]=alliono[3];
			}else if(*it=="O2+")
			{
				resultat[*it]=alliono[4];
			}else if(*it=="NO+")
			{
				resultat[*it]=alliono[5];
			}else
			{
				Error err("Earth::atmomodel","Your species does not exists in the binary file","The specie "+(*it)+" is not valid for the earth msis model");
				throw err;

			}
		}

	}else
	{

		Error err("Earth::atmomodel","not defined","your type  "+ntostr(vType)+" is not defined");
		throw err;
	}
	return resultat;

}
ublas::vector<double> Earth::ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType)
{

	ublas::vector<double> resu;
	switch(vType)
	{
		case 1:{

			       if(!mbIsModelLoaded)
			       {
				       LoadEarth();
			       }
			       return (mEarth->ReturnIono(vAltGridKm))[0];
			       break;
		       };
		default:{
				Error err("Earth::Electron density","bad type","Your type for the electron density is not well defined");
				throw err;
			};


	};
	return resu;


}
ublas::vector<double> Earth::ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{

	ublas::vector<double> resu;
	switch(vType)
	{
		case 1:{
			       if(!mbIsModelLoaded)
			       {
				       LoadEarth();
			       }
			       return (mEarth->ReturnETemp(vAltGridKm));
			       break;
		       };
		default:{
				Error err("Earth::Electron density","bad type","Your type for the electron density is not well defined");
				throw err;
			};


	};
	return resu;
}
ublas::vector<double> Earth::IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	ublas::vector<double> resu;
	switch(vType)
	{
		case 1:{

			       if(!mbIsModelLoaded)
			       {
				       LoadEarth();
			       }

			       return (mEarth->ReturnITemp(vAltGridKm));
			       break;
		       };
		default:{
				Error err("Earth::Electron density","bad type","Your type for the electron density is not well defined");
				throw err;
			};
	};

	return resu;
}
		


void Earth::LoadEarth()
{

	int modatmo=0;
	double f107=0.;
	double f107bar=0.;
	int year=0;
	int doy=0;
	double ut=0.;
	double tinf=-500;
	double ap=0.;


	mpParameter->ExistsOrDie("/aero_main/planet/model_atmo","You should define the model_atmo value to work with msis");
	mpParameter->ExistsOrDie("/aero_main/planet/year","You should define the year to work with the earth");
	mpParameter->ExistsOrDie("/aero_main/planet/day_of_year","You should define the day of year to work with the earth");
	mpParameter->ExistsOrDie("/aero_main/planet/Ap","You should define the Ap to work with the earth");
	if(mpParameter->Exists("/aero_main/planet/Tinf"))
	{
		mpParameter->GetValue("/aero_main/planet/Tinf",tinf);
	}
	mpParameter->GetValue("/aero_main/planet/model_atmo",modatmo);
	mpParameter->GetValue("/aero_main/planet/year",year);
	mpParameter->GetValue("/aero_main/planet/day_of_year",doy);

	// Ut is required only if we work with planet longitude and latitude
	if(mIsLoadCoordError)
	{
		mpParameter->ExistsOrDie("/aero_main/planet/UT","You should define the universal time to work with the earth");
		mpParameter->GetValue("/aero_main/planet/UT",ut);
	}
	mpParameter->GetValue("/aero_main/planet/Ap",ap);


	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You have to put the f107 parameter to use the msis model");
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107av","You have to put the f107av (average of f107 over three monthes) parameter to use the msis model");

	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	mpParameter->GetValue("/aero_main/sun/model/f107av",f107bar);

	mEarth=new EarthIriMsis(modatmo,f107,f107bar,year,doy,ut,tinf,ap);
	if(!mIsLoadCoordError)
	{
		mEarth->InitSubsolarCoords(mLatDegree,mLoDegree);
	}else
	{
		mpParameter->ExistsOrDie("/aero_main/planet/planet_lat","You should define lat (subsolar), or planet_lat (local) to use earth local latitude");
		mpParameter->ExistsOrDie("/aero_main/planet/planet_long","You should define long (subsolar), or planet_long (local) to use earth local longitude");
		double planet_lat=0;
		double planet_long=0;
		mpParameter->GetValue("/aero_main/planet/planet_lat",planet_lat);
		mpParameter->GetValue("/aero_main/planet/planet_long",planet_long);
		planet_long+=15.*ut;// We face the sun at 12
		mEarth->InitLocalCoords(planet_lat,planet_long);
		mLatDegree=mEarth->ReturnLatDeg();
		mLoDegree=mEarth->ReturnLongDeg();
		mSZADegree=mEarth->ReturnSZADeg();

	}

//	Log::SetPriority(Log::CONFIG);
	Log::mL<<"SZA : "<<mSZADegree<<endl;

	mIsLoadCoordError=false;

	mbIsModelLoaded=true;
//	bool mbIsModelLoaded=true;
}

ublas::vector<double> Earth::ReturndB_B(const ublas::vector<double>& vAltGridKm)
{
	ublas::vector<double> db_b;
	//int resu;
	//resu = 
	int type = ReaddB_B(vAltGridKm, db_b);
	if(1 == type)
	{
		Log::mI<<"We compute the db_b parameter"<<endl;
		mpParameter->ExistsOrDie("/aero_main/planet/dB_B/alt_init","You have to define the initial altitude for dB_B");
		mpParameter->ExistsOrDie("/aero_main/planet/dB_B/dip_angle","You have to define the dip angle for dB_B");
		double alti=0;
		double dipa=0;
		mpParameter->GetValue("/aero_main/planet/dB_B/dip_angle",dipa);
		mpParameter->GetValue("/aero_main/planet/dB_B/alt_init",alti);
		double sd = sin(dipa * PI / 180.);
		for(size_t i = 0; i < db_b.size(); ++i)
		{
			db_b[i] = -1.5 * sd / (6378. + alti + vAltGridKm[i]);
		}
	}
	return db_b;
}
