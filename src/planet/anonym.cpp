/**
 * \file anonym.cpp
 * \brief Implements the anonym planet
 * Copyright G Gronoff Nov2009
 * Last Modification $Id: anonym.cpp 1362 2011-11-19 15:29:56Z gronoff $
 */


#include "anonym.hpp"
using namespace std;

AnonymousPlanet::AnonymousPlanet(XmlParameters* pParam):Planete(pParam)
{
	Log::mL<<"Anonymous planet..."<<endl;
	mName="Anonymous";

	mGms_2=0;
	mRKm=0;
	mUA=0;

	mpParameter->ExistsOrDie("/aero_main/planet/go","You must define the g parameter of the planet");
	mpParameter->GetValue("/aero_main/planet/go",mGms_2);
	mpParameter->ExistsOrDie("/aero_main/planet/radius","You must define the radius of the planet");
	mpParameter->GetValue("/aero_main/planet/radius",mRKm);
	mpParameter->ExistsOrDie("/aero_main/planet/UA","You must define the distance of the planet");
	mpParameter->GetValue("/aero_main/planet/UA",mUA);

	if(!LoadCoords())
	{
		cout<<"Error in the coordinates you entered"<<endl;
		Error err("AnonymousPlanet::AnonymousPlanet","Coordinates with errors","The coordinates are not ok");
		throw err;
	}


}

AnonymousPlanet::~AnonymousPlanet()
{
	Log::mL<<"Bye bye anonymous planet"<<endl;
}


AnonymousPlanet::AnonymousPlanet(XmlParameters* pParam,double vUA):Planete(pParam,vUA)
{
	Log::mL<<"Anonymous planet..."<<endl;
	mName="Anonymous";

	mGms_2=0;
	mRKm=0;

	mpParameter->ExistsOrDie("/aero_main/planet/go","You must define the g parameter of the planet");
	mpParameter->GetValue("/aero_main/planet/go",mGms_2);
	mpParameter->ExistsOrDie("/aero_main/planet/radius","You must define the radius of the planet");
	mpParameter->GetValue("/aero_main/planet/radius",mRKm);

	if(!LoadCoords())
	{
		cout<<"Error in the coordinates you entered"<<endl;
		Error err("AnonymousPlanet::AnonymousPlanet","Coordinates with errors","The coordinates are not ok");
		throw err;
	}


}



std::map< std::string, ublas::vector<double> > AnonymousPlanet::AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType)
{

	map<string, ublas::vector<double> > resultat;
#ifdef DEBUG 
	// To avoid the non-used method
	Log::mL<<vAltitudeGridKm.size()<<" "<<vSpNames.size()<<" "<<vType<<endl;
#endif	
	Error err("AnonymousPlanet::AtmoModel","anonymous planet does not have atmosphere model!"," solution : define an atmospheric model inside the xml file");
	throw err;
	return resultat;
}




ublas::vector<double> AnonymousPlanet::ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType)
{

	ublas::vector<double> edens;
#ifdef DEBUG
	// To avoid the non-used method
	Log::mL<<vAltGridKm.size()<<" "<<vType<<endl;
#endif
	Error err("AnonymousPlanet::ElectronDensity","anonymous planet does not have an electron density model!"," solution : define an electron density model inside the xml file");
	throw err;
	return edens;
}
ublas::vector<double> AnonymousPlanet::ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	ublas::vector<double> etemp;

#ifdef DEBUG
	// To avoid the non-used method
	Log::mL<<vAltGridKm.size()<<" "<<vType<<endl;
#endif
	Error err("AnonymousPlanet::ElectronTemperature","anonymous planet does not have an electron temperature model!"," solution : define an electron temperature model inside the xml file");
	throw err;
	return etemp;
}
ublas::vector<double> AnonymousPlanet::IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{

#ifdef DEBUG
	// To avoid the non-used method
	Log::mL<<vAltGridKm.size()<<" "<<vType<<endl;
#endif
	ublas::vector<double> itemp;
	Error err("AnonymousPlanet::IonTemperature","anonymous planet does not have an ion temperature model!"," solution : define an ion temperature model inside the xml file");
	throw err;
	return itemp;
}


ublas::vector<double> AnonymousPlanet::ReturndB_B(const ublas::vector<double>& vAltGridKm)
{
	ublas::vector<double> db_b;
	//int resu;
	//resu = ReaddB_B(vAltGridKm, db_b);
	ReaddB_B(vAltGridKm, db_b);
	return db_b;
}


