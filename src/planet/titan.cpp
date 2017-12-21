#include "titan.hpp"
using namespace std;

Titan::Titan(XmlParameters* pParam):Planete(pParam)
{
	Log::mL<<"Bonjour Titan"<<endl;
	mName="Titan";
	mGms_2=1.35;// g in m/s2
	mRKm=2576;
	mUA=9.55;

	if(!LoadCoords())
	{
		cout<<"Error in the coordinates you entered"<<endl;
		Error err("Titan::Titan","Coordinates with errors","The coordinates are not ok");
		throw err;
	}

}

Titan::~Titan()
{
	Log::mL<<"Bye Bye Titan"<<endl;
}

Titan::Titan(XmlParameters* pParam,double vUA):Planete(pParam,vUA)
{
	Log::mL<<"Bonjour Titan"<<endl;
	mName="Titan";
	mGms_2=1.35;// g in m/s2
	mRKm=2576;

	if(!LoadCoords())
	{
		cout<<"Error in the coordinates you entered"<<endl;
		Error err("Titan::Titan","Coordinates with errors","The coordinates are not ok");
		throw err;
	}

}



std::map< std::string, ublas::vector<double> > Titan::AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType)
{
	map<string, ublas::vector<double> > resultat;

	deque< ublas::vector<double> > atmo;
	switch(vType)
	{
		case 1:	
			{
				Log::mI<<"Wodarg 2000 neutral atmosphere"<<endl;
				atmo=TitanData::WodargNeutralAtmo(vAltitudeGridKm);

				deque<string>::iterator it;
				for(it=vSpNames.begin();it!=vSpNames.end();++it)
				{
					if(*it=="N2")
					{
						resultat[*it]=(atmo.at(0));
					}else if(*it=="CH4")
					{
						resultat[*it]=(atmo.at(1));
					}else if(*it=="H2")
					{ // From VMR data by Michel Dobrijevic
						resultat[*it]=(atmo.at(2));
					}else
					{
					//	Log::SetPriority(Log::ERROR);
						Log::mE<<"Error: specie not found!!!"<<endl;
						Log::mE<<"Unfortunately, your specie "<<*it<<" is not taken into account in the VTS3 model (extended with O2)"<<endl;
						Error err("Venus::Venus"," specie not found"," solution : please choose a valid species for the neutral atmosphere");
						throw err;
					}
				}
				mIsTemperatureDefined=true; // Vts3 defined a temperature model
				mTemperatureModelGridK=TitanData::WodargTn(vAltitudeGridKm);
			}
			break;

		case 2:	
			{
				Log::mL<<"Waite 2004 neutral atmosphere"<<endl;
				atmo=TitanData::WaiteNeutralAtmo(vAltitudeGridKm);

				deque<string>::iterator it;
				for(it=vSpNames.begin();it!=vSpNames.end();++it)
				{	
					if(*it=="N2")
					{
						resultat[*it]=(atmo.at(0));
					}else if(*it=="CH4")
					{
						resultat[*it]=(atmo.at(1));
					}else if(*it=="H2")
					{ // From VMR data by Michel Dobrijevic
						resultat[*it]=(atmo.at(2));
					}else
					{
						//Log::SetPriority(Log::ERROR);
						Log::mE<<"Error: specie not found!!!"<<endl;
						Log::mE<<"Unfortunately, your specie "<<*it<<" is not taken into account in the VTS3 model (extended with O2)"<<endl;
						Error err("Venus::Venus"," specie not found"," solution : please choose a valid species for the neutral atmosphere");
						throw err;
					}
				}
				mIsTemperatureDefined=true; // Vts3 defined a temperature model
				mTemperatureModelGridK=TitanData::WaiteTn(vAltitudeGridKm);
			}
			break;
		case 3:	{

				Log::mL<<"Cui 2009 neutral atmosphere"<<endl;
				atmo=TitanData::CuiNeutralAtmo(vAltitudeGridKm);

				deque<string>::iterator it;
				for(it=vSpNames.begin();it!=vSpNames.end();++it)
				{	
					if(*it=="N2")
					{
						resultat[*it]=(atmo.at(0));
					}else if(*it=="CH4")
					{
						resultat[*it]=(atmo.at(1));
					}else if(*it=="H2")
					{ // From VMR data by Michel Dobrijevic
						resultat[*it]=(atmo.at(2));
					}else
					{
					//	Log::SetPriority(Log::ERROR);
						Log::mE<<"Error: specie not found!!!"<<endl;
						Log::mE<<"Unfortunately, your specie "<<*it<<" is not taken into account in the VTS3 model (extended with O2)"<<endl;
						Error err("Venus::Venus"," specie not found"," solution : please choose a valid species for the neutral atmosphere");
						throw err;
					}
				}
				mIsTemperatureDefined=true; // Vts3 defined a temperature model
				mTemperatureModelGridK=TitanData::CuiTn(vAltitudeGridKm);
			}
			break;
		default:
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Atmosphere type  not defined"<<endl;
			Error err("Titan::Titan"," Atmosphere type not found"," solution : please choose a valid type for the neutral atmosphere");
			throw err;
	};

	return resultat;
}


ublas::vector<double> Titan::ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	ublas::vector<double> edens;
	switch(vType)
	{
		case 1:
			edens=TitanData::WodargDensE(vAltGridKm);
			break;
		default:
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Electron density type  not defined"<<endl;
			Error err("Titan::Titan"," Electron density type not found"," solution : please choose a valid type for the electron density");
			throw err;

	};
	return edens;
}
ublas::vector<double> Titan::ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	ublas::vector<double> etemp;
	switch(vType)
	{
		case 1:
			etemp=TitanData::WodargTempE(vAltGridKm);
			break;
		default:
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Electron temperature type  not defined"<<endl;
			Error err("Titan::Titan"," Electron temperature type not found"," solution : please choose a valid type for the electron temperature");
			throw err;

	};
	return etemp;
}

ublas::vector<double> Titan::IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	ublas::vector<double> itemp;
	switch(vType)
	{
		case 1:
			itemp=TitanData::WodargTempI(vAltGridKm);
			break;
		default:
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Electron temperature type  not defined"<<endl;
			Error err("Titan::Titan"," Electron temperature type not found"," solution : please choose a valid type for the electron temperature");
			throw err;

	};
	return itemp;
}


ublas::vector<double> Titan::ReturndB_B(const ublas::vector<double>& vAltGridKm)
{
	ublas::vector<double> db_b;
	//int resu;
	//resu =
	ReaddB_B(vAltGridKm, db_b);
	return db_b;
}
