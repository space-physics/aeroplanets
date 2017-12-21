/**
 * \file mars.cpp
 * \brief Implements Mars
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: mars.cpp 1362 2011-11-19 15:29:56Z gronoff $
 */



#include "mars.hpp"
using namespace std;







        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo--------------------------+
         |                                         |
         |Maintenant, on va definir la planete Mars|
         |                                         |
         +----------------------------------------*/


Mars::Mars(XmlParameters* pParam) : Planete(pParam)
{
	mName="Mars";
	mUA=1.52; // Mean value
	mGms_2=3.71; // g in ms-2
	mRKm=3396.2; // Radius in km
	mLsDegree=0.;
	mbIsBinaryLoaded=false;
	mbIsMartimLoaded=false;

	if(pParam->Exists("/aero_main/planet/Ls"))
	{// If the Ls is defined, it is taken into account
		pParam->GetValue("/aero_main/planet/Ls",mLsDegree);
		LsInit(mLsDegree);
		//Log::SetPriority(Log::CONFIG);
		Log::mI<<" Ls = "<<mLsDegree<<" allows to compute distance in UA : "<<mUA<<endl;
	}
	if(pParam->Exists("/aero_main/planet/UA"))
	{
	
		pParam->GetValue("/aero_main/planet/UA",mUA);
	}
	if(!LoadCoords())
	{
	//	Log::SetPriority(Log::DEBUGG);
		Log::mI<<"Error in the coordinates you entered"<<endl;
		Log::mI<<"This is not a problem if local coordinates have been also defined"<<endl;
		Log::mI<<"Will raise an error if martim neutral atmosphere is not used"<<endl;
	//	exit(1);
	}


}

Mars::Mars(XmlParameters* pParam,double vUA) : Planete(pParam,vUA)
{
	mName="Mars";
	mGms_2=3.71; // g in ms-2
	mRKm=3396.2; // Radius in km
	mLsDegree=0.;
	mbIsBinaryLoaded=false;
	mbIsMartimLoaded=false;
	if(!LoadCoords())
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mL<<"Error in the coordinates you entered"<<endl;
		Log::mL<<"This is not a problem if local coordinates have been also defined"<<endl;
		Log::mL<<"Will raise an error if martim neutral atmosphere is not used"<<endl;
	//	exit(1);
	}
}


void Mars::LsInit(double vLs)
{
	mUA=1.5237*7.73396*((1.-0.934*0.934))/(1.+0.0934*cos((vLs+109.)*M_PI/180.));



}


void Mars::LoadMarTim()
{
	mpParameter->ExistsOrDie("/aero_main/planet/martim_temp_file","Solution : edit the xml file to define the file, or modify the option not to use it");
	mpParameter->ExistsOrDie("/aero_main/planet/martim_compo_file","Solution : edit the xml file to define the file, or modify the option not to use it");

	string compofile=mpParameter->SubFileName("/aero_main/planet/martim_compo_file");
	string tempfile=mpParameter->SubFileName("/aero_main/planet/martim_temp_file");
	
	mMarTim=new MarsAtmoTim(compofile,tempfile);

	if(!mIsLoadCoordError)
	{
		mMarTim->InitSubsolarCoords(mLsDegree,mLatDegree,mLoDegree);
	}else
	{
		mpParameter->ExistsOrDie("/aero_main/planet/planet_lat","You should define lat (subsolar), or planet_lat (local) to use martim");
		mpParameter->ExistsOrDie("/aero_main/planet/planet_long","You should define long (subsolar), or planet_long (local) to use martim");
		double planet_lat=0;
		double planet_long=0;
		mpParameter->GetValue("/aero_main/planet/planet_lat",planet_lat);
		mpParameter->GetValue("/aero_main/planet/planet_long",planet_long);
		mMarTim->InitLocalCoords(mLsDegree,planet_lat,planet_long);
		mLatDegree=mMarTim->ReturnLatDeg();
		mLoDegree=mMarTim->ReturnLongDeg();
		mSZADegree=mMarTim->ReturnSZADeg();
	}
	mIsLoadCoordError=false;
	mbIsMartimLoaded=true;
}

void Mars::LoadDump()
{
	mpParameter->ExistsOrDie("/aero_main/planet/mars_binary_file","Solution : edit the xml file to define the file, or modify the option not to use it");

	string binfile=mpParameter->GetSubFileName("/aero_main/planet/mars_binary_file");
	if(!mpParameter->Exists("/aero_main/planet/mars_binary_file_byte"))
	{
		mBin=new DumpData(binfile);
	}else
	{
		int byte;
		mpParameter->GetValue("/aero_main/planet/mars_binary_file_byte",byte);
		Log::mL<<"You have a "<<byte<<" byte file!"<<endl;
		mBin=new DumpData(binfile,byte);
	}

	mbIsBinaryLoaded=true;
	mBin->ReadTableau();


}

Mars::~Mars()
{
	cout<<"Bye Bye mars"<<endl;
	if(mbIsMartimLoaded)
		delete mMarTim;
	if(mbIsBinaryLoaded)
		delete mBin;
	Log::mL<<"End Mars destructor"<<endl;
}



std::map< std::string,ublas::vector<double> > Mars::AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType)
{

	map<string, ublas::vector<double> > resultat;

//	Log::SetPriority(Log::CONFIG);
	if(vType==1)
	{
		Log::mI<<"Binary file neutral atmosphere "<<endl;

		if(!mbIsBinaryLoaded)
		{
			LoadDump();
		}

		Log::mD<<"return neutral atmo"<<endl;
		deque< ublas::vector<double> > allatmo=mBin->ReturnNeutralAtmo(vAltitudeGridKm);
		deque< ublas::vector<double> > alliono=mBin->ReturnIono(vAltitudeGridKm);
		mIsTemperatureDefined=true;
		Log::mD<<"temperature"<<endl;
		mTemperatureModelGridK=mBin->ReturnNTemp(vAltitudeGridKm);
		Log::mD<<"Ntemp fini"<<endl;
		deque<string>::iterator it;
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			if(*it=="H")
			{
				resultat[*it]=allatmo[0];
			}else if(*it=="O")
			{
				resultat[*it]=allatmo[1];
			}else if(*it=="N2")
			{
				resultat[*it]=allatmo[2];
			}else if(*it=="O2")
			{
				resultat[*it]=allatmo[3];
			}else if(*it=="CO2")
			{
				resultat[*it]=allatmo[4];
			}else if(*it=="H+")
			{
				resultat[*it]=alliono[0];
			}else if(*it=="O+")
			{
				resultat[*it]=alliono[1];
			}else if(*it=="O2+")
			{
				resultat[*it]=alliono[2];
			}else if(*it=="CO2+")
			{
				resultat[*it]=alliono[3];
			}else
			{
				Error err("Mars::atmomodel","Your species does not exists in the binary file","The specie "+(*it)+" is not valid for the martian binary file");
				throw err;
			}
		}
		Log::mD<<"Fin mise en place atmo"<<endl;

	}
	if(vType==2)
	{
		Log::mI<<"Martim neutral atmosphere "<<endl;

		if(!mbIsMartimLoaded)
		{
			LoadMarTim();
		}
		mIsTemperatureDefined=true;
		mTemperatureModelGridK=mMarTim->ReturnNTemp(vAltitudeGridKm);
		deque< ublas::vector<double> > allatmo=mMarTim->ReturnNeutralAtmo(vAltitudeGridKm);
		deque< ublas::vector<double> > alliono=mMarTim->ReturnIono(vAltitudeGridKm);
		deque<string>::iterator it;
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			if(*it=="N2")
			{
				resultat[*it]=allatmo[0];
			}else if(*it=="O")
			{
				resultat[*it]=allatmo[1];
			}else if(*it=="CO2")
			{
				resultat[*it]=allatmo[2];
			}else if(*it=="Ar")
			{
				resultat[*it]=allatmo[3];
			}else if(*it=="CO")
			{
				resultat[*it]=allatmo[4];
			}else if(*it=="O2")
			{
				resultat[*it]=allatmo[5];
			}else if(*it=="NO")
			{
				resultat[*it]=allatmo[6];
			}else if(*it=="CO2+")
			{
				resultat[*it]=alliono[0];
			}else if(*it=="N2+")
			{
				resultat[*it]=alliono[1];
			}else if(*it=="O+")
			{
				resultat[*it]=alliono[2];
			}else if(*it=="O2+")
			{
				resultat[*it]=alliono[3];
			}else if(*it=="CO+")
			{
				resultat[*it]=alliono[4];
			}else
			{
				Error err("Mars::atmomodel","Your species does not exists in Martim atmosphere","The specie "+(*it)+" is not valid for the martian binary file");
				throw err;
			}
		}

	}
	if(vType==3 || vType==4)
	{


		deque< ublas::vector<double> > allatmo;
		mIsTemperatureDefined=true;
		if(vType==4)
		{
			Log::mI<<"Mariner neutral atmosphere "<<endl;
			allatmo=MarsData::MarinerNeutralAtmo(vAltitudeGridKm);
			mTemperatureModelGridK=MarsData::MarinerNTemp(vAltitudeGridKm);
		}else
		{
			Log::mI<<"Viking neutral atmosphere "<<endl;
			allatmo=MarsData::VikingNeutralAtmo(vAltitudeGridKm);
			mTemperatureModelGridK=MarsData::VikingNTemp(vAltitudeGridKm);
		}	


		deque<string>::iterator it;
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			if(*it=="CO2")
			{
				resultat[*it]=allatmo[0];
			}else if(*it=="N2")
			{
				resultat[*it]=allatmo[1];
			}else if(*it=="CO")
			{
				resultat[*it]=allatmo[2];
			}else if(*it=="O")
			{
				resultat[*it]=allatmo[3];
			}else if(*it=="O2")
			{
				resultat[*it]=allatmo[4];
			}else if(*it=="H")
			{
				resultat[*it]=allatmo[5];
			}else
			{
				Error err("Mars::atmomodel","Your species does not exists in the Mariner/Viking neutral atmosphere","The specie "+(*it)+" is not valid for the martian binary file");
				throw err;
			}
		}

	}

	if(vType<0||vType>4)
	{
		Error err("Mars::atmomodel","not defined","your type  "+ntostr(vType)+" is not defined");
		throw err;
	}


	if(mIsLoadCoordError)
	{
		Error err("Mars::LoadCoord","Error in the coordinates you entered","Please provide real coordinates, or use the Martim atmosphere");
		throw err;
	}
	return resultat;
}



ublas::vector<double> Mars::ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	Log::mD<<"Electron density"<<endl;
	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();



	if(vType==1)
	{
		if(!mbIsBinaryLoaded)
		{
			LoadDump();
		}
		deque< ublas::vector<double> > iono=mBin->ReturnIono(vAltGridKm);
		for(unsigned sp=0;sp<iono.size();++sp)
		{
			for(unsigned i=0;i<vAltGridKm.size();++i)
			{
				machin[i]+=iono[sp][i];
			}
		}
	}


	if(vType==2)
	{
		if(!mbIsMartimLoaded)
		{
			LoadMarTim();
		}
		deque< ublas::vector<double> > iono=mMarTim->ReturnIono(vAltGridKm);

		machin=iono[MARTIM_POSITION_EDENS];
	}

	if(vType<0||vType>2)
	{
		Error err("Mars::electron density","not defined","your type  "+ntostr(vType)+" is not defined");
		throw err;
	}
	return machin;
}

ublas::vector<double> Mars::ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	Log::mD<<"Electron temperature"<<endl;
	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();

	if(vType==1)
	{
		if(!mbIsBinaryLoaded)
		{
			LoadDump();
		}
		machin=mBin->ReturnETemp(vAltGridKm);
	}

	if(vType==2)
	{
		if(!mbIsMartimLoaded)
		{
			LoadMarTim();
		}
	//	Log::SetPriority(Log::WARNING,"Martim::ElectronTemperature");
		Log::mW<<"The electron temp of martim is not defined : we use the neutral temperature instead"<<endl;
		machin=mMarTim->ReturnNTemp(vAltGridKm);

	}

	if(vType<0||vType>2)
	{
		Error err("Mars::electron temperature","not defined","your type  "+ntostr(vType)+" is not defined");
		throw err;
	}

	return machin;
}

ublas::vector<double> Mars::IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	Log::mL<<"Ion temperature"<<endl;
	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();

	if(vType==1)
	{
		if(!mbIsBinaryLoaded)
		{
			LoadDump();
		}
		machin=mBin->ReturnITemp(vAltGridKm);
	}

	if(vType==2)
	{
		if(!mbIsMartimLoaded)
		{
			LoadMarTim();
		}
	//	Log::SetPriority(Log::WARNING,"Martim::IonTemperature");
		Log::mW<<"The ion temp of martim is not defined : we use the neutral temperature instead"<<endl;
		machin=mMarTim->ReturnNTemp(vAltGridKm);

	}

	if(vType<0||vType>2)
	{
		Error err("Mars::Iontemperature","not defined","your type  "+ntostr(vType)+" is not defined");
		throw err;
	}

	return machin;
}


ublas::vector<double> Mars::ReturndB_B(const ublas::vector<double>& vAltGridKm)
{
	ublas::vector<double> db_b;
	//int resu;
	//resu =
	ReaddB_B(vAltGridKm, db_b);
	return db_b;
}
