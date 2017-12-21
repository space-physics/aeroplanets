/**
 * \file venus.cpp
 * \brief Implements the Venus planet
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: venus.cpp 1362 2011-11-19 15:29:56Z gronoff $
 */

#include "venus.hpp"
using namespace std;
using namespace PhysTime;
        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo---------------------------+
         |                                          |
         |Maintenant, on va definir la planete Venus|
         |                                          |
         +-----------------------------------------*/

Venus::Venus(XmlParameters* pParam):Planete(pParam)
{
	Log::mL<<"Bonjour Venus"<<endl;
	mName="Venus";
	mGms_2=8.87; // g in m/s2
	mRKm=6051.8;
	mUA=0.723;
	if(pParam->Exists("/aero_main/planet/UA"))
	{
		pParam->GetValue("/aero_main/planet/UA",mUA);
	}
	if(!LoadCoords())
	{
		Log::mE<<"Error in the coordinates you entered"<<endl;
		Error err("Venus::Venus","Coordinates with errors","The coordinates are not ok");
		throw err;
		
	}
}

Venus::~Venus()
{
	Log::mL<<"Bye bye Venus"<<endl;
}


Venus::Venus(XmlParameters* pParam,double vUA):Planete(pParam,vUA)
{
	mName="Venus";
	mGms_2=8.87; //g in m/s2
	mRKm=6051.8;
	if(!LoadCoords())
	{
		Log::mE<<"Error in the coordinates you entered"<<endl;
		Error err("Venus::Venus","Coordinates with errors","The coordinates are not ok");
		throw err;
	}
}



std::map< std::string,ublas::vector<double> > Venus::AtmoModel(const ublas::vector<double>& vAltitudeGridKm, std::deque< std::string >  vSpNames,int vType)
{

	deque< ublas::vector<double> > densites;



	float f107,f107av;

	map<string, ublas::vector<double> > resultat;
	if(vType==1)
	{

		mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You have to put the f107 parameter to use the VTS3 - extended to O2 - model");
		mpParameter->ExistsOrDie("/aero_main/sun/model/f107av","You have to put the f107av (average of f107 over three monthes) parameter to use the VTS3 - extended to O2 - model");

		mpParameter->GetValue("/aero_main/sun/model/f107",f107);
		mpParameter->GetValue("/aero_main/sun/model/f107av",f107av);
		mIsTemperatureDefined=true; // Vts3 defined a temperature model

		//mTemperatureModelGridK.erase(mTemperatureModelGridK.begin(),mTemperatureModelGridK.end());
		mTemperatureModelGridK.clear();
		mTemperatureModelGridK.resize(vAltitudeGridKm.size());
		for(unsigned i=0;i<vAltitudeGridKm.size();++i)
		{
			vector<float> densi(7);
			vector<float> temp(2);

			float alti=(float)vAltitudeGridKm[i];
			float xlat=(float)mLatDegree;
			float xloc=(float)mHrLoc;
			int mas=48;
			float yrd=0;

			vts3_(&alti,&yrd,&xlat,&xloc,&f107av,&f107,&mas,&densi[0],&temp[0]);


			//		for(unsigned j=0;j<7;++j)
			//			cout<<"Test VTS3 alti, densi[1] : "<<alti<<" - "<<densi[j]<<" -- "<<temp[0]<<" T "<<temp[1]<<" "<<xlat<<" "<<xloc<<" "<<f107av<<" "<<f107<<" "<<mas<<" ---"<<j<<endl;
			ublas::vector<double> dens(densi.size());
			copy(densi.begin(),densi.end(),dens.begin());
			densites.push_back(dens);
			mTemperatureModelGridK[i]= static_cast<double>(temp[1]);

		}

		deque<string>::iterator it;
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			// Just for the test
			//	vector<double> machin;
			if(*it=="CO2")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][1]);
				}
			}else if(*it=="O")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][2]);
				}
			}else if(*it=="CO")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][3]);
				}
			}else if(*it== "He")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][4]);
				}
			}else if(*it== "N")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][5]);
				}
			}else if(*it=="N2")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					resultat[*it][i]=(densites[i][6]);
				}
			}else if(*it=="O2")
			{
				resultat[*it].resize(densites.size());
				for(unsigned i=0;i<densites.size();++i)
				{
					// To compute the density of O2, we assume that it is the density of CO2 divided by 1E3
					//resultat[*it][i]=densites[i][1]/1E3;
					resultat[*it][i]=(densites[i][1]/1E3);
				}
			}else
			{
			//	Log::SetPriority(Log::ERROR);
				Log::mE<<"Error: specie not found!!!"<<endl;
				Log::mE<<"Unfortunately, your specie "<<*it<<" is not taken into account in the VTS3 model (extended with O2)"<<endl;
				Error err("Venus::Venus"," specie not found"," solution : please choose a valid species for the neutral atmosphere");
					throw err;
			}
		}
	}else if(vType == 2)
	{
		mIsTemperatureDefined=true; // Vts3 defined a temperature model

		//mTemperatureModelGridK.erase(mTemperatureModelGridK.begin(),mTemperatureModelGridK.end());
		//mTemperatureModelGridK.clear();
		ublas::matrix<double> atmocm_3;
		LaunchVenGram( vAltitudeGridKm, atmocm_3 , mTemperatureModelGridK);

		deque<string>::iterator it;
		size_t nalt = vAltitudeGridKm.size();
		for(it=vSpNames.begin();it!=vSpNames.end();++it)
		{
			// Just for the test
			//	vector<double> machin;
			if(*it=="CO2")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,0);
				}
			}else if(*it=="N2")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,1);
				}
			}else if(*it=="CO")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,2);
				}
			}else if(*it=="O")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,3);
				}
			}else if(*it=="O2")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,4);
				}
			}else if(*it=="He")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,5);
				}
			}else if(*it=="N")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,6);
				}
			}else if(*it=="H")
			{
				resultat[*it].resize(nalt);
				for(size_t i=0;i<nalt;++i)
				{
					resultat[*it][i]=atmocm_3(i,7);
				}
			}else
			{
				Log::mE<<"Error: specie not found!!!"<<endl;
				Log::mE<<"Unfortunately, your specie "<<*it<<" is not taken into account in the VTS3 model (extended with O2)"<<endl;
				Error err("Venus::Venus"," specie not found"," solution : please choose a valid species for the neutral atmosphere");
					throw err;
			}
		}

	}else
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Atmosphere type  not defined"<<endl;
		Error err("Venus::Venus"," Atmosphere type not found"," solution : please choose a valid type for the neutral atmosphere");
		throw err;
	}

	return resultat;
}


ublas::vector<double> Venus::ElectronDensity(const ublas::vector<double> &vAltGridKm,const int& vType)
{

	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();
	
	if(vType==1)
	{
		
		return TheisModelNecm_3(vAltGridKm,mSZADegree);
	


	}else
	{

		Error err("Venus::ElectronDensity","Model not defined","The model n "+ntostr(vType)+" is not defined for the electron density");
		//cout<<"Electron Density  Model for Venus NOT DEFINED"<<endl;
		//cout<<"Especially for the type"<<vType<<endl;
		//exit(1);
		throw err;
	}
	return machin;
}

ublas::vector<double> Venus::ElectronTemperature(const ublas::vector<double> & vAltGridKm,const int& vType)
{
	if(vType==1)
	{
		return TheisModelTeK(vAltGridKm,mSZADegree);

	}else
	{
		Error err("Venus::ElectronTemperature","Model not defined","The model n "+ntostr(vType)+" is not defined for the electron temperature");
		throw err;
	}

	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();
//	acout<<"Electron temperature  Model for Venus NOT DEFINED"<<endl;
//	cout<<"Especially for the type"<<vType<<endl;
//	exit(1);
	return machin;
}

ublas::vector<double> Venus::IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)
{
	if(vType==1)
	{

		double f107;
		mpParameter->GetValue("/aero_main/sun/model/f107",f107);
		return IonKTempScannedModel(vAltGridKm,f107);
	}else
	{
		Error err("Venus::IonTemperature","Model not defined","The model n "+ntostr(vType)+" is not defined for the ion temperature");
		throw err;

	}

	ublas::vector<double> machin(vAltGridKm.size());
	machin.clear();
	return machin;
}

double Venus::FsModNecm_3(double vAltKm,double vSZADegree)
{// GG translated this from fortran...

	// The parameter array, 
	double p[26]={5.334022,252.5587,120.4444,9.1163822E-02,-6.195264,37.38898,3.331818,245.5529,0.7587564,0.1453792,-0.9461438,127.1485,-0.3391595,2.4936652E-02,0.1879423,-5.3126566E-02,2.9439656E-02,1.654282,415.7859,2.6281483E-02,145.3951,2.045347,0.4421963,3.036170,9.1009791E-04,141.5824	};

	double sza=vSZADegree*PI/180;
	double s2=sin(sza)*sin(sza);
	double szap=sza*(1-p[22]*s2*exp(-p[23]*s2-p[24]*(vAltKm-p[25])*(vAltKm-p[25])));
	double sinz=sin(p[17]+szap);
	double s=p[0]+p[1]/(vAltKm-p[2])+sinz*p[12]+p[18]*exp(-p[19]*((vAltKm-p[20])*(vAltKm-p[20]))-p[21]*szap*szap);
	double o=p[3]+p[4]/(vAltKm-p[5])+sinz*p[13];
	double c=p[6]+p[7]/(vAltKm-p[8])+sinz*p[14];
	double m=p[9]+p[10]/(vAltKm-p[11])+sinz*p[15];
	double ss=s*cos(p[16]+szap);
	double calc=(1.77245*ss*erfc(-ss)+exp(-ss*ss)+o);
	if(calc<=0)
		calc=0.00001;
	return pow(10,c+m*log(calc));
}


ublas::vector<double> Venus::TheisModelNecm_3(ublas::vector<double>  vAltGridKm,double vSZADegree)
{
	ublas::vector<double> resu(vAltGridKm.size());
	resu.clear();
	double corrfactor=-1;
	unsigned size=vAltGridKm.size();
	double invert=false;
	if(vAltGridKm[0]<vAltGridKm[size-1])
	{
		invert=true;
		std::reverse(vAltGridKm.begin(),vAltGridKm.end());
	}

	double f107;
	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	ublas::vector<double> autre_dene=NeDensityScannedModel(vAltGridKm,f107);

	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		if(vAltGridKm[i]>149)
		{
			resu[i]=(FsModNecm_3(vAltGridKm[i],vSZADegree));
		}else
		{
			if(corrfactor<0)
			{
				corrfactor=resu[i-1]/autre_dene[i];
			}
			resu[i]=(autre_dene[i]*corrfactor);
		}

	}
	if(invert)
	{
		std::reverse(resu.begin(),resu.end());
	}
	return resu;
}


ublas::vector<double> Venus::NeDensityScannedModel(ublas::vector<double> vAltGridKm,double vF107)
{

	unsigned model_nbpt=100;
	// Altitude in km for our model
	double alt_km_model[]={90.0299988,90.1200027,90.2799988,90.5,90.7799988,
		91.1200027,91.5199966,91.9800034,92.5100021,93.0999985,93.75,
		94.4599991,95.2399979,96.0800018,96.9800034,97.9400024,98.9599991,
		100.040001,101.190002,102.400002,103.669998,105.,106.400002,
		107.860001,109.379997,110.959999,112.599998,114.300003,116.07,
		117.900002,119.790001,121.739998,123.760002,125.839996,127.970001,
		130.179993,132.440002,134.759995,137.149994,139.600006,142.110001,
		144.679993,147.320007,150.020004,152.770004,155.600006,158.479996,
		161.419998,164.429993,167.5,170.630005,173.820007,177.080002,
		180.399994,183.770004,187.220001,190.720001,194.279999,197.910004,
		201.600006,205.350006,209.160004,213.039993,216.979996,220.979996,
		225.039993,229.160004,233.339996,237.589996,241.899994,246.270004,
		250.699997,255.199997,259.76001,264.380005,269.059998,273.799988,
		278.600006,283.470001,288.399994,293.390015,298.440002,303.559998,
		308.73999,313.980011,319.279999,324.640015,330.059998,335.549988,
		341.100006,346.709991,352.380005,358.119995,363.920013,369.769989,
		375.700012,381.679993,387.720001,393.829987,400.};
	//
	//       Valeurs à f107=60 de la densité électronique
	//	Electron density at f107=60
	//
	double dene_cm_3_f60[]={239.689682,276.354126,341.530884,431.151642,545.213318,
		683.718933,846.662415,1034.05298,1249.95581,1485.04285,1610.35559,
		1747.23511,1897.60974,2002.8877,2019.71985,2037.67419,2056.75049,
		2058.01465,2031.20898,2004.44812,2005.91455,2007.4502,2009.06665,
		2045.64075,2169.96606,2480.18384,2980.75977,3854.54956,4466.104,
		5704.2002,11532.8789,22192.293,39516.9961,67227.9531,80571.7109,
		93437.5859,130449.094,189950.062,259237.047,295263.844,294497.688,
		274202.625,246228.359,219949.078,197509.094,171762.016,156127.266,
		137219.625,118699.352,99213.7812,86476.0938,74495.1406,64475.8945,
		60401.5469,56394.4219,52723.293,49567.3633,48401.6836,45432.4453,
		43052.6445,40888.7422,38608.5781,36749.4922,34867.0664,33338.1836,
		31986.7246,30137.9395,28535.3652,27778.75,25966.8652,25687.3672,
		24874.8418,23837.7012,23175.5137,21569.0664,20470.8926,19637.3496,
		19050.8301,18482.5859,17577.9492,16847.2773,16432.4707,15846.9766,
		15195.1924,14836.2598,14008.8281,13696.3223,13286.6357,12871.6572,
		12452.1416,12028.0928,11599.5068,11165.6318,10727.2197,10285.0303,
		9836.79102,9384.77539,8928.22266,8466.37988,8000};

	//
	//        Valeurs à f107=200 de la densité électronique
	//	Electron density at f107=200
	//
	double dene_cm_3_f200[]={4842.4873,4876.02734,4935.65039,5017.63428,5121.97656,
		5248.67969,5397.73877,5569.16113,5742.22803,5901.1001,6075.12988,
		6211.65869,6361.64795,6497.10059,6500.23193,6503.57227,6431.69775,
		6332.81299,6229.77393,6122.44482,6252.12598,6445.55615,6655.40479,
		6927.86084,7355.23779,9065.67285,11328.1621,12776.9756,18822.1348,
		30715.1699,55772.1055,102649.094,133499.938,152851.922,172669.219,
		229127.078,294381.312,364434.812,422354.188,459075.188,434730.781,
		398497.125,358729.5,336370.625,309828.094,277876.094,256043.,
		246721.922,230178.906,209651.562,192373.578,177213.391,163630.141,
		163921.141,164216.594,164519.75,164828.375,165142.844,165463.578,
		165274.547,162185.312,157946.078,149795.703,143109.422,141926.203,
		137225.703,130529.289,126209.805,121186.797,116892.195,114062.398,
		109193.188,103861.984,100626.75,97031.3125,94119.0547,89461.2734,
		85806.0156,83119.8828,80010.3359,77945.7188,74716.5156,73011.2578,
		69454.1562,67309.9297,65299.5586,63156.9727,61173.2305,59633.6016,
		57316.9258,56141.6016,55312.0664,54643.1094,53728.,53191.4258,
		51266.9375,50096.7109,48173.8906,46627.5312,45800.0469};

	ublas::vector<double> tmp_dene(model_nbpt);
	for(unsigned i=0;i<model_nbpt;++i)
	{
		tmp_dene[i]=( (dene_cm_3_f200[i]-dene_cm_3_f60[i])*(vF107-60)/140.+dene_cm_3_f60[i]);
	}

	ublas::vector<double> tmp_gridalt(model_nbpt);
	std::copy(alt_km_model,alt_km_model+model_nbpt,tmp_gridalt.begin());

	return MathFunction::IntLin(tmp_gridalt,tmp_dene,vAltGridKm);
}




ublas::vector<double> Venus::IonKTempScannedModel(ublas::vector<double> vAltGridKm,double vF107)
{
	unsigned model_nbpt=100;
	double alt_km_model[]={90.0299988,90.1200027,90.2799988,90.5,90.7799988,
		91.1200027,91.5199966,91.9800034,92.5100021,93.0999985,93.75,
		94.4599991,95.2399979,96.0800018,96.9800034,97.9400024,98.9599991,
		100.040001,101.190002,102.400002,103.669998,105.,106.400002,
		107.860001,109.379997,110.959999,112.599998,114.300003,116.07,
		117.900002,119.790001,121.739998,123.760002,125.839996,127.970001,
		130.179993,132.440002,134.759995,137.149994,139.600006,142.110001,
		144.679993,147.320007,150.020004,152.770004,155.600006,158.479996,
		161.419998,164.429993,167.5,170.630005,173.820007,177.080002,
		180.399994,183.770004,187.220001,190.720001,194.279999,197.910004,
		201.600006,205.350006,209.160004,213.039993,216.979996,220.979996,
		225.039993,229.160004,233.339996,237.589996,241.899994,246.270004,
		250.699997,255.199997,259.76001,264.380005,269.059998,273.799988,
		278.600006,283.470001,288.399994,293.390015,298.440002,303.559998,
		308.73999,313.980011,319.279999,324.640015,330.059998,335.549988,
		341.100006,346.709991,352.380005,358.119995,363.920013,369.769989,
		375.700012,381.679993,387.720001,393.829987,400.};


	//
	//       Valeurs à f107=200 de la température ionique
	//
	double tempion_K_f200[]={182.106506,181.997391,181.805573,181.546326,
		181.223938,180.843826,180.412628,179.938049,179.419586,178.878067,
		178.324951,177.772873,177.22905,176.716797,176.252441,175.853424,
		175.538437,175.32724,175.240692,175.303665,175.539658,175.973343,
		177.005295,177.316116,177.270935,177.289536,177.110794,177.485519,
		178.722443,179.289429,179.042023,179.43277,180.443573,181.44873,
		182.216202,185.874451,187.115372,192.196655,201.638031,218.049484,
		241.726151,260.743927,279.599548,297.560394,318.388489,344.022552,
		377.39444,406.576447,435.600433,472.859924,503.11322,526.084595,
		547.882263,565.959595,580.194763,589.640503,594.209778,599.620605,
		606.877136,616.478271,626.989136,637.964172,652.084839,669.211975,
		690.725098,713.542114,737.134033,764.707825,795.465759,829.413147,
		861.844421,902.776062,947.751282,996.602478,1057.05969,1111.10547,
		1177.901,1240.34741,1317.65161,1396.10889,1471.76135,1552.39954,
		1623.93726,1692.22717,1763.4729,1851.45435,1923.80298,1980.58374,
		2050.46509,2104.69824,2155.93286,2217.69971,2260.6543,2299.1936,
		2336.91187,2352.87769,2394.93164,2405.04004,2429.14893,2456.85278};

	//
	//       Valeurs à f107=60 de la température ionique
	//
	double tempion_K_f60[]={182.106506,181.997391,181.805573,181.546326,
		181.223938,180.843826,180.412628,179.938049,179.419586,178.878067,
		178.324951,177.772873,177.22905,176.716797,176.252441,175.853424,
		175.538437,175.32724,175.240692,175.303665,175.539658,175.973343,
		177.005295,177.316116,177.270935,177.289536,177.110794,177.485519,
		178.722443,179.289429,179.042023,179.43277,180.443573,181.44873,
		182.56636,183.552017,186.433945,189.494186,192.388931,204.819168,
		219.199951,236.670425,250.339813,262.923859,281.547974,307.437012,
		330.276337,363.756317,391.89389,416.240601,447.924957,475.590973,
		495.797485,511.716095,524.326538,533.465881,540.090881,544.105469,
		551.546753,560.627991,570.721008,578.659485,597.257385,614.083679,
		632.228088,657.211975,692.260315,709.510498,736.523804,767.636475,
		800.514832,841.900269,889.486511,930.218811,986.945557,1043.67834,
		1115.18091,1193.29138,1256.77197,1349.54285,1405.18811,1484.16895,
		1596.71851,1689.04932,1757.09473,1837.99561,1898.59058,1975.09253,
		2034.73523,2089.42261,2148.11938,2202.53564,2238.29248,2280.54224,
		2329.13892,2361.54077,2384.47021,2424.51538,2427.33838,2453.3811};

	ublas::vector<double> tmp_iont(model_nbpt);
	for(unsigned i=0;i<model_nbpt;++i)
	{
		tmp_iont[i]=( (tempion_K_f200[i]-tempion_K_f60[i])*(vF107-60)/140.+tempion_K_f60[i]);
	}

	ublas::vector<double> tmp_gridalt(model_nbpt);
	std::copy(alt_km_model,alt_km_model+model_nbpt,tmp_gridalt.begin());

	return MathFunction::IntLin(tmp_gridalt,tmp_iont,vAltGridKm);


}

ublas::vector<double> Venus::ElectronKTempScannedModel(ublas::vector<double> vAltGridKm)
{
	unsigned model_nbpt=100;
	double alt_km_model[]={90.0299988,90.1200027,90.2799988,90.5,90.7799988,
		91.1200027,91.5199966,91.9800034,92.5100021,93.0999985,93.75,
		94.4599991,95.2399979,96.0800018,96.9800034,97.9400024,98.9599991,
		100.040001,101.190002,102.400002,103.669998,105.,106.400002,
		107.860001,109.379997,110.959999,112.599998,114.300003,116.07,
		117.900002,119.790001,121.739998,123.760002,125.839996,127.970001,
		130.179993,132.440002,134.759995,137.149994,139.600006,142.110001,
		144.679993,147.320007,150.020004,152.770004,155.600006,158.479996,
		161.419998,164.429993,167.5,170.630005,173.820007,177.080002,
		180.399994,183.770004,187.220001,190.720001,194.279999,197.910004,
		201.600006,205.350006,209.160004,213.039993,216.979996,220.979996,
		225.039993,229.160004,233.339996,237.589996,241.899994,246.270004,
		250.699997,255.199997,259.76001,264.380005,269.059998,273.799988,
		278.600006,283.470001,288.399994,293.390015,298.440002,303.559998,
		308.73999,313.980011,319.279999,324.640015,330.059998,335.549988,
		341.100006,346.709991,352.380005,358.119995,363.920013,369.769989,
		375.700012,381.679993,387.720001,393.829987,400.};

	//
	//        Valeurs à  de la température électronique (indépendante de f107)
	//
	double    tempelec_K[]={175.274826,175.275879,175.277756,175.280334,
		175.283615,175.287613,175.292297,175.297699,175.303909,175.310822,
		175.318451,175.326782,175.335938,175.345779,175.356339,175.367615,
		175.379578,175.392258,175.405762,175.419968,175.434875,175.4505,
		175.533133,175.544159,175.666431,177.7508,177.836481,177.960022,
		179.380737,180.476334,180.700897,187.138885,198.658142,208.32019,
		219.43988,232.678757,246.937195,262.876221,284.674438,301.03775,
		320.66391,344.493011,376.844208,421.958923,491.436554,600.764343,
		702.263367,867.160278,1075.64966,1262.33008,1522.93359,1803.76794,
		2003.547,2137.32275,2197.5127,2255.44922,2323.49854,2401.46948,
		2475.95801,2548.43506,2610.60083,2675.51685,2744.01929,2810.14331,
		2874.75073,2933.01514,2989.53369,3042.89429,3109.9729,3166.65747,
		3204.31909,3246.17676,3289.5415,3339.35864,3386.65771,3429.49512,
		3449.18677,3498.88794,3540.25073,3563.63745,3590.35107,3617.58887,
		3646.67725,3706.52588,3740.08032,3769.86499,3800.23047,3831.17578,
		3862.78931,3895.00806,3928.44434,3962.92847,3960.82007,3986.23657,
		4052.17334,4048.96021,4084.40039,4127.27344,4134.08008,4135.77637};


	ublas::vector<double> tmp_et(model_nbpt);
	std::copy(tempelec_K,tempelec_K+model_nbpt,tmp_et.begin());

	ublas::vector<double> tmp_gridalt(model_nbpt);
	std::copy(alt_km_model,alt_km_model+model_nbpt,tmp_gridalt.begin());

	return MathFunction::IntLin(tmp_gridalt,tmp_et,vAltGridKm);
}


ublas::vector<double> Venus::TheisModelTeK(ublas::vector<double> vAltGridKm,double vSZADegree)
{
	ublas::vector<double> resu(vAltGridKm.size());
	resu.clear();
	double corrfactor=-1;
	unsigned size=vAltGridKm.size();
	double invert=false;
	if(vAltGridKm[0]<vAltGridKm[size-1])
	{
		invert=true;
		std::reverse(vAltGridKm.begin(),vAltGridKm.end());
	}

	ublas::vector<double> autre_tempe=ElectronKTempScannedModel(vAltGridKm);

	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		if(vAltGridKm[i]>149)
		{
			resu[i]=(FsModTeK(vAltGridKm[i],vSZADegree));
		}else
		{
			if(corrfactor<0)
			{
				corrfactor=resu[i-1]/autre_tempe[i];
			}

			resu[i]=(autre_tempe[i]*corrfactor);
		}
	}
	if(invert)
	{
		std::reverse(resu.begin(),resu.end());
	}
	return resu;



}

double Venus::FsModTeK(double vAltKm,double vSZADegree)
{
	double p[28]={
		3.567296,1.3836078E-02,1.5086544E-02,6.0729804E-03,4.0306677E-03,
		-1.3217842E-02,7.8089219E-03,-2775.023,-123.7310,8.828382,
		-96.94563,15.84634,138.9812,-139.7511,-84.03277,
		0.2649297,-1.444215,0.6347786,-2.192531,-0.5396207,
		0.6054179,1.0569388E-04,-8.0908003E-06,1.3877957E-05,-1.2327257E-05,
		-1.1256760E-05,1.3830228E-05,-1.4350858E-05};
	double sza=vSZADegree*PI/180;

	double s=sin(sza);
	double c=cos(sza);
	double s2=sin(sza*2);
	double c2=cos(sza*2);
	double s3=sin(sza*3);
	double c3=cos(sza*3);
	int k=0;

	double pb[4]={0.,0.,0.,0.};
	for(int j=1;j<5;++j)
	{
		k=j*7-1;// It makes the translation easier from the fortran
      		pb[j-1]=p[k-6]+p[k-5]*s+p[k-4]*c+p[k-3]*s2+p[k-2]*c2+p[k-1]*s3+p[k]*c3;
	}
	return pow(10.,pb[0]+pb[1]/((vAltKm+pb[2])*(vAltKm+pb[2]))+pb[3]*vAltKm);
}



void Venus::LaunchVenGram(ublas::vector<double> vAltGridKm, ublas::matrix<double>& rAtmocm_3 , ublas::vector<double>& rTempK)
{
	rTempK.resize(vAltGridKm.size());
	rTempK.clear();
	rAtmocm_3.resize(vAltGridKm.size(),8);
	rAtmocm_3.clear();


	double clat = mLatDegree;
	double clon = mLoDegree;

	int year = 2011;
	int month = 9;
	int day = 9;
	double ut = 12.42;

	mpParameter->ExistsOrDie("/aero_main/planet/year","You should define the year to work with the earth");
	mpParameter->ExistsOrDie("/aero_main/planet/month","You should define the year to work with the earth");
	mpParameter->ExistsOrDie("/aero_main/planet/day","You should define the year to work with the earth");
	mpParameter->ExistsOrDie("/aero_main/planet/UT","You should define the year to work with the earth");
	
	mpParameter->GetValue("/aero_main/planet/UT",ut);
	mpParameter->GetValue("/aero_main/planet/year",year);
	mpParameter->GetValue("/aero_main/planet/month",month);
	mpParameter->GetValue("/aero_main/planet/day",day);

	double day0 = PhysTime::CalToJul(year, month, day, ut);

	string sdatadir;

	mpParameter->ExistsOrDie("/aero_main/planet/VengramDir","You need to define /aero_main/planet/VengramDir, the directory where the Vengram model data are");
	sdatadir = mpParameter->GetSubFileName("/aero_main/planet/VengramDir");

//	Log::mL<<"THE VENGRAM MODEL IS USED : "<<endl;
//	Log::mL<<"THE TEST IS : "<< mpParameter->SubFileName("/aero_main/planet/VengramDir")<<endl;
//	Log::mL<<"THE DIR IS : "<< sdatadir <<endl;


	sdatadir+="/ ";
//	char* DATADIR = new char[sdatadir.size()];
//	strcpy(DATADIR,sdatadir.c_str());
//	Log::mL<<"THE DATADIR IS : "<< DATADIR <<endl;
	for(size_t j = 0; j<vAltGridKm.size(); ++j)
	{
		// Parameters modified with altitude
		double chgt = vAltGridKm[j];
		// Input constants
		int i=0;
		double csec=0;
		double RH0d=0; // For the perturbation of the density : put 0
		double RHOu=0; // For the perturbation of the density : put 0
		double RHOv=0; // For the perturbation of the density : put 0
		int eof=0; // Put 0, is 1 if a problems appears with the files
		double DELHGT=0; // Modification in altitude : put 0
		double DELLAT=0;// Modification in latitude : put 0
		double DELLON=0;// Modification in longitude : put 0
		double DELTIME=0;// Modification in time : put 0
		double dsunlat=0; // Forces the latitude of the Sun if dradau>0
		double dsunlon=0; // Forces the longitude of the Sun if dradau>0
		double dsunLs=0; // Forces the Ls of the Sun if dradau>0
		double dradau=0; // Forces the distance to the Sun, and the 3 previous parameters
		double dowlt=0; // Venus-earth one way light time (minutes) if dradau>0
		int LonEast=1; // 1 for east longitudes positives
		int IERT=0; // 1 for time input as Earth-receive time, or 0 (recommended) Venus-event time
		int IUTC=1; // 1 for time input as UTC, or 0 for Terrestrial dynamical time

		double pertstep = 0; // perturbation step (0)
		double corlmin = 0; // minimum relative step size for perturbation (0)
		int iupdate = 0; // Storage for warnings
		double profnear =1; // Lat-lon radius within which weight for auxiliary profile is 1.0
		double proffar = 0; // Lat-lon radius beyond which weight for auxiliary profile is 0.0
		int nprof = 1; // Number of profiles (1)


		// The outputs
		double TEMP=0;  // The temperature OUTPUT
		double PRES=0; // The pressure OUTPUT
		double DENSLO=0; // OUTPUT
		double DENS=0;// The density OUTPUT
		double DENSHI=0;// OUTPUT
		double DENSP=0; // OUTPUT
		double EWWIND=0; // OUTPUT
		double EWpert=0; // OUTPUT
		double NSWIND=0; //OUTPUT
		double NSpert=0; // OUTPUT
		double Hrho=0; // OUTPUT
		double HSCALE=0; // OUTPUT
		double corlim=0; //OUTPUT
		double fmolCO2=0; // Molar fraction CO2 OUTPUT
		double fmolN2=0;// Molar fraction N2 OUTPUT
		double fmolO=0;// Molar fraction O OUTPUT
		double fmolCO=0;// Molar fraction CO OUTPUT
		double fmolHe=0;// Molar fraction He OUTPUT
		double fmolN=0;// Molar fraction N OUTPUT
		double fmolH=0;// Molar fraction H OUTPUT
		double ALS=0; // OUTPUT
		double SZA=0; // OUTPUT
		double owlt=0; //OUTPUT 
		double sunlat=0; // OUTPUT
		double sunlon=0; // OUTPUT
		double VenusAU=0; // OUTPUT
		double TLOCAL=0; // OUTPUT		
		double AMz=0; // OUTPUT		
		double DENSTOT=0; //OUTPUT

		datastep_v05_(&i, // Position: used if there is a change. 0 here
				&chgt, // altitude
				&clat, // latitude
				&clon, // longitude
				&csec, // variation in time: put 0
				&day0, // the julian day - use the function
				&RH0d, // For the perturbation of the density : put 0
				&RHOu, // For the perturbation of the density : put 0
				&RHOv, // For the perturbation of the density : put 0
				&eof, // Put 0, is 1 if a problems appears with the files
				&DELHGT, // Modification in altitude : put 0
				&DELLAT,// Modification in latitude : put 0
				&DELLON,// Modification in longitude : put 0
				&DELTIME,// Modification in time : put 0
				&TEMP,  // The temperature OUTPUT
				&PRES, // The pressure OUTPUT
				&DENSLO, // OUTPUT
				&DENS,// The density OUTPUT
				&DENSHI,// OUTPUT
				&DENSP, // OUTPUT
				&EWWIND, // OUTPUT
				&EWpert, // OUTPUT
				&NSWIND, //OUTPUT
				&NSpert, // OUTPUT
				&Hrho, // OUTPUT
				&HSCALE, // OUTPUT
				&dsunlat, // Forces the latitude of the Sun if dradau>0
				&dsunlon, // Forces the longitude of the Sun if dradau>0
				&dsunLs, // Forces the Ls of the Sun if dradau>0
				&dradau, // Forces the distance to the Sun, and the 3 previous parameters
				&dowlt, // Venus-earth one way light time (minutes) if dradau>0
				&LonEast, // 1 for east longitudes positives
				&corlim, //OUTPUT
				&DENSTOT, //OUTPUT
				&IERT, // 1 for time input as Earth-receive time, or 0 (recommended) Venus-event time
				&IUTC, // 1 for time input as UTC, or 0 for Terrestrial dynamical time
				&fmolCO2, // Molar fraction CO2 OUTPUT (! in %)
				&fmolN2,// Molar fraction N2 OUTPUT (! in %)
				&fmolO,// Molar fraction O OUTPUT (! in %)
				&fmolCO,// Molar fraction CO OUTPUT (! in %)
				&fmolHe,// Molar fraction He OUTPUT (! in %)
				&fmolN,// Molar fraction N OUTPUT (! in %)
				&fmolH,// Molar fraction H OUTPUT (! in %)
				&AMz, // OUTPUT
				&pertstep, // perturbation step (0)
				&corlmin, // minimum relative step size for perturbation (0)
				&iupdate, // Storage for warnings
				&ALS, // OUTPUT
				&SZA, // OUTPUT
				&owlt, //OUTPUT 
				&sunlat, // OUTPUT
				&sunlon, // OUTPUT
				&VenusAU, // OUTPUT
				&TLOCAL, // OUTPUT
				&profnear, // Lat-lon radius within which weight for auxiliary profile is 1.0
				&proffar, // Lat-lon radius beyond which weight for auxiliary profile is 0.0
				&nprof, // Number of profiles (1)
				const_cast<char*>(sdatadir.c_str())); // Array of 300 char for the directory where the data are!
		//		DATADIR); // Array of 300 char for the directory where the data are!
		double densf = DENSTOT  / AMz  / 1E5 * AVOGADRO; // 6.02214179E-17; //AVOGADRO ;	
		// DENSTOT is in kg/m**3, AMz in kg/k-mole therefore DENSTOT/AMz * 1E3 is in mole/m**3; /1E6 to be in /cm**3; /100 because the molar fractions are in percentage and AVOGADRO for the nb of particles
		rAtmocm_3(j,0) = densf * fmolCO2;
		rAtmocm_3(j,1) = densf * fmolN2;
		rAtmocm_3(j,2) = densf * fmolCO;
		rAtmocm_3(j,3) = densf * fmolO;
		rAtmocm_3(j,4) = densf * fmolCO2 / 1E3; // We apply the Fox technique to estiimate O2
		rAtmocm_3(j,5) = densf * fmolHe;
		rAtmocm_3(j,6) = densf * fmolN;
		rAtmocm_3(j,7) = densf * fmolH;
		rTempK[j] = TEMP;
	}
	//delete[] DATADIR;
}

ublas::vector<double> Venus::ReturndB_B(const ublas::vector<double>& vAltGridKm)
{
	ublas::vector<double> db_b;
	//int resu;
	//resu = 
	ReaddB_B(vAltGridKm, db_b);
	return db_b;
}
