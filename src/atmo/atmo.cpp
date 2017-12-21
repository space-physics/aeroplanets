/** 
 * \file atmo.cpp
 * \brief Implements the Atmo class : allows to define the grids and start the different computation for ionization and excitation
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: atmo.cpp 1921 2014-02-26 05:37:14Z gronoff $
 *
 */

#include "atmo.hpp"
using namespace std;
using namespace boost::assign;

Atmo::Atmo(std::string vParameterFile)
{
	mbIsChem=false;
	mIsBendingActivated=false;
	mIsPlanet=false;
	mParameterFile=vParameterFile;
	mpParameter=new XmlParameters(mParameterFile);

	// The multiplicators are declared as 0 here
	// this is not a mistake: these x will be
	// updated after.
	// 
	// They are declared as 0 because if someone wants to use them 
	// before their real initialisation, the error will be obvious!
	mElecDensMult=0;
	mElecTempMult=0;
	mTempMult=0;


	if(mpParameter->Exists("/aero_main/Log"))
	{
		if(mpParameter->KeyExists("/aero_main/Log","verbosity"))
		{
			string verbo=trim(mpParameter->GetKey("/aero_main/Log","verbosity"));
			if("Config"==verbo)
			{
				Log::msCoutPriority=Log::CONFIG;
			}
			if("Info"==verbo)
			{
				Log::msCoutPriority=Log::INFO;
			}
			if("Section"==verbo)
				Log::msCoutPriority=Log::SECTION;
			if("Warning"==verbo)
			{
				Log::msCoutPriority=Log::WARNING;
			}

			if("ERROR"==verbo)
			{
				Log::msCoutPriority=Log::ERROR;
				Log::msFilePriority=Log::INFO;
			}

		}
	}

	if(mpParameter->Exists("/aero_main/ExtraInfo"))
	{
		try{
		string extrainfo = mpParameter->Elem("/aero_main/ExtraInfo");
		Log::msMessageLog += "# " + extrainfo + "\n";
		}
		catch(...)
		{
			Log::mW<<"The ExtraInfo markup is probably empty."<<endl;
		}
	}



	// We want to use the MonteCarlo system if it is defined
	if(mpParameter->Exists("/aero_main/SetMonteCarloActive"))
	{
		Log::mI<<"Monte Carlo analysis activated!!!"<<endl;
		mpParameter->SetMonteCarloActive();

	}


	mIsMainParameter=true;
	mIsModelAtmosphere=false;
	mIsPhotoionizationModel=false;
	mIsElectroionizationModel=false;
	mIsGaussianAngle=false;
	mbIsPhotoFluxDefined=false;
	mbIsCosmoFluxDefined=false;
	mbIsProtonEFluxDefined=false;
	mbIsElecPrecipFluxDefined=false;
	mbIsElecModelDefined=false;
	mbIsPhotoFluxSZAs=false;
	mbIsCosmoFluxSZAs=false;
	mbIsProtonFluxSZAs=false;
	mbIsCosmoModel=false;
	mbIsProtonModel=false;
	InitGrids();
	InitPlanet();
	InitNeutralAtmo();
	InitNeutralTemp();
	InitIonosphere();
	mpModelAtmosphere->ComputeVerticalColumnDensity(mNeutralTemperature);

	mMultiplePhoto=false;
	mMultipleCosmo=false;
	mMultipleProton=false;
	mMultipleElec=false;
	mMultipleChem=false;
	mMultipleEmit=false;
}

Atmo::~Atmo()
{
	if(mIsMainParameter)
	{
		// On supprime le parametre de fichier
		delete mpParameter;
	}
	if(mIsModelAtmosphere)
	{
		delete mpModelAtmosphere;
	}

	if(mIsPlanet)
	{
		delete mpPlanet;
	}
	if(mIsGaussianAngle)
	{
		delete mpGAngle;
	}
	if(mIsPhotoionizationModel)
	{
		Log::mD<<"Delete mpPhotoModel"<<endl;
		delete mpPhotomodel;
		Log::mD<<"End delete"<<endl;
	}

	if(mIsElectroionizationModel)
	{
		Log::mD<<"Delete mpElecModel"<<endl;
		delete mpElecModel;
		Log::mD<<"End delete"<<endl;
	}
	// We delete the vectors of species.
	for(unsigned i=0;i<mPhotoionizationResu.size();++i)
	{
		delete mPhotoionizationResu[i];
	}
	for(unsigned i=0;i<mCosmoionizationResu.size();++i)
	{
		delete mCosmoionizationResu[i];
	}
	for(unsigned i=0;i<mProtoionizationResu.size();++i)
	{
		delete mProtoionizationResu[i];
	}


	for(unsigned i=0;i<mTotalResu.size();++i)
	{
		delete mTotalResu[i];
	}

	for(unsigned i=0;i<mElectronImpactResu.size();++i)
	{
		delete mElectronImpactResu[i];
	}


	if(mbIsPhotoFluxDefined)
	{
		delete mpPhotoFlux;
	}

	if(mbIsCosmoFluxDefined)
	{
		delete mpCosmoFlux;
	}
	if(mbIsCosmoModel)
	{
		delete mpCosmoModel;
	}

	if(mbIsProtonEFluxDefined)
	{
		delete mpProtonEFlux;
	}
	if(mbIsProtonModel)
	{
		delete mpProtonModel;
	}
	if(mbIsElecPrecipFluxDefined)
	{
		delete mpTotalFlux;
		delete mpElecPrecipFlux;
	}
	if(mbIsElecModelDefined)
	{
		delete mpElecModel;
	}


	for(unsigned i=0;i<mPhotoFluxSZAs.size();++i)
	{
		delete mPhotoFluxSZAs[i];
	}
	for(unsigned i=0;i<mResuPhotoSZAs.size();++i)
	{
		for(unsigned j=0;j<mResuPhotoSZAs[i].size();++j)
		{
			delete mResuPhotoSZAs[i][j];
		}
	}

	for(unsigned i=0;i<mCosmoFluxSZAs.size();++i)
	{
		delete mCosmoFluxSZAs[i];
	}
	if(mResuCosmoSZAs.size()>0)
	{
		for(unsigned i=0;i<mResuCosmoSZAs[0].size();++i)
		{
			delete mResuCosmoSZAs[0][i];
		}
	}

	for(unsigned i=0;i<mProtonEFluxSZAs.size();++i)
	{
		delete mProtonEFluxSZAs[i];
	}
	if(mResuProtonSZAs.size()>0)
	{
		for(unsigned i=0;i<mResuProtonSZAs[0].size();++i)
		{
			delete mResuProtonSZAs[0][i];
		}
	}
	for(size_t i=0;i<mAtmoBentSpecies.size();++i)
	{
		delete mAtmoBentSpecies[i];
	}
	for(size_t i=0;i<mAtmoSpeciesSZAs.size();++i)
	{
		for(size_t j=0;j<mAtmoSpeciesSZAs[i].size();++j)
		{
			delete mAtmoSpeciesSZAs[i][j];
		}

	}
}


void Atmo::InitChem()
{

	//Log::SetPriority(Log::DEBUGG);
	// we init the altitude grid
	Log::mS<<"Initialization of the chemistry"<<endl;


	mChem.reset(	new Chem(mpParameter,
				mpPlanet,
				mAltGridKm,
				mpModelAtmosphere->mAtmoSpecies,
				mTotalResu,
				mPhotoionizationResu,
				mElectronImpactResu,
				mProtoionizationResu,
				mCosmoionizationResu,
				mNeutralTemperature,
				mElectronTemperature,
				mIonTemperature,
				mElectronDensity,
				mElecDensMult)
		   );
}

void Atmo::InitGrids()
{
	//Log::SetPriority(Log::DEBUGG);
	// we init the altitude grid
	Log::mS<<"Initialization of the altitude grid"<<endl;
	InitAltGrid();
	// we init the electron grids
	Log::mS<<"Initialization of the electron grid"<<endl;
	InitElecGrids();
	// we init the photon grid
	Log::mS<<"Initialization of the photons grid"<<endl;
	InitPhotGrid();
	// We init the proton grid
	Log::mS<<"Initialization of the proton grid"<<endl;
	InitProtGrid();

	Log::mS<<"End of grid initialization"<<endl;

}

void Atmo::InitPhotGrid()
{
	int photgridparam=mpParameter->NbParams("/aero_main/sun/grid/use_model","type");
	if(photgridparam!=1)
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Error: you do not define, or you multiply defined /sun/grid/use_model in the file "<<mParameterFile<<endl;
		Error err("InitElecGrid","Electrons grid in error","You do not define, or you multiply defined /sun/grid/use_model in the file "+mParameterFile);
		throw err;
	}
	int modeltype=0;
	mpParameter->GetNKey("/aero_main/sun/grid/use_model","type",modeltype);
	if(modeltype==0)
	{// Toor et Toor initialization. There is no more overlapping in that grid

		std::vector<double> PhotonGrideVmax,PhotonGrideVmin;

		PhotonGrideVmax+=652.548421,413.280667,247.9684,123.9842,82.6561333,61.9921,49.59368,48.3846391,48.3746391,43.6433627,41.3280667,40.8870565,40.8770565,40.8238126,40.8138126,35.4240571,33.6949512,33.6849512,30.99605,27.5520444,26.6606599,26.6506599,24.79684,22.5425818,22.3748827,22.3648827,21.2281815,21.2181815,20.6640333,20.3332787,19.6984697,19.6884697,19.0744923,17.7120286,17.6386701,17.6286702,16.5312267,16.2139077,16.2039077,16.103275,16.093275,15.7169271,15.7069271,15.498025,14.5863765,13.7760222,13.0509684,12.7000371,12.6900371,12.39842,12.0975288,12.0875288,12.0250207,12.0150207;
		PhotonGrideVmin+= 413.280667,247.9684,123.9842,82.6561333,61.9921,49.59368,48.3846391,48.3746391,43.6333627,41.3280667,40.8870565,40.8770565,40.8238126,40.8138126,35.4240571,33.6949512,33.6849512,30.99605,27.5520444,26.6606599,26.6506599,24.79684,22.5425818,22.3748827,22.3648827,21.2281815,21.2181815,20.6640333,20.3332787,19.6984697,19.6884697,19.0744923,17.7120286,17.6386701,17.6286702,16.5312267,16.2139077,16.2039077,16.103275,16.093275,15.7169271,15.7069271,15.498025,14.5863765,13.7760222,13.0509684,12.7000371,12.6900371,12.39842,12.0975288,12.0875288,12.0250207,12.0150207,11.8080191;


		mPhotonGrideVmax.resize(PhotonGrideVmax.size());
		mPhotonGrideVmin.resize(PhotonGrideVmin.size());
		std::copy(PhotonGrideVmax.begin(),PhotonGrideVmax.end(),mPhotonGrideVmax.begin());
		std::copy(PhotonGrideVmin.begin(),PhotonGrideVmin.end(),mPhotonGrideVmin.begin());
		mPhotonGrideV.resize(mPhotonGrideVmin.size());
		for(unsigned i=0;i<mPhotonGrideVmin.size();++i)
		{
			mPhotonGrideV[i]=(mPhotonGrideVmin[i]+(mPhotonGrideVmax[i]-mPhotonGrideVmin[i])*0.5);
		}
		assert(mPhotonGrideV.size()==mPhotonGrideVmin.size());
		return;
	}else if(modeltype==1)
	{
		int type=0;
		int number=0;
		double min=0.;
		double max=0.;

		//cout<<"You are using a standard grid for the photons..."<<endl;
		mpParameter->GetNKey("/aero_main/sun/grid/st_grid","type",type);
		mpParameter->GetValue("/aero_main/sun/grid/st_grid/emin",min);
		mpParameter->GetValue("/aero_main/sun/grid/st_grid/emax",max);
		mpParameter->GetValue("/aero_main/sun/grid/st_grid/number",number);
		//		StandardGrid
		if(type<0 and type > 2)
		{
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Photon grid type not found"<<endl;
			Error err("InitPhotGrid","Photon grid type not found","The photon grid has a bad type");
			throw err;
		}
		if(type==0)
		{
			mPhotonGrideVmin=MathGrid::GridExp(min,max,number);
			std::reverse(mPhotonGrideVmin.begin(),mPhotonGrideVmin.begin());
		}
		if(type==1)
		{
			mPhotonGrideVmin=MathGrid::GridPow(min,max,number);
			std::reverse(mPhotonGrideVmin.begin(),mPhotonGrideVmin.begin());
		}
		if(type==2)
		{
			mPhotonGrideVmin=MathGrid::GridCst(min,max,number);
			std::reverse(mPhotonGrideVmin.begin(),mPhotonGrideVmin.begin());
		}
		if(number>2)
		{// Si on n'a pas assez de points...
			mPhotonGrideVmax.resize(number);
			mPhotonGrideVmax[0]=(2*mPhotonGrideVmin[0]-mPhotonGrideVmin[1]);
			for(int i=1;i<number;++i)
			{
				mPhotonGrideVmax[i]=(mPhotonGrideVmin[i-1]);
			}

			SolarGridLines lines(mpParameter);
			Log::mD<<"Before SL : "<<mPhotonGrideVmin.size()<<endl;
			lines.AddSolarLines(mPhotonGrideVmin,mPhotonGrideVmax);
			Log::mD<<"After SL : "<<mPhotonGrideVmin.size()<<endl;



			assert(mPhotonGrideVmax.size()==mPhotonGrideVmin.size());



			mPhotonGrideV.resize(mPhotonGrideVmax.size());
			for(unsigned i=0;i<mPhotonGrideVmin.size();++i)
			{
				mPhotonGrideV[i]=(mPhotonGrideVmin[i]+(mPhotonGrideVmax[i]-mPhotonGrideVmin[i])*0.5);
			}
			assert(mPhotonGrideV.size()==mPhotonGrideVmin.size());

		}

		Log::mD<<"Your photon grid"<<endl;
		Log::mD<<mPhotonGrideV<<endl;
		//MathString::print1d(mPhotonGrideV);
		return;
	}else if(modeltype==2)
	{// Toor et Toor initialization, with the OW data There is no more overlapping in that grid

		std::vector<double> PhotonGrideVmax,PhotonGrideVmin;

	//	PhotonGrideVmax+=652.548421,413.280667,247.9684,123.9842,82.6561333,61.9921,49.59368,48.3846391,48.3746391,43.6433627,41.3280667,40.8870565,40.8770565,40.8238126,40.8138126,35.4240571,33.6949512,33.6849512,30.99605,27.5520444,26.6606599,26.6506599,24.79684,22.5425818,22.3748827,22.3648827,21.2281815,21.2181815,20.6640333,20.3332787,19.6984697,19.6884697,19.0744923,17.7120286,17.6386701,17.6286702,16.5312267,16.2139077,16.2039077,16.103275,16.093275,15.7169271,15.7069271,15.498025,14.5863765,13.7760222,13.0509684,12.7000371,12.6900371,12.39842,12.0975288,12.0875288,12.0250207,12.0150207,11.808019047619048,10.332016666666668,10.208836855396614,10.188836855396614,9.0499416058394164,8.9056635288715427,8.8856635288715431,8.8560142857142861,8.8485266294545806,8.828526629454581,8.5506344827586211,8.2656133333333344,8.0183323106337081,8.0050089310471577,7.998980645161291,7.9983323106337085,7.9850089310471581,7.9526137091607954,7.9326137091607958,7.749012500000001,7.5141939393939401,7.4915471880279991,7.4715471880279996,7.2931882352941182;
	//	PhotonGrideVmin+=413.280667,247.9684,123.9842,82.6561333,61.9921,49.59368,48.3846391,48.3746391,43.6333627,41.3280667,40.8870565,40.8770565,40.8238126,40.8138126,35.4240571,33.6949512,33.6849512,30.99605,27.5520444,26.6606599,26.6506599,24.79684,22.5425818,22.3748827,22.3648827,21.2281815,21.2181815,20.6640333,20.3332787,19.6984697,19.6884697,19.0744923,17.7120286,17.6386701,17.6286702,16.5312267,16.2139077,16.2039077,16.103275,16.093275,15.7169271,15.7069271,15.498025,14.5863765,13.7760222,13.0509684,12.7000371,12.6900371,12.39842,12.0975288,12.0875288,12.0250207,12.0150207,11.8080191,10.332016666666668,10.208836855396614,10.188836855396614,9.0499416058394164,8.9056635288715427,8.8856635288715431,8.8560142857142861,8.8485266294545806,8.828526629454581,8.5506344827586211,8.2656133333333344,8.0183323106337081,8.0050089310471577,7.998980645161291,7.9983323106337085,7.9850089310471581,7.9526137091607954,7.9326137091607958,7.749012500000001,7.5141939393939401,7.4915471880279991,7.4715471880279996,7.2931882352941182,7.0848114285714292;


		PhotonGrideVmax+=652.54842105,413.28066667,247.9684,123.9842,82.65613333,61.9921,49.59368,48.37563909,48.37363909,43.63436266,43.63236266,41.32806667,40.87805648,40.87605648,40.81481263,40.81281263,35.42405714,33.68595123,33.68395123,30.99605,27.55204444,26.6516599,26.6496599,24.79684,22.54258182,22.36588266,22.36388266,21.21918151,21.21718151,20.66403333,20.33427867,20.33227867,19.68946966,19.68746966,19.07449231,17.71202857,17.62967015,17.62767015,16.53122667,16.20490773,16.20290773,16.09427501,16.09227501,15.70792713,15.70592713,15.498025,14.58637647,13.77602222,13.05096842,12.69103705,12.68903705,12.39842,12.08852876,12.08652876,12.01602069,12.01402069,11.80801905,10.33201667,10.19983686,10.19783686,9.04994161,8.89666353,8.89466353,8.85601429,8.83952663,8.83752663,8.55063448,8.26561333,8.00933231,8.00733231,7.99898065,7.99600893,7.99400893,7.94361371,7.94161371,7.7490125,7.51419394,7.48254719,7.48054719,7.29318824;
		PhotonGrideVmin+=413.28066667,247.9684,123.9842,82.65613333,61.9921,49.59368,48.37563909,48.37363909,43.63436266,43.63236266,41.32806667,40.87805648,40.87605648,40.81481263,40.81281263,35.42405714,33.68595123,33.68395123,30.99605,27.55204444,26.6516599,26.6496599,24.79684,22.54258182,22.36588266,22.36388266,21.21918151,21.21718151,20.66403333,20.33427867,20.33227867,19.68946966,19.68746966,19.07449231,17.71202857,17.62967015,17.62767015,16.53122667,16.20490773,16.20290773,16.09427501,16.09227501,15.70792713,15.70592713,15.498025,14.58637647,13.77602222,13.05096842,12.69103705,12.68903705,12.39842,12.08852876,12.08652876,12.01602069,12.01402069,11.80801905,10.33201667,10.19983686,10.19783686,9.04994161,8.89666353,8.89466353,8.85601429,8.83952663,8.83752663,8.55063448,8.26561333,8.00933231,8.00733231,7.99898065,7.99600893,7.99400893,7.94361371,7.94161371,7.7490125,7.51419394,7.48254719,7.48054719,7.29318824,7.08481143;
		mPhotonGrideVmax.resize(PhotonGrideVmax.size());
		mPhotonGrideVmin.resize(PhotonGrideVmin.size());
		std::copy(PhotonGrideVmax.begin(),PhotonGrideVmax.end(),mPhotonGrideVmax.begin());
		std::copy(PhotonGrideVmin.begin(),PhotonGrideVmin.end(),mPhotonGrideVmin.begin());
		mPhotonGrideV.resize(mPhotonGrideVmin.size());
		for(unsigned i=0;i<mPhotonGrideVmin.size();++i)
		{
			mPhotonGrideV[i]=(mPhotonGrideVmin[i]+(mPhotonGrideVmax[i]-mPhotonGrideVmin[i])*0.5);
		}
		assert(mPhotonGrideV.size()==mPhotonGrideVmin.size());
		return;
	}




	//Log::SetPriority(Log::ERROR);
	Log::mE<<"Error : sun model undefined"<<endl;
	Error err("InitPhotGrid","sun model undefined","Define your sun model");
	throw err;

}

void Atmo::InitAltGrid()
{
	int altgridparam=mpParameter->NbParams("/aero_main/atmosphere/alt_grid/use_model","type");
	if(altgridparam!=1)
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Error: you do not define, or you multiply defined /atmosphere/alt_grid/use_model in the file "<<mParameterFile<<endl;
		Error err("InitElecGrid","Electrons grid in error","You do not define, or you multiply defined /atmosphere/alt_grid/use_model in the file "+mParameterFile);
		throw err;
	}
	int modeltype=0;

	mpParameter->GetNKey("/aero_main/atmosphere/alt_grid/use_model","type",modeltype);

	if(modeltype==0)
	{// We use the standard model
		int type=0;
		int number=0;
		double min=0.;
		double max=0.;

		mpParameter->GetNKey("/aero_main/atmosphere/alt_grid/st_grid","type",type);
		mpParameter->GetValue("/aero_main/atmosphere/alt_grid/st_grid/altmin",min);
		mpParameter->GetValue("/aero_main/atmosphere/alt_grid/st_grid/altmax",max);
		mpParameter->GetValue("/aero_main/atmosphere/alt_grid/st_grid/number",number);
		//		StandardGrid
		if(type<0 and type > 2)
		{

			//Log::SetPriority(Log::ERROR);
			Log::mE<<"Altitude grid type not found"<<endl;
			Error err("InitAltGrid","Altitude grid type not found","The altitude grid has a bad type");
			throw err;
		}

		if(type==0)
		{
			mAltGridKm=MathGrid::GridExp(min,max,number);
		}
		if(type==1)
		{
			mAltGridKm=MathGrid::GridPow(min,max,number);
		}
		if(type==2)
		{
			mAltGridKm=MathGrid::GridCst(min,max,number);
		}
	}else
	{
		MathString::LitToutString(mpParameter->Elem("/aero_main/atmosphere/alt_grid/altdata"),mAltGridKm);
	}


}

void Atmo::InitElecGrids()
{

	mpParameter->ExistsOrDie("/aero_main/electron/nb_angles","You have to define the number of angle for the multistream computation. Even if you do not compute this....");
	mpParameter->GetValue("/aero_main/electron/nb_angles",mNbAngle);
	mpGAngle=new MathFunction::GaussianAngle(mNbAngle);
	mIsGaussianAngle=true;
	
	int elecgridparam=mpParameter->NbParams("/aero_main/electron/grid/st_grid","type");
	if(elecgridparam!=1)
	{

		///Log::SetPriority(Log::ERROR);
		Log::mE<<"Error: you do not define, or you multiply defined /aero_main/electron/grid/st_grid in the file "<<mParameterFile<<endl;
		Error err("InitElecGrid","Electrons grid in error","You do not define, or you multiply defined /electron/grid/st_grid in the file "+mParameterFile);
		throw err;
	}

	int type=0;
	int number=0;
	double min=0.;
	double max=0.;

	mpParameter->GetNKey("/aero_main/electron/grid/st_grid","type",type);
	mpParameter->GetValue("/aero_main/electron/grid/st_grid/emin",min);
	mpParameter->GetValue("/aero_main/electron/grid/st_grid/emax",max);
	mpParameter->GetValue("/aero_main/electron/grid/st_grid/number",number);
	//		StandardGrid
	if(type<0 and type > 1)
	{

		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Electron grid type not found"<<endl;
		Error err("InitElecGrid","Electron grid type not found","The electron grid has a bad type");
		throw err;
	}
	if(type==0)
	{
		mElecCentEeV=MathGrid::GridExp(min,max,number);
		mSpfactor=-1;
		mElecDdengeV=MathGrid::WidthGrid(mElecCentEeV,1);
	}
	if(type==1)
	{
		//	cout<<"GridPolo"<<endl;
		if(!MathGrid::GridPolo(number,min,max,mElecCentEeV,mElecDdengeV,mSpfactor))
		{
			Error err("InitElecGrids","Gridpolo","Impossible to initialize your energy grid");
			throw err;
		}
		//cout<<"Endgridpolo"<<endl;
		std::reverse(mElecCentEeV.begin(),mElecCentEeV.end());
		std::reverse(mElecDdengeV.begin(),mElecDdengeV.end());
	}
	// we compute the bottom grid
	unsigned esize=mElecCentEeV.size();
	double ebot=0.001;
	mElecBotEeV.resize(esize,0.);
	mElecBotEeV[esize-1]=mElecCentEeV[esize-1]-mElecDdengeV[esize-1]/2.;
	//cout<<"Compute bot and centE"<<endl;
	for(int ien=esize-2;ien>-1;--ien)
	{
		mElecBotEeV[ien]=mElecBotEeV[ien+1]+mElecDdengeV[ien+1];
	}
	if (mElecBotEeV[esize-1]<ebot)
	{
		mElecBotEeV[esize-1]=ebot;
		mElecDdengeV[esize-1]-=ebot;
	}

	if (min>1.)
	{
		for(unsigned ien=0;ien<esize;++ien)
		{
			mElecBotEeV[ien]+=min-1.;
			mElecCentEeV[ien]+=min-1.;
		}
	}

	// We check that the grid is decreasing
	assert(mElecCentEeV[0]>*(mElecCentEeV.end()-1));

}

void Atmo::InitProtGrid()
{
	if(!mpParameter->Exists("/aero_main/proton/use_proton"))
	{
		Log::mS<<"You do not use protons"<<endl;
		return;
	}

	mpParameter->ExistsOrDie("/aero_main/proton/nb_angles","You have to define the number of angle for the multistream computation. Even if you do not compute this....");
	mpParameter->GetValue("/aero_main/proton/nb_angles",mNbProtonAngle);
	boost::shared_ptr<MathFunction::GaussianAngle> tmpga(new MathFunction::GaussianAngle(mNbProtonAngle));
	mpProtonGaussianAngle = tmpga;
	int protgridparam=mpParameter->NbParams("/aero_main/proton/grid/st_grid","type");
	if(protgridparam!=1)
	{

		///Log::SetPriority(Log::ERROR);
		Log::mE<<"Error: you do not define, or you multiply defined /aero_main/proton/grid/st_grid in the file "<<mParameterFile<<endl;
		Error err("InitProtGrid","Proton grid in error","You do not define, or you multiply defined /proton/grid/st_grid in the file "+mParameterFile);
		throw err;
	}

	int type=0;
	int number=0;
	double min=0.;
	double max=0.;

	mpParameter->GetNKey("/aero_main/proton/grid/st_grid","type",type);
	mpParameter->GetValue("/aero_main/proton/grid/st_grid/emin",min);
	mpParameter->GetValue("/aero_main/proton/grid/st_grid/emax",max);
	mpParameter->GetValue("/aero_main/proton/grid/st_grid/number",number);
	//		StandardGrid
	if(type<0 and type > 2)
	{

		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Proton grid type not found"<<endl;
		Error err("InitProtGrid","Proton grid type not found","The proton grid has a bad type");
		throw err;
	}
	if(type==0)
	{
		mProtonGrideV=MathGrid::GridExp(min,max,number);
	//	mElecCentEeV=MathGrid::GridExp(min,max,number);
	//	mSpfactor=-1;
		mProtonWidthGrideV=MathGrid::WidthGrid(mProtonGrideV,1);
	}
	if(type==1)
	{
		//	cout<<"GridPolo"<<endl;
		//ublas::vector<double> dummy;
		double factor;
		if(!MathGrid::GridPolo(number,min,max,mProtonGrideV,mProtonWidthGrideV,factor))
		{
			Error err("InitElecGrids","Gridpolo","Impossible to initialize your energy grid");
			throw err;
		}
		std::reverse(mProtonGrideV.begin(),mProtonGrideV.end());
		std::reverse(mProtonWidthGrideV.begin(),mProtonWidthGrideV.end());
	}

	if(2==type)
	{
		mpParameter->Get1DArray("/aero_main/proton/grid/st_grid/data",mProtonGrideV);
		if (mProtonGrideV[0]<*(mProtonGrideV.end()-1))
		{
			std::reverse(mProtonGrideV.begin(),mProtonGrideV.end());
		}
		mProtonWidthGrideV=MathGrid::WidthGrid(mProtonGrideV,1);
		assert(mProtonGrideV.size() > 0);
	/*	mProtonWidthGrideV[0] = 0;
		for(size_t e = 1; e < mProtonGrideV.size(); ++e)
		{
			mProtonWidthGrideV[e] = mProtonGrideV[e] - mProtonGrideV[e - 1];

		}*/
	}

	// We check that the grid is decreasing
	assert(mProtonGrideV[0]>*(mProtonGrideV.end()-1));

}

void Atmo::InitPlanet()
{

	//Log::SetPriority(Log::INFO);


	if(!mpParameter->Exists("/aero_main/planet/name"))
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Error: you must define a valid planet"<<endl;
		Log::mE<<"ie Titan, Mars or Venus for now"<<endl;
		Error err("Atmo::InitPlanet","exists /planet/name","Error: you must define a valid planet (mars or venus now...)");
		throw err;
	}
	string planetname=trim(mpParameter->Elem("/aero_main/planet/name"));
	Log::AddVersionInfo("# Planet:"+planetname);

	if(planetname=="Venus")
	{
		Log::mS<<"You are working on Venus: take care of temperature and pressure"<<endl;
		Log::mS<<"It is a bit hot here"<<endl;
		mpPlanet=new Venus(mpParameter);	
		mIsPlanet=true;
		return;
	}


	if(planetname=="Mars")
	{
		Log::mS<<"You are working on Mars: it can be cold here. But it is easy to boil water!"<<endl;
		mpPlanet=new Mars(mpParameter);
		mIsPlanet=true;
		return;
	}

	if(planetname=="Titan")
	{
		Log::mS<<"Titan, I whish you can admire the Saturn's rings! It's a bit cold here, but pressure is Ok! Don't read too much Baxter books"<<endl;
		mpPlanet=new Titan(mpParameter);
		mIsPlanet=true;
		return;
	}

	if(planetname=="Earth" || planetname=="Terre")// hehehe option secrète!
	{
		Log::mS<<"Earth... This is really an original choice, did you manage to move to work today?"<<endl;
		mpPlanet=new Earth(mpParameter);
		mIsPlanet=true;
		return;
	}

	if(planetname=="Anonymous" || planetname=="Pwett")// re option secrete
	{
		Log::mS<<"Anonymous! I can't imagine the tons of discoveries; add me in the co-author list ;-) "<<endl;
		mpPlanet=new AnonymousPlanet(mpParameter);
		mIsPlanet=true;
		return;
	 }

	Error err("Atmo::InitPlanet","No existing planet detected","Houston, we'we got a problem: the planet does not exist!");
	throw err;

}





void Atmo::InitNeutralAtmo()
{

	Log::mS<<"Initialisation of the neutral atmosphere"<<endl;
	mpModelAtmosphere=new NeutralAtmo(mpParameter,mAltGridKm,mpPlanet,&mPhotonGrideV,&mElecCentEeV,&mElecDdengeV,&mProtonGrideV);
	mdB_B = mpPlanet->ReturndB_B(mAltGridKm);
	mIsModelAtmosphere=true;

}

void Atmo::InitNeutralTemp()
{
	/// By default now, but can evolve 
	mNeutralTemperature=mpModelAtmosphere->Temperature();

	if(!mpParameter->Exists("/aero_main/atmosphere/neutral/temperature/multiplicator"))
	{
		mTempMult=1.;
	}else
	{
		mpParameter->GetValue("/aero_main/atmosphere/neutral/temperature/multiplicator",mTempMult);
		mNeutralTemperature*=mTempMult;
	}
}


void Atmo::InitIonosphere()
{
        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-+
         |                |
         |Electron density|
         |                |
         +---------------*/
	//Log::SetPriority(Log::DEBUGG);
	Log::mS<<"Initialization of the ionosphere"<<endl;
	Log::mD<<"Init electron density"<<endl;
	mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/density/use_model","You have to define the electron density model");

	int type=0;
	mpParameter->GetNKey("/aero_main/atmosphere/iono/electron/density/use_model","type",type);

	if(type==0)
	{
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/density/alt","You have to define the altitude, because you use the data model");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/density/ne","You have to define the electron density, because you use the data model");
		ublas::vector<double> alt,ne;
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/electron/density/alt",alt);
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/electron/density/ne",ne);
		if(mpParameter->KeyExists("/aero_main/atmosphere/iono/electron/density/ne","DeleteZero"))
		{
			Log::mD<<"The zero values in your electron densities are set to 1E-40"<<endl;
			for(unsigned i=0; i< ne.size(); ++i)
			{
				if(ne[i] < 1E-40)
					ne[i] = 1E-40;
			}
		}	
		
		mElectronDensity=MathFunction::IntLog(alt,ne,mAltGridKm);
	}else if(type==-1)
	{
		TiXmlNode* node=mpParameter->GetNode("/aero_main/atmosphere/iono/electron/density");

		int model=0;

		mpParameter->GetNKey(node,"//Model","type",model);
		switch(model)
		{
			case 0:
				{
					// Not really physical here, but interesting to check
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					double T0=0;
					mpParameter->GetValue(node,"//Model/To",T0);
					double shape=0;
					mpParameter->GetValue(node,"//Model/Shape",shape);
					//unsigned pos=CloseInVector(alt0,mAltGridKm);
					mElectronDensity=MathFunction::BatesWalkerProfile(mAltGridKm,alt0,dens0,T0,Texo,shape,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);

				}
				break;

			case 1:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double SH=0;
					mpParameter->GetValue(node,"//Model/SH",SH);
					mElectronDensity=MathFunction::ChapmanProfile(mAltGridKm,alt0,dens0,SH);
				}
				break;
			case 2:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double stddev=0;
					mpParameter->GetValue(node,"//Model/Dev",stddev);
					mElectronDensity=MathFunction::GaussianProfile(mAltGridKm,alt0,dens0,stddev);
				}
				break;
			case 3:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					mElectronDensity=MathFunction::ExpProfile(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);

				}
				break;
			case 4:
				{
					ublas::vector<double> alts,logval;
					mpParameter->Get1DArray(node,"//Model/altitudes",alts);
					mpParameter->Get1DArray(node,"//Model/logvalues",logval);
					if(alts.size()!=logval.size())
					{
						Log::mE<<"Size mismatch: The size of Model/altitudes and Model/logvalues are not compatible. please check"<<endl;
						Error err("ReadParameters: add_parametrized_species case4","Size mismatch","The size of Model//altitudes and Model/logvalues are not compatible. please check");
						throw err;
					}
					mElectronDensity=MathFunction::SplineInterpExp(alts,logval,mAltGridKm);
				}
				break;
			case 5:
				{
					// Not really physical here, but interesting to check
					double alt0=0;
					mpParameter->GetValue(node,"//Model/AltT",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/DensT",dens0);
					double alt1=0;
					mpParameter->GetValue(node,"//Model/AltM",alt1);
					double dens1=0;
					mpParameter->GetValue(node,"//Model/DensM",dens1);
					if(dens0<0)
						dens0=1E-50;
					if(dens1<0)
						dens1=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					double Tmeso=0;
					mpParameter->GetValue(node,"//Model/Tmeso",Tmeso);
					double mixmassamu=0;
					mpParameter->GetValue(node,"//Model/MixMassAmu",mixmassamu);
					//unsigned pos=CloseInVector(alt0,mAltGridKm);
					mElectronDensity=MathFunction::DoubleExpProfile(mAltGridKm,alt0,alt1,dens0,dens1,Texo,Tmeso,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.,mixmassamu);

				}
				break;
			case 6:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					double SZA =0;
					mpParameter->GetValue(node,"//Model/SZA",SZA);
					mElectronDensity=MathFunction::ChapmanCos(mAltGridKm,alt0,dens0,Texo,SZA,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);

				}
				break;
			case 7:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					double C =0;
					mpParameter->GetValue(node,"//Model/C",C);
					mElectronDensity=MathFunction::ChapmanVar(mAltGridKm,alt0,dens0,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);

				}
				break;
			case 8:
				{
					double alt0=0;
					mpParameter->GetValue(node,"//Model/alt",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/dens",dens0);
					if(dens0<0)
						dens0=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					mElectronDensity=MathFunction::Epstein(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);

				}
				break;
			case 9:
				{
					// Not really physical here, but interesting to check
					double alt0=0;
					mpParameter->GetValue(node,"//Model/AltT",alt0);
					double dens0=0;
					mpParameter->GetValue(node,"//Model/DensT",dens0);
					double alt1=0;
					mpParameter->GetValue(node,"//Model/AltM",alt1);
					double dens1=0;
					mpParameter->GetValue(node,"//Model/DensM",dens1);
					if(dens0<0)
						dens0=1E-50;
					if(dens1<0)
						dens1=1E-50;
					double Texo=0;
					mpParameter->GetValue(node,"//Model/Texo",Texo);
					double Tmeso=0;
					mpParameter->GetValue(node,"//Model/Tmeso",Tmeso);
					double mixmassamu=0;
					mpParameter->GetValue(node,"//Model/MixMassAmu",mixmassamu);
					double C = 0;
					mpParameter->GetValue(node,"//Model/C",C);

					//unsigned pos=CloseInVector(alt0,mAltGridKm);
					mElectronDensity=MathFunction::ExpHyperbolaProfile(mAltGridKm,alt1,alt0,dens1,dens0,Tmeso,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.,mixmassamu);

				}
				break;
			default:
				Error err("ReadParameters: add_parametrized_species","Model error"," The model "+ntostr(model)+" is not a valid model for the neutral density profile");
		}

	}
	else
	{
		mElectronDensity=mpPlanet->ElectronDensity(mAltGridKm,type);
	}

	if(!mpParameter->Exists("/aero_main/atmosphere/iono/electron/density/multiplicator"))
	{
		mElecDensMult=1.;
	}else
	{
		mpParameter->GetValue("/aero_main/atmosphere/iono/electron/density/multiplicator",mElecDensMult);
		mElectronDensity*=mElecDensMult;
	}

        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-----+
         |                    |
         |Electron Temperature|
         |                    |
         +-------------------*/

	Log::mD<<"Init electron temperature"<<endl;
	mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/temperature/use_model","You have to define the electron temperature model");
	mpParameter->GetNKey("/aero_main/atmosphere/iono/electron/temperature/use_model","type",type);

	if(type==0)
	{
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/temperature/alt","You have to define the altitude, because you use the data model");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/electron/temperature/T","You have to define the electron temperature, because you use the data model");
		ublas::vector<double> alt,T;
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/electron/temperature/alt",alt);
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/electron/temperature/T",T);
		mElectronTemperature=MathFunction::IntLog(alt,T,mAltGridKm);
	}else
	{
		mElectronTemperature=mpPlanet->ElectronTemperature(mAltGridKm,type);
	}

	if(!mpParameter->Exists("/aero_main/atmosphere/iono/electron/temperature/multiplicator"))
	{
		mElecTempMult=1.;
	}else
	{
		mpParameter->GetValue("/aero_main/atmosphere/iono/electron/temperature/multiplicator",mElecTempMult);
		mElectronTemperature*=mElecTempMult;
	}


        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-+
         |                |
         |Ion temperature |
         |                |
         +---------------*/


	Log::mD<<"Init ion temperature"<<endl;
	mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/ions/temperature/use_model","You have to define the ion temperature model");
	mpParameter->GetNKey("/aero_main/atmosphere/iono/ions/temperature/use_model","type",type);

	if(type==0)
	{
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/ions/temperature/alt","You have to define the altitude, because you use the data model");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/iono/ions/temperature/T","You have to define the ion temperature, because you use the data model");
		ublas::vector<double> alt,T;
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/ions/temperature/alt",alt);
		mpParameter->Get1DArray("/aero_main/atmosphere/iono/ions/temperature/T",T);
		mIonTemperature=MathFunction::IntLog(alt,T,mAltGridKm);
	}else
	{
		mIonTemperature=mpPlanet->IonTemperature(mAltGridKm,type);
	}

	if(!mpParameter->Exists("/aero_main/atmosphere/iono/ions/temperature/multiplicator"))
	{
		mIonTempMult=1.;
	}else
	{
		mpParameter->GetValue("/aero_main/atmosphere/iono/ions/temperature/multiplicator",mIonTempMult);
		mIonTemperature*=mIonTempMult;
	}
}


void Atmo::PrintSpecies(std::string vFilenamePrefix)
{
	for(unsigned sp=0;sp<mpModelAtmosphere->mAtmoSpecies.size();++sp)
	{
		mpModelAtmosphere->mAtmoSpecies[sp]->PrintCrossSections(vFilenamePrefix);

	}
}


void Atmo::PrintBentIonosphere(std::string vFilenamePrefix)
{

	assert(mIsBendingActivated);
	string filename_prefix=trim(vFilenamePrefix);

	string densfname=filename_prefix+"_density.dat";
	string tempfname=filename_prefix+"_temp.dat";
        
	/**********************/
        /*                    */
        /* Temperature output */
        /*                    */
        /**********************/

	
	if(FileExists(tempfname))
	{
		Log::mD<<"Low level warning : We overwrite"<<tempfname<<endl;
	}
	ofstream of(tempfname.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Atmosphere temperatures "<<endl;
	of<<"# Temperature in K"<<endl;
	of<<"# Length in km (bent atmosphere)"<<endl;
	of<<"# Neutral T"<<endl;
	of<<"# Electron T"<<endl;
	of<<"# Ion T"<<endl;


	assert(mAtmoPath.mLengthKm.size()==mBentNeutralTemperature.size());
	assert(mAtmoPath.mLengthKm.size()==mBentElectronTemperature.size());
	assert(mAtmoPath.mLengthKm.size()==mBentIonTemperature.size());
//	assert(mAltGridKm.size()==mElectronTemperature.size());
//	assert(mAltGridKm.size()==mIonTemperature.size());

	for(unsigned i=0;i<mAtmoPath.mLengthKm.size();++i)
	{
		of<<mAtmoPath.mLengthKm[i]<<"\t"<<mBentNeutralTemperature[i]<<"\t"<<mBentElectronTemperature[i]<<"\t"<<mBentIonTemperature[i]<<endl;
	}
	of.close();

        /******************/
        /*                */
        /* Density output */
        /*                */
        /******************/


	if(FileExists(densfname))
	{
		Log::mD<<"Low level warning : We overwrite"<<densfname<<endl;
	}
	ofstream ofd(densfname.c_str());
	ofd.precision(9);
	ofd.setf(ios::scientific);
	ofd<<"# Atmosphere temperatures "<<endl;
	ofd<<"# Density in cm-3"<<endl;
	ofd<<"# Length in km (bent atmosphere)"<<endl;
	ofd<<"# Electron"<<endl;

	assert(mAtmoPath.mLengthKm.size()==mBentElectronDensity.size());
	for(unsigned i=0;i<mAtmoPath.mLengthKm.size();++i)
	{
		ofd<<mAtmoPath.mLengthKm[i]<<"\t"<<mBentElectronDensity[i]<<endl;
	}
	ofd.close();


}


void  Atmo::PrintIonosphere(std::string vFilenamePrefix)
{

	string filename_prefix=trim(vFilenamePrefix);

	string densfname=filename_prefix+"_density.dat";
	string tempfname=filename_prefix+"_temp.dat";
        
	/**********************/
        /*                    */
        /* Temperature output */
        /*                    */
        /**********************/

	
	if(FileExists(tempfname))
	{
		Log::mD<<"Low level warning : We overwrite"<<tempfname<<endl;
	}
	ofstream of(tempfname.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Atmosphere temperatures "<<endl;
	of<<"# Temperature in K"<<endl;
	of<<"# Altitude in km"<<endl;
	of<<"# Neutral T"<<endl;
	of<<"# Electron T"<<endl;
	of<<"# Ion T"<<endl;


	assert(mAltGridKm.size()==mNeutralTemperature.size());
	assert(mAltGridKm.size()==mElectronTemperature.size());
	assert(mAltGridKm.size()==mIonTemperature.size());

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t"<<mNeutralTemperature[i]<<"\t"<<mElectronTemperature[i]<<"\t"<<mIonTemperature[i]<<endl;
	}

	of<<Log::msMessageLog<<endl;
	of.close();

        /******************/
        /*                */
        /* Density output */
        /*                */
        /******************/


	if(FileExists(densfname))
	{
		Log::mD<<"Low level warning : We overwrite"<<densfname<<endl;
	}
	ofstream ofd(densfname.c_str());
	ofd.precision(9);
	ofd.setf(ios::scientific);
	ofd<<"# Atmosphere temperatures "<<endl;
	ofd<<"# Density in cm-3"<<endl;
	ofd<<"# Altitude in km"<<endl;
	ofd<<"# Electron (TEC "<<MathFunction::TrapzInt(mAltGridKm,mElectronDensity)*1E5<<" cm-2)"<<endl;

	assert(mAltGridKm.size()==mElectronDensity.size());
	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		ofd<<mAltGridKm[i]<<"\t"<<mElectronDensity[i]<<endl;
	}

	of<<Log::msMessageLog<<endl;
	ofd.close();


}

void Atmo::PrintNeutralAtmo(std::string vFile)
{
	mpModelAtmosphere->PrintNeutralAtmo(vFile);
}


void Atmo::PhotoIonize()
{
	mpPhotomodel=new Photoionization(mpParameter);

	// Init (std::vector< double > *photonGrid, std::vector< double > *photonGridmin, std::vector< double > *photonGridmax, std::vector< double > *mElecBotEeV, std::vector< double > *mElecCentEeV, std::vector< double > *mElecDdengeV, std::vector< double > *altitude, double o_UA, double o_R, double o_SZA)
	mpPhotomodel->Init(&mPhotonGrideV,&mPhotonGrideVmin,&mPhotonGrideVmax,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV,&mAltGridKm,mpPlanet->mUA,mpPlanet->mRKm,mpPlanet->mSZADegree);
	
	mpPhotoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);

	mbIsPhotoFluxDefined=true;

	mpPhotomodel->ComputePhotoionization((mpModelAtmosphere->mAtmoSpecies), mPhotoionizationResu, *mpPhotoFlux);

	delete mpPhotomodel;
}


void Atmo::PhotoAbsorb()
{
	mpPhotomodel=new Photoionization(mpParameter);

	// The default SZA is set at the planet one
	double sza = mpPlanet->mSZADegree;

	if(mpParameter->Exists("/aero_main/sun/transmission_sza"))
	{
		mpParameter->GetValue("/aero_main/sun/transmission_sza",sza);
	}
	mpParameter->ExistsOrDie("/aero_main/sun/transmission_file","You have to define transmission_file for writing your results");
	
	std::string outputfile = mpParameter->Elem("/aero_main/sun/transmission_file");
	
	mpPhotomodel->Init(&mPhotonGrideV,&mPhotonGrideVmin,&mPhotonGrideVmax,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV,&mAltGridKm,mpPlanet->mUA,mpPlanet->mRKm,sza);
	
//	mpPhotoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
//	mbIsPhotoFluxDefined=true;

	mpPhotomodel->ComputeTransmission((mpModelAtmosphere->mAtmoSpecies), outputfile);

	delete mpPhotomodel;
}


void Atmo::PhotoIonizeBent()
{
	assert(mComputeSZA.size()==mAtmoSpeciesSZAs.size());
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		for(unsigned j=0;j<mAtmoSpeciesSZAs[i].size();++j)
		{
			mAtmoSpeciesSZAs[i][j]->ClearProd();
		}
		Photoionization* photomodel=new Photoionization(mpParameter);
		photomodel->Init(&mPhotonGrideV,&mPhotonGrideVmin,&mPhotonGrideVmax,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV,&mAltGridKm,mpPlanet->mUA,mpPlanet->mRKm,mComputeSZA[i]);
		EFlux* tmpflux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
		std::deque<Specie*> resuphoto;
		photomodel->ComputePhotoionization(mAtmoSpeciesSZAs[i],resuphoto,*tmpflux);
		mPhotoFluxSZAs.push_back(tmpflux);
		mResuPhotoSZAs.push_back(resuphoto);
		delete photomodel;
	}

	mbIsPhotoFluxDefined=true;
	Log::mD<<"bending the photoelectron flux"<<endl;
	mpPhotoFlux = new EFlux(mPhotoFluxSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	Log::mD<<"bending the photoionization species result"<<endl;
	assert(mResuPhotoSZAs.size()==mComputeSZA.size());
	mPhotoionizationResu=SpecieUtils::BendAtmosphere(mResuPhotoSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	Log::mD<<"End photoionization bending"<<endl;
	mbIsPhotoFluxSZAs=true;
	
}


void Atmo::CosmoIonize()
{
	mpCosmoModel=new ProtoCosIonization(mpParameter,&mAltGridKm);

	mpCosmoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);

	mbIsCosmoFluxDefined=true;

	mpCosmoModel->ComputeCosmoionization((mpModelAtmosphere->mAtmoSpecies), mCosmoionizationResu, *mpCosmoFlux);
	delete mpCosmoModel;

}

void Atmo::CosmoIonizeBent()
{

	mpCosmoModel=new ProtoCosIonization(mpParameter,&mAltGridKm);
	EFlux* cosmoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		for(unsigned j=0;j<mAtmoSpeciesSZAs[i].size();++j)
		{
			mAtmoSpeciesSZAs[i][j]->ClearProd();
		}
	}
	std::deque<Specie*> cosmoioni;
	mpCosmoModel->ComputeCosmoionization(mAtmoSpeciesSZAs[0], cosmoioni, *cosmoFlux);
	delete mpCosmoModel;
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		mCosmoFluxSZAs.push_back(cosmoFlux);
		mResuCosmoSZAs.push_back(cosmoioni);
	}
	mpCosmoFlux = new EFlux(mCosmoFluxSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	mCosmoionizationResu=SpecieUtils::BendAtmosphere(mResuCosmoSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	mbIsCosmoFluxDefined=true;
	mbIsCosmoFluxSZAs=true;
}

void Atmo::ProtonIonize()
{
	Log::mI<<"We start the proton model in the non-bent case"<<endl;
	mpProtonModel=new ProtonHydrogenTransport(mpParameter,&mAltGridKm, mpProtonGaussianAngle, &mdB_B);
	mpProtonEFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	
	boost::shared_ptr<PHFlux> protonflux(new PHFlux(mpParameter,mpProtonGaussianAngle,&mProtonGrideV,&mProtonWidthGrideV,"proton"));
	boost::shared_ptr<PHFlux> hydrogenflux(new PHFlux(mpParameter,mpProtonGaussianAngle,&mProtonGrideV,&mProtonWidthGrideV,"hydrogen"));
	protonflux->InitVoid(mAltGridKm.size()); // The proton flux will be filled by the transport function
	hydrogenflux->InitVoid(mAltGridKm.size()); // The hydrogen flux will be filled by the transport function
	
	protonflux->ReadPrecipitation(mAltGridKm, mpModelAtmosphere->mAtmoSpecies);
	hydrogenflux->ReadPrecipitation(mAltGridKm, mpModelAtmosphere->mAtmoSpecies);
	mbIsProtonEFluxDefined=true;

	mpProtonModel->ComputeProtonImpact((mpModelAtmosphere->mAtmoSpecies), protonflux, hydrogenflux, mProtoionizationResu, *mpProtonEFlux);
	mbIsProtonModel = true;
	//delete mpProtonModel;
	//delete protonflux;
	//delete hydrogenflux;
}
void Atmo::ProtonIonizeBent()
{
	Log::mI<<"We start the proton model in the bent case"<<endl;
	mpProtonModel=new ProtonHydrogenTransport(mpParameter,&mAltGridKm,mpProtonGaussianAngle, &mdB_B);
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		for(unsigned j=0;j<mAtmoSpeciesSZAs[i].size();++j)
		{
			mAtmoSpeciesSZAs[i][j]->ClearProd();
		}
	}
	std::deque<Specie*> protioni;
	//EFlux* testflux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	boost::shared_ptr<PHFlux> protonflux(new PHFlux(mpParameter,mpProtonGaussianAngle,&mProtonGrideV,&mProtonWidthGrideV,"proton"));
	boost::shared_ptr<PHFlux> hydrogenflux(new PHFlux(mpParameter,mpProtonGaussianAngle,&mProtonGrideV,&mProtonWidthGrideV,"hydrogen"));
	protonflux->InitVoid(mAltGridKm.size()); // The proton flux will be filled by the transport function
	hydrogenflux->InitVoid(mAltGridKm.size()); // The hydrogen flux will be filled by the transport function
	
	EFlux* protonEFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	
	protonflux->ReadPrecipitation(mAtmoPath.mLengthKm, mAtmoBentSpecies);
	hydrogenflux->ReadPrecipitation(mAtmoPath.mLengthKm, mAtmoBentSpecies);
	mpProtonModel->ComputeProtonImpact(mAtmoSpeciesSZAs[0],protonflux, hydrogenflux, protioni, *protonEFlux);
	delete mpProtonModel;
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		mProtonEFluxSZAs.push_back(protonEFlux);
		mResuProtonSZAs.push_back(protioni);
	}
	mpProtonEFlux = new EFlux(mCosmoFluxSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	mProtoionizationResu=SpecieUtils::BendAtmosphere(mResuProtonSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	mbIsProtonEFluxDefined=true;
	mbIsProtonFluxSZAs=true;
	//delete testflux;
	delete protonEFlux;
}

void Atmo::ElectroIonizeBent()
{
	Log::mI<<"We add the different electron fluxes before the computation of the bent atmosphere"<<endl;
	// We load the electron precipitation
	mpElecPrecipFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	mpElecPrecipFlux->InitVoid(mAtmoPath.mAltitudeKm.size());
	mbIsElecPrecipFluxDefined=true;
	Log::mD<<"We read the precipitation!"<<endl;
	mpElecPrecipFlux->ReadPrecipitation(mAtmoPath.mLengthKm,mAtmoBentSpecies);

	mpTotalFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	*mpTotalFlux=*mpElecPrecipFlux;
	if(mbIsPhotoFluxDefined)
	{
		Log::mL<<"Total photoelectron flux (bent) "<<mpPhotoFlux->FluxEnergy(mAtmoPath.mLengthKm)<<endl;
		(*mpTotalFlux)+=*mpPhotoFlux;
	}

	if(mbIsCosmoFluxDefined)
	{
		Log::mL<<"Total cosmoelectron flux (bent) "<<mpCosmoFlux->FluxEnergy(mAtmoPath.mLengthKm)<<endl;
		(*mpTotalFlux)+=*mpCosmoFlux;
	}
	
	if(mbIsProtonEFluxDefined)
	{
		Log::mL<<"Total protoelectron flux (bent) "<<mpProtonEFlux->FluxEnergy(mAtmoPath.mLengthKm)<<endl;
		(*mpTotalFlux)+=*mpProtonEFlux;
	}
	Log::mD<<"We load the bent electron impact object"<<endl;
	Log::mD<<" Alt "<<mAtmoPath.mAltitudeKm<<endl;
	Log::mD<<" ET : "<<mBentElectronTemperature<<endl;
	Log::mD<<" ED : "<<mBentElectronDensity<<endl;
	mpElecModel=new ElectronImpactIonization(mpParameter,&mBentElectronDensity,&mBentElectronTemperature,&mAtmoPath.mLengthKm);
	Log::mL<<"We lauch the electron impact computation"<<endl;
	mpElecModel->ComputeElectronImpact(mAtmoBentSpecies,*mpTotalFlux,mElectronImpactResu);
	
	mbIsElecModelDefined=true;
}


void Atmo::ElectroIonize()
{
	Log::mI<<"We add the different electron fluxes before the computation"<<endl;
	// We load the electron precipitation
	mpElecPrecipFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
	mpElecPrecipFlux->InitVoid(mAltGridKm.size());
	mbIsElecPrecipFluxDefined=true;
	Log::mD<<"We read the precipitation!"<<endl;
	mpElecPrecipFlux->ReadPrecipitation(mAltGridKm,(mpModelAtmosphere->mAtmoSpecies));
	

//	mpElecPrecipFlux->PrintPrecip("mpelecprecipflux");

	Log::mD<<"We add the different electron fluxes"<<endl;
	// We add the electron precipitation and the photoelectron flux
	// please note that if you want to add another process, it is here ;-)
	


	mpTotalFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
//	EFlux electrons=*mpElecPrecipFlux;
	*mpTotalFlux=*mpElecPrecipFlux;//+*mpPhotoFlux;
	if(mbIsPhotoFluxDefined)
	{
		Log::mL<<"Total photoelectron flux "<<mpPhotoFlux->FluxEnergy(mAltGridKm)<<endl;
//		mpPhotoFlux->PrintElecProfile("elecprofile", mAltGridKm);
		(*mpTotalFlux)+=*mpPhotoFlux;
	//	electrons+=*mpPhotoFlux;
	}

	if(mbIsCosmoFluxDefined)
	{
		Log::mD<<"Total cosmoelectron flux "<<mpCosmoFlux->FluxEnergy(mAltGridKm)<<endl;
		(*mpTotalFlux)+=*mpCosmoFlux;
	}
	
	if(mbIsProtonEFluxDefined)
	{
		Log::mL<<"Total protoelectron flux "<<mpProtonEFlux->FluxEnergy(mAltGridKm)<<endl;
		(*mpTotalFlux)+=*mpProtonEFlux;
	}
//	electrons.PrintPrecip("electronsprecipflux");
//	electrons.PrintEnergyFlux("out_theenergyflux",mAltGridKm);
//	EFlux electrons=*mpPhotoFlux;
//	EFlux electrons=*mpElecPrecipFlux;
/*	cout<<"Address meleccente"<<&mElecCentEeV<<endl;
	cout<<"Address electrons melec..."<<endl;
	electrons.printelecadress();
	cout<<"Autre :"<<endl;
	*mpElecPrecipFlux->printelecadress();
*/
	Log::mD<<"We load the electron impact object"<<endl;
	mpElecModel=new ElectronImpactIonization(mpParameter,&mElectronDensity,&mElectronTemperature,&mAltGridKm);
	Log::mD<<"We lauch the electron impact computation"<<endl;
	mpElecModel->ComputeElectronImpact((mpModelAtmosphere->mAtmoSpecies),*mpTotalFlux,mElectronImpactResu);
//	mpElecModel->ComputeElectronImpact((mpModelAtmosphere->mAtmoSpecies),electrons,mElectronImpactResu);
	
	mbIsElecModelDefined=true;
//	delete mpElecModel;
}


void Atmo::BendMyAtmosphere()
{
	Log::mS<<"Initialize atmosphere lists for bending the atmosphere"<<endl;
	mAtmoSpeciesSZAs.resize(0);
	mNeutralTempSZAs.resize(0);
	mElecTempSZAs.resize(0);
	mIonTempSZAs.resize(0);
	mElectronDensitySZAs.resize(0);
	for(unsigned i=0;i<mComputeSZA.size();++i)
	{
		mAtmoSpeciesSZAs.push_back(SpecieUtils::CopyAtmo(mpModelAtmosphere->mAtmoSpecies));
		mNeutralTempSZAs.push_back(&mNeutralTemperature);
		mElecTempSZAs.push_back(&mElectronTemperature);
		mIonTempSZAs.push_back(&mIonTemperature);
		mElectronDensitySZAs.push_back(&mElectronDensity);
	}

	Log::mL<<"Bend the atmosphere"<<endl;
	mAtmoBentSpecies=SpecieUtils::BendAtmosphere(mAtmoSpeciesSZAs,mAltGridKm,mComputeSZA,mAtmoPath);
	
	for(unsigned i=0;i<mAtmoBentSpecies.size();++i)
	{
		Log::mI<<"Species : "<<mAtmoBentSpecies[i]->mName<<endl;
		Log::mI<<"Density : "<<mAtmoBentSpecies[i]->mTotDensitycm_3<<endl;
		Log::mI<<"ScaleH : "<<mAtmoBentSpecies[i]->mScaleHcm<<endl;
		Log::mI<<"ColDens : "<<mAtmoBentSpecies[i]->mColDenscm_2<<endl;
	}
	
	
	Log::mD<<"Bend the neutral temperature"<<endl;
	mBentNeutralTemperature=MathFunction::IntLogPath(mAltGridKm,mComputeSZA,mNeutralTempSZAs,mAtmoPath.mAltitudeKm,mAtmoPath.mSZADegree);
	Log::mD<<"Bend the ion temperature"<<endl;
	mBentIonTemperature=MathFunction::IntLogPath(mAltGridKm,mComputeSZA,mIonTempSZAs,mAtmoPath.mAltitudeKm,mAtmoPath.mSZADegree);
	Log::mD<<"Bend the electron temperature"<<endl;
	mBentElectronTemperature=MathFunction::IntLogPath(mAltGridKm,mComputeSZA,mElecTempSZAs,mAtmoPath.mAltitudeKm,mAtmoPath.mSZADegree);
	Log::mD<<"Bend the electron density"<<endl;
	mBentElectronDensity=MathFunction::IntLogPath(mAltGridKm,mComputeSZA,mElectronDensitySZAs,mAtmoPath.mAltitudeKm,mAtmoPath.mSZADegree);

	MinValue(mBentElectronDensity,1E-42);

}






void Atmo::ComputeBent()
{
	Log::mS<<"BENT ATMOSPHERE"<<endl;
	mIsBendingActivated=true;
	Log::AddVersionInfo("# BENT ATMOSPHERE");
//	std::ublas<double> computeSZA;
	mComputeSZA.push_back(mpPlanet->mSZADegree);
	if(mpParameter->Exists("/aero_main/atmosphere/bent/compute_SZA"))
	{
		ublas::vector<double> computeSZAu;
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/compute_SZA",computeSZAu);
		
		
		for(unsigned i=0;i<computeSZAu.size();++i)
		{
			mComputeSZA.push_back(computeSZAu[i]);
		}
		sort(mComputeSZA.begin(),mComputeSZA.end());
		
	}
	Log::mD<<"SZA array computed"<<endl;

	bool isSZA=mpParameter->Exists("/aero_main/atmosphere/bent/path_SZA");
	bool isXYZ=mpParameter->Exists("/aero_main/atmosphere/bent/path_XYZ");
	if(isSZA==isXYZ)
	{ // none of them are defined, or the two are defined... == equivalent to !Xor
		Error err("Atmo::ComputeBent","Bad definition in you XML parameter file"," You define both/aero_main/atmosphere/bent/path_SZA and aero_main/atmosphere/bent/path_SZA it is not possible, one should be deleted");
		throw err;
	}
	Log::mD<<"We search for the path"<<endl;
	if(isSZA)
	{
		Log::mD<<" SZA path"<<endl;
		mpParameter->ExistsOrDie("/aero_main/atmosphere/bent/path_SZA/altitude","You should define the altitude parameter");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/bent/path_SZA/length","You should define the length parameter");
		ublas::vector<double> alts,lengths;
		
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_SZA/altitude",alts);
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_SZA/length",lengths);

		// We define the Path!
		if(mpParameter->Exists("/aero_main/atmosphere/bent/path_SZA/SZA"))
		{
			ublas::vector<double> SZAs;
			mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_SZA/SZA",SZAs);
			mAtmoPath.InitAltLenSZA(alts,lengths,SZAs);
		}else
		{
			mAtmoPath.InitAltLen(alts,lengths,mpPlanet->mSZADegree);
		}
	}
	if(isXYZ)
	{
		Log::mD<<" XYZ path"<<endl;
		mpParameter->ExistsOrDie("/aero_main/atmosphere/bent/path_XYZ/X","You should define the X parameter");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/bent/path_XYZ/Y","You should define the Y parameter");
		mpParameter->ExistsOrDie("/aero_main/atmosphere/bent/path_XYZ/Z","You should define the Z parameter");
		ublas::vector<double> x,y,z;
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_XYZ/X",x);
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_XYZ/Y",y);
		mpParameter->Get1DArray("/aero_main/atmosphere/bent/path_XYZ/Z",z);
		mAtmoPath.CartesianToPath(x,y,z,mpPlanet->mRKm);
	}

	// Here, we should create some code to have an atmosphere - ionosphere dependant on the SZA.
	

	// bent neutral atmosphere
	// bent electrons temp and density
	// bent ion temp
	Log::mD<<"We bend the neutral atmosphere"<<endl;
	BendMyAtmosphere();



	// It should implements all the arrays of the atmosphere for the other computations
	// It should use the SZA in computeSZA to compute the atmospheric parameters.


	if(mpParameter->Exists("/aero_main/sun/use_sun"))
	{
		Log::AddVersionInfo("# Photoionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mI<<"Start photoionization (and excitation...)"<<endl;
		PhotoIonizeBent();
	}
	Log::mD<<"Fin photoionization"<<endl;

	if(mpParameter->Exists("/aero_main/proton_cosmic"))
	{
		Log::AddVersionInfo("# Cosmic ray ionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start cosmic ray ionization"<<endl;
		CosmoIonizeBent();
	}

	if(mpParameter->Exists("/aero_main/proton/use_proton"))
	{
		Log::AddVersionInfo("# Proton/Hydrogen ionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start proton/hydrogen ionization"<<endl;
		ProtonIonizeBent();
	}
	if(mpParameter->Exists("/aero_main/electron/use_electron"))
	{
		Log::AddVersionInfo("# Electron impact ionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start electron impact ionization (and excitation...)"<<endl;
		ElectroIonizeBent();
	}
	
	// These resu contains now the result on the BENT ATMOSPHERE
	// then we can add them, and use them in the big result...
	//mTotalResu=SpecieUtils::MergeResult(mPhotoionizationResu,mElectronImpactResu);
	//mTotalResu=SpecieUtils::MergeResult(mCosmoionizationResu,mTotalResu);

	std::deque<Specie* > totalResutmp2=SpecieUtils::MergeResult(mPhotoionizationResu,mElectronImpactResu);
	std::deque<Specie* > totalResutmp=SpecieUtils::MergeResult(totalResutmp2,mProtoionizationResu);
	mTotalResu=SpecieUtils::MergeResult(mCosmoionizationResu,totalResutmp); // the overload is problematic for deleting the species... maybe I should use a smart pointer now
	for(size_t i=0;i<totalResutmp.size();++i)
		delete totalResutmp[i];
	for(size_t i=0;i<totalResutmp2.size();++i)
		delete totalResutmp2[i];
}

void Atmo::ComputePhotodissociation(std::string suffix)
{
	mpParameter->ExistsOrDie("/aero_main/sun/use_sun", "You need to define the SUN");
	Log::AddVersionInfo("# Photoionization used");
	mpModelAtmosphere->ClearProductions();
	Log::mS<<"Start photoionization (and excitation...)"<<endl;
	mpPhotomodel=new Photoionization(mpParameter);
	mpPhotomodel->Init(&mPhotonGrideV,&mPhotonGrideVmin,&mPhotonGrideVmax,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV,&mAltGridKm,mpPlanet->mUA,mpPlanet->mRKm,mpPlanet->mSZADegree);
	mpPhotomodel->ComputePhotodissociationBranching((mpModelAtmosphere->mAtmoSpecies), suffix);
	delete mpPhotomodel;
}

void Atmo::Compute()
{
	if(mpParameter->Exists("/aero_main/atmosphere/use_bent_atmosphere"))
	{// Here, we work with the bent atmosphere
		ComputeBent();
		return;
	}
	
	if(mpParameter->Exists("/aero_main/sun/use_sun"))
	{
		Log::AddVersionInfo("# Photoionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start photoionization (and excitation...)"<<endl;
		PhotoIonize();
	}
	if(mpParameter->Exists("/aero_main/sun/use_transmission"))
	{
		Log::AddVersionInfo("# PhotoTransmission used");
		Log::mS<<"Start phototransmission"<<endl;
		PhotoAbsorb();
	}

	Log::mD<<"Fin photoionization"<<endl;
	
	if(mpParameter->Exists("/aero_main/proton_cosmic"))
	{
		Log::AddVersionInfo("# Cosmic ray ionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start cosmic ray ionization"<<endl;
		CosmoIonize();

	}
	if(mpParameter->Exists("/aero_main/proton/use_proton"))
	{
		Log::AddVersionInfo("# Proton/Hydrogen ionization used");
		mpModelAtmosphere->ClearProductions();
		Log::mS<<"Start proton/hydrogen ionization"<<endl;
		ProtonIonize();

	}

	// Very important, we want a new computation, so we need to clear the states!
	mpModelAtmosphere->ClearProductions();
	if(mpParameter->Exists("/aero_main/electron/use_electron"))
	{
		Log::AddVersionInfo("# Electron impact ionization used");
		Log::mS<<"Start electron impact ionization (and excitation...)"<<endl;
		ElectroIonize();
	}



	std::deque<Specie* > totalResutmp2 = SpecieUtils::MergeResult(mPhotoionizationResu,mElectronImpactResu);
	std::deque<Specie* > totalResutmp = SpecieUtils::MergeResult(totalResutmp2,mProtoionizationResu);
	mTotalResu=SpecieUtils::MergeResult(mCosmoionizationResu,totalResutmp2); // the overload is problematic for deleting the species... maybe I should use a smart pointer now
	for(size_t i=0;i<totalResutmp.size();++i)
		delete totalResutmp[i];
	for(size_t i=0;i<totalResutmp2.size();++i)
		delete totalResutmp2[i];

}
void Atmo::ReadProductions()
{

	mPhotoionizationResu.erase(mPhotoionizationResu.begin(),mPhotoionizationResu.end());
	mCosmoionizationResu.erase(mCosmoionizationResu.begin(),mCosmoionizationResu.end());
	mProtoionizationResu.erase(mProtoionizationResu.begin(),mProtoionizationResu.end());
	mElectronImpactResu.erase(mElectronImpactResu.begin(),mElectronImpactResu.end());
	mTotalResu.erase(mTotalResu.begin(),mTotalResu.end());
	mpModelAtmosphere->ClearProductions();


	// we read the production in the xml file
	vector<TiXmlNode*> prods = mpParameter->GetNodes("/aero_main/atmosphere/production/specie");
	
	//size_t s = prods.size();
	//We check the number of processes to read
	//
	//And we modify the 1st species (in mpModelAtmosphere->mAtmoSpecies) accordingly
	//this species will embeed all the productions

	if(mpModelAtmosphere->mAtmoSpecies.size()==0)
		return;
	Specie* sp = mpModelAtmosphere->mAtmoSpecies[0];

		//std::deque<std::string> mProcessNames;
	std::deque<std::string> processNames;
		//std::deque< std::deque< SpecieId > > mCreatedSpecies;
	std::deque< std::deque< SpecieId > > createdSpecies;
		//std::deque< std::deque<double> > mMultiplicativeFactor;
	std::deque< std::deque<double> > mmultiplicativeFactor;
		//std::deque< ublas::vector<double> > mSpeciesProductioncm_3s_1;
	std::deque< ublas::vector<double> > mspeciesProductioncm_3s_1;

	// Loop on the different species
	for(vector<TiXmlNode*>::iterator it = prods.begin();it!=prods.end();++it)
	{
		Log::mL<<"hello"<<endl;
		std::deque< SpecieId > tmpcspecies;
		std::deque<double> tmpmultfact;
		//Read the species
		string spname=mpParameter->GetKey(*it,"/","name");
		string spstate=mpParameter->GetKey(*it,"/","state");
		SpecieId tmpsp(spname,spstate);
		tmpcspecies.push_back(tmpsp);
		tmpmultfact.push_back(1.);
		ublas::vector<double> prod, alt, uncert;

		//Read the altitudes corresponding to the production
		mpParameter->Get1DArray(*it,"//alt",alt);
		//Read the production
		mpParameter->Get1DArray(*it,"//prod",prod);
		assert(alt.size()==prod.size());

		if(mpParameter->Exists(*it,"//use_mc"))
		{
			//Read the uncertainties
			mpParameter->Get1DArray(*it,"//uncert",uncert);
			//perform the monte carlo modification if necessary
			assert(alt.size()==uncert.size());
			if(mpParameter->Exists(*it,"//use_lognormal"))
			{
				for(size_t i = 0;i<alt.size();++i)
				{
					if(prod[i]<=0)
					{
						prod[i]=0;
					}else{
						double lp = log(prod[i]); 
						mpParameter->ApplyMCP(lp,-1,uncert[i]/prod[i]);
						prod[i] = exp(lp);
					}
				}
			}else{

				for(size_t i = 0;i<alt.size();++i)
				{
					mpParameter->ApplyMCP(prod[i],-1,uncert[i]);
					if(prod[i]<0)
						prod[i]=0;
				}
			}
		}

		if(mpParameter->Exists(*it,"//use_mc_fact"))
		{
			//Read the uncertainties
			mpParameter->Get1DArray(*it,"//uncert",uncert);
			//perform the monte carlo modification if necessary
			assert(alt.size()==uncert.size());
			double factor = 0;
			mpParameter->ApplyMCP(factor,-1,1.);// We use a normal law, centered at 0, with sigma = 1
			if(mpParameter->Exists(*it,"//use_lognormal"))
			{
				for(size_t i = 0;i<alt.size();++i)
				{
					if(prod[i]<=0)
					{
						prod[i]=0;
					}else{
						double sp = uncert[i] / prod[i];
						prod[i] *= exp(factor * sp);
					}
				}
			}else{
				for(size_t i = 0;i<alt.size();++i)
				{
					//mpParameter->ApplyMC(prod[i],-1,uncert[i]);
					prod[i]+= uncert[i] * factor;
					if(prod[i]<0)
						prod[i]=0;
				}
			}
		}
		//Interpolate on the current grid
		//
		ublas::vector<double> newprod = MathFunction::IntLin(alt,prod,mAltGridKm);

		createdSpecies.push_back(tmpcspecies);
		mmultiplicativeFactor.push_back(tmpmultfact);
		mspeciesProductioncm_3s_1.push_back(newprod);
	}
	//Write in the mAtmoSpecies  (we just put in the first one....

	sp->mProcessNames = processNames;
	sp->mCreatedSpecies = createdSpecies;
	sp->mMultiplicativeFactor = mmultiplicativeFactor;
	sp->mSpeciesProductioncm_3s_1 = mspeciesProductioncm_3s_1;


	// The productions are now in mpModelAtmosphere->mAtmoSpecies

	SpecieUtils::SpeciesToResu(mpModelAtmosphere->mAtmoSpecies,mTotalResu);
	mpModelAtmosphere->ClearProductions();
}

void Atmo::ProceedEmissions(std::string suffix)
{
	Log::mD<<"Search for emissions"<<endl;
	if(mpParameter->Exists("/aero_main/chem/use_chem"))
	{
		mbIsChem=true;
		Log::mD<<"Proceed emissions"<<endl;
		InitChem();
		//InitEmission();
		mEmit.reset(new Emission(mpParameter,mChem,mpPlanet,mAltGridKm));

		mEmit->ComputeEmissions(suffix);
		Log::mD<<"End emissions"<<endl;
	}
	
}


void Atmo::ProceedOutputs(std::string suffix)
{
	if(mbIsChem && mpParameter->Exists("/aero_main/chem/chem_atmo"))
	{
		string file=mpParameter->Elem("/aero_main/chem/chem_atmo");

		if(suffix!="")
		{
			file+=suffix;
		}
		mChem->PrintChemDensModels(file);

	}



	if(mpParameter->Exists("/aero_main/output/neutral_atmo"))
	{
		string file=mpParameter->Elem("/aero_main/output/neutral_atmo");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mIsBendingActivated)
		{
			mpModelAtmosphere->PrintBentNeutralAtmo(mAtmoBentSpecies,mAtmoPath.mLengthKm,file);
		}else
		{
			mpModelAtmosphere->PrintNeutralAtmo(file);
		}
	}

	if(mpParameter->Exists("/aero_main/output/neutral_atmo_colden"))
	{
		string file=mpParameter->Elem("/aero_main/output/neutral_atmo_colden");

		if(suffix!="")
		{
			file+=suffix;
		}

		if(mIsBendingActivated)
		{
			mpModelAtmosphere->PrintBentNeutralColdens(mAtmoBentSpecies,mAtmoPath.mLengthKm,file);
		}else
		{
			mpModelAtmosphere->PrintNeutralColdens(file);
		}
	}

	if(mpParameter->Exists("/aero_main/output/ionosphere"))
	{
		string file=mpParameter->Elem("/aero_main/output/ionosphere");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mIsBendingActivated)
		{
			PrintBentIonosphere(file);
		}else
		{
			PrintIonosphere(file);
		}
	}
	if(mpParameter->Exists("/aero_main/output/cross_sections"))
	{
		string file=mpParameter->Elem("/aero_main/output/cross_sections");
	//	Log::SetPriority(Log::DEBUGG);
	//	Log::mL<<"We print the cross sections"<<endl;
		if(suffix!="")
		{
			file+=suffix;
		}

		PrintSpecies(file);
	}


	if(mpParameter->Exists("/aero_main/output/PWOM") && mbIsElecModelDefined)
	{
		string file=mpParameter->Elem("/aero_main/output/PWOM");
		if(suffix!="")
		{
			file+=suffix;
		}
		mpElecModel->PrintPWOM(file);

	}

	if(mpParameter->Exists("/aero_main/output/fluxes") && mbIsElecModelDefined)
	{
		vector<TiXmlNode*> fluxes=mpParameter->GetNodes("/aero_main/output/fluxes");
		vector<TiXmlNode*>::iterator sit;
		for(sit=fluxes.begin();sit!=fluxes.end();++sit)
		{
			unsigned key=0;
			mpParameter->GetNKey(*sit,"/","option",key);
			mpParameter->ExistsOrDie(*sit,"//fluxes_filename","you must give a filename");
			string file="";
			file=mpParameter->Elem(*sit,"//fluxes_filename");

			if(suffix!="")
			{
				file+=suffix;
			}
			if(!mIsBendingActivated)
			{

				ublas::vector<double> mesaltitudes;
				if( mpParameter->Exists(*sit,"//altitudes"))//,"You should give the altitudes where you want to see the fluxes");
				{
					mpParameter->Get1DArray(*sit,"//altitudes",mesaltitudes);
				}
				mpElecModel->PrintFluxes(file,mesaltitudes,key);
			}else
			{
				ublas::vector<double> mesaltitudes,meslengths;
				if(mpParameter->Exists(*sit,"//altitudes"))
				{
					Log::mD<<"Loading altitude for the flux printing"<<endl;
					Log::mD<<"Be careful, due to the possible non unicity of the altitudes in the case of a bent atmosphere, you get only the first check of this altitude in the Path"<<endl;
					mpParameter->Get1DArray(*sit,"//altitudes",mesaltitudes);
				}
				if(mpParameter->Exists(*sit,"//lengths"))
				{
					Log::mD<<"Loading lengths for the flux printing"<<endl;
					mpParameter->Get1DArray(*sit,"//lengths",meslengths);
				}
				mpElecModel->PrintBentFluxes(file,mAtmoPath,mesaltitudes,meslengths,key);
			}
		}
	}

	if(mpParameter->Exists("/aero_main/output/proton_fluxes") && mbIsProtonModel)
	{
		vector<TiXmlNode*> fluxes=mpParameter->GetNodes("/aero_main/output/proton_fluxes");
		vector<TiXmlNode*>::iterator sit;
		for(sit=fluxes.begin();sit!=fluxes.end();++sit)
		{
			unsigned key=0;
			mpParameter->GetNKey(*sit,"/","option",key);
			mpParameter->ExistsOrDie(*sit,"//fluxes_filename","you must give a filename");
			string file="";
			file=mpParameter->Elem(*sit,"//fluxes_filename");
			bool is_proton = true;
			if(mpParameter->Exists(*sit,"//H"))
			{
				is_proton = false;
			}


			if(suffix!="")
			{
				file+=suffix;
			}
			if(!mIsBendingActivated)
			{

				ublas::vector<double> mesaltitudes;
				if( mpParameter->Exists(*sit,"//altitudes"))//,"You should give the altitudes where you want to see the fluxes");
				{
					mpParameter->Get1DArray(*sit,"//altitudes",mesaltitudes);
				}
				mpProtonModel->PrintFluxes(file,mesaltitudes,key, is_proton);
			}else
			{
				ublas::vector<double> mesaltitudes,meslengths;
				if(mpParameter->Exists(*sit,"//altitudes"))
				{
					Log::mD<<"Loading altitude for the flux printing"<<endl;
					Log::mD<<"Be careful, due to the possible non unicity of the altitudes in the case of a bent atmosphere, you get only the first check of this altitude in the Path"<<endl;
					mpParameter->Get1DArray(*sit,"//altitudes",mesaltitudes);
				}
				if(mpParameter->Exists(*sit,"//lengths"))
				{
					Log::mD<<"Loading lengths for the flux printing"<<endl;
					mpParameter->Get1DArray(*sit,"//lengths",meslengths);
				}
				//	mpElecModel->PrintBentFluxes(file,mAtmoPath,mesaltitudes,meslengths,key);
				Error err("Atmo::ProceedOutputs() / proton fluxes","Bent atmo used","The bent atmosphere is not implemented for the proton computation  yet");
				throw err;
			}
		}
	}
	Log::mI<<"Output production"<<endl;
	//	if(mpParameter->Exists("/output/production_selected_species"))
	if(mpParameter->Exists("/aero_main/output/selection"))
	{
		vector<TiXmlNode*> selections=mpParameter->GetNodes("/aero_main/output/selection");
		vector<TiXmlNode*>::iterator sit;
		for(sit=selections.begin();sit!=selections.end();++sit)
		{

			if(!mpParameter->Exists(*sit,"//production_selected_species"))
			{
				Log::mW<<"There is a void /output/selection !!!!"<<endl;
			}

			string file=mpParameter->Elem(*sit,"//production_selected_species");

			if(suffix!="")
			{
				file+=suffix;
			}

			deque< SpecieId > mes_id;

			if(mpParameter->Exists(*sit,"//select_species"))
			{
				vector<TiXmlNode*> species=mpParameter->GetNodes(*sit,"//select_species/Specie");
				vector<TiXmlNode*>::iterator jt;
				//			Log::mL<<"Number of nodes :"<<species.size()<<endl;
				for(jt=species.begin();jt!=species.end();++jt)
				{

					string spname=mpParameter->GetKey(*jt,"/","name");
					string spstate=mpParameter->GetKey(*jt,"/","state");
					SpecieId tmpsp(spname,spstate);
					mes_id.push_back(tmpsp);
				}



			}else
			{
				SpecieId CO2p("CO2+","X");
				mes_id.push_back(CO2p);
				SpecieId COp("CO+","X");
				mes_id.push_back(COp);
				SpecieId Op("O+","X");
				mes_id.push_back(Op);
				SpecieId Cp("C+","X");
				mes_id.push_back(Cp);
				SpecieId Cpp("C++","X");
				mes_id.push_back(Cpp);
				SpecieId CO2pp("CO2++","X");
				mes_id.push_back(CO2pp);
				// Try an unknown specie...
				SpecieId CO4mm("CO4--","X");
				mes_id.push_back(CO4mm);
			}

			if(!mIsBendingActivated)
			{
				string photfile=file+"-photon.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mPhotoionizationResu,mAltGridKm,photfile);
				string cosmofile=file+"-cosmo.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mCosmoionizationResu,mAltGridKm,cosmofile);
				string protofile=file+"-proton.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mProtoionizationResu,mAltGridKm,protofile);
				string elecfile=file+"-electron.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mElectronImpactResu,mAltGridKm,elecfile);

				SpecieUtils::SelectedPrintProduction(mes_id,mTotalResu,mAltGridKm,file);
			}else
			{
				string photfile=file+"-photon.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mPhotoionizationResu,mAtmoPath.mLengthKm,photfile,true);
				string cosmofile=file+"-cosmo.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mCosmoionizationResu,mAtmoPath.mLengthKm,cosmofile,true);
				string protofile=file+"-proton.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mProtoionizationResu,mAtmoPath.mLengthKm,protofile,true);
				string elecfile=file+"-electron.dat";
				SpecieUtils::SelectedPrintProduction(mes_id,mElectronImpactResu,mAtmoPath.mLengthKm,elecfile,true);

				SpecieUtils::SelectedPrintProduction(mes_id,mTotalResu,mAtmoPath.mLengthKm,file,true);
			}
		}
	}

	if(mpParameter->Exists("/aero_main/output/electron_production"))
	{
		string file=mpParameter->Elem("/aero_main/output/electron_production");

		if(suffix!="")
		{
			file+=suffix;
		}
			if(!mIsBendingActivated)
			{
				SpecieUtils::PrintElectronProduction(mpModelAtmosphere->mAtmoSpecies,mAltGridKm,file);
			}else
			{
				//Log::SetPriority(Log::WARNING);
				Log::mW<<"The electron production printing does not work with a bent atmosphere"<<endl;
			}
	}


	if(mpParameter->Exists("/aero_main/output/production"))
	{
		string file=mpParameter->Elem("/aero_main/output/production");

		if(suffix!="")
		{
			file+=suffix;
		}
		bool excitations=false;
		if(mpParameter->NbParams("/aero_main/output/production","setExcitation"))
		{
			string val=mpParameter->GetKey("/aero_main/output/production","setExcitation");
			if(val=="active")
			{
				Log::mI<<"Excitations allowed for the output of the productions"<<endl;
				excitations=true;
			}

		}

		if(!mIsBendingActivated)
		{
			string photfile=file+"-photon.dat";
			SpecieUtils::PrintProduction(mPhotoionizationResu,mAltGridKm,photfile,excitations);
			string cosmofile=file+"-cosmo.dat";
			SpecieUtils::PrintProduction(mCosmoionizationResu,mAltGridKm,cosmofile,excitations);
			string protofile=file+"-proton.dat";
			SpecieUtils::PrintProduction(mProtoionizationResu,mAltGridKm,protofile,excitations);
			string elecfile=file+"-electron.dat";
			SpecieUtils::PrintProduction(mElectronImpactResu,mAltGridKm,elecfile,excitations);

			SpecieUtils::PrintProduction(mTotalResu,mAltGridKm,file,excitations);
		}else
		{

			string photfile=file+"-photon.dat";
			Log::mL<<"Print photon"<<endl;
			SpecieUtils::PrintProduction(mPhotoionizationResu,mAtmoPath.mLengthKm,photfile,excitations,true);
			string cosmofile=file+"-cosmo.dat";
			Log::mL<<"Print cosmo"<<endl;
			SpecieUtils::PrintProduction(mCosmoionizationResu,mAtmoPath.mLengthKm,cosmofile,excitations,true);
			string protofile=file+"-proton.dat";
			Log::mL<<"Print proto"<<endl;
			SpecieUtils::PrintProduction(mProtoionizationResu,mAtmoPath.mLengthKm,protofile,excitations,true);
			string elecfile=file+"-electron.dat";
			Log::mL<<"Print electron"<<endl;
			SpecieUtils::PrintProduction(mElectronImpactResu,mAtmoPath.mLengthKm,elecfile,excitations,true);

			Log::mL<<"Print total"<<endl;
			SpecieUtils::PrintProduction(mTotalResu,mAtmoPath.mLengthKm,file,excitations,true);
		}
	}


	Log::mI<<"Output energy"<<endl;
	if(mpParameter->Exists("/aero_main/output/photoelectron_energy_flux"))
	{
		string file=mpParameter->Elem("/aero_main/output/photoelectron_energy_flux");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mbIsPhotoFluxDefined)
		{
			if(!mIsBendingActivated)
			{
				mpPhotoFlux->PrintEnergyFlux(file,mAltGridKm);
			}else
			{
				mpPhotoFlux->PrintEnergyFlux(file,mAtmoPath.mLengthKm,true);
			}
		}
	}
	if(mpParameter->Exists("/aero_main/output/proton_energy_flux") && mbIsProtonModel)
	{
		string file=mpParameter->Elem("/aero_main/output/proton_energy_flux");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mbIsProtonModel)
		{
			if(!mIsBendingActivated)
			{
				mpProtonModel->PrintEnergyFluxes(file, mAltGridKm, true, true);
			}else
			{
				mpProtonModel->PrintEnergyFluxes(file, mAltGridKm, true);
			}
		}
	}
	if(mpParameter->Exists("/aero_main/output/H_energy_flux") && mbIsProtonModel)
	{
		string file=mpParameter->Elem("/aero_main/output/H_energy_flux");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mbIsProtonModel)
		{
			if(!mIsBendingActivated)
			{
				mpProtonModel->PrintEnergyFluxes(file, mAltGridKm, false, true);
			}else
			{
				mpProtonModel->PrintEnergyFluxes(file, mAltGridKm, false);
			}
		}
	}
	if(mpParameter->Exists("/aero_main/output/pair_production"))
	{
		string file=mpParameter->Elem("/aero_main/output/pair_production");

		if(suffix!="")
		{
			file+=suffix;
		}
		mpElecModel->PrintEnergyPairs(file);
	}
	if(mpParameter->Exists("/aero_main/output/production_mean_energy"))
	{
		string file=mpParameter->Elem("/aero_main/output/production_mean_energy");

		if(suffix!="")
		{
			file+=suffix;
		}
		mpElecModel->PrintEnergyPerProduction(mElectronImpactResu, file);
	}
	
	if(mpParameter->Exists("/aero_main/output/energy_conservation"))
	{
		string file=mpParameter->Elem("/aero_main/output/energy_conservation");

		if(suffix!="")
		{
			file+=suffix;
		}
		mpElecModel->PrintEnergyConservation(file);
		if(mbIsProtonModel)
		{
			mpProtonModel->PrintEnergyConservation(file + "_proton.xml");
		}
	}
		
	if(mpParameter->Exists("/aero_main/output/electron_heating"))
	{
		string file=mpParameter->Elem("/aero_main/output/electron_heating");

		if(suffix!="")
		{
			file+=suffix;
		}
		mpElecModel->PrintElectronHeating(file);
	}
	

	Log::mI<<"Output chemistry"<<endl;
	if(mpParameter->Exists("/aero_main/output/chem_list"))
	{
		string file=mpParameter->Elem("/aero_main/output/chem_list");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mbIsChem)
		{
			mChem->PrintChemList(file);
		}

	}
	if(mpParameter->Exists("/aero_main/output/chem_atmo"))
	{
		string file=mpParameter->Elem("/aero_main/output/chem_atmo");

		if(suffix!="")
		{
			file+=suffix;
		}
		if(mbIsChem)
		{
			mChem->PrintChemDensModels(file);
		}

	}
	Log::mI<<"Output computed density"<<endl;
	//	if(mpParameter->Exists("/output/production_selected_species"))
	if(mpParameter->Exists("/aero_main/output/chem_selection"))
	{
		vector<TiXmlNode*> selections=mpParameter->GetNodes("/aero_main/output/chem_selection");
		vector<TiXmlNode*>::iterator sit;
		for(sit=selections.begin();sit!=selections.end();++sit)
		{
			if(!mpParameter->Exists(*sit,"//production_selected_species"))
			{
				Log::mW<<"There is a void /output/chem_selection !!!!"<<endl;
			}
			string file=mpParameter->Elem(*sit,"//production_selected_species");
			if(suffix!="")
			{
				file+=suffix;
			}
			deque< SpecieId > mes_id;
			if(mpParameter->Exists(*sit,"//select_species"))
			{
				vector<TiXmlNode*> species=mpParameter->GetNodes(*sit,"//select_species/Specie");
				vector<TiXmlNode*>::iterator jt;
				for(jt=species.begin();jt!=species.end();++jt)
				{
					string spname=mpParameter->GetKey(*jt,"/","name");
					string spstate=mpParameter->GetKey(*jt,"/","state");
					SpecieId tmpsp(spname,spstate);
					mes_id.push_back(tmpsp);
				}
			}else
			{
				SpecieId CO2pp("CO2++","X");
				mes_id.push_back(CO2pp);
				SpecieId Opp("O++","X");
				mes_id.push_back(Opp);
				SpecieId N2pp("N2++","X");
				mes_id.push_back(N2pp);
			}

			mChem->SelectedPrintDensities(mes_id, mAltGridKm, file);
		}
	}
	Log::mI<<"Fin output"<<endl;

}


void Atmo::InitMultipleComp()
{
	//Log::SetPriority(Log::INFO);
	Log::mS<<"We initialize the multiple computation"<<endl;

	if(mpParameter->Exists("/aero_main/atmosphere/use_bent_atmosphere"))
	{// Here, we work with the bent atmosphere
		Error err("InitMultipleComp()","Bent atmo used","The bent atmosphere is not valid for multiple computation yet");
		throw err;
	}

	if(mpParameter->Exists("/aero_main/sun/use_sun"))
	{
		mMultiplePhoto=true;
		mpPhotomodel=new Photoionization(mpParameter);
		mpPhotomodel->Init(&mPhotonGrideV,&mPhotonGrideVmin,&mPhotonGrideVmax,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV,&mAltGridKm,mpPlanet->mUA,mpPlanet->mRKm,mpPlanet->mSZADegree);

		mpPhotoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);

		mbIsPhotoFluxDefined=true;
		mIsPhotoionizationModel=true;
	}

	if(mpParameter->Exists("/aero_main/proton_cosmic"))
	{
		mMultipleCosmo=true;
		mpCosmoModel=new ProtoCosIonization(mpParameter,&mAltGridKm);

		mpCosmoFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);

		mbIsCosmoFluxDefined=true;

	}
	
	if(mpParameter->Exists("/aero_main/proton/use_proton"))
	{
		mMultipleProton=true;
		mpProtonModel=new ProtonHydrogenTransport(mpParameter,&mAltGridKm, mpProtonGaussianAngle, &mdB_B);

		mpProtonEFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);

		mbIsProtonEFluxDefined=true;

	}


	if(mpParameter->Exists("/aero_main/electron/use_electron"))
	{
		mMultipleElec=true;
		mpElecPrecipFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
		mpElecPrecipFlux->InitVoid(mAltGridKm.size());
		mbIsElecPrecipFluxDefined=true;
		Log::mD<<"We read the precipitation!"<<endl;
		mpElecPrecipFlux->ReadPrecipitation(mAltGridKm,(mpModelAtmosphere->mAtmoSpecies));
		mpTotalFlux=new EFlux(mpParameter,mpGAngle,&mElecBotEeV,&mElecCentEeV,&mElecDdengeV);
		mpElecModel=new ElectronImpactIonization(mpParameter,&mElectronDensity,&mElectronTemperature,&mAltGridKm);
		mbIsElecModelDefined=true;
	}
	if(mpParameter->Exists("/aero_main/chem/use_chem"))
	{
		mMultipleChem=true;
		mbIsChem=true;
		InitChem(); // We initalize the chemistry with currents parameters
		// when electron and species densities are changed, we have to reset it in the chemistry with the adapted parameters!
		// We have to check if the modifications of the chemistry modifies the emissions... (probably, because shared pointer

		if(mpParameter->Exists("/aero_main/emissions/use_emissions"))
		{
			mMultipleEmit=true;
			mEmit.reset(new Emission(mpParameter,mChem,mpPlanet,mAltGridKm));
			mEmit->InitMultipleEmit();
		}
	}
}

void Atmo::FinishMultipleComp()
{
	//Log::SetPriority(Log::INFO);
	Log::mS<<"We proceed to the uninitialization of the multiple computation"<<endl;
	// a lot of the defined objects are automatically deleted during the destruction of Atmo

}





void Atmo::MCompute()
{
	/*
	deque<Specie*>::iterator it;
	for(it=mTotalResu.begin();it!=mTotalResu.end();++it)
	{
		(*it)->ClearProd();
	}

	for(it=mPhotoionizationResu.begin();it!=mPhotoionizationResu.end();++it)
	{
		(*it)->ClearProd();
	}
	for(it=mCosmoionizationResu.begin();it!=mCosmoionizationResu.end();++it)
	{
		(*it)->ClearProd();
	}
	for(it=mElectronImpactResu.begin();it!=mElectronImpactResu.end();++it)
	{
		(*it)->ClearProd();
	}*/
	mPhotoionizationResu.erase(mPhotoionizationResu.begin(),mPhotoionizationResu.end());
	mCosmoionizationResu.erase(mCosmoionizationResu.begin(),mCosmoionizationResu.end());
	mProtoionizationResu.erase(mProtoionizationResu.begin(),mProtoionizationResu.end());
	mElectronImpactResu.erase(mElectronImpactResu.begin(),mElectronImpactResu.end());
	mTotalResu.erase(mTotalResu.begin(),mTotalResu.end());


	if(mMultiplePhoto)
	{
		Log::mS<<"Compute photoionization"<<endl;
		mpModelAtmosphere->ClearProductions();
		mpPhotomodel->ComputePhotoionization((mpModelAtmosphere->mAtmoSpecies), mPhotoionizationResu, *mpPhotoFlux);
	}
	if(mMultipleCosmo)
	{
		Log::mS<<"Compute cosmoionization"<<endl;
		mpModelAtmosphere->ClearProductions();
		mpCosmoModel->ComputeCosmoionization((mpModelAtmosphere->mAtmoSpecies), mCosmoionizationResu, *mpCosmoFlux);
	}
	if(mMultipleProton)
	{
		Log::mS<<"Compute protoionization"<<endl;
		mpModelAtmosphere->ClearProductions();
		//mpProtonModel->ComputeCosmoionization((mpModelAtmosphere->mAtmoSpecies), mCosmoionizationResu, *mpCosmoFlux);
	}
	if(mMultipleElec)
	{
		Log::mS<<"Compute elecionization"<<endl;
		mpModelAtmosphere->ClearProductions();
		*mpTotalFlux=*mpElecPrecipFlux;//+*mpPhotoFlux;
		if(mbIsPhotoFluxDefined)
			(*mpTotalFlux)+=*mpPhotoFlux;
		if(mbIsCosmoFluxDefined)
			(*mpTotalFlux)+=*mpCosmoFlux;
		if(mbIsProtonEFluxDefined)
			(*mpTotalFlux)+=*mpProtonEFlux;
		mpElecModel->ComputeElectronImpact((mpModelAtmosphere->mAtmoSpecies),*mpTotalFlux,mElectronImpactResu);
	}

	Log::mD<<"Merge species"<<endl;
	std::deque<Specie* > totalResutmp2=SpecieUtils::MergeResult(mPhotoionizationResu,mElectronImpactResu);
	std::deque<Specie* > totalResutmp=SpecieUtils::MergeResult(totalResutmp2,mProtoionizationResu);
	mTotalResu=SpecieUtils::MergeResult(mCosmoionizationResu,totalResutmp); // the overload is problematic for deleting the species... maybe I should use a smart pointer now
	for(size_t i=0;i<totalResutmp.size();++i)
		delete totalResutmp[i];
	for(size_t i=0;i<totalResutmp2.size();++i)
		delete totalResutmp2[i];


	if(mMultipleChem)
	{
		Log::mD<<"Use of chemistry"<<endl;
		// Reset atmosphere and electron density for the chemistry
	}
	if(mMultipleEmit)
	{
		Log::mD<<"Consider emissions"<<endl;
		//ProceedEmissions("");
	//	InitChem();
	//	mChem->ResetDensity(mpModelAtmosphere->mAtmoSpecies);
	/*	mChem->ResetProduction(
				mTotalResu,
				mPhotoionizationResu,
				mElectronImpactResu,
				mProtoionizationResu,
				mCosmoionizationResu
				);

		ublas::matrix<double> subResults;
		std::deque<std::string> subResultsNames;
		std::string warnings;
		SpecieId o1s("O","1S");
		ublas::vector<double> dens=mChem->GetDensity(o1s,subResults,subResultsNames,warnings);
		Log::mL<<"GGTEST DENSITY "<<dens<<endl;*/
		InitChem();
		mEmit.reset(new Emission(mpParameter,mChem,mpPlanet,mAltGridKm));
		mEmit->InitMultipleEmit();
	//	mEmit->ComputeEmissions("");
		mEmit->CompMultipleEmit();
		//mEmit.reset(new Emission(mpParameter,mChem,mpPlanet,mAltGridKm));
		//mEmit->ComputeEmissions();
	}
	Log::mS<<"End of multiple computation"<<endl;
}



/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |Reset e- precip |
 |                |
 +---------------*/


void Atmo::ResetElecPrecip(int vModel, std::deque<double> vParams, std::deque<double> vAddParams)
{
	if(mbIsElecPrecipFluxDefined)
	{
		mpElecPrecipFlux->ResetElecPrecip(vModel, vParams, vAddParams);
	}else
	{
		Error err("ResetElecPrecip","Electron precipitation object not defined","The initialisation of the electron precipitation went wrong: you probably did not select the electron impact model");
		throw err;
	}
}



/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |Reset atmosphere|
 |                |
 +---------------*/


void Atmo::ResetSpecieDens(std::string vName,int vModel, std::deque<double> vParams, std::deque<double> vAddParams)
{
	if(mIsModelAtmosphere)
	{
		mpModelAtmosphere->ResetSpecieDens(vName,vModel,vParams,vAddParams,mNeutralTemperature);
	}else
	{
		Error err("ResetSpecieDens","Model atmosphere not defined","The initialisation of the atmosphere went probably wrong");
		throw err;
	}
}


void Atmo::ResetSpecieDensInterp(std::string vName, ublas::vector<double> vAltKm, ublas::vector<double> vDenscm_3)
{
	if(mIsModelAtmosphere)
	{
		ublas::vector<double> newdens=MathFunction::IntLin(vAltKm,vDenscm_3,mAltGridKm);
		mpModelAtmosphere->ResetSpecieDensInterp(vName,newdens,mNeutralTemperature);
	}else
	{
		Error err("ResetSpecieDens","Model atmosphere not defined","The initialisation of the atmosphere went probably wrong");
		throw err;
	}
}

void Atmo::ResetElectronDens(int vModel, std::deque<double> vParams, std::deque<double> vAddParams)
{
	switch(vModel)
	{
		case 0:
			{
				assert(vParams.size()==4);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(1);
				double T0=vParams.at(2);
				double shape=vParams.at(3);
				//unsigned pos=CloseInVector(alt0,mAltGridKm);
				mElectronDensity=MathFunction::BatesWalkerProfile(mAltGridKm,alt0,dens0,T0,Texo,shape,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000);

			}
			break;

		case 1:
			{
				assert(vParams.size()==3);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				mElectronDensity=MathFunction::ChapmanProfile(mAltGridKm,vParams.at(0),dens0,vParams.at(2));
			}
			break;
		case 2:
			{
				assert(vParams.size()==3);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				mElectronDensity=MathFunction::GaussianProfile(mAltGridKm,vParams.at(0),dens0,vParams.at(2));
			}
			break;
		case 3:
			{
				assert(vParams.size()==2);
				double alt0=vAddParams.at(0);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(1);
				mElectronDensity=MathFunction::ExpProfile(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000);
				
			}
			break;
		case 4:
			{
				unsigned siz=static_cast<unsigned>(vAddParams.at(0));
				ublas::vector<double> alt(siz),logval(siz);
				for(unsigned k=0;k<siz;++k)
				{
					alt[k]=vAddParams.at(k+1);
					logval[k]=vParams.at(k);
				}
				mElectronDensity = MathFunction::SplineInterpExp(alt,logval,mAltGridKm);
			}
			break;
		case 5:
			{
				assert(vParams.size()==4);
				assert(vAddParams.size()==3);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double alt1=vAddParams.at(1);
				double mixmassamu=vAddParams.at(2);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double dens1=vParams.at(1);
				if(dens1<0)
					dens1=1E-50;

				double Texo=vParams.at(2);
				double Tmeso=vParams.at(3);
				mElectronDensity=MathFunction::DoubleExpProfile(mAltGridKm,alt0,alt1,dens0,dens1,Texo,Tmeso,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.,mixmassamu);
			}
			break;
		case 6:
			{
				assert(vParams.size()==4);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				double SZA=vParams.at(3);
				mElectronDensity=MathFunction::ChapmanCos(mAltGridKm,alt0,dens0,Texo,SZA,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);
			}
			break;
		case 7:
			{
				assert(vParams.size()==3);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				double C=vAddParams.at(0);
				mElectronDensity=MathFunction::ChapmanVar(mAltGridKm,alt0,dens0,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);
			}
			break;
		case 8:
			{
				assert(vParams.size()==3);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				mElectronDensity=MathFunction::Epstein(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.);
			}
			break;
		case 9:
			{
				assert(vParams.size()==5);
				assert(vAddParams.size()==3);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double alt1=vAddParams.at(1);
				double mixmassamu=vAddParams.at(2);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double dens1=vParams.at(1);
				if(dens1<0)
					dens1=1E-50;

				double Texo=vParams.at(2);
				double Tmeso=vParams.at(3);
				double C = vParams.at(4);
				mElectronDensity=MathFunction::ExpHyperbolaProfile(mAltGridKm,alt1,alt0,dens1,dens0,Tmeso,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,1/2000.,mixmassamu);
			}
			break;
		default:
			Error err("Reset electron density","Model error"," The model "+ntostr(vModel)+" is not a valid model for the electron density profile");
	}
	if(mbIsChem)
	{
		mChem->ResetEdens(mElectronDensity);
	}
}


void Atmo::ResetElectronDensInterp(ublas::vector<double> vAltKm, ublas::vector<double> vDenscm_3)
{
	mElectronDensity=MathFunction::IntLin(vAltKm,vDenscm_3,mAltGridKm);
	if(mbIsChem)
	{
		mChem->ResetEdens(mElectronDensity);
	}
}

ublas::vector<double> Atmo::RetrieveObservable(SpecieId vId, std::deque<int> vParams, ublas::vector<double> vComparisonScale )
{
	if(!mbIsChem)
	{
		Error err("RetrieveObservable","Chemistry must be activated to retrieve observable","");
		throw err;
	}

	switch(vParams.at(0))
	{
		case 0:
			{ // production
				ublas::vector<double> data=mChem->GetTotProd(vId.mName,vId.mState);
				if(vComparisonScale.size()==0)
					return data;
				return MathFunction::IntLin(mAltGridKm,data,vComparisonScale);
			}
			break;
		case 1:
			{// VER
				ublas::vector<double> data=mEmit->GetVER(vId,vParams.at(1));
				if(vComparisonScale.size()==0)
					return data;
				return MathFunction::IntLin(mAltGridKm,data,vComparisonScale);
			}
			break;
		case 2:
			{// Limb
				ublas::vector<double> data=mEmit->GetLimb(vId,vParams.at(1));
				if(vComparisonScale.size()==0)
					return data;
				return MathFunction::IntLin(mEmit->GetLimbParameter(),data,vComparisonScale);
			}
			break;
		case 3:
			{//Density
				ublas::matrix<double> subResults;
				std::deque<std::string> subResultsNames;
				std::string warnings;
				ublas::vector<double> dens=mChem->GetDensity(vId,subResults,subResultsNames,warnings);
				if(0==dens.size())
				{
					dens=mChem->GetDens(vId.mName,vId.mState);
				}
				return MathFunction::IntLin(mAltGridKm,dens,vComparisonScale);
			}
			break;


		default:
			Error err("Retrieve observable","impossible to find your type of observable","");
			throw err;
	}

	ublas::vector<double> resu;
	assert(vParams.size()>0);
	assert(vId.StandardName()!="");
	return resu;
}










/*
void Atmo::StandardGrid(int number,int type,double min,double max,std::vector<double>& cent,std::vector<double> ddeng,double& spfac)
{
	if(type==0)
	{
		cent=MathGrid::gridexp(min,max,number);
		spfac=-1;
		ddeng=MathGrid::WidthGrid(cent,1);
		return;
	}
	if(type==1)
	{
		MathGrid::gridpolo(number,min,max,cent,ddeng,spfac);

	}

	if(type==2)
	{

	}


	cout<<"Nous ne trouvons pas votre type de grille"<<endl;
	cout<<"you grid type is not found"<<endl;
	exit(1);

}
*/
