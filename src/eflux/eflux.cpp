/** 
 * \file eflux.cpp
 * \brief Implements the EFlux class 
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: eflux.cpp 1682 2013-02-18 20:24:31Z gronoff $
 *
 */



#include "eflux.hpp"
using namespace std;
typedef std::deque< ublas::vector<double>* > DequeUblas;
typedef std::deque< ublas::matrix<double>* > DequeUblasMat;
typedef std::deque< ublas::vector< ublas::matrix<double> >* > DequeVecMat;

EFlux::EFlux(XmlParameters* pParam,MathFunction::GaussianAngle* pAngle,ublas::vector<double>* pEBotE,ublas::vector<double>* pECentE,ublas::vector<double>* pEDdeng)
{
	mpParameter=pParam;
	mpGAngle=pAngle;
	mpElecBotEeV=pEBotE;
	mpElecCentEeV=pECentE;
	mpElecDdengeV=pEDdeng;
	mPrecipitationDefined=false;
	mWarningMissingCrs=false;


}


void EFlux::InitVoid(unsigned vNbAlt)
{
//	vector<double> zeroenergy(mpElecBotEeV->size(),0.);
//	vector< vector<double> > zeroaltener(vNbAlt,zeroenergy);
//	mElecEcm_3s_1eV_1=zeroaltener;
	mElecEcm_3s_1eV_1.resize(vNbAlt,mpElecBotEeV->size());
	mElecEcm_3s_1eV_1.clear();

//	vector<double> zeroalt(vNbAlt,0.);
//	mEleccm_3s_1=zeroalt;
	mEleccm_3s_1.resize(vNbAlt);
	mEleccm_3s_1.clear();	
	
//	mFluxcm_2s_1eV_1=zeroaltener;
	mFluxcm_2s_1eV_1.resize(vNbAlt,mpElecBotEeV->size());
	mFluxcm_2s_1eV_1.clear();


//	vector<double> zeroa(mpGAngle->mNbAngles,0.);
//	vector< vector<double> > zeroenergya(mpElecBotEeV->size(),zeroa);
//	vector< vector< vector<double> > > zeroaltenera(vNbAlt,zeroenergya);
//	mFluxAcm_2s_1eV_1sr_1=zeroaltenera;
	
	mFluxAcm_2s_1eV_1sr_1.resize(vNbAlt);
	for(unsigned i=0;i<vNbAlt;++i)
	{
		mFluxAcm_2s_1eV_1sr_1[i].resize(mpElecBotEeV->size(),mpGAngle->mNbAngles);
		mFluxAcm_2s_1eV_1sr_1[i].clear();
	}



}


void EFlux::InitIsotropic(ublas::matrix<double>  vElecE,ublas::vector<double> vElec,std::deque<Specie*> vSp)
{

	unsigned nbalt=vElec.size();
	assert(nbalt>0);// It is too complicated else
	unsigned nbener=vElecE.size2();
	unsigned nbang=mpGAngle->mNbAngles;
	// If the specie has no defined cross section, we do not throw an error, but we put a warning flag in the class!
	mElecEcm_3s_1eV_1=(vElecE);
	mEleccm_3s_1=(vElec);
	
	
//	vector<double> tmpener(nbener,0);
	// Densig ialt, ien
	//vector< vector<double> > densig(nbalt,tmpener);
	ublas::matrix<double> densig(nbalt,nbener);
	densig.clear();
	// We resize the mFlux, the same way
//	mFluxcm_2s_1eV_1.resize(nbalt,tmpener);
	mFluxcm_2s_1eV_1.resize(nbalt,nbener);
	mFluxcm_2s_1eV_1.clear();


//	vector<double> tmpa(nbang,0.);
//	vector< vector<double> > tmpe(nbener,tmpa);
	// We resize the other mFlux
//	mFluxAcm_2s_1eV_1sr_1.resize(nbalt,tmpe);

	mFluxAcm_2s_1eV_1sr_1.resize(nbalt);
	for(unsigned i=0;i<nbalt;++i)
	{
		mFluxAcm_2s_1eV_1sr_1[i].resize(nbener,nbang);
		mFluxAcm_2s_1eV_1sr_1[i].clear();
	}



	for(unsigned i=0;i<vSp.size();++i)
	{
		if(!vSp.at(i)->CheckElecCrs())
		{
			//Log::SetPriority(Log::WARNING);

			mWarningMissingCrs=true;
			mMissingCrsMsg.push_back("Lack of electrons cross section for your specie "+vSp[i]->mName);
			Log::mW<<"WARNING CROSS SECTION"<<endl;
			Log::mW<<"Lack of electrons cross section for your specie "+vSp[i]->mName<<endl;

		}else
		{
			for(unsigned j=0;j<nbalt;++j)
			{
				for(unsigned k=0;k<nbener;++k)
				{
					densig(j,k)+=((*(vSp[i]->mpElecCrs)).mCrsCincm2[k] + (*(vSp[i]->mpElecCrs)).mElasticCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
					/*if(isnan(densig(j,k)))
					{
					Log::mL<<j<<" "<<k<<endl;
					if(isnan((*vSp[i]).mTotDensitycm_3[j]))
					Log::mL<<"bien ce que je pensais"<<endl;
					}*/
				}
			}
		}
	}// Ok, we finish  the computation of the densig -> absorption/cm


	for(unsigned j=0;j<nbalt;++j)
	{
		for(unsigned k=0;k<nbener;++k)
		{
			if(densig(j,k)!=0)
			{
				mFluxcm_2s_1eV_1(j,k)=mElecEcm_3s_1eV_1(j,k)/densig(j,k);
				if(isnan(mFluxcm_2s_1eV_1(j,k)))
				{
				/*	Log::mL<<j<<" "<<k<<endl;*/
					Log::mD<< mFluxcm_2s_1eV_1(j,k)<<"  =  "<<mElecEcm_3s_1eV_1(j,k)<<" /  --->"<<densig(j,k);
					Error err("Overflow","Init isotropic","autre flux en nan");
					throw err;
				}
			}else
			{// really easier here
				//Log::SetPriority(Log::WARNING);
				Log::mW<<"WARNING, densig equal to 0 here"<<endl;
				mFluxcm_2s_1eV_1(j,k)=0;
				if(!mWarningMissingCrs)
				{
					Log::mW<<"WARNING, densig equal to 0 here"<<endl;
					mWarningMissingCrs=true;
					mMissingCrsMsg.push_back("Warning, densig equal to 0, check your cross sections for the species");
				}
			}
//			Log::mL<<"altitude n "<<j<<" energy :"<<(*mpElecCentEeV)[k]<<endl;
			for(unsigned ang=0;ang<nbang;++ang)
			{// We divide by 4Pi because we have a production -> production/sr
				mFluxAcm_2s_1eV_1sr_1(j)(k,ang)=mFluxcm_2s_1eV_1(j,k)/(4.*PI);
		
				if(isnan(mFluxAcm_2s_1eV_1sr_1(j)(k,ang)))
				{
					Log::mD<< mFluxAcm_2s_1eV_1sr_1[j](k,ang)<<" "<<mFluxcm_2s_1eV_1(j,k)/(4.*PI)<<endl;
					Error err("Overflow","Init isotropic","mflux en nan");
					throw err;
				
				}


			}//Nb if we have had a production by angle, (N electrons produced between angle a and a+theta)
			// Then we should have divided by 2 Pi
		}
	}



}


void EFlux::ResetElecPrecip(int vModel, std::deque<double> vParams, std::deque<double> vAddParams)
{
	mPrecipitationDefined = true;
	unsigned nang = mpGAngle->mNbAngles;
	mPrecipFluxDowncm_2s_1sr_1.clear();
	mPrecipFluxUpcm_2s_1sr_1.clear();

/*
	unsigned isotro=0, powlaw=0;
	double entot=0.,E0=0.;
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/E0",E0);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/entot",entot);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/powlaw",powlaw);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/isotro",isotro);
	bool b_power_law=false;
	if(powlaw>0)
		b_power_law=true;
*/
	switch(vModel)
	{
		case 0: {
				Log::mL<<"Null precipitation"<<endl;
			}
			break;
		case 1:
			{
				Log::mL<<"Maxwell precipitation"<<endl;
				double entot = vParams[0];
				double E0 = vParams[1];
				unsigned powlaw = static_cast<unsigned>(vAddParams[0]);
				unsigned isotro = static_cast<unsigned>(vAddParams[1]);
				bool b_power_law=false;
				if(powlaw>0)
					b_power_law=true;
				MathFlux::InMaxw(entot, E0, b_power_law, nang, *mpElecCentEeV, isotro, mpGAngle->mXmu, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1);
				MathFlux::NormFlux(entot, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1, mpGAngle->mWeight, mpGAngle->mXmu, *mpElecCentEeV, *mpElecDdengeV);
			}
			break;
		case 2:
			{
				Log::mL<<"Gaussian precipitation"<<endl;
				double entot = vParams[0];
				double E0 = vParams[1];
				unsigned powlaw = static_cast<unsigned>(vAddParams[0]);
				unsigned isotro = static_cast<unsigned>(vAddParams[1]);
				bool b_power_law=false;
				if(powlaw>0)
					b_power_law=true;
				MathFlux::InGauss(entot, E0, b_power_law, nang, *mpElecCentEeV, isotro, mpGAngle->mXmu, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1);
				MathFlux::NormFlux(entot, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1, mpGAngle->mWeight, mpGAngle->mXmu, *mpElecCentEeV, *mpElecDdengeV);
			}
			break;
		case 3:
			{
				Log::mL<<"Dirac precipitation"<<endl;
				double entot = vParams[0];
				double E0 = vParams[1];
				unsigned powlaw = static_cast<unsigned>(vAddParams[0]);
				unsigned isotro = static_cast<unsigned>(vAddParams[1]);
				bool b_power_law=false;
				if(powlaw>0)
					b_power_law=true;
				MathFlux::InDirac(entot, E0, b_power_law, nang, *mpElecCentEeV, isotro, mpGAngle->mXmu, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1);
				Log::mD<<"Normalisation"<<endl;
				MathFlux::NormFlux(entot, mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1, mpGAngle->mWeight, mpGAngle->mXmu, *mpElecCentEeV, *mpElecDdengeV);
				Log::mD<<"Return"<<endl;
			}
			break;
		case 4:
			{
				Log::mL<<"Spline precipitation"<<endl;
				ublas::vector<double> enes,logval, flux;
				enes = StdToUblas(vAddParams);
				logval = StdToUblas(vParams);
				flux = MathFunction::SplineInterpExp(enes, logval, *mpElecCentEeV);
				MathFlux::FillFlux(mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1, flux);
			}
			break;
		default:
			Log::mE<<"Impossible to find your precipitation model"<<endl;

			Error err("EFlux::ResetElecPrecip","switch type","Impossible to find your precipitation model, type should be between 0 and 4");
			throw err;
	};
}


void EFlux::ReadPrecipitation(const ublas::vector<double> & vAltitudes,std::deque<Specie*> vSp)
{
	mPrecipitationDefined=true;

	unsigned nang=mpGAngle->mNbAngles;
	unsigned nang_2=nang/2;//mpGAngle->mNbAngles/2;
	unsigned nener=(*mpElecCentEeV).size();
	//	vector<double> tmpvec(nang_2,0.);
	mPrecipFluxDowncm_2s_1sr_1.resize(nener,nang_2);
	mPrecipFluxUpcm_2s_1sr_1.resize(nener,nang_2);
	mPrecipFluxDowncm_2s_1sr_1.clear();
	mPrecipFluxUpcm_2s_1sr_1.clear();


	//Log::SetPriority(Log::CONFIG);
	if(!mpParameter->Exists("/aero_main/electron/use_precipitation"))
	{
		Log::mL<<"No precipitations"<<endl;
		return;// No precipitations
	}
	mpParameter->ExistsOrDie("/aero_main/electron/precipitation","You want to use precipitations, YOU HAVE TO DEFINE THEM");

	if(mpParameter->Exists("/aero_main/electron/precipitation/measurement"))
	{
		Log::mL<<"We use electron precipitation measurement"<<endl;

		ReadMeasuredPrecipitation();
		return;
	}

	mpParameter->ExistsOrDie("/aero_main/electron/precipitation/use_model","You don't want to use measured precipitation, it's your right, but then you have to give me a precipitation model... or not to use precipitations.");
	unsigned type=0;
	mpParameter->GetNKey("/aero_main/electron/precipitation/use_model","type",type);
	unsigned isotro=0, powlaw=0;
	double entot=0.,E0=0.;
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/E0",E0);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/entot",entot);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/powlaw",powlaw);
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/isotro",isotro);
	bool b_power_law=false;
	if(powlaw>0)
		b_power_law=true;
	//	bool b_isotrope=false;
	//	if(isotro>0)
	//		b_isotrope=true;
//	Log::SetPriority(Log::CONFIG);
	switch(type)
	{
		case 0: Log::mL<<"Null precipitation"<<endl;
			break;
		case 1:
			Log::mL<<"Maxwell precipitation"<<endl;
			MathFlux::InMaxw(entot,E0,b_power_law,nang,*mpElecCentEeV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpElecCentEeV,*mpElecDdengeV);
			break;
		case 2:
			Log::mL<<"Gaussian precipitation"<<endl;
			MathFlux::InGauss(entot,E0,b_power_law,nang,*mpElecCentEeV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpElecCentEeV,*mpElecDdengeV);
			break;
		case 3:
			Log::mL<<"Dirac precipitation"<<endl;
			MathFlux::InDirac(entot,E0,b_power_law,nang, *mpElecCentEeV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			Log::mD<<"Normalisation"<<endl;
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpElecCentEeV,*mpElecDdengeV);
			Log::mD<<"Return"<<endl;
			break;
		case 4:
			{
				Log::mL<<"Spline precipitation"<<endl;
				ublas::vector<double> enes,logval, flux;
				mpParameter->Get1DArray("/aero_main/electron/precipitation/use_model/energies", enes);
				mpParameter->Get1DArray("/aero_main/electron/precipitation/use_model/logflux", logval);
				flux = MathFunction::SplineInterpExp(enes, logval, *mpElecCentEeV);
				MathFlux::FillFlux(mPrecipFluxDowncm_2s_1sr_1, mPrecipFluxUpcm_2s_1sr_1, flux);
			}
			break;
		default:
			Log::mE<<"Impossible to find your precipitation model"<<endl;

			Error err("EFlux::ReadPrecipitation","switch type","Impossible to find your precipitation model, type should be between 0 and 4");
			throw err;

	};
	if(mpParameter->Exists("/aero_main/electron/precipitation/use_model/altitude"))
	{
		double altprecip=vAltitudes[0];
	mpParameter->GetValue("/aero_main/electron/precipitation/use_model/altitude",altprecip);
		Log::mD<<"You want your precipitations at the position "<<altprecip<<" km "<<endl;
		double distance=nabs(vAltitudes[0]-altprecip);
		unsigned position=0;
		for(unsigned i=1;i<vAltitudes.size();++i)
		{
			if(nabs(vAltitudes[i]-altprecip)<distance)
			{
				position=i;
				distance=nabs(vAltitudes[i]-altprecip);
			}
		}
		Log::mD<<"The closes position is : "<<vAltitudes[position]<<" km"<<endl;



		if(mFluxAcm_2s_1eV_1sr_1.size()!=vAltitudes.size())
		{
			Error err("uninitialized Flux","","");

			throw err; 
		}
/*
		ublas::matrix<double> elecE(vAltitudes.size(),nener);
		ublas::vector<double> elec(vAltitudes.size());
		elecE.clear();
		elec.clear();

		for(unsigned e=0;e<nener;++e)
		{
			elecE(position,e)=mPrecipFluxDowncm_2s_1sr_1(e,0)*4*PI;

		}
*/

		
		Log::mD<<"Remplissage du vecteur"<<endl;
		ublas::vector< ublas::matrix<double> > myflux(nener);
		for(unsigned e=0;e<nener;++e)
		{
			myflux(e).resize(vAltitudes.size(),nang);
			myflux(e).clear();// init at 0

				for(unsigned k=0;k<nang;++k)
				{
					if(k<nang_2)
					{
						myflux(e)(position,k)=mPrecipFluxDowncm_2s_1sr_1(e,k);
					}else
					{
						myflux(e)(position,k)=mPrecipFluxUpcm_2s_1sr_1(e,k-nang_2);
					}
				}
		}

		/*for(unsigned i=0;i<nang;++i)
		{
			ublas::matrix_column< ublas::matrix<double> > ro(mFluxAcm_2s_1eV_1sr_1[position],i);
			if(i<nang_2)
			{
				ublas::matrix_column< ublas::matrix<double> > rdwn(mPrecipFluxDowncm_2s_1sr_1,i);
				ro=rdwn;
			}else
			{
				ublas::matrix_column< ublas::matrix<double> > rdwn(mPrecipFluxUpcm_2s_1sr_1,i-nang_2);
				ro=rdwn;

			}
		}*/
		//Log::mL<<"Your matrix : "<<mFluxAcm_2s_1eV_1sr_1[position]<<endl;

		//The end, we put the fluxes at 0, and we delete the precipitation definition (it is inside the flux vector!)
		mPrecipitationDefined=false;
		mPrecipFluxDowncm_2s_1sr_1.resize(0,0);
		mPrecipFluxUpcm_2s_1sr_1.resize(0,0);
		Log::mD<<"C'est partit pour l'anisotropique"<<endl;
		AddAnisotropicFlux(myflux, vSp);
	}
}


void EFlux::ReadMeasuredPrecipitation()
{
	// If no downward...
	mpParameter->ExistsOrDie("/aero_main/electron/precipitation/measurement/downward/mflux","You have to provide at least one downward flux for your measured precipitation");

	mpParameter->ExistsOrDie("/aero_main/electron/precipitation/measurement/Egrid","You have to defined the energetic grid for your measured precipitations");

	ublas::vector<double> elec_grid;
	mpParameter->Get1DArray("/aero_main/electron/precipitation/measurement/Egrid",elec_grid);

	// 1) We read the number of upward and downward values:
	unsigned dwn_number=mpParameter->Numbers("/aero_main/electron/precipitation/measurement/downward/mflux");
	unsigned upw_number=mpParameter->Numbers("/aero_main/electron/precipitation/measurement/upward/mflux");
	// 2) We get the values
	vector<TiXmlNode*> dwn_node=mpParameter->GetNodes("/aero_main/electron/precipitation/measurement/downward/mflux");
	vector<TiXmlNode*> upw_node=mpParameter->GetNodes("/aero_main/electron/precipitation/measurement/upward/mflux");

	deque<double> upw_angle_degree;// The upward angles in degree
	deque<double> dwn_angle_degree;// The downward angles in degree

	deque< ublas::vector<double> > upw_fluxes ;// The upward fluxes
	deque< ublas::vector<double> > dwn_fluxes ;// The downward fluxes

	assert(dwn_node.size()==dwn_number);
	assert(upw_node.size()==upw_number);
	Log::mD<<"For loop dwn"<<endl;
	for(unsigned i=0;i<dwn_number;++i)
	{
		double tmp_angle=0;
		ublas::vector<double> tmp_flux;

		Log::mD<<"Read angle for i="<<i<<endl;
		try
		{
		string str=mpParameter->GetKey(dwn_node[i],"/","angle");
		}
		catch(Error & err)
		{
			Log::mE<<"MEGAPROBLEME"<<endl;
			err.Affiche();
		}
		catch(...)
		{
			Log::mE<<"GIGAPROBLEME!!!!"<<endl;
			Error err("Pb","Pb","Pb");
			throw err;
		}

		mpParameter->GetNKey(dwn_node[i],"/","angle",tmp_angle);
		Log::mD<<"tmp angle : "<<tmp_angle<<endl;
		dwn_angle_degree.push_back(tmp_angle);
		mpParameter->Get1DArray(dwn_node[i],"/",tmp_flux);
		dwn_fluxes.push_back(tmp_flux);
	}

	Log::mD<<"For loop upw"<<endl;
	for(unsigned i=0;i<upw_number;++i)
	{
		double tmp_angle=0;
		ublas::vector<double> tmp_flux;
		mpParameter->GetNKey(upw_node[i],"/","angle",tmp_angle);
		upw_angle_degree.push_back(tmp_angle);
		mpParameter->Get1DArray(upw_node[i],"/",tmp_flux);
		upw_fluxes.push_back(tmp_flux);
	}
//mPrecipFluxDowncm_2s_1sr_1


	Log::mL<<"Redistribution downward flux"<<endl;
	// 3) Redistribution of the downward flux
	RedistributeFlux(elec_grid,dwn_angle_degree,dwn_fluxes,mPrecipFluxDowncm_2s_1sr_1);

	Log::mL<<"Redistribution upward flux"<<endl;
	// 3) Redistribution of the upward flux
	if(upw_number>0)
		RedistributeFlux(elec_grid,upw_angle_degree,upw_fluxes,mPrecipFluxUpcm_2s_1sr_1,false);
	Log::mD<<"Bingo : end of the flux reading"<<endl;
}

void EFlux::RedistributeFlux(ublas::vector<double> vElecGrideV,
				std::deque<double> vAnglesDegree,
				const std::deque< ublas::vector<double> >& vInputFluxcm_2sr_1s_1,
				ublas::matrix<double>& rFluxcm_2sr_1s_1,bool vbDown)
{// We redistribute from the old grid to the new one...


	//	unsigned nb_input_angle=vAnglesDegree.size();
	//	assert(vElecGrideV.size()==vAnglesDegree.size());
	//	cout<<rFluxcm_2sr_1s_1.size()<<endl;
	unsigned nang2=mpGAngle->mAngzb.size()/2;

	if(vAnglesDegree.size()==1)
	{// We have one flux -> isotrope!!!
		//Log::SetPriority(Log::INFO);
		Log::mI<<"One measured flux: isotrope!"<<endl;
		ublas::vector<double> interp_flux=MathFunction::IntLog(vElecGrideV,vInputFluxcm_2sr_1s_1[0],*mpElecCentEeV);

		// The flux has the good size	
		//	rFluxcm_2sr_1s_1.resize(mpGAngle->mNbAngles/2);
		assert(interp_flux.size()==mpElecCentEeV->size());
		assert(rFluxcm_2sr_1s_1.size1()==interp_flux.size());
		
		for(unsigned i=0;i<rFluxcm_2sr_1s_1.size1();++i)
		{
			for(unsigned j=0;j<rFluxcm_2sr_1s_1.size2();++j)
			{
				rFluxcm_2sr_1s_1(i,j)=interp_flux[i];
			}
		}
		return;
	}




	// The map, when we use an iterator, is sorted!!!
	map< double,ublas::vector<double> > flux_redist;

	for(unsigned i=0;i<vInputFluxcm_2sr_1s_1.size();++i)
	{ 
		flux_redist[vAnglesDegree[i]]=MathFunction::IntLog(vElecGrideV,vInputFluxcm_2sr_1s_1[i],*mpElecCentEeV);
	}


	std::deque<double> tangle;
	std::deque<double> gangle;
	unsigned itmp=0;
	for(map< double,ublas::vector<double> >::iterator it=flux_redist.begin();it!=flux_redist.end();++it)
	{
		double angle=(*it).first;
		if(vbDown)
		{
			angle*=-1;
			angle+=90;
		}else
		{
			angle+=90;
		}

		Log::mD<<"Resulting angle :"<<angle<<" ["<<itmp<<"]"<<endl;
		tangle.push_back(angle);
		itmp++;
	}
	ublas::vector<double> tangle2(tangle.size());
	std::copy(tangle.begin(),tangle.end(),tangle2.begin());
	for(unsigned i=0;i<nang2;++i)
	{

		if(vbDown)
		{
			gangle.push_back(mpGAngle->mAngzb[i]);
		}else
		{
			gangle.push_back(mpGAngle->mAngzb[i+nang2]);
		}
		Log::mD<<"gangle : "<<gangle[i]<<endl;
	}
	ublas::vector<double> gangle2(gangle.size());
	std::copy(gangle.begin(),gangle.end(),gangle2.begin());
	//assert(itmp==nang2);

	for(unsigned i=0;i<(*mpElecCentEeV).size();++i)
	{
		std::deque<double> tmp;
		for(map< double,ublas::vector<double> >::iterator it=flux_redist.begin();it!=flux_redist.end();++it)
		{
			tmp.push_back(((*it).second)[i]);
		}
		ublas::vector<double> tmp2(tmp.size());
		std::copy(tmp.begin(),tmp.end(),tmp2.begin());
		/*
		   for(unsigned it=0;it<mpGAngle->mAngzb.size();++it)
		   {
		   Log::mL<<"angzb : "<<mpGAngle->mAngzb[it]<<endl;
		   }
		   */

		ublas::vector<double> fresu=MathFunction::IntLin(tangle2,tmp2,gangle2); // Log interpolation : some values must be 0

		assert(fresu.size()==rFluxcm_2sr_1s_1.size2());

		for(unsigned j=0;j<fresu.size();++j)
		{
			if(fresu[j]<0)
				fresu[j]=0.;
			rFluxcm_2sr_1s_1(i,j)=fresu[j];
		}
	}
}

EFlux operator+(const EFlux& vA , const EFlux&vB)
{// The grids are supposed to be the same

	Log::mD<<"We do a flux addition!!!"<<endl;
	// We do strong assert: not only it is the same grid, but the same pointer
	assert(vA.mpElecBotEeV==vB.mpElecBotEeV);
	assert(vA.mpElecCentEeV==vB.mpElecCentEeV);
	assert(vA.mpElecDdengeV==vB.mpElecDdengeV);
	// It works because the function is friend
	assert(vA.mpGAngle==vB.mpGAngle);


	EFlux resu(vA.mpParameter,vA.mpGAngle,vA.mpElecBotEeV,vA.mpElecCentEeV,vA.mpElecDdengeV);

	// We add the electron production by altitude
	resu.mEleccm_3s_1=AddOrFullVec(vA.mEleccm_3s_1,vB.mEleccm_3s_1);

	// We add the electron production by altitude and energy
	resu.mElecEcm_3s_1eV_1=AddOrFullVec2(vA.mElecEcm_3s_1eV_1,vB.mElecEcm_3s_1eV_1);
	
	// We add the electron flux function of altitude and energy
	resu.mFluxcm_2s_1eV_1=AddOrFullVec2(vB.mFluxcm_2s_1eV_1,vB.mFluxcm_2s_1eV_1);
	// We add the electron flux function of altitude, energy, and angle
	resu.mFluxAcm_2s_1eV_1sr_1=AddOrFullVec3(vB.mFluxAcm_2s_1eV_1sr_1,vB.mFluxAcm_2s_1eV_1sr_1);

	resu.mPrecipitationDefined=false; // Default

	if(vA.mPrecipitationDefined&&!vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined=true;
		resu.mPrecipFluxDowncm_2s_1sr_1=vA.mPrecipFluxDowncm_2s_1sr_1;
		resu.mPrecipFluxUpcm_2s_1sr_1=vA.mPrecipFluxUpcm_2s_1sr_1;
	}

	if(!vA.mPrecipitationDefined&&vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined=true;
		resu.mPrecipFluxDowncm_2s_1sr_1=vB.mPrecipFluxDowncm_2s_1sr_1;
		resu.mPrecipFluxUpcm_2s_1sr_1=vB.mPrecipFluxUpcm_2s_1sr_1;
	}
	if(vA.mPrecipitationDefined&&vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined=true;
		resu.mPrecipFluxDowncm_2s_1sr_1=(vA.mPrecipFluxDowncm_2s_1sr_1+vB.mPrecipFluxDowncm_2s_1sr_1);
		resu.mPrecipFluxUpcm_2s_1sr_1=(vA.mPrecipFluxUpcm_2s_1sr_1+vB.mPrecipFluxUpcm_2s_1sr_1);
	}

	return resu;
}


EFlux& EFlux::operator+=(const EFlux& vB)
{
// The grids are supposed to be the same

	Log::mD<<"We do a flux internal addition (+=)!!!"<<endl;
	// We do strong assert: not only it is the same grid, but the same pointer
	assert(mpElecBotEeV==vB.mpElecBotEeV);
	assert(mpElecCentEeV==vB.mpElecCentEeV);
	assert(mpElecDdengeV==vB.mpElecDdengeV);
	assert(mpGAngle==vB.mpGAngle);



	// We add the electron production by altitude
	mEleccm_3s_1=AddOrFullVec(mEleccm_3s_1,vB.mEleccm_3s_1);

	// We add the electron production by altitude and energy
	mElecEcm_3s_1eV_1=AddOrFullVec2(mElecEcm_3s_1eV_1,vB.mElecEcm_3s_1eV_1);
	
	// We add the electron flux function of altitude and energy
	mFluxcm_2s_1eV_1=AddOrFullVec2(mFluxcm_2s_1eV_1,vB.mFluxcm_2s_1eV_1);
	// We add the electron flux function of altitude, energy, and angle
	mFluxAcm_2s_1eV_1sr_1=AddOrFullVec3(mFluxAcm_2s_1eV_1sr_1,vB.mFluxAcm_2s_1eV_1sr_1);

//	resu.mPrecipitationDefined=false; // Default

	// No modifications!
//	if(mPrecipitationDefined&&!vB.mPrecipitationDefined)
//	{
//		mPrecipitationDefined=true;
//		mPrecipFluxDowncm_2s_1sr_1=vA.mPrecipFluxDowncm_2s_1sr_1;
//		mPrecipFluxUpcm_2s_1sr_1=vA.mPrecipFluxUpcm_2s_1sr_1;
//	}

	if(!mPrecipitationDefined&&vB.mPrecipitationDefined)
	{
		mPrecipitationDefined=true;
		mPrecipFluxDowncm_2s_1sr_1=vB.mPrecipFluxDowncm_2s_1sr_1;
		mPrecipFluxUpcm_2s_1sr_1=vB.mPrecipFluxUpcm_2s_1sr_1;
	}
	if(mPrecipitationDefined&&vB.mPrecipitationDefined)
	{
		mPrecipitationDefined=true;
		mPrecipFluxDowncm_2s_1sr_1=(mPrecipFluxDowncm_2s_1sr_1+vB.mPrecipFluxDowncm_2s_1sr_1);
		mPrecipFluxUpcm_2s_1sr_1=(mPrecipFluxUpcm_2s_1sr_1+vB.mPrecipFluxUpcm_2s_1sr_1);
	}

	return *this;


}


EFlux::EFlux(std::deque< EFlux* > vEFlux,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath)
{
	if(vEFlux.size()==0)
	{
		return;
	}
	EFlux* f1=vEFlux[0];
	mpGAngle=f1->mpGAngle;
	mpParameter=f1->mpParameter;
	mpElecBotEeV=f1->mpElecBotEeV;
	mpElecCentEeV=f1->mpElecCentEeV;
	mpElecDdengeV=f1->mpElecDdengeV;

	size_t nbf=vEFlux.size();

	// Initialization bend
	//	mElecEcm_3s_1eV_1 (matrix)
	if(vEFlux[0]->mElecEcm_3s_1eV_1.size1()>0)
	{
		Log::mD<<"Bending elecE"<<endl;
		//
		size_t size1=vEFlux[0]->mElecEcm_3s_1eV_1.size1();
		size_t size2=vEFlux[0]->mElecEcm_3s_1eV_1.size2();
		DequeUblasMat tmpeeev;
		for(unsigned i=0;i<nbf;++i)
		{
			assert(size1==vEFlux[i]->mElecEcm_3s_1eV_1.size1());
			assert(size2==vEFlux[i]->mElecEcm_3s_1eV_1.size2());
			tmpeeev.push_back(&(vEFlux[i]->mElecEcm_3s_1eV_1)); 
		}
		mElecEcm_3s_1eV_1=MathFunction::MatIntLogPath(vAltitudesKm,vSZADegree,tmpeeev,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	//	mEleccm_3s_1 (vector)
	if(vEFlux[0]->mEleccm_3s_1.size()>0)
	{
		//
		Log::mD<<"Bending elec"<<endl;
		DequeUblas tmpee;
		for(unsigned i=0;i<nbf;++i)
		{
			MinValue((vEFlux[i]->mEleccm_3s_1),1E-42);
			tmpee.push_back(&(vEFlux[i]->mEleccm_3s_1)); 
		}
		mEleccm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmpee,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	//	mFluxcm_2s_1eV_1 (matrix)

	if(vEFlux[0]->mFluxcm_2s_1eV_1.size1()>0)
	{
		Log::mD<<"Bending fluxcm2"<<endl;
		DequeUblasMat tmpfev;
		for(unsigned i=0;i<nbf;++i)
		{
			tmpfev.push_back(&(vEFlux[i]->mFluxcm_2s_1eV_1)); 
		}
		mFluxcm_2s_1eV_1=MathFunction::MatIntLogPath(vAltitudesKm,vSZADegree,tmpfev,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	if(vEFlux[0]->mFluxAcm_2s_1eV_1sr_1.size()>0)
	{
		//

		//	mFluxAcm_2s_1eV_1sr_1 (vectormatrix)
		Log::mD<<"Bending fluxAcm2"<<endl;
		DequeVecMat tmpfaev;
		for(unsigned i=0;i<nbf;++i)
		{
			tmpfaev.push_back(&(vEFlux[i]->mFluxAcm_2s_1eV_1sr_1)); 
		}
		mFluxAcm_2s_1eV_1sr_1=MathFunction::VecMatIntLogPath(vAltitudesKm,vSZADegree,tmpfaev,vPath.mAltitudeKm,vPath.mSZADegree);
	}


	// The precipitations. Theoretically, not here, but we copy f1
	mPrecipitationDefined=f1->mPrecipitationDefined;
	mPrecipFluxDowncm_2s_1sr_1=f1->mPrecipFluxDowncm_2s_1sr_1;
	mPrecipFluxUpcm_2s_1sr_1=f1->mPrecipFluxUpcm_2s_1sr_1;
}




double EFlux::FluxEnergy(const ublas::vector<double> & vAltGridKm)
{
	double totenergy=0;// Energy in eV in one cm2 column

	unsigned nbenerg=mpElecCentEeV->size();

	if(mElecEcm_3s_1eV_1.size1()!=0)
	{// We compute the energy over the different altitudes
		unsigned nalt=mElecEcm_3s_1eV_1.size1();
		assert(nalt==vAltGridKm.size());
		assert(nbenerg==mElecEcm_3s_1eV_1.size2());
		ublas::vector<double> altenergy(nalt);
		altenergy.clear();
		for(unsigned i=0;i<nalt;++i)
		{
			for(unsigned j=0;j<nbenerg;++j)
			{
				altenergy[i]+=mElecEcm_3s_1eV_1(i,j)*(*mpElecCentEeV)[j]*(*mpElecDdengeV)[j];
			}
		}

		totenergy=MathFunction::TrapzInt(vAltGridKm,altenergy)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
	}

	
	if(mPrecipFluxDowncm_2s_1sr_1.size1()!=0)
	{
		//cout<<"Energie totale : "<<endl;
		//cout<<mPrecipFluxDowncm_2s_1sr_1<<endl;
		/*
		for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
		{
			for(unsigned j=0;j<nbenerg;++j)
			{
				Log::mL<<"precip\t"<<j<<"\t"<<k<<"\t"<<mPrecipFluxDowncm_2s_1sr_1(j,k)<<endl;
			}
		}*/
		for(unsigned j=0;j<nbenerg;++j)
		{
			for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			{// We did not inverted the angle orientation -> I need to verify this
				//totenergy+=mPrecipFluxDowncm_2s_1sr_1(j,k)*( mpGAngle->mXmu[k]*mpGAngle->mAngzb[k] )*(*mpElecCentEeV)[j];//*(*mpElecDdengeV)[j];
				totenergy+=mPrecipFluxDowncm_2s_1sr_1(j,k)*( mpGAngle->mXmu[k]*mpGAngle->mWeight[k] )*(*mpElecCentEeV)[j]*(*mpElecDdengeV)[j]*2*PI;
			}
		}
	}
	return totenergy;
}


void EFlux::PrintPrecip(std::string vFilename)
{

	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Precipitating flux "<<endl;
	of<<"# Energy in eV "<<endl;
	for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
		of<<"# Flux depending on the angle number"<<k<<endl;
	
	for(unsigned i=0;i<mPrecipFluxDowncm_2s_1sr_1.size1();++i)
	{
		for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			of<<mPrecipFluxDowncm_2s_1sr_1(i,k)<<"\t";
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
	

}



void EFlux::PrintEnergyFlux(std::string vFilename,const ublas::vector<double>& vAltGridKm,bool vIsLength)
{
	unsigned nen=mpElecDdengeV->size();
	unsigned nbang=mpGAngle->mNbAngles;

	if(FileExists(vFilename))
	{
	//	Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	if(vIsLength)
	{
		of<<"# Energy flux in function of length (bent atmosphere) "<<endl;
		of<<"# Length in km "<<endl;
	}else
	{
		of<<"# Energy flux in function of altitude "<<endl;
		of<<"# Altitude in km "<<endl;
	}
	of<<"# Energy Flux in eV/cm2/s at each altitude"<<endl;

	for(unsigned i=0;i<vAltGridKm.size();++i)
	{
		double energy=0.;
		for(unsigned j=0;j<nen;++j)
		{
			for(unsigned k=0;k<nbang;++k)
			{ // correction 23 november 2011 : 2 Pi
				energy+= 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(i)(j,k)* (*mpElecDdengeV)[j];
			}
		}	

		of<<vAltGridKm[i]<<"\t"<<energy<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}



void EFlux::PrintElecProfile(std::string vFilename,const ublas::vector<double>& vAltGridKm)
{
	unsigned nen=mpElecDdengeV->size();

	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Energy flux in function of altitude "<<endl;
	of<<"# Altitude in km "<<endl;
	for(unsigned j=0;j<nen;++j)
	{
		of<<"# Number of electrons in the "<<(*mpElecCentEeV)[j]<< " eV box"<<endl;
	}
	for(unsigned i=0;i<vAltGridKm.size();++i)
	{	of<<vAltGridKm[i]<<"\t";
		for(unsigned j=0;j<nen;++j)
		{
			of<<mElecEcm_3s_1eV_1(i,j)<<"\t";
		}	

		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
	vFilename+="primelec.dat";
	ofstream ofs(vFilename.c_str());
	ofs.precision(9);
	ofs.setf(ios::scientific);

	ofs<<"# Energy flux in function of altitude "<<endl;
	ofs<<"# Altitude in km "<<endl;
		ofs<<"# Number of electrons "<<endl;
	for(unsigned i=0;i<vAltGridKm.size();++i)
	{	ofs<<vAltGridKm[i]<<"\t";
		ofs<<mEleccm_3s_1(i)<<"\t";
		ofs<<endl;
	}

	/// Very important: display the  credits!
	ofs<<Log::msMessageLog<<endl;
	of.close();




}



void EFlux::AddAnisotropicFlux(ublas::vector< ublas::matrix<double> > vFluxcm_3s_1eV_1sr_1,std::deque<Specie*> vSp)
{
	unsigned nbe=vFluxcm_3s_1eV_1sr_1.size();
	if(nbe==0)
		return;
	unsigned nbalt=vFluxcm_3s_1eV_1sr_1[0].size1();
	unsigned nbang=vFluxcm_3s_1eV_1sr_1[0].size2();
	ublas::matrix<double> neleccm_3s_1eV_1(nbalt,nbe);
	ublas::vector<double> neleccm_3s_1(nbalt);
	neleccm_3s_1eV_1.clear();
	neleccm_3s_1.clear();
	
	//Log::SetPriority(Log::DEBUGG);
	Log::mD<<"Initialization new mElec"<<endl;
	for(unsigned i=0;i<nbalt;++i)
	{
		for(unsigned e=0;e<nbe;++e)
		{
			for(unsigned k=0;k<nbang;++k)
			{
				neleccm_3s_1eV_1(i,e)+=vFluxcm_3s_1eV_1sr_1(e)(i,k)*mpGAngle->mWeight[k]*2*PI;// 2 PI: sr to integrated value* weight
			}
			//neleccm_3s_1(i)+=neleccm_3s_1eV_1(i,e)*(*mpElecCentEeV)[e]*(*mpElecDdengeV)[e];
			neleccm_3s_1(i)+=neleccm_3s_1eV_1(i,e)*(*mpElecDdengeV)[e];
		}
	}
	Log::mD<<"We implement them"<<endl;
	mElecEcm_3s_1eV_1=AddOrFullVec2(mElecEcm_3s_1eV_1,neleccm_3s_1eV_1);
	mEleccm_3s_1=AddOrFullVec(mEleccm_3s_1,neleccm_3s_1);

	ublas::matrix<double> densig(nbalt,nbe);
	densig.clear();


	Log::mD<<"Densig computation"<<endl;
	for(unsigned i=0;i<vSp.size();++i)
	{
		if(!vSp.at(i)->CheckElecCrs())
		{
			//Log::SetPriority(Log::WARNING);

			mWarningMissingCrs=true;
			mMissingCrsMsg.push_back("Lack of electrons cross section for your specie "+vSp[i]->mName);
			Log::mW<<"WARNING CROSS SECTION"<<endl;
			Log::mW<<"Lack of electrons cross section for your specie "+vSp[i]->mName<<endl;

		}else
		{
			for(unsigned j=0;j<nbalt;++j)
			{
				for(unsigned k=0;k<nbe;++k)
				{
					densig(j,k)+=((*(vSp[i]->mpElecCrs)).mCrsCincm2[k] + (*(vSp[i]->mpElecCrs)).mElasticCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
				}
			}
		}
	}// Ok, we finish  the computation of the densig -> absorption/cm


	
	ublas::matrix<double> newFluxcm_2s_1eV_1(nbalt,nbe);
	newFluxcm_2s_1eV_1.clear();

	ublas::vector< ublas::matrix<double> >  newFluxAcm_2s_1eV_1sr_1(nbalt);

	Log::mD<<"Loop for the fluxes cm_2"<<endl;
	for(unsigned j=0;j<nbalt;++j)
	{
		newFluxAcm_2s_1eV_1sr_1(j).resize(nbe,nbang);
		newFluxAcm_2s_1eV_1sr_1(j).clear();

		for(unsigned k=0;k<nbe;++k)
		{
			if(densig(j,k)!=0)
			{
				newFluxcm_2s_1eV_1(j,k)=neleccm_3s_1eV_1(j,k)/densig(j,k);
				if(isnan(newFluxcm_2s_1eV_1(j,k)))
				{
					Log::mD<< newFluxcm_2s_1eV_1(j,k)<<"  =  "<<neleccm_3s_1eV_1(j,k)<<" /  --->"<<densig(j,k);
					Error err("Overflow","Init isotropic","autre flux en nan");
					throw err;
				}
			}else
			{
				//Log::SetPriority(Log::WARNING);
				Log::mW<<"WARNING, densig equal to 0 here"<<endl;
				newFluxcm_2s_1eV_1(j,k)=0;
				if(!mWarningMissingCrs)
				{
					Log::mW<<"WARNING, densig equal to 0 here"<<endl;
					mWarningMissingCrs=true;
					mMissingCrsMsg.push_back("Warning, densig equal to 0, check your cross sections for the species");
				}
			}
			//			Log::mL<<"altitude n "<<j<<" energy :"<<(*mpElecCentEeV)[k]<<endl;
			for(unsigned ang=0;ang<nbang;++ang)
			{// We divide by 4Pi because we have a production -> production/sr

				if(densig(j,k)!=0)
				{
					newFluxAcm_2s_1eV_1sr_1(j)(k,ang)= vFluxcm_3s_1eV_1sr_1(k)(j,ang) /densig(j,k);// the equivalent to the 4pi division as already be done
					if(isnan(newFluxAcm_2s_1eV_1sr_1(j)(k,ang)))
					{
						Error err("Overflow","Init isotropic","mflux en nan");
						throw err;

					}
				}else
				{ // warning already done!
					newFluxAcm_2s_1eV_1sr_1(j)(k,ang)=0;// vFluxcm_3s_1eV_1sr_1(k)(j,ang) /densig(j,k);
				}


			}//Nb if we have had a production by angle, (N electrons produced between angle a and a+theta)
			// Then we should have divided by 2 Pi
		}
	}

	Log::mD<<"Add fluxes in the main flux"<<endl;

	mFluxcm_2s_1eV_1=AddOrFullVec2(mFluxcm_2s_1eV_1,newFluxcm_2s_1eV_1);
	mFluxAcm_2s_1eV_1sr_1=AddOrFullVec3(mFluxAcm_2s_1eV_1sr_1,newFluxAcm_2s_1eV_1sr_1);


}
