/** 
 * \file phflux.cpp
 * \brief Implements the PHFlux class 
 * Copyright G Gronoff Nov 2011
 * Last Modification $Id$
 */




#include "phflux.hpp"
using namespace std;

typedef std::deque< ublas::vector<double>* > DequeUblas;
typedef std::deque< ublas::matrix<double>* > DequeUblasMat;
typedef std::deque< ublas::vector< ublas::matrix<double> >* > DequeVecMat;

PHFlux::PHFlux(XmlParameters* pParam, boost::shared_ptr<MathFunction::GaussianAngle> pAngle, ublas::vector<double>* vpGrideV, ublas::vector<double>* vpWidthGrideV, std::string vParticleName): mpGAngle(pAngle), mpParameter(pParam),  mParticleName(vParticleName), mpGrideV(vpGrideV), mpWidthGrideV(vpWidthGrideV)
{
	mPrecipitationDefined=false;
	mWarningMissingCrs=false;
	mParticlePosition = "";

	if(mParticleName == "proton")
	{
		mParticlePosition = "/aero_main/proton/proton_precip";

	}
	if(mParticleName == "hydrogen")
	{
		mParticlePosition = "/aero_main/proton/hydrogen_precip";
	}
}


void PHFlux::InitVoid(unsigned vNbAlt)
{
	mProdcm_3s_1eV_1.resize(vNbAlt,mpGrideV->size());
	mProdcm_3s_1eV_1.clear();

	mProdcm_3s_1.resize(vNbAlt);
	mProdcm_3s_1.clear();	
	
	mFluxcm_2s_1eV_1.resize(vNbAlt,mpGrideV->size());
	mFluxcm_2s_1eV_1.clear();
	
	mFluxAcm_2s_1eV_1sr_1.resize(vNbAlt);
	for(unsigned i=0;i<vNbAlt;++i)
	{
		mFluxAcm_2s_1eV_1sr_1[i].resize(mpGrideV->size(),mpGAngle->mNbAngles);
		mFluxAcm_2s_1eV_1sr_1[i].clear();
	}
}


void PHFlux::InitIsotropic(ublas::matrix<double> vProdcm_3s_1eV_1, ublas::vector<double> vProdcm_3s_1, std::deque<Specie*> vSp)
{

	unsigned nbalt=vProdcm_3s_1.size();
	assert(nbalt>0);
	unsigned nbener=vProdcm_3s_1eV_1.size2();
	unsigned nbang=mpGAngle->mNbAngles;
	
	// If the specie has no defined cross section, we do not throw an error, but we put a warning flag in the class!
	mProdcm_3s_1eV_1 = vProdcm_3s_1eV_1;
	mProdcm_3s_1 = vProdcm_3s_1;
	
	ublas::matrix<double> densig(nbalt,nbener);
	densig.clear();
	
	// We resize the mFlux, the same way
	mFluxcm_2s_1eV_1.resize(nbalt,nbener);
	mFluxcm_2s_1eV_1.clear();

	// We resize the other mFlux
	mFluxAcm_2s_1eV_1sr_1.resize(nbalt);
	for(unsigned i=0;i<nbalt;++i)
	{
		mFluxAcm_2s_1eV_1sr_1[i].resize(nbener,nbang);
		mFluxAcm_2s_1eV_1sr_1[i].clear();
	}
	
	for(unsigned i=0;i<vSp.size();++i)
	{
		
		if( (mParticleName == "proton" || mParticleName == "hydrogen") && vSp.at(i)->CheckProtCrs())
		{
			if(mParticleName == "proton")
			{
				for(unsigned j=0;j<nbalt;++j)
				{
					for(unsigned k=0;k<nbener;++k)
					{
						densig(j,k)+= ((*(vSp[i]->mpProtCrs)).mTotalCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
					}
				}
			}
			if(mParticleName == "hydrogen")
			{
				for(unsigned j=0;j<nbalt;++j)
				{
					for(unsigned k=0;k<nbener;++k)
					{
						densig(j,k)+= ((*(vSp[i]->mpHCrs)).mTotalCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
					}
				}
			}
		}else
		{

			mWarningMissingCrs=true;
			mMissingCrsMsg.push_back("Lack of " + mParticleName + " cross section for your specie "+vSp[i]->mName);
			Log::mW<<"WARNING CROSS SECTION"<<endl;
			Log::mW<<"Lack of  " + mParticleName + " cross section for your specie "+vSp[i]->mName<<endl;
		}
	}// Ok, we finish  the computation of the densig -> absorption/cm


	for(unsigned j=0;j<nbalt;++j)
	{
		for(unsigned k=0;k<nbener;++k)
		{
			if(densig(j,k)!=0)
			{
				mFluxcm_2s_1eV_1(j,k)=mProdcm_3s_1eV_1(j,k)/densig(j,k);
				if(isnan(mFluxcm_2s_1eV_1(j,k)))
				{
					Error err("Overflow","PHFLUX Init isotropic","mFluxcm_2s_1eV_1(j,k) is nan");
					throw err;
				}
			}else
			{// really easier here
				Log::mW<<"WARNING, densig equal to 0 here"<<endl;
				mFluxcm_2s_1eV_1(j,k)=0;
				if(!mWarningMissingCrs)
				{
					Log::mW<<"WARNING, densig equal to 0 here"<<endl;
					mWarningMissingCrs=true;
					mMissingCrsMsg.push_back("Warning, densig equal to 0, check your cross sections for the species");
				}
			}
			for(unsigned ang=0;ang<nbang;++ang)
			{// We divide by 4Pi because we have a production -> production/sr
				mFluxAcm_2s_1eV_1sr_1(j)(k,ang)=mFluxcm_2s_1eV_1(j,k)/(4.*PI);
		
				if(isnan(mFluxAcm_2s_1eV_1sr_1(j)(k,ang)))
				{
					Log::mD<< mFluxAcm_2s_1eV_1sr_1[j](k,ang)<<" "<<mFluxcm_2s_1eV_1(j,k)/(4.*PI)<<endl;
					Error err("Overflow","Init isotropic","mFluxAcm_2s_1eV_1sr_1(j)(k,ang) is nan");
					throw err;
				
				}


			}//Nb if we have had a production by angle, (N electrons produced between angle a and a+theta)
			// Then we should have divided by 2 Pi
		}
	}
}



void PHFlux::ReadPrecipitation(const ublas::vector<double> & vAltitudes, std::deque<Specie*> vSp)
{
	mPrecipitationDefined=true;

	unsigned nang=mpGAngle->mNbAngles;
	unsigned nang_2=nang/2;//mpGAngle->mNbAngles/2;
	unsigned nener = (*mpGrideV).size();
	//	vector<double> tmpvec(nang_2,0.);
	mPrecipFluxDowncm_2s_1sr_1.resize(nener,nang_2);
	mPrecipFluxUpcm_2s_1sr_1.resize(nener,nang_2);
	mPrecipFluxDowncm_2s_1sr_1.clear();
	mPrecipFluxUpcm_2s_1sr_1.clear();


	//Log::SetPriority(Log::CONFIG);
	if(!mpParameter->Exists(mParticlePosition + "/use_precipitation"))
	{
		Log::mL<<"No precipitations for "<<mParticleName<<endl;
		return;// No precipitations
	}
	mpParameter->ExistsOrDie(mParticlePosition + "/precipitation","You want to use precipitations, YOU HAVE TO DEFINE THEM");

	if(mpParameter->Exists(mParticlePosition + "/precipitation/measurement"))
	{
		Log::mL<<"We use "<<mParticleName <<" precipitation measurement"<<endl;

		ReadMeasuredPrecipitation();
		return;
	}

	mpParameter->ExistsOrDie(mParticlePosition + "/precipitation/use_model","You don't want to use measured precipitation, it's your right, but then you have to give me a precipitation model... or not to use precipitations.");
	unsigned type=0;

	mpParameter->GetNKey(mParticlePosition + "/precipitation/use_model","type",type);
	unsigned isotro=0, powlaw=0;
	double entot=0.,E0=0.;
	mpParameter->GetValue(mParticlePosition + "/precipitation/use_model/E0",E0);
	mpParameter->GetValue(mParticlePosition + "/precipitation/use_model/entot",entot);
	mpParameter->GetValue(mParticlePosition + "/precipitation/use_model/powlaw",powlaw);
	mpParameter->GetValue(mParticlePosition + "/precipitation/use_model/isotro",isotro);
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
			MathFlux::InMaxw(entot,E0,b_power_law,nang,*mpGrideV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpGrideV,*mpWidthGrideV);
			break;
		case 2:
			Log::mL<<"Gaussian precipitation"<<endl;
			MathFlux::InGauss(entot,E0,b_power_law,nang,*mpGrideV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpGrideV,*mpWidthGrideV);
			break;
		case 3:
			Log::mL<<"Dirac precipitation"<<endl;
			MathFlux::InDirac(entot,E0,b_power_law,nang,*mpGrideV,isotro,mpGAngle->mXmu,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1);
			Log::mD<<"Normalisation"<<endl;
			MathFlux::NormFlux(entot,mPrecipFluxDowncm_2s_1sr_1,mPrecipFluxUpcm_2s_1sr_1,mpGAngle->mWeight,mpGAngle->mXmu,*mpGrideV,*mpWidthGrideV);
			Log::mD<<"Return"<<endl;
			break;
		default:
			Log::mE<<"Impossible to find your precipitation model"<<endl;

			Error err("PHFlux::ReadPrecipitation","switch type","Impossible to find your precipitation model, type should be between 0 and 4");
			throw err;

	};
	if(mpParameter->Exists(mParticlePosition + "/precipitation/use_model/altitude"))
	{
		double altprecip=vAltitudes[0];
	mpParameter->GetValue(mParticlePosition + "/precipitation/use_model/altitude",altprecip);
		Log::mD<<"You want your precipitations at the position "<<altprecip<<" km "<<endl;
		double distance=nabs(vAltitudes[0] - altprecip);
		unsigned position=0;
		for(unsigned i=1;i<vAltitudes.size();++i)
		{
			if(nabs(vAltitudes[i] - altprecip) < distance)
			{
				position = i;
				distance = nabs(vAltitudes[i] - altprecip);
			}
		}
		Log::mD<<"The closes position is : "<<vAltitudes[position]<<" km"<<endl;



		if(mFluxAcm_2s_1eV_1sr_1.size() != vAltitudes.size())
		{
			Error err("uninitialized Flux", "", "");
			throw err; 
		}		
		
		Log::mD<<"Remplissage du vecteur"<<endl;
		ublas::vector< ublas::matrix<double> > myflux(nener);
		for(unsigned e=0; e < nener; ++e)
		{
			myflux(e).resize(vAltitudes.size(), nang);
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

		//The end, we put the fluxes at 0, and we delete the precipitation definition (it is inside the flux vector!)
		mPrecipitationDefined=false;
		mPrecipFluxDowncm_2s_1sr_1.resize(0,0);
		mPrecipFluxUpcm_2s_1sr_1.resize(0,0);
		Log::mD<<"C'est partit pour l'anisotropique"<<endl;
		AddAnisotropicFlux(myflux,vSp);

	}


}


void PHFlux::ReadMeasuredPrecipitation()
{
	// If no downward...
	mpParameter->ExistsOrDie(mParticlePosition + "/precipitation/measurement/downward/mflux","You have to provide at least one downward flux for your measured precipitation");

	mpParameter->ExistsOrDie(mParticlePosition + "/precipitation/measurement/Egrid","You have to defined the energetic grid for your measured precipitations");

	ublas::vector<double> elec_grid;
	mpParameter->Get1DArray(mParticlePosition + "/precipitation/measurement/Egrid",elec_grid);

	// 1) We read the number of upward and downward values:
	unsigned dwn_number=mpParameter->Numbers(mParticlePosition + "/precipitation/measurement/downward/mflux");
	unsigned upw_number=mpParameter->Numbers(mParticlePosition + "/precipitation/measurement/upward/mflux");
	// 2) We get the values
	vector<TiXmlNode*> dwn_node=mpParameter->GetNodes(mParticlePosition + "/precipitation/measurement/downward/mflux");
	vector<TiXmlNode*> upw_node=mpParameter->GetNodes(mParticlePosition + "/precipitation/measurement/upward/mflux");

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
	//	Log::mD<<"tmp angle : "<<tmp_angle<<endl;
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

void PHFlux::RedistributeFlux(ublas::vector<double> vGrideV,
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
		ublas::vector<double> interp_flux=MathFunction::IntLog(vGrideV, vInputFluxcm_2sr_1s_1[0], *mpGrideV);

		// The flux has the good size	
		//	rFluxcm_2sr_1s_1.resize(mpGAngle->mNbAngles/2);
		assert(interp_flux.size() == mpGrideV->size());
		assert(rFluxcm_2sr_1s_1.size1() == interp_flux.size());
		
		for(unsigned i=0; i < rFluxcm_2sr_1s_1.size1(); ++i)
		{
			for(unsigned j=0; j < rFluxcm_2sr_1s_1.size2(); ++j)
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
		flux_redist[vAnglesDegree[i]]=MathFunction::IntLog(vGrideV,vInputFluxcm_2sr_1s_1[i],*mpGrideV);
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

	for(unsigned i=0;i<(*mpGrideV).size();++i)
	{
		std::deque<double> tmp;
		for(map< double,ublas::vector<double> >::iterator it=flux_redist.begin();it!=flux_redist.end();++it)
		{
			tmp.push_back(((*it).second)[i]);
		}
		ublas::vector<double> tmp2(tmp.size());
		std::copy(tmp.begin(),tmp.end(),tmp2.begin());
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

PHFlux operator+(const PHFlux& vA , const PHFlux& vB)
{// The grids are supposed to be the same

	Log::mD<<"We do a PHflux addition!!!"<<endl;
	// We do strong assert: not only it is the same grid, but the same pointer
	assert(vA.mpGrideV==vB.mpGrideV);
	assert(vA.mpWidthGrideV==vB.mpWidthGrideV);
	// It works because the function is friend
	assert(vA.mpGAngle==vB.mpGAngle);

	assert(vA.mParticleName == vB.mParticleName);
	
	assert(vA.mParticlePosition == vB.mParticlePosition);


	PHFlux resu(vA.mpParameter, vA.mpGAngle, vA.mpGrideV, vA.mpWidthGrideV, vA.mParticleName);

	// We add the particle production by altitude
	resu.mProdcm_3s_1 = AddOrFullVec(vA.mProdcm_3s_1, vB.mProdcm_3s_1);

	// We add the particle production by altitude and energy
	resu.mProdcm_3s_1eV_1 = AddOrFullVec2(vA.mProdcm_3s_1eV_1, vB.mProdcm_3s_1eV_1);
	
	// We add the particle flux function of altitude and energy
	resu.mFluxcm_2s_1eV_1 = AddOrFullVec2(vB.mFluxcm_2s_1eV_1, vB.mFluxcm_2s_1eV_1);
	// We add the particle flux function of altitude, energy, and angle
	resu.mFluxAcm_2s_1eV_1sr_1 = AddOrFullVec3(vB.mFluxAcm_2s_1eV_1sr_1, vB.mFluxAcm_2s_1eV_1sr_1);

	resu.mPrecipitationDefined = false; // Default

	if(vA.mPrecipitationDefined && !vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined = true;
		resu.mPrecipFluxDowncm_2s_1sr_1 = vA.mPrecipFluxDowncm_2s_1sr_1;
		resu.mPrecipFluxUpcm_2s_1sr_1 = vA.mPrecipFluxUpcm_2s_1sr_1;
	}

	if(!vA.mPrecipitationDefined && vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined = true;
		resu.mPrecipFluxDowncm_2s_1sr_1 = vB.mPrecipFluxDowncm_2s_1sr_1;
		resu.mPrecipFluxUpcm_2s_1sr_1 = vB.mPrecipFluxUpcm_2s_1sr_1;
	}
	if(vA.mPrecipitationDefined && vB.mPrecipitationDefined)
	{
		resu.mPrecipitationDefined = true;
		resu.mPrecipFluxDowncm_2s_1sr_1 = (vA.mPrecipFluxDowncm_2s_1sr_1+vB.mPrecipFluxDowncm_2s_1sr_1);
		resu.mPrecipFluxUpcm_2s_1sr_1 = (vA.mPrecipFluxUpcm_2s_1sr_1+vB.mPrecipFluxUpcm_2s_1sr_1);
	}
	return resu;
}


PHFlux& PHFlux::operator+=(const PHFlux& vB)
{
// The grids are supposed to be the same

	Log::mD<<"We do a flux internal addition (+=)!!!"<<endl;
	// We do strong assert: not only it is the same grid, but the same pointer
	assert(mParticleName==vB.mParticleName);
	assert(mParticlePosition == vB.mParticlePosition);
	assert(mpGrideV==vB.mpGrideV);
	assert(mpWidthGrideV==vB.mpWidthGrideV);
	assert(mpGAngle==vB.mpGAngle);



	// We add the particle production by altitude
	mProdcm_3s_1 = AddOrFullVec(mProdcm_3s_1, vB.mProdcm_3s_1);

	// We add the particle production by altitude and energy
	mProdcm_3s_1eV_1 = AddOrFullVec2(mProdcm_3s_1eV_1, vB.mProdcm_3s_1eV_1);
	
	// We add the particle flux function of altitude and energy
	mFluxcm_2s_1eV_1 = AddOrFullVec2(mFluxcm_2s_1eV_1, vB.mFluxcm_2s_1eV_1);
	// We add the particle flux function of altitude, energy, and angle
	mFluxAcm_2s_1eV_1sr_1 = AddOrFullVec3(mFluxAcm_2s_1eV_1sr_1, vB.mFluxAcm_2s_1eV_1sr_1);

	if(!mPrecipitationDefined && vB.mPrecipitationDefined)
	{
		mPrecipitationDefined = true;
		mPrecipFluxDowncm_2s_1sr_1 = vB.mPrecipFluxDowncm_2s_1sr_1;
		mPrecipFluxUpcm_2s_1sr_1 = vB.mPrecipFluxUpcm_2s_1sr_1;
	}
	if(mPrecipitationDefined && vB.mPrecipitationDefined)
	{
		mPrecipitationDefined = true;
		mPrecipFluxDowncm_2s_1sr_1 = (mPrecipFluxDowncm_2s_1sr_1+vB.mPrecipFluxDowncm_2s_1sr_1);
		mPrecipFluxUpcm_2s_1sr_1 = (mPrecipFluxUpcm_2s_1sr_1+vB.mPrecipFluxUpcm_2s_1sr_1);
	}

	return *this;
}


PHFlux::PHFlux(std::deque< PHFlux* > vPHFlux,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath)
{
	if(vPHFlux.size()==0)
	{
		return;
	}
	PHFlux* f1=vPHFlux[0];
	mpGAngle=f1->mpGAngle;
	mpParameter=f1->mpParameter;
	mpGrideV=f1->mpGrideV;
	mpWidthGrideV=f1->mpWidthGrideV;
	mParticleName=f1->mParticleName;
	mParticlePosition=f1->mParticlePosition;
//	mpElecDdengeV=f1->mpElecDdengeV;

	size_t nbf=vPHFlux.size();

	// Initialization bend
	//	mProdcm_3s_1eV_1 (matrix)
	if(vPHFlux[0]->mProdcm_3s_1eV_1.size1()>0)
	{
		Log::mD<<"Bending Prod matrix"<<endl;
		//
		size_t size1=vPHFlux[0]->mProdcm_3s_1eV_1.size1();
		size_t size2=vPHFlux[0]->mProdcm_3s_1eV_1.size2();
		DequeUblasMat tmpeeev;
		for(unsigned i=0;i<nbf;++i)
		{
			assert(size1==vPHFlux[i]->mProdcm_3s_1eV_1.size1());
			assert(size2==vPHFlux[i]->mProdcm_3s_1eV_1.size2());
			tmpeeev.push_back(&(vPHFlux[i]->mProdcm_3s_1eV_1)); 
		}
		mProdcm_3s_1eV_1 = MathFunction::MatIntLogPath(vAltitudesKm,vSZADegree,tmpeeev,vPath.mAltitudeKm,vPath.mSZADegree);
	}
	//	mEleccm_3s_1 (vector)
	if(vPHFlux[0]->mProdcm_3s_1.size()>0)
	{
		//
		Log::mD<<"Bending prod"<<endl;
		DequeUblas tmpee;
		for(unsigned i=0;i<nbf;++i)
		{
			MinValue((vPHFlux[i]->mProdcm_3s_1),1E-42);
			tmpee.push_back(&(vPHFlux[i]->mProdcm_3s_1)); 
		}
		mProdcm_3s_1=MathFunction::IntLogPath(vAltitudesKm,vSZADegree,tmpee,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	//	mFluxcm_2s_1eV_1 (matrix)

	if(vPHFlux[0]->mFluxcm_2s_1eV_1.size1()>0)
	{
		Log::mD<<"Bending fluxcm2"<<endl;
		DequeUblasMat tmpfev;
		for(unsigned i=0;i<nbf;++i)
		{
			tmpfev.push_back(&(vPHFlux[i]->mFluxcm_2s_1eV_1)); 
		}
		mFluxcm_2s_1eV_1=MathFunction::MatIntLogPath(vAltitudesKm,vSZADegree,tmpfev,vPath.mAltitudeKm,vPath.mSZADegree);
	}

	if(vPHFlux[0]->mFluxAcm_2s_1eV_1sr_1.size()>0)
	{
		//

		//	mFluxAcm_2s_1eV_1sr_1 (vectormatrix)
		Log::mD<<"Bending fluxAcm2"<<endl;
		DequeVecMat tmpfaev;
		for(unsigned i=0;i<nbf;++i)
		{
			tmpfaev.push_back(&(vPHFlux[i]->mFluxAcm_2s_1eV_1sr_1)); 
		}
		mFluxAcm_2s_1eV_1sr_1=MathFunction::VecMatIntLogPath(vAltitudesKm,vSZADegree,tmpfaev,vPath.mAltitudeKm,vPath.mSZADegree);
	}


	// The precipitations. Theoretically, not here, but we copy f1
	mPrecipitationDefined=f1->mPrecipitationDefined;
	mPrecipFluxDowncm_2s_1sr_1=f1->mPrecipFluxDowncm_2s_1sr_1;
	mPrecipFluxUpcm_2s_1sr_1=f1->mPrecipFluxUpcm_2s_1sr_1;
}




double PHFlux::FluxEnergy(const ublas::vector<double> & vAltGridKm)
{
	double totenergy=0;// Energy in eV in one cm2 column

	unsigned nbenerg=mpGrideV->size();

	if(mProdcm_3s_1eV_1.size1()!=0)
	{// We compute the energy over the different altitudes
		unsigned nalt=mProdcm_3s_1eV_1.size1();
		assert(nalt==vAltGridKm.size());
		assert(nbenerg==mProdcm_3s_1eV_1.size2());
		ublas::vector<double> altenergy(nalt);
		altenergy.clear();
		for(unsigned i=0;i<nalt;++i)
		{
			for(unsigned j=0;j<nbenerg;++j)
			{
				altenergy[i]+=mProdcm_3s_1eV_1(i,j)*(*mpGrideV)[j]*(*mpWidthGrideV)[j];
			}
		}

		totenergy=MathFunction::TrapzInt(vAltGridKm,altenergy)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
	}

	
	if(mPrecipFluxDowncm_2s_1sr_1.size1()!=0)
	{
		for(unsigned j=0;j<nbenerg;++j)
		{
			for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			{// We did not inverted the angle orientation -> I need to verify this
				totenergy+=mPrecipFluxDowncm_2s_1sr_1(j,k)*( mpGAngle->mXmu[k]*mpGAngle->mWeight[k] )*(*mpGrideV)[j]*(*mpWidthGrideV)[j]*2*PI;
			}
		}
	}
	return totenergy;
}

double PHFlux::FluxEnergyDown()
{
	double totenergy=0;// Energy in eV in one cm2 column
	unsigned nbenerg=mpGrideV->size();
	if(mPrecipFluxDowncm_2s_1sr_1.size1()!=0)
	{
		for(unsigned j=0;j<nbenerg;++j)
		{
			for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			{
				totenergy+=mPrecipFluxDowncm_2s_1sr_1(j,k)*( mpGAngle->mXmu[k]*mpGAngle->mWeight[k] )*(*mpGrideV)[j]*(*mpWidthGrideV)[j]*2*PI;
			}
		}
	}
	return totenergy;
}

double PHFlux::FluxEnergyUp()
{
	double totenergy=0;// Energy in eV in one cm2 column
	unsigned nbenerg=mpGrideV->size();
	if(mPrecipFluxDowncm_2s_1sr_1.size1()!=0)
	{
		for(unsigned j=0;j<nbenerg;++j)
		{

			for(unsigned k=0;k<static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			{
				totenergy+=mPrecipFluxUpcm_2s_1sr_1(j,k)*( mpGAngle->mXmu[k]*mpGAngle->mWeight[k] )*(*mpGrideV)[j]*(*mpWidthGrideV)[j]*2*PI;
			}
		}
	}
	return totenergy;
}
double PHFlux::FluxEnergyDown2()
{
	double totenergy=0;// Energy in eV in one cm2 column
	unsigned nbenerg=mpGrideV->size();
		for(unsigned j=0;j<nbenerg;++j)
		{
			for(unsigned k=  0; k < static_cast<unsigned>(mpGAngle->mNbAngles/2);++k)
			{
				totenergy += mFluxAcm_2s_1eV_1sr_1[0](j,k)*( nabs(mpGAngle->mXmu[k]*mpGAngle->mWeight[k]) )*(*mpGrideV)[j]*(*mpWidthGrideV)[j] * 2 * PI;
			}
		}
	return totenergy;
}

double PHFlux::FluxEnergyUp2()
{
	double totenergy=0;// Energy in eV in one cm2 column
	unsigned nbenerg=mpGrideV->size();
	for(unsigned j=0;j<nbenerg;++j)
	{
		for(unsigned k= static_cast<unsigned>(mpGAngle->mNbAngles/2); k < static_cast<unsigned>(mpGAngle->mNbAngles); ++k)
		{
			totenergy += mFluxAcm_2s_1eV_1sr_1[0](j,k)*( nabs(mpGAngle->mXmu[k]*mpGAngle->mWeight[k]) )*(*mpGrideV)[j]*(*mpWidthGrideV)[j] * 2 * PI;
		}
	}
	return totenergy;
}

void PHFlux::PrintPrecip(std::string vFilename)
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



void PHFlux::PrintEnergyFlux(std::string vFilename,const ublas::vector<double>& vAltGridKm,bool vIsLength)
{
	unsigned nen=mpGrideV->size();
	unsigned nbang=mpGAngle->mNbAngles;

	if(FileExists(vFilename))
	{
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
			{// Correction 23 november 2011: add 2Pi
				energy+= 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(i)(j,k)* (*mpWidthGrideV)[j];
			}
		}	
		of<<vAltGridKm[i]<<"\t"<<energy<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}

double PHFlux::FluxInt(unsigned vAlt, unsigned vEne, unsigned vOption)
{
	double resu = 0;
	assert(vEne < mpGrideV->size());
	assert(vAlt < mFluxAcm_2s_1eV_1sr_1.size());
	unsigned nbang=mpGAngle->mNbAngles;
	switch(vOption)
	{
		case 0:
			//				precision="down";
			for(unsigned k=0;k<nbang/2;++k)
			{
				resu += 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(vAlt)(vEne,k);
			}
			break;
		case 1:
			//				precision="up";
			for(unsigned k=nbang/2;k<nbang;++k)
			{
				resu += 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(vAlt)(vEne,k);
			}
			break;
		case 2:
			//				precision="total";
			for(unsigned k=0;k<nbang;++k)
			{
				resu += 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(vAlt)(vEne,k);
			}
			break;

		case 3:
			//precision="net";
			for(unsigned k=0;k<nbang / 2;++k)
			{
				resu += 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(vAlt)(vEne,k);
			}
			for(unsigned k=nbang/2;k<nbang;++k)
			{
				resu -= 2*PI*mpGAngle->mWeight[k]*mFluxAcm_2s_1eV_1sr_1(vAlt)(vEne,k);
			}
			break;
		default:
			Error err("PHFlux::FluxInt","Wrong option","The option for fluxes is between 0 and 3.");
			throw err;
	}


	return resu;
}

void PHFlux::PrintProfile(std::string vFilename,const ublas::vector<double>& vAltGridKm)
{
	unsigned nen=mpGrideV->size();

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
		of<<"# Number of " + mParticleName + " in the "<<(*mpGrideV)[j]<< " eV box"<<endl;
	}
	for(unsigned i=0;i<vAltGridKm.size();++i)
	{	of<<vAltGridKm[i]<<"\t";
		for(unsigned j=0;j<nen;++j)
		{
			of<<mProdcm_3s_1eV_1(i,j)<<"\t";
		}	
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
	vFilename+="prim_" + mParticleName + "_.dat";
	ofstream ofs(vFilename.c_str());
	ofs.precision(9);
	ofs.setf(ios::scientific);

	ofs<<"# Energy flux in function of altitude "<<endl;
	ofs<<"# Altitude in km "<<endl;
		ofs<<"# Number of  particle"<<endl;
	for(unsigned i=0;i<vAltGridKm.size();++i)
	{	ofs<<vAltGridKm[i]<<"\t";
		ofs<<mProdcm_3s_1(i)<<"\t";
		ofs<<endl;
	}

	/// Very important: display the  credits!
	ofs<<Log::msMessageLog<<endl;
	ofs.close();
}



void PHFlux::AddAnisotropicFlux(ublas::vector< ublas::matrix<double> > vFluxcm_3s_1eV_1sr_1,std::deque<Specie*> vSp)
{
	unsigned nbe=vFluxcm_3s_1eV_1sr_1.size();
	if(nbe==0)
		return;
	unsigned nbalt=vFluxcm_3s_1eV_1sr_1[0].size1();
	unsigned nbang=vFluxcm_3s_1eV_1sr_1[0].size2();
	ublas::matrix<double> npartcm_3s_1eV_1(nbalt,nbe);
	ublas::vector<double> npartcm_3s_1(nbalt);
	npartcm_3s_1eV_1.clear();
	npartcm_3s_1.clear();
	
	//Log::SetPriority(Log::DEBUGG);
	Log::mD<<"Initialization new mpart"<<endl;
	for(unsigned i=0;i<nbalt;++i)
	{
		for(unsigned e=0;e<nbe;++e)
		{
			for(unsigned k=0;k<nbang;++k)
			{
				npartcm_3s_1eV_1(i,e)+=vFluxcm_3s_1eV_1sr_1(e)(i,k)*mpGAngle->mWeight[k]*2*PI;// 2 PI: sr to integrated value* weight
			}
			npartcm_3s_1(i)+=npartcm_3s_1eV_1(i,e)*(*mpWidthGrideV)[e];
		}
	}
	Log::mD<<"We implement them"<<endl;
	mProdcm_3s_1eV_1=AddOrFullVec2(mProdcm_3s_1eV_1,npartcm_3s_1eV_1);
	mProdcm_3s_1=AddOrFullVec(mProdcm_3s_1,npartcm_3s_1);

	ublas::matrix<double> densig(nbalt,nbe);
	densig.clear();


	Log::mD<<"Densig computation"<<endl;
	for(unsigned i=0;i<vSp.size();++i)
	{
		if( (mParticleName == "proton" || mParticleName == "hydrogen") && vSp.at(i)->CheckProtCrs())
		{
			if(mParticleName == "proton")
			{
				for(unsigned j=0;j<nbalt;++j)
				{
					for(unsigned k=0;k < nbe;++k)
					{
						densig(j,k)+= ((*(vSp[i]->mpProtCrs)).mTotalCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
					}
				}
			}
			if(mParticleName == "hydrogen")
			{
				for(unsigned j=0;j<nbalt;++j)
				{
					for(unsigned k=0;k < nbe;++k)
					{
						densig(j,k)+= ((*(vSp[i]->mpHCrs)).mTotalCrscm2[k]) * (*vSp[i]).mTotDensitycm_3[j];
					}
				}
			}
		}else
		{

			mWarningMissingCrs=true;
			mMissingCrsMsg.push_back("Lack of " + mParticleName + " cross section for your specie "+vSp[i]->mName);
			Log::mW<<"WARNING CROSS SECTION"<<endl;
			Log::mW<<"Lack of  " + mParticleName + " cross section for your specie "+vSp[i]->mName<<endl;
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
				newFluxcm_2s_1eV_1(j,k)=npartcm_3s_1eV_1(j,k)/densig(j,k);
				if(isnan(newFluxcm_2s_1eV_1(j,k)))
				{
					Log::mD<< newFluxcm_2s_1eV_1(j,k)<<"  =  "<<npartcm_3s_1eV_1(j,k)<<" /  --->"<<densig(j,k);
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


			}//Nb if we have had a production by angle, (N particle produced between angle a and a+theta)
			// Then we should have divided by 2 Pi
		}
	}

	Log::mD<<"Add fluxes in the main flux"<<endl;

	mFluxcm_2s_1eV_1=AddOrFullVec2(mFluxcm_2s_1eV_1,newFluxcm_2s_1eV_1);
	mFluxAcm_2s_1eV_1sr_1=AddOrFullVec3(mFluxAcm_2s_1eV_1sr_1,newFluxAcm_2s_1eV_1sr_1);


}

void PHFlux::AnisotropicFluxToAverage()
{// 
	unsigned nalt = mFluxAcm_2s_1eV_1sr_1.size();
	mFluxcm_2s_1eV_1.resize(nalt,mpGrideV->size());
	mFluxcm_2s_1eV_1.clear();
	unsigned nbang=mpGAngle->mNbAngles;

	for(unsigned i = 0 ; i < nalt ; ++i)
	{
		for(unsigned e = 0 ; e < mpGrideV->size() ; ++e)
		{

			for(unsigned ang=0;ang<nbang;++ang)
			{
				mFluxcm_2s_1eV_1(i,e) += 2.*PI* mFluxAcm_2s_1eV_1sr_1(i)(e,ang) * mpGAngle->mWeight[ang];
			}
		}
	}
}


void PHFlux::PrintFluxInt(std::string vFilename, const ublas::vector<double>& vAltGridKm)
{
	unsigned nalt = mFluxcm_2s_1eV_1.size1();
	ofstream ofs(vFilename.c_str());
	ofs.precision(9);
	ofs.setf(ios::scientific);

	ofs<<"# Energy flux in function of altitude "<<endl;
	ofs<<"# Flux"<<endl;
	ofs<<"# Altitude in km "<<endl;
	for(unsigned i=0;i<   mpGrideV->size();++i)
	{	ofs<<"# "<<(*mpGrideV)[i]<<endl;
	}

	
	
	for(unsigned i = 0 ; i < nalt ; ++i)
	{
		ofs<<vAltGridKm[i]<<"\t";
		for(unsigned e = 0 ; e < mpGrideV->size() ; ++e)
		{

				ofs<<mFluxcm_2s_1eV_1(i,e)<<"\t";
		}
		ofs<<endl;
	}
	/// Very important: display the  credits!
	ofs<<Log::msMessageLog<<endl;
	ofs.close();

}






