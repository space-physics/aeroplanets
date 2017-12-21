/**
 * \file photoionization.cpp
 * \brief The implementation of the photoionization class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: photoionization.cpp 1922 2014-02-26 22:38:52Z gronoff $
 */


#include "photoionization.hpp"
using namespace std;


Photoionization::Photoionization(XmlParameters* pParam):mpParameter(pParam)
{
		mbIsSolfluxDefined=false;
}


Photoionization::~Photoionization()
{
	if(mbIsSolfluxDefined)
	{
		delete mpSolflux;
	}
}


void Photoionization::Init(ublas::vector<double>* pPhotonGrideV,ublas::vector<double>* pPhotonGridmineV,ublas::vector<double>* pPhotonGridmaxeV,ublas::vector<double>* pElecBotEeV,ublas::vector<double>* pElecCentEeV,ublas::vector<double>* pElecDdengeV,ublas::vector<double>* vpAltGridKm,double vUA,double vRKm,double vSZADegree)
{
	mpPhotonGreV=pPhotonGrideV;
	mpPhotonGrMineV=pPhotonGridmineV;
	mpPhotonGrMaxeV=pPhotonGridmaxeV;
	mpEBotEeV=pElecBotEeV;
	mpECentEeV=pElecCentEeV;
	mpEDdengeV=pElecDdengeV;
	mpAltGridKm=vpAltGridKm;
	mUA=vUA;
	mSZADegree=vSZADegree;
	mRKm=vRKm;


	mpParameter->ExistsOrDie("/aero_main/sun/model/flux_model","You need to define the model used");
	int type;

	mpParameter->GetNKey("/aero_main/sun/model/flux_model","type",type);

	switch(type)
	{

		case 0:
			Log::mI<<"You use the originial 39 boxes solar flux model"<<endl;
			mpSolflux=new  Solar39Boxes(mpParameter);
			mpSolflux->RetrieveFlux(mUA,*mpPhotonGreV,*mpPhotonGrMineV,*mpPhotonGrMaxeV,mFluxPhcm_2s_1);
			mbIsSolfluxDefined=true;
			break;
		case 1:
			Log::mI<<"You use the user defined solar flux model"<<endl;
			mpSolflux=new  SolarUserDefined(mpParameter);
			mpSolflux->RetrieveFlux(mUA,*mpPhotonGreV,*mpPhotonGrMineV,*mpPhotonGrMaxeV,mFluxPhcm_2s_1);
			mbIsSolfluxDefined=true;
			break;
		case 2:
			Log::mI<<"You use the Trans* solar flux model (39 boxes extended to low energies)"<<endl;
			mpSolflux=new  SolarTransBoxes(mpParameter);
			mpSolflux->RetrieveFlux(mUA,*mpPhotonGreV,*mpPhotonGrMineV,*mpPhotonGrMaxeV,mFluxPhcm_2s_1);
			mbIsSolfluxDefined=true;
			break;
		case 3:
			Log::mI<<"You use the HEUVAC flux model (1.5- 105.5nm; bin 1 nm)"<<endl;
			mpSolflux=new  HeuvacFlux(mpParameter);
			mpSolflux->RetrieveFlux(mUA,*mpPhotonGreV,*mpPhotonGrMineV,*mpPhotonGrMaxeV,mFluxPhcm_2s_1);
			mbIsSolfluxDefined=true;
			break;

		case 4:
			Log::mI<<"You use the HEUVAC flux model (1.5- 165.7nm) GridBoxes does not work"<<endl;
			mpSolflux=new  HeuvacFluxLow(mpParameter);
			mpSolflux->RetrieveFlux(mUA,*mpPhotonGreV,*mpPhotonGrMineV,*mpPhotonGrMaxeV,mFluxPhcm_2s_1);
			mbIsSolfluxDefined=true;
			break;

		default:
			//Log::SetPriority(Log::ERROR);
			Log::mE<<"The sun model type defined, "<<type<<", was not implemented. Please reconsider another flux model"<<endl;
			Error err("Photoionization::Init","sun model not defined","The sun model number"+ntostr(type)+" is not defined");
			throw err;
	}

	if(mpParameter->GetMonteCarlo()&&mpParameter->Exists("/aero_main/sun/model/flux_uncertainty"))
	{
		double uncert=0;
		mpParameter->GetValue("/aero_main/sun/model/flux_uncertainty",uncert);
		Log::mI<<"The uncertainty on your solar flux is set to "<<uncert<<"%, and is active"<<endl;
		for(size_t i=0;i<mFluxPhcm_2s_1.size();++i)
		{
			double sigma=mFluxPhcm_2s_1[i]*uncert/100;
			mFluxPhcm_2s_1[i]=nmax(MathRandom::GetNormal(mFluxPhcm_2s_1[i],sigma),0.);
		}

		if(mpParameter->Exists("/aero_main/sun/model/flux_uncertainty_factor"))
		{
			Log::mW<<"Warning, you have a solar flux uncertainty and a solar flux uncertainty factor. Both will be taken into account!"<<endl;
		}
	}
	if(mpParameter->GetMonteCarlo()&&mpParameter->Exists("/aero_main/sun/model/flux_uncertainty_factor"))
	{
		double uncert=0;
		mpParameter->GetValue("/aero_main/sun/model/flux_uncertainty_factor",uncert);
		Log::mI<<"The uncertainty on your solar flux factor is set to "<<uncert<<"%, and is active"<<endl;
		double uncertfact = MathRandom::GetNormal(1,uncert/100.);
		for(size_t i=0;i<mFluxPhcm_2s_1.size();++i)
		{
			mFluxPhcm_2s_1[i]=nmax(mFluxPhcm_2s_1[i]*uncertfact,0.);
		}

	}
	if(mpParameter->Exists("/aero_main/sun/model/flux_factor"))
	{
		double ffact=0;
		mpParameter->GetValue("/aero_main/sun/model/flux_factor",ffact);
		Log::mI<<"The factor on your solar flux factor is set to "<<ffact<<" and is active"<<endl;
		for(size_t i=0;i<mFluxPhcm_2s_1.size();++i)
		{
			mFluxPhcm_2s_1[i]=nmax(mFluxPhcm_2s_1[i]*ffact,0.);
		}
	}
}


void Photoionization::ChapInit(std::deque<Specie*>& vpSp)
{
	std::deque<Specie*>::iterator it;

	for(it=vpSp.begin();it!=vpSp.end();++it)
	{
		mSpecieToChap[(*it)]=Chap(mSZADegree,*mpAltGridKm,(*it)->mScaleHcm);
	}

}


double Photoionization::Xchap(int iAlt,Specie* vpSp)
{
	return mSpecieToChap[vpSp][iAlt];
}


double Photoionization::Sperfc(double d)
{
        if(d<=8.){
          return (1.0606963+0.55643831*d) /(1.0619896+1.7245609*d+d*d);
	}else
	{
          return 0.56498823/(0.06651874+d);
	} 
	//Log::SetPriority(Log::ERROR);
	Log::mE<<"Fatal error in the function sperfc : check if the parameter is a number != Nan||inf"<<endl;
	Error err("Photoionization::Sperfc","pb","pb impossible!");
	throw err;
}


ublas::vector<double> Photoionization::Chap(double vChiDeg,const ublas::vector<double>& vAltKm,const ublas::vector<double>& vScaleHcm)
{
	ublas::vector<double> chapsp(vScaleHcm.size());

	assert(vAltKm.size()==vScaleHcm.size());
	for(unsigned i=0;i<vScaleHcm.size();++i)
	{
		double hg=(mRKm+vAltKm[i])*1E5/vScaleHcm[i]; // 1E5 : altitude from km to cm
		double coschi=cos(vChiDeg*PI/180);
		double hf=0.5*hg*(coschi*coschi);

		double a = sqrt(0.5*PI*hg);
		double d=Sperfc(sqrt(hf));
		if(vChiDeg<90)
		{
			chapsp[i]=(a*d);
		}else
		{
			double sinchi=sin(vChiDeg*PI/180);
			double c=sqrt(sinchi*exp(hg*(1.-sinchi)));
			double b=c-0.5*d;
			chapsp[i]=(a*b);
		}


	}
	return chapsp;
}




void Photoionization::ComputeTransmission(std::deque<Specie*>& vSp,std::string vOutFile)
{
	if((vSp).size()==0)
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"You must have a defined number of species to photoionize-excite"<<endl;
		Error err("ComputeTransmission","Void species","You must have a defined number of species to photoabsorb");
		throw err;
		return;
	}

	// We initialize the specie_to_density
	// And we initialize each species for photoionization
	
	deque<Specie*>::iterator isp;
	// Initialization of the chapman function
	ChapInit(vSp);
	// Production of electrons
	//ublas::vector<double> proelec(mpAltGridKm->size());
	//proelec.clear();
	//ublas::matrix<double> proelecE(mpAltGridKm->size(),mpECentEeV->size());
	//proelecE.clear();
	ublas::matrix<double> transmission(mpAltGridKm->size(),mFluxPhcm_2s_1.size());
	transmission.clear();
	for(unsigned i=0;i<mFluxPhcm_2s_1.size();++i)
	{
		for(unsigned j=0;j<mpAltGridKm->size();++j)
		{
			double exa=0;
			//double fl=mFluxPhcm_2s_1[i]; // The flux: not needed here!
			for(unsigned k=0;k<vSp.size();++k)
			{
				exa+=Xchap(j,(vSp)[k])*( (vSp)[k]->mpPhotoCrs->mTotalCrscm2[i])* (vSp)[k]->mColDenscm_2[j];
			}

			// Now the Absorption is equal to exp(-exa)
			// Therefore, the Transmission is 1 - exp(- exa)
			transmission(j,i) = exp(-exa);
		}
	}
	if(FileExists(vOutFile))
	{
		Log::mD<<"Low level warning : We overwrite"<<vOutFile<<endl;
	}
	ofstream of(vOutFile.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Transmission for SZA = "<<mSZADegree<<" Degree"<<endl;
	of<<"# Transmission"<<endl;
	of<<"# Altitude (km)"<<endl;

	for(unsigned i=0;i<mFluxPhcm_2s_1.size();++i)
	{
		of<<"# Box "<<(*mpPhotonGrMineV)[i]<<" - "<<(*mpPhotonGrMaxeV)[i]<<endl;
	}
	for(unsigned j=0;j<mpAltGridKm->size();++j)
	{
		of<<(*mpAltGridKm)[j]<<" ";
		for(unsigned i=0;i<mFluxPhcm_2s_1.size();++i)
		{
			of<<transmission(j,i)<<" ";
		}
		of<<endl;
	}

	of<<Log::msMessageLog<<endl;
	of.close();

}

void Photoionization::PhotodissocieSpecies(std::deque<Specie*>& vSp,double vEnergy, unsigned iEv, double vFl)
{
	std::deque<Specie*>::iterator it;
	for(it=vSp.begin();it!=vSp.end();++it)
	{	
		for(unsigned i=0;i<(*it)->mProcessNames.size();++i)
		{
			double threshold=(*it)->mpPhotoCrs->mThresholdseV[i];
			if(vEnergy>threshold)
			{
				double produ= vFl * (*it)->mpPhotoCrs->mCrscm2.at(i)[iEv];
				(*it)->mSpeciesProductionProbabilitys_1.at(i)+=produ;
				double elnumber=((*it)->mpPhotoCrs->mNumberOfElectrons.at(i));
				(*it)->mPhotoElecProductionProbabilitys_1+=produ*elnumber;
			}
		}
		(*it)->mPhotoDissociationProbabilitys_1 += vFl * (*it)->mpPhotoCrs->mTotalCrscm2[iEv];
	}
}

void Photoionization::IonizeSpecies(std::deque<Specie*>& vSp,double vEnergy,unsigned iEv,unsigned iAlt,double vFl,double vExa,ublas::vector<double>& rProelec,ublas::matrix<double> & rProelecE)
{

	std::deque<Specie*>::iterator it;
	double depr=vFl*exp(-vExa);

//	cout<<"Photoionize energy,ialt "<<energy<<", "<<ialt<<endl;
	for(it=vSp.begin();it!=vSp.end();++it)
	{


		//	double prodtot=depr*(*it)->PhotoCrs->TotalCrs[i]*specie_to_density[*it];
		double prod=depr*(*(mSpecieToDensity[*it]))[iAlt];

		for(unsigned i=0;i<(*it)->mProcessNames.size();++i)
		{

			double threshold=(*it)->mpPhotoCrs->mThresholdseV[i];

			if(vEnergy>threshold)
			{
				double produ=prod*(*it)->mpPhotoCrs->mCrscm2.at(i)[iEv];
				(*it)->mSpeciesProductioncm_3s_1.at(i)[iAlt]+=produ;
				double elnumber=((*it)->mpPhotoCrs->mNumberOfElectrons.at(i));
				double elenergy=vEnergy-threshold;
				rProelec[iAlt]+=produ*elnumber;
				(*it)->mPhotoElecProductioncm_3s_1[iAlt]+=produ*elnumber; // The species electron production

				unsigned energypos=Search(elenergy);
				rProelecE(iAlt,energypos)+=produ/(*mpEDdengeV)[energypos];// The energy is considered as ONE electron!
				if((*it)->mpPhotoCrs->mIsAuger.at(i))
				{// If the auger process is defined
					for(unsigned j=0;j<(*it)->mpPhotoCrs->mAugerEnergy.at(i).size();++j)
					{
						unsigned energypos2=Search((*it)->mpPhotoCrs->mAugerEnergy.at(i).at(j));//  We search the position for the auger electron
						rProelecE(iAlt,energypos2)+=produ/(*mpEDdengeV)[energypos2]*(*it)->mpPhotoCrs->mAugerEfficiency.at(i).at(j) ;
						(*it)->mPhotoElecProductioncm_3s_1[iAlt]+=produ/(*mpEDdengeV)[energypos2]*(*it)->mpPhotoCrs->mAugerEfficiency.at(i).at(j);
					}
				}
			}
		}
	}
}

void Photoionization::ComputePhotoionization(std::deque<Specie*>& vSp,
					     std::deque<Specie*>& rResult,
					     EFlux& rResultFlux)
{
	if((vSp).size()==0)
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"You must have a defined number of species to photoionize-excite"<<endl;
		Error err("ComputePhotoionization","Void species","You must have a defined number of species to photoionize-excite");
		throw err;
		return;
	}

	// We initialize the specie_to_density
	// And we initialize each species for photoionization
	
	deque<Specie*>::iterator isp;
	for(isp=vSp.begin();isp!=vSp.end();++isp)
	{
		if(!((*isp)->mpPhotoCrs->mIsDefinedTotalCrs))
		{
			Error err("Photoionization",(*isp)->mName,"Total cross section for this specie not defined");
			throw err;
		}
		mSpecieToDensity[*isp]=&((*isp)->mTotDensitycm_3);
		(*isp)->InitPhotoIonization(mpAltGridKm->size());
	}
	// Initialization of the chapman function
	ChapInit(vSp);

	// Production of electrons
//	std::vector<double> proelec(mpAltGridKm->size(),0);
	ublas::vector<double> proelec(mpAltGridKm->size());
	proelec.clear();
//	std::vector<double> tmp(mpECentEeV->size(),0);
//	std::vector< std::vector<double> > proelecE(mpAltGridKm->size(),tmp);
	ublas::matrix<double> proelecE(mpAltGridKm->size(),mpECentEeV->size());
	proelecE.clear();


	for(unsigned i=0;i<mFluxPhcm_2s_1.size();++i)
	{
		double deV=(*mpPhotonGrMaxeV)[i]-(*mpPhotonGrMineV)[i];

		if(deV<1.)
		{// On a line : we consider the whole ionization
			for(unsigned j=0;j<mpAltGridKm->size();++j)
			{
				double exa=0;
				double fl=mFluxPhcm_2s_1[i]; // The flux
				for(unsigned k=0;k<vSp.size();++k)
				{
					exa+=Xchap(j,(vSp)[k])*( (vSp)[k]->mpPhotoCrs->mTotalCrscm2[i])* (vSp)[k]->mColDenscm_2[j];
				}

				if(exa<35)
				{// evaluation of the different ionizations.

					IonizeSpecies(vSp,(*mpPhotonGreV)[i],i,j,fl,exa,proelec,proelecE);
				}

			}


		}else
		{// On the continuum: we cut the box into several small intervals.
			// The number of intervals
			unsigned ninter= (int)(deV)-1;
			// The width of the intervals
			double deltaE=deV/((double)ninter);
			double fl=mFluxPhcm_2s_1[i]/((double)ninter); // The flux

			for(unsigned u=0;u<ninter;++u)
			{
				double ener=(*mpPhotonGrMineV)[i]+deltaE/2.+double(u)*deltaE;
				for(unsigned j=0;j<mpAltGridKm->size();++j)
				{
					double exa=0;
					for(unsigned k=0;k<vSp.size();++k)
					{
						exa+=Xchap(j,(vSp)[k])*( (vSp)[k]->mpPhotoCrs->mTotalCrscm2[i])* (vSp)[k]->mColDenscm_2[j];
					}

					if(exa<35)
					{// evaluation of the different ionizations.
						IonizeSpecies(vSp,ener,i,j,fl,exa,proelec,proelecE);
					}

				}


			}


		}


	} // Qu'elle est grande cette boucle! This is a big looop

	// We put the electron fluxes into the appropriate object
	rResultFlux.InitIsotropic(proelecE,proelec,vSp);

	// We now transfer the species between the input and resu vector
	// For that, we use a Specie function 

	SpecieUtils::SpeciesToResu(vSp,rResult);

	//(*sp)[0]->print_productions("CO2_photoionization.dat",*alt_grid);

//	cout<<"Result.size : a faire"<<Result->size()<<endl;


}

unsigned Photoionization::Search(double vEnergy)
{
	unsigned size=(*mpEBotEeV).size();
	if(vEnergy<(*mpEBotEeV)[size-1])
	{
		return size-1;
	}
	if(vEnergy>(*mpEBotEeV)[0])
	{
		return 0;
	}
	
	/*for(unsigned i=0;i<size;++i)
	{
		double inf=(*mpEBotEeV)[i];
		double sup=inf+(*mpEDdengeV)[i];

		if(vEnergy>inf && vEnergy<sup)
		{
			return i;
		}
	}

	Log::SetPriority(Log::ERROR,"Photoionization::search");
	Error err("Photoionization search error","Photoionization::search","pb in the energy position search...");
	Log::mL<<"Photoionization search error : does not match"<<endl;
	throw err;
	return 0;*/
	// The divide and conquer algorithm should be more efficient!
	// The idea: the energy grid is strictly decreasing!
		
	unsigned pmin=0;
	unsigned pmax=size-1;
	unsigned pmiddle=0;
//	unsigned pmiddle2=vPosMin+1;
	while(pmax-pmin>0)
	{
/*		if(pmiddle2==pmiddle)
		{
			cout<<pmiddle<<endl<<pmin<<pmax<<endl;
			cout<<"Gros probleme"<<endl;
		}
		pmiddle2=pmiddle;
*/
		if((pmax-pmin)==1)
		{
			if((*mpEBotEeV)[pmax]<=vEnergy)
			{
				return pmax;
			}else
			{
				return pmin;
			}


		}else
		{
			pmiddle=(pmax-pmin)/2+pmin;
			if((*mpEBotEeV)[pmiddle]>vEnergy)
			{
				pmin=pmiddle;
			}else
			{
				pmax=pmiddle;
			}
		}
	}// pmax==pmin
	return pmin;

}



void Photoionization::ComputePhotodissociationBranching(std::deque<Specie*>& vSp, std::string vSuffix)
{
	if((vSp).size()==0)
	{
		Log::mE<<"You must have a defined number of species to photoionize-excite"<<endl;
		Error err("ComputePhotoionization","Void species","You must have a defined number of species to photoionize-excite");
		throw err;
		return;
	}
	// And we initialize each species for photoionization
	deque<Specie*>::iterator isp;
	for(isp=vSp.begin();isp!=vSp.end();++isp)
	{
		if(!((*isp)->mpPhotoCrs->mIsDefinedTotalCrs))
		{
			Error err("Photoionization",(*isp)->mName,"Total cross section for this specie not defined");
			throw err;
		}
		//	mSpecieToDensity[*isp]=&((*isp)->mTotDensitycm_3);
		(*isp)->InitPhotoIonization(mpAltGridKm->size());
	}
	for(unsigned i=0;i<mFluxPhcm_2s_1.size();++i)
	{

		double deV=(*mpPhotonGrMaxeV)[i]-(*mpPhotonGrMineV)[i];
		if(deV<1.)
		{// We go a line
			double fl=mFluxPhcm_2s_1[i]; // The flux
			PhotodissocieSpecies(vSp,(*mpPhotonGreV)[i], i, fl);
		}else
		{
			unsigned ninter= (int)(deV)-1;
			// The width of the intervals
			double deltaE=deV/((double)ninter);
			double fl=mFluxPhcm_2s_1[i]/((double)ninter); // The flux
			for(unsigned u=0;u<ninter;++u)
			{
				double ener=(*mpPhotonGrMineV)[i]+deltaE/2.+double(u)*deltaE;
				// PHOTODISSSPECIES(vSP,ener, fl)
				PhotodissocieSpecies(vSp, ener, i, fl);
			}
		}
	}
	
	
	mpParameter->ExistsOrDie("/aero_main/sun/probability_file","You have to define probability_file for writing your results");
	std::string outFile = mpParameter->Elem("/aero_main/sun/probability_file") + vSuffix ;
	if(FileExists(outFile))
	{
		Log::mD<<"Low level warning : We overwrite"<<outFile<<endl;
	}
	ofstream of(outFile.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"<PhotoProbability>"<<endl;
	for(isp=vSp.begin();isp!=vSp.end();++isp)
	{
		of << (*isp)->ReturnPhotoProbabilityInfo() <<endl;
	}
	of<<"<!-- "<<Log::msMessageLog<<"-->"<<endl;
	of<<"</PhotoProbability>"<<endl;
	of.close();
}









/*void Photoionization::compute_photoionization(std::vector<Specie*>*sp,
					     std::vector<double> altitude,
					     std::map< std::string, vector<double> > sp_density,
					     std::vector<Specie*>* Result,
					     EFlux* result_flux
					     
					     );

*/
