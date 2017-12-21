 /** 
 * \file eleccrs.cpp
 * \brief Implements the electron cross section class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: eleccrs.cpp 1213 2011-01-28 17:56:54Z gronoff $
 *
 */

#include "eleccrs.hpp"
using namespace std;

void ElecCrossSection::PrintCrs(std::string vFileName,std::string vRedisFilename)
{
	// We print the standard filename
	//mIsDefinedTotalCrs=false;
	CrossSection::PrintCrs(vFileName);

	string filename2=vFileName+"-elec-specific.dat";
	if(FileExists(filename2))
	{
		Log::mD<<"Low level warning : We overwrite"<<filename2<<endl;
	}
	ofstream of(filename2.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Cross section for "<<mSpecies<<endl;
	of<<"# Energy in eV"<<endl;
	of<<"# Cross section in cm2"<<endl;
	of<<"# inelastic cross section "<<endl;
	of<<"# elastic cross section "<<endl;



	unsigned nb_e=mCrsCincm2.size();

	for(unsigned i=0;i<nb_e;++i)
	{
		of<<mGrideV[i]<<"\t";
		of<<mCrsCincm2[i]<<"\t"<<mElasticCrscm2[i]<<endl;;
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();



	// We do the same processus for the redistribution
	Log::mL<<"Redistribution file name: "<<vRedisFilename<<endl;
	if(FileExists(vRedisFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vRedisFilename<<endl;
	}

	ofstream orf(vRedisFilename.c_str());
	orf.precision(9);
	orf.setf(ios::scientific);
	orf<<"# Cross section for "<<mSpecies<<endl;
	orf<<"# Energy in eV"<<endl;
	orf<<"# Cross section in cm2"<<endl;
	for(unsigned j=0;j<nb_e;++j)
	{
		orf<<"#\t"<<mGrideV[j]<<" eV"<<endl;
	}
	for(unsigned i=0;i<nb_e;++i)
	{
		orf<<mGrideV[i]<<"\t";
		for(unsigned j=0;j<nb_e;++j)
		{
			orf<<mOutSecondarycm2(i,j)<<"\t";
		}
		orf<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();



}



bool ElecCrossSection::LoadCrs(std::string vFilename,std::string vName,bool bIsMonteCarlo,bool bLogInterp, double vReesAlpha, double vReesBeta, double vReesGamma, bool bIsOpal, double vOpalEbareV)
{
	
	//XmlParameters*param=new XmlParameters(vFilename);
	mpParam=XmlParameters::AttachFileParameter(vFilename);
	mParamLoaded=true;
	mSpecies=vName;
	bool resu=LoadCrs(mpParam,vName,bIsMonteCarlo,bLogInterp, vReesAlpha, vReesBeta, vReesGamma, bIsOpal, vOpalEbareV);
//	delete param;
	return resu;
}
bool ElecCrossSection::LoadCrs(XmlParameters* pParams,std::string vName,bool bIsMonteCarlo, bool bLogInterp, double vReesAlpha, double vReesBeta, double vReesGamma, bool bIsOpal, double vOpalEbareV)
{// We initialize the cross sections
//	Log::SetPriority(Log::INFO);

	mReesAlpha = vReesAlpha;
	mReesBeta = vReesBeta;
	mReesGamma = vReesGamma;
	mbIsOpal = bIsOpal;
	mEbareV = vOpalEbareV;

	Log::mL<<"Load electron impact cross section for"<<vName<<endl;
	string name=trim(vName);
	if(!	CrossSection::LoadCrs(pParams,name,bIsMonteCarlo,bLogInterp))
	{
		return false;
	}
	if(pParams->Exists("/crs/"+name+"/ZeroCrs"))
	{
		mIsDefinedCrs=false;
		return true;
	}
	// We check for the inelastic cross section
	if(!pParams->Exists("/crs/"+name+"/ElasticCrs"))
	{
		Error err("ElecCrossSection::LoadCrs","ElasticCrs","Your elastic cross section has not been defined");
		throw err;
	}
	mbIsDefinedElastic=true;

	bool extrapolate_elastic_crs=true;
	TiXmlNode* elasticnode=pParams->GetNode("/crs/"+name+"/ElasticCrs");
	if(pParams->Exists("/crs/"+name+"/ElasticCrs/NoStdExtrapolate"))
	{
		extrapolate_elastic_crs=false;
	}
	if(bLogInterp)
	{// nb no threshold for elastic cross section
		mElasticCrscm2=ExtractCrsLog(pParams,elasticnode,0.,mMaxDefinedEnergyElasticCrseV);

	}else
	{
		mElasticCrscm2=ExtractCrs(pParams,elasticnode,0.,mMaxDefinedEnergyElasticCrseV);
	}
	// We check if we have excitation cross sections

	if(mbTotalSumIneel)
	{
		mIsDefinedTotalCrs=true;
		mTotalCrscm2.resize(mElasticCrscm2.size());
		mTotalCrscm2.clear();

	}



	//
	// the excitation cross section is like
	// <ExcitationCrs threshold="12.5>
	//  .... (normal crs)
	//  </ExcitationCrs>
	//
	//
	if(pParams->Exists("/crs/"+name+"/ExcitationCrs"))
	{
		TiXmlNode* excnode=pParams->GetNode("/crs/"+name+"/ExcitationCrs");
		double threshold;
		pParams->GetNKey("/crs/"+name+"/ExcitationCrs","threshold",threshold);
		double tmpdef;
		if(bLogInterp)
		{// nb no threshold for elastic cross section
			mExcitationUniqueCrscm2=ExtractCrsLog(pParams,excnode,threshold,tmpdef);

		}else
		{
			mExcitationUniqueCrscm2=ExtractCrs(pParams,excnode,threshold,tmpdef);
		}

		mExcitationCrscm2.erase(mExcitationCrscm2.begin(),mExcitationCrscm2.end());
		mExcitationCrscm2.push_back(&mExcitationUniqueCrscm2);
		if(mbTotalSumIneel)
		{
			mTotalCrscm2+=mExcitationUniqueCrscm2;
		}
		mExcitationThresholdeV.erase(mExcitationThresholdeV.begin(),mExcitationThresholdeV.end());
		mExcitationThresholdeV.push_back(threshold);

	}else
	{
		assert(mExcitationThresholdeV.size()==mExcitationCrsPosition.size());//mExcitationCrscm2.size());

		for(unsigned j=0;j<mExcitationCrsPosition.size();++j)
		{
			unsigned k=mExcitationCrsPosition[j];
			//		Log::mL<<"Super important pour le check de k "<<k<<endl;
			mExcitationCrscm2.push_back(&mCrscm2[k]);
			if(mbTotalSumIneel)
			{
				mTotalCrscm2+=mCrscm2[k];
			}
		}
		if(mExcitationCrscm2.size()==0)
		{// Nothing was defined : error
			Error err("ElecCrossSection::LoadCrs","No excitation cross section","The code was unable to find excitation cross section (even null). It is considered as an error");
			throw err;
		}
	}


	// We check if we have ionization cross sections


	//
	// the ionization cross section is like
	// <IonizationCrs threshold="12.5>
	//  .... (normal crs)
	//  </IonizationCrs>
	//
	//
	if(pParams->Exists("/crs/"+name+"/IonizationCrs"))
	{
		TiXmlNode* excnode=pParams->GetNode("/crs/"+name+"/IonizationCrs");
		double threshold;
		pParams->GetNKey("/crs/"+name+"/IonizationCrs","threshold",threshold);
		double tmpdef;
		if(bLogInterp)
		{// nb no threshold for elastic cross section
			mIonizationUniqueCrscm2=ExtractCrsLog(pParams,excnode,threshold,tmpdef);

		}else
		{
			mIonizationUniqueCrscm2=ExtractCrs(pParams,excnode,threshold,tmpdef);
		}
		if(mbTotalSumIneel)
		{
			mTotalCrscm2+=mIonizationUniqueCrscm2;
		}
		mIonizationCrscm2.erase(mIonizationCrscm2.begin(),mIonizationCrscm2.end());
		mIonizationCrscm2.push_back(&mIonizationUniqueCrscm2);
		mIonizationThresholdeV.erase(mIonizationThresholdeV.begin(),mIonizationThresholdeV.end());
		mIonizationThresholdeV.push_back(threshold);

	}else
	{
		mIonizationPosition=true;
		Log::mD<<"C est partit pour une bonne ionization"<<endl;

		assert(mIonizationThresholdeV.size()==mIonizationCrsPosition.size());

		for(unsigned j=0;j<mIonizationCrsPosition.size();++j)
		{
			unsigned k=mIonizationCrsPosition[j];
			mIonizationCrscm2.push_back(&(mCrscm2[k]));
			if(mbTotalSumIneel)
			{
				mTotalCrscm2+=mCrscm2[k];
			}
		}


		if(mIonizationCrscm2.size()==0)
		{// Nothing was defined : error
			Error err("ElecCrossSection::LoadCrs","No ionization cross section","The code was unable to find ionization cross section (even null). It is considered as an error");
			throw err;
		}
	}

	// Extrapolations
	//
	// For elastic cross sections



	if(extrapolate_elastic_crs)
	{
		ExtrapolateElastic(mMaxDefinedEnergyElasticCrseV,mElasticCrscm2);
	}else
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"The elastic extrapolation was not performed"<<endl;
	}
	if(mbTotalSumIneel)
	{
		mTotalCrscm2+=mElasticCrscm2;
	}
	// We consider that the extrapolation for inelastic cross section is already done
	// in fact, it is in the consideration of the logarithmic interpolation


	// Now, we do the serious game : the redistribution
	//
	//
	//-> Initialization of mOmDegradcm2
	Log::mD<<"Initialization of Degrad matrix"<<endl;
	unsigned n=mGrideV.size();
	assert(n==mGrEngddeV.size());
	mOmDegradcm2.resize(n,n);
	mOmDegradcm2.clear();
	mOutSecondarycm2.resize(n,n);
	mOutSecondarycm2.clear();
	mCrsCincm2.resize(n);
	mCrsCincm2.clear();
	//-> Redistribution over the ionization processes
	Log::mD<<"Redistribute ionization"<<endl;
	// For the test
	/*
	   (*mIonizationCrscm2[0])[0]=8.31090067E-17;
	   (*mIonizationCrscm2[0])[1]=9.37847849E-17;
	   (*mIonizationCrscm2[0])[2]=1.05830315E-16;
	   (*mIonizationCrscm2[0])[3]=1.19425227E-16;
	   (*mIonizationCrscm2[0])[4]=1.34770128E-16;
	   (*mIonizationCrscm2[0])[5]=1.52091913E-16;
	   */
	//	Log::mL<<mIonizationCrscm2.size()<<" processes"<<endl;
	//	Log::mL<<" position of ionizations "<<StdToUblas(mIonizationCrsPosition)<<"/"<<mCrscm_2.size()<<endl;
	try
	{
		RedistributeIonization(mIonizationCrscm2,mIonizationThresholdeV);
	}
	catch(Error& err)
	{
		//Log::SetPriority(Log::ERROR);
		Log::mE<<"Error in the ionization redistribution"<<endl;
		err.Affiche();
		throw err; // We throw the error again!!!
		return false;
	}
	catch(...)
	{
		Log::mE<<"Unknown error in ionization redistribution"<<endl;
		return false;
	}
	//	Log::mL<<"Fin redistribution ionization"<<endl;
	//-> Redistribution over the excitation processes
	//	Log::mL<<"Redistribute excitation"<<endl;
	//	Log::mL<<mExcitationCrscm2.size()<<" processes"<<endl;
	//	Log::mL<<" position of excitations "<<StdToUblas(mExcitationCrsPosition)<<"/"<<mCrscm_2.size()<<endl;
	RedistributeExcitation(mExcitationCrscm2,mExcitationThresholdeV);
	//	cout<<"Fin excitation"<<endl;
	Log::mD<<"End of redistribution"<<endl;
	return true;
}


void ElecCrossSection::ExtrapolateElastic(double vEnerMaxeV,ublas::vector<double>& rCrscm2)
{// Warning, the extrapolation at special energies is not yet done
	unsigned grid_size=mGrideV.size();
	bool defined_gt400=false;
	bool defined_gt2200=false;
	double fac=0;
	unsigned i400=0;
	unsigned i2200=0;
	for(int i=(grid_size-1);i>-1;--i)
	{
		double tmp_e=mGrideV[i];
		if(tmp_e>vEnerMaxeV)
		{
			if(tmp_e>400&&tmp_e<2000)
			{// HERE, the cross section is proportional
				// to E^-0.65
				if(!defined_gt400)
				{
					defined_gt400=true;
					fac=pow(tmp_e,0.65);
					i400=i;
				}else
				{
					rCrscm2[i]=rCrscm2[i400]*fac/pow(tmp_e,0.65);
				}


			}else if(tmp_e>2000)
			{	
				// And here, it is proportional to E^-1

				if(!defined_gt2200)
				{
					// We assure the continuity
					rCrscm2[i]=rCrscm2[i400]*fac/pow(tmp_e,0.65);
					defined_gt2200=true;
					fac=tmp_e;
					i2200=i;
				}else
				{
					rCrscm2[i]=rCrscm2[i2200]*fac/tmp_e;
				}
			}
		}

#ifdef DEBUG
		if(isnan(rCrscm2[i]))
		{
			Error err("Extrapolate elastic","nan","nan in extrapolate elatic");
			throw err;
		}
#endif
	}


}




ElecCrossSection::ElecCrossSection(ublas::vector<double>vGreV,ublas::vector<double> vGrEngddeV):CrossSection(vGreV),mGrEngddeV(vGrEngddeV)
{

	mIonizationPosition=false;
	mIsDefinedTotalCrs=false;
	mbIsOpal=false;
	mEbareV = 0; // Technically, if Opal is selected, the parameter should be defined there
	mReesAlpha = 1/31.5; // The original parameters (Rees 1969)
	mReesBeta = 339.; // put by default here, 
	mReesGamma = 1/2.49; // but can be modified, or perturbated
	// Actually, they are overloaded in the function loadcrs

	mNben=mGrideV.size();
	mCenterGrideV.resize(mNben);
	for(unsigned n=0;n<mNben;++n)
	{
		mCenterGrideV[n]=mGrideV[n]-mGrEngddeV[n]/2.;
	}
}






void ElecCrossSection::RedistributeIonization(std::deque< ublas::vector<double>* >  vIonProcessCrscm2, std::deque<double> vIonProcessThresholdeV)
{//mGrideV
	// mGrEngddeV
	// mOmDegradcm2

	assert(vIonProcessCrscm2.size()==vIonProcessThresholdeV.size());

	unsigned nen=mGrideV.size();
	unsigned nproc=vIonProcessCrscm2.size();
	// In this function, we assume that lsec=true. (It was
	// a parameter in the fortran version)
	// Thus, it is a Rees formulae
	for(unsigned p=0;p<nproc;++p)
	{
		mCrsCincm2[nen-1]+=(*(vIonProcessCrscm2[p]))[nen-1];
	}
	for(unsigned n=0;n<nen-1;++n)
	{// you read n<nen-1 and not n<nen due to boundary conditions!
		double e=mGrideV[n];
		double crosp=0.;
		double crosec=0.;
		//		cout<<"Nombre energy : "<<nen<<" valeur :"<<n<<endl;
		for(unsigned p=0;p<nproc;++p)
		{
			//	cout<<"Nombre process : "<<nproc<<" valeur :"<<p<<endl;
			double del=vIonProcessThresholdeV[p]; // energy loss
	//		Log::mL<<" Process p="<<p<<" threshold = "<<vIonProcessThresholdeV[p]<<"  energy ="<<e<<" position = "<<n<<endl;

			if(del<e)
			{
				unsigned nk=PosEner(e-del,n+1);
				// In trans, there is something strange here with nk, giving an ek which was rewrote into
				double ek=e-nmax(mGrEngddeV[n]/2.,del);
				ek=nmax(0.,ek);
				double es=nmax(0.,e-del);
				double wi=e-ek; // Min energy loss of primary!
				//	cout<<"Initialization cros"<<endl;
				//	cout<<"taille : "<<(*(vIonProcessCrscm2[p])).size()<<endl;
				double cros=(*(vIonProcessCrscm2[p]))[n];//*del/wi;
//#ifdef OLDCODE
//				double cros=(*(vIonProcessCrscm2[p]))[n]*del/wi;
//				Log::mW<<" You are using the old code technique!"<<endl;
//#endif

				//Log::mL<<" CROS P N *del/wi :"<<(*(vIonProcessCrscm2[p]))[n]<<"------>"<<cros<<" del :"<<del<<" ek,wi "<<ek<<" "<<wi<< " Process = "<<p<< "energy pos "<<n<<endl;
				crosp+=cros;
				crosec=(*(vIonProcessCrscm2[p]))[n];
				//			double fac=FNorm;
				//			double facs=FNorm;
				double ffac=0.;
				double ffacs=0.;
				//	cout<<"Initialization ffac"<<endl;
				for(unsigned o=nen-1;o>=nk;--o)
				{
					//				cout<<" ffac, ffacs "<<ffac<<" ====== "<<ffacs<<endl;
					ffac+=AverageEnergyPrimary(o,n,wi,ek)*mGrEngddeV[o];
					ffacs+=AverageEnergySecondary(o,n,del,es)*mGrEngddeV[o];
				}
				double fac=ffac;
				double facs=ffacs;
				//			cout<<"FAC,FACS,nk"<<fac<<" "<<facs<<" and nk"<<nk<<endl;
				//			cout<<"Ek Es"<<ek<<" "<<es<<endl;
				// ???????????????
				//Log::mL<<"Process number "<<p<<" of "<<nproc<<endl;
				if(fac<=0. or facs<=0.)
				{
					if(nen-1-nk<1)
					{
						Log::mD<<"Process number "<<p<<" of "<<nproc<<" for the species "<<mSpecies<<endl;
						Log::mD<<"The energy grid reached a threshold for this process... Passed"<<endl;
						continue;
					}
					Log::mD<<"FAC,FACS,nen,nk"<<fac<<" "<<facs<<" "<<nen<<" "<<nk<<endl;
					Log::mD<<"Processus number : "<<p<<" for species "<<mSpecies<<" "<<mProcessNames[p]<<endl;
					Log::mD<<"Threshold of the process : "<<del<<endl;
					Log::mD<<"Value of e - threshold : "<<es<<endl;
					Log::mD<<"Minimum energy loss of primary : (wi)"<<wi<<endl;
					Log::mD<<"e-threshold or width grid/2"<<ek<<endl;
					Log::mD<<"Nen-1: "<<nen-1<<" nk "<<nk<<" (o=nen-1;o>=nk)"<<endl;
					double a1=AverageEnergyPrimary(nen-1,n,wi,ek);
					Log::mD<<"Average Primary  (nen-1) : "<<a1<<endl;
					a1=AverageEnergySecondary(nen-1,n,del,es);
					Log::mD<<"Average Secondary  (nen-1) : "<<a1<<endl;
					a1=AverageEnergyPrimary(nen-2,n,wi,ek);
					Log::mD<<"Average Primary  (nen-2) : "<<a1<<endl;
					a1=AverageEnergySecondary(nen-2,n,del,es);
					Log::mD<<"Average Secondary  (nen-2) : "<<a1<<endl;

					Error err("Redistribute ionization","fac<=0 or facs<0.","Transition impossible for you inelastic phenomenum. Please check cross sections and threshold");
					throw err;
				}
				fac=cros/fac; // fac and facs corresponds to the sigma(Ep)/N(Ep) factor in Rees 1969, eq 6
				facs=crosec/facs;
				if(mIonizationPosition)
				{
					
					if(mIsAuger.at( mIonizationCrsPosition[p]))
					{
						for(unsigned j=0;j<mAugerEnergy.at( mIonizationCrsPosition[p]).size();++j)
						{
							unsigned energypos2=PosEner(mAugerEnergy.at(mIonizationCrsPosition[p]).at(j),0);//  We search the position for the auger electron
							//rProelecE(iAlt,energypos2)+=produ/(*mpEDdengeV)[energypos2]*(*it)->mpPhotoCrs->mAugerEfficiency.at(i).at(j) ;
							//double old=mOmDegradcm2(energypos2,n);
							double tmp=cros*mAugerEfficiency.at(mIonizationCrsPosition[p]).at(j);
							mOmDegradcm2(energypos2,n)+=tmp;// We have to divide by the grid width, because it is multiplied by it later!
//#ifdef OLDCODE
//							mOmDegradcm2(energypos2,n)+=tmp/mGrEngddeV[energypos2];// We have to divide by the grid width, because it is multiplied by it later!
//#endif
							mOutSecondarycm2(energypos2,n)+=tmp;
							//mOmDegradcm2(n,energypos2)+=tmp;
							//mOutSecondarycm2(n,energypos2)+=tmp;
						}
					}
				}
				//	cout<<"Remplissage Degrad"<<endl;
				for(unsigned ns=nen-1;ns>=nk;--ns)
				{
					// Secondaries
					//
			/// TEST GUIGUI 21 JAN
					double sectmp=AverageEnergySecondary(ns,n,del,es)*facs; // Double differential cross section
					//double sectmp=AverageEnergySecondary(ns,n,del,es)*facs*mGrEngddeV[ns];
					mOmDegradcm2(ns,n)+=sectmp*mGrEngddeV[n];// We multiply by the size of the final grid, to have a value in cm2
					// Actually, there is a problem here: I think the explanation follows:
					// In theory, we should multiply by mGrEngddeV[ns], to have a value in cm2
					// but, to account for the different sizes of the grids, we multiply by
					// mGrEngddeV[n] to have the total number of electron in the flux
					// and then we divide by mGrEngddeV[ns] to have a value in term of eV-1
					// the net result is just a multiplication by mGrEngddeV.
//#ifdef OLDCODE
// 					mOmDegradcm2(ns,n)+=sectmp*mGrEngddeV[n];
//#endif
					mOutSecondarycm2(ns,n)+=sectmp; 
					// Actually, the value should be in cm2/eV, so outsecondary in that form should have the good unit!

#ifdef DEBUG
					if(isnan(mOmDegradcm2(ns,n)))
					{
						Log::mD<<"average secondaries "<<AverageEnergySecondary(ns,n,del,es)<<" facs"<<facs<<endl;
						Error err("mOmDegrad","nan","secondaries");
						throw err;
					}
#endif
					// Primaries
			/// TEST GUIGUI 21 JAN
					mOmDegradcm2(n,ns)+=AverageEnergyPrimary(ns,n,wi,ek)*fac*mGrEngddeV[n];
//#ifdef OLDCODE
//					mOmDegradcm2(n,ns)+=AverageEnergyPrimary(ns,n,wi,ek)*fac*mGrEngddeV[n];
//#endif

#ifdef DEBUG
					if(isnan(mOmDegradcm2(n,ns)))
					{
						Log::mD<<"average primaries "<<AverageEnergyPrimary(ns,n,wi,ek)<<" facs"<<facs<<endl;
						Error err("mOmDegrad","nan","primaries");
						throw err;
					}
#endif
				}

			}



		}
		assert(!isnan(crosp));
		assert(!isinf(crosp));
		mCrsCincm2[n]+=crosp;

	}


}



void ElecCrossSection::RedistributeExcitation(std::deque< ublas::vector<double>* > vExcProcessCrscm2,std::deque<double> vExcProcessThresholdeV)
{
	//        Degrade primaries due to exitation. For reference see Swartz
	//        "Optimization of discrete energy degradation", JGR 1985, p.6587
	// We work with increasing energies
	unsigned nen=mGrideV.size();// Number of energies
	assert(mGrEngddeV.size()==nen);
	assert(mOmDegradcm2.size1()==nen);
	assert(mOmDegradcm2.size2()==nen);
	unsigned nproc= vExcProcessCrscm2.size();// Number of processes 
	assert(mGrideV[0]>*(mGrideV.end()-1)); // Check that the energy grid is decreasing

	for(unsigned p=0;p<nproc;++p)
	{
		mCrsCincm2[nen-1]+=(*(vExcProcessCrscm2[p]))[nen-1];
	}
	for(int n=nen-2;n>-1;--n) // Increasing in energies
	{
		double dde=mGrEngddeV[n];
		double dde2=dde/2.;
		double e=mGrideV[n];
		double cros=0;
		double crosp=0;
		for(unsigned p=0;p<nproc;++p)
		{

			double del=vExcProcessThresholdeV[p]; // energy loss
			if(e>del)
			{
#ifdef TESTNK				
				unsigned nk=PosEner(e-del,n+1) ;//-> the position of the new energy; =================== NOT USED HERE!!! IS THERE A PROBLEM IN THE OLD DEGRAD
				// but used in the corrected function
				
				double ew=nmax(e-mGrideV[n+1],e-mGrideV[nk]);
				assert(ew>0);
				assert(e-mGrideV[n+1]>0);
				assert(e-mGrideV[nk]>0);
#else
				double ew = nmax(del,dde);
#endif
			//	Modification with Jean 6 april 2010 | Questionnement de guillaume 8 april 2010
			//	double ew=del;
			//	GG 18 jan 2011, maybe put del instead of ew will make the work here and assure the conservation in energy
			//	emin and emax are the boundaries for the energy of the degraded electron
				double emin=e-ew-dde2;// Spred primaries over sink : 
				double emax=e-ew+dde2;//cells so that sum of sink (??? GG 20 sept 09)
				// Upper energy position (top, but lower in our decreasing grid)
				unsigned ntop=PosEner(emax,n+1); //-> the top position of energy, corrected
				// Lower energy position(bottom, but higher in our decreasing grid)
				unsigned nbot=PosEner(emin,ntop);//-> the bottom position of energy;
				//nbot=nmax(nmin(nbot,ntop+1),0); -> correct this for inversed energies. -> that gives
				nbot=nmin(nmax(nbot,ntop+1),nen-1);
				//			cout<<"taille : "<<(*(vExcProcessCrscm2[p])).size()<<" Position"<<p<<" energy "<<n<<endl;
				cros=(*(vExcProcessCrscm2[p]))[n]*del/ew;
				crosp+=cros;
				//for(unsigned i=ntop+1;i>nbot;--i)// correct this for inversed
				//do i=nbot+1,ntop-1,1
				//for(unsigned i=ntop+1;i<nbot;++i)
				for(unsigned i=nbot-1;i>=ntop+1;--i)
				{
			/// TEST GUIGUI 21 JAN
					mOmDegradcm2(n,i)+=cros;
//#ifdef OLDCODE
//					mOmDegradcm2(n,i)+=cros/dde;
//#endif

#ifdef DEBUG
					if(isnan(mOmDegradcm2(n,i)))
					{
						Log::mD<<"redist excitattion +cros/dde "<<cros<<" /"<<dde<<endl;
						Error err("mOmDegrad","nan","secondaries");
						throw err;
					}
#endif
				}
				double fac=(emax-mGrideV[ntop]+mGrEngddeV[ntop]/2)/mGrEngddeV[ntop];
				if(fac<0)
				{
					Log::mD<<"fac : "<<fac<<"put to 0"<<endl;
					fac=0.;
				//	Error err("fac < 0.","momdeg","pb");
				//	throw err;
				}
			/// TEST GUIGUI 21 JAN
				mOmDegradcm2(n,ntop)+=cros*fac;
//#ifdef OLDCODE
//				mOmDegradcm2(n,ntop)+=cros*fac/dde;
//#endif

#ifdef DEBUG
				if(isnan(mOmDegradcm2(n,ntop)))
				{
					Log::mD<<"redist excitattion ntop +cros/dde*fac "<<cros<<" /"<<dde<<" * "<<fac<<endl;
					Error err("mOmDegrad","nan","secondaries");
					throw err;
				}
#endif
				fac=(-emin+mGrideV[nbot]+mGrEngddeV[nbot]/2)/mGrEngddeV[nbot];
				if(fac<0)
				{
					Log::mD<<"fac : "<<fac<<"put to 0"<<endl;
					fac=0.;
				//	Log::mL<<"fac : "<<fac<<endl;
				//	Error err("fac < 0.","momdeg","pb");
				//	throw err;
				}
				mOmDegradcm2(n,nbot)+=cros*fac;
//#ifdef OLDCODE
//				mOmDegradcm2(n,nbot)+=cros*fac/dde;
//#endif

#ifdef DEBUG
				if(isnan(mOmDegradcm2(n,nbot)))
				{
					Log::mD<<"redist excitattion nbot +cros/dde*fac "<<cros<<" /"<<dde<<" * "<<fac<<endl;
					Error err("mOmDegrad","nan","secondaries");
					throw err;
				}
#endif
			}
		}
		if(crosp<0)
		{
			Error err("crossp < 0.","momdeg","pb");
			throw err;
		}
		mCrsCincm2[n]+=crosp;
	}
}


unsigned ElecCrossSection::PosEner(double vEnereV,unsigned vPosMin)
{ // was called nlev in Trans*
//	unsigned nen=mGrideV.size();
	double enemin=mGrideV[mNben-1];
	//if(vEnereV<0 or vPosMin>nen)
	if(vEnereV<enemin or vPosMin>mNben-1)
	{
		return mNben-1;
	}
	if(vEnereV>mCenterGrideV[vPosMin])
		return vPosMin;
// divide and conquer algorithm
	
	unsigned pmin=vPosMin;
	unsigned pmax=mNben-1;
	unsigned pmiddle=vPosMin;
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
			if(mCenterGrideV[pmax]<=vEnereV)
			{
				return pmax;
			}else
			{
				return pmin;
			}


		}else
		{
			pmiddle=(pmax-pmin)/2+pmin;
			if(mCenterGrideV[pmiddle]>vEnereV)
			{
				pmin=pmiddle;
			}else
			{
				pmax=pmiddle;
			}
		}
	}// pmax==pmin
	return pmin;
	/*
	for(unsigned n=vPosMin;n<mNben;++n)
	{
	//	double de=mGrEngddeV[n]/2.;
	//	if(vEnereV>mGrideV[n]-de)
		if(vEnereV>mCenterGrideV[n])
		{// If our energy is greater than the 
		// bottom energy of the grid -> Ok
			cout<<"val "<<n<<" "<<mNben<<endl;
			return n;
		}
	}

	Error err("ElecCrossSection::PosEner"," Impossible to find the position","It was impossible to find the position of your energy in the position. Please try to put the electron grid in lower values for the minimum of energy");
        throw err;	
	return 0;
	*/
}





double ElecCrossSection::AverageEnergySecondary(unsigned vN, unsigned vNp, double vW, double vEmaxeV)
{ // Was avf
//	cout<<mGrideV[vN]<<"            "<<mGrideV[vNp]<<endl;
	double dde=mGrEngddeV[vN]/2.;
	double e1=mGrideV[vN]-dde;
	double e2=nmin(mGrideV[vN]+dde,vEmaxeV);
	if(e2<=e1)
	{
		return 0.;
	}
//	cout<<"Normalised : "<<e1<<" "<<e2<<" "<<mGrideV[vNp]<<" "<<vW<<endl;
	return NormalizedSecondaryElectron(e1,e2,mGrideV[vNp],vW)/(e2-e1);

}

double ElecCrossSection::AverageEnergyPrimary(unsigned vN, unsigned vNp, double vW, double vEmaxeV)
{	
	double dde=mGrEngddeV[vN]/2.;
	double e1=mGrideV[vN]-dde;
	double e2=nmin(mGrideV[vN]+dde,vEmaxeV);
	if(e2<=e1)
	{
		return 0.;
	}
//	cout<<"NORMALIZED: "<<(mGrideV[vNp]-e2-vW)<<" "<<(mGrideV[vNp]-e1-vW)<<" "<<mGrideV[vNp]<<" "<<vW<<" "<<(e2-e1)<<endl;;
	return NormalizedSecondaryElectron(
			mGrideV[vNp]-e2-vW,
			mGrideV[vNp]-e1-vW,
			mGrideV[vNp],vW)/(e2-e1);
}



double ElecCrossSection::NormalizedSecondaryElectron(double vEmineV,double vEmaxeV,double vYy,double vZ)
{
	//ww and xx are not from -17;17, but 0:35
	
	double ww[35]={1.36435219e-07, 1.64996788e-06, 1.40066322e-05, 8.77183047e-05,
      4.23141464e-04, 1.63230451e-03, 5.20187337e-03, 1.40839620e-02,
       3.31747346e-02, 6.93263561e-02, 1.30518615e-01, 2.23881721e-01,
       3.52479190e-01, 5.11343718e-01, 6.84302330e-01, 8.44149709e-01,
       9.58392859e-01, 1.00000000e+00, 9.58392859e-01, 8.44149709e-01,
       6.84302330e-01, 5.11343718e-01, 3.52479190e-01, 2.23881721e-01,
       1.30518615e-01, 6.93263561e-02, 3.31747346e-02, 1.40839620e-02,
       5.20187337e-03, 1.63230451e-03, 4.23141464e-04, 8.77183047e-05,
       1.40066322e-05, 1.64996788e-06, 1.36435219e-07};
	double xx[35]={-1.00000000e+00,-9.99999821e-01,-9.99998450e-01,-9.99988973e-01,
      -9.99938786e-01,-9.99728441e-01,-9.99006689e-01,-9.96921241e-01,
      -9.91719127e-01,-9.80284095e-01,-9.57760513e-01,-9.17497337e-01,
      -8.51593614e-01,-7.52291977e-01,-6.14237905e-01,-4.37139213e-01,
      -2.27766961e-01, 0.00000000e+00, 2.27766961e-01, 4.37139213e-01,
      6.14237905e-01, 7.52291977e-01, 8.51593614e-01, 9.17497337e-01,
      9.57760513e-01, 9.80284095e-01, 9.91719127e-01, 9.96921241e-01,
      9.99006689e-01, 9.99728441e-01, 9.99938786e-01, 9.99988973e-01,
      9.99998450e-01, 9.99999821e-01, 1.00000000e+00};

	if(vEmineV>vEmaxeV)
	{// In trans*, this is a warning. Here, it is an error, at least for the tests
		Error err("ElecCrossSection::NormalizedSecondaryElectron","vEmineV>vEmaxeV","Error in the parameters...");
		throw err;
	}
	double fnorm=0;
	const  double delmax=100;
	const double h=2.30999470e-01;
	unsigned nn=static_cast<unsigned>((vEmaxeV-vEmineV)/delmax);
	double d=(vEmaxeV-vEmineV)/(double(nn)+1);

	for(unsigned i=0;i<nn+1;++i)
	{
		double aa=vEmineV+i*d;
		double bb=vEmineV+(i+1)*d;
//		cout<<" aa et bb"<<aa<<" et "<<bb<<" vmin, vmax"<<vEmineV<<" "<<vEmaxeV<<endl;
//		cout<<"nn et d : "<<nn<<" "<<d<<endl;
		double fno=0;
		for(unsigned j=0;j<35;++j)
		{
			double tmp=((bb-aa)*xx[j]+aa+bb)/2;
			fno+=ww[j]*SecondaryElectronFunction(tmp,vYy,vZ);
		}
		fno*=h*d/2.;
		fnorm+=fno;
	}
	return fnorm;
}


double ElecCrossSection::SecondaryElectronFunction(double vX,double vY,double vZ)
{
	const double e_limit=5600;

	if( vX <= 0 || vX >= (vY-vZ) )
		return 0.;
	if(mbIsOpal)
	{
		double f=1.+pow((vX/mEbareV),2.1);
		f=1./f;
		if(vX>((vY-vZ)/2.))
			return 0;
		return f;
	}else
	{ // Rees double differential function

		if( (vX+vZ)/2.49 <= e_limit )
			return exp(-(vX+vZ) * mReesAlpha - mReesBeta *exp(-(vX+vZ)*mReesGamma))/(vX+vZ)
			       *log(
						(sqrt(vY)+sqrt(nmax(0.,vY-vX-vZ)))
						/(sqrt(vY)-sqrt(nmax(0.,vY-vX-vZ)))
				    );
	}
	return 0.;
}

void ElecCrossSection::TestPrint()
{
	//Log::SetPriority(Log::INFO);
	Log::mI<<"We print the inelastic cross section :"<<endl;
	Log::mI<<(mElasticCrscm2)<<endl;



	Log::mI<<"We test the Secondary electronfunction : SecondaryElectronFunction(2,14,6) :"<<SecondaryElectronFunction(2.,14.,6.)<<endl;


	double tmp=1.;
	for(unsigned i=0;i<50;++i)
	{
		tmp+=1.;
		Log::mI<<" SecondaryElectronFunction("<<tmp<<",14,6) :"<<SecondaryElectronFunction(tmp,14.,6.)<<endl;
	}


	Log::mI<<"We test the electron normalization function NormalizedSecondaryElectron(5,50,14,14)"<<NormalizedSecondaryElectron(5,50,14,14)<<endl;

	tmp=1.;
	for(unsigned i=0;i<50;++i)
	{
		tmp+=1.;
		Log::mI<<i<<" NormalizedSecondaryElectron(5,50,"<<tmp<<",14) :"<<NormalizedSecondaryElectron(5.,50.,tmp,14.)<<endl;
	}


	Log::mI<<"test average energy Primary"<<endl;
	tmp=1.;
	const unsigned i1=13;
	const unsigned i2=10;
	for(unsigned i=0;i<50;++i)
	{
		tmp+=1.;
		Log::mI<<i<<" => "<<mGrideV[i1]<<" "<<mGrideV[i2]<<" AverageEnergySecondary(20,10,2.,"<<tmp<<") :"<<AverageEnergySecondary(i1,i2,2.,tmp)<<endl;
	}



	Log::mI<<"Maintenant : section efficace ionization"<<endl;
	for(unsigned p=0;p<mIonizationCrscm2.size();++p)
	{
		Log::mI<<" Processus ionization numero "<<p<<endl;
		//MathString::print1d(*mIonizationCrscm2[p]);
		Log::mI<<(*mIonizationCrscm2[p])<<endl;;

	}

	Log::mI<<"Maintenant : section efficace excitation"<<endl;
	for(unsigned p=0;p<mExcitationCrscm2.size();++p)
	{
		Log::mI<<" Processus excitation numero "<<p<<endl;
		//MathString::print1d(*mExcitationCrscm2[p]);
		Log::mI<<*mExcitationCrscm2[p]<<endl;
	}
}





