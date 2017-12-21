/**
 * \file electronionization.cpp
 * \brief To compute the electron impact ionization and excitation.
 *        It solves the transport of the electrons.
 *        Depends on some Fortran functions NOW
 *        Copyright G Gronoff Sept 2009, 
 *        see the other files for Fortran-truc copyright
 *        Last Modification : $Id: electronionization.cpp 1900 2014-01-16 22:30:20Z gronoff $
 *
 */

#include "electronionization.hpp"
using namespace std;
using boost::timer;
using boost::progress_timer;
using boost::progress_display;


ElectronImpactIonization::ElectronImpactIonization(XmlParameters* pParam,ublas::vector<double>* vpElectronDensitycm_3,ublas::vector<double>* vpElectronTemperatureK,ublas::vector<double>* vpAltGridKm): mpParameter(pParam),mpElectronDensitycm_3(vpElectronDensitycm_3),mpElectronTemperatureK(vpElectronTemperatureK),mpAltGridKm(vpAltGridKm)
{
	mIsDisortChecking=1; // disort is checking by default...
	mError=0.;
	mProdTot=0.;
	mProdIonTot=0.;
	mTotalEnergy=0.;
	mReflectedEnergy=0.;
	mTransmittedEnergy=0.;
	mAbsorbedEnergy=0.;
	mInelasticEnergy=0.;
	mResultTotalEnergy=0.;
	mElos1=0.;
	mElos2=0.;
	mElos3=0.;
	mIlos1=0.;
	mIlos2=0.;
	mIlos3=0.;
	mSecondarySumeV=0.;
	mDegradedSumeV=0.;
	mCoulombianSumeV=0.;
	mFactUncertLE=1.;
	mNbalt=vpAltGridKm->size();
}




void ElectronImpactIonization::Mstream(unsigned vPosEner,const ublas::matrix<double>& rQInt,ublas::vector< ublas::matrix<double> >& rQIntensity )
{
//	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	int nalt=static_cast<int>(mNbalt);
	int ndbl=2*nalt-1;
	int nblayer=nalt-1;
	int nang=static_cast<int>(mNbang);
	int nango2=static_cast<int>(mNbang2);

	//ublas::zero_matrix<double> tmp_src(nang,3*(nalt-1)); // The photoelectron source! <---------------
	ublas::matrix<double> src(mNbang,3*(nalt-1)); // The photoelectron source! <---------------
	src.clear();


	// GG: jan 2011, comments following the writing of the uncertainty paper
	// In this loop, we add the Qphoto to the Qredist
	// It corresponds to the sum of Qphoto and the redistribution term 
	// in equation 3.3b in Lummerzheim PHD thesis
	//
	// Actually, Qphoto is the flux coming initially from the photoelectrons
	// since we added other sources, it is mpElecflux->mFluxAcm_2s_1sr_1
	// To that, we add the electrons coming from the redistribution in energy
	// in rQint. This rQint part is computed in QMstream


	for(unsigned i=0;i<static_cast<unsigned>(nblayer);++i)
	{
//		Log::mL<<"Loop for layers for src : "<<i<<endl;
		for(unsigned k=0;k<mNbang2;++k)
		{
			unsigned k2=k+nango2+1;
			unsigned km=nango2-k-1;

			unsigned fk=k+nango2;
			unsigned fmk=nango2-k-1;

			//Log::SetPriority(Log::INFO);
			double quelle=rQInt(2*i,k2)+mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i](vPosEner,fk);

#ifdef DEBUG
			if(isnan(quelle))
			{
				Log::mD<<"quelle : "<<quelle<<" rqint : "<<rQInt(2*i,k2)<<"  flux :"<<mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i](vPosEner,fk)<<endl;
				Log::mD<<" k2 "<<k2<<" i "<<i<<" vPosEner :"<<vPosEner<<endl;
				
				Error err("nan quelle","","");
				throw(err);
			}
#endif 
			src(fk,3*i)=nmax(quelle,0.);

			quelle=rQInt(2*i,km)+mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i](vPosEner,fmk);
#ifdef DEBUG
			if(isnan(quelle))
			{
				Log::mD<<"quelle : "<<quelle<<" rqint : "<<rQInt(2*i,km)<<"  flux :"<<mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i](vPosEner,fmk)<<endl;
				Log::mD<<" fmk "<<fmk<<" i "<<i<<" vPosEner :"<<vPosEner<<endl;
				Error err("nan quelle","","");
				throw(err);
			}
#endif
			src(fmk,3*i)=nmax(quelle,0.);


			quelle=rQInt(2*i+2,km)+mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i+1](vPosEner,fmk);
#ifdef DEBUG
			if(isnan(quelle))
			{
				Log::mD<<"quelle : "<<quelle<<" rqint : "<<rQInt(2*i+2,km)<<"  flux :"<<mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i+1](vPosEner,fmk)<<endl;
				Log::mD<<" fmk "<<fmk<<" i "<<i+1<<" vPosEner :"<<vPosEner<<endl;
				Error err("nan quelle","","");
				throw(err);
			}
#endif
			src(fmk,3*i+2)=nmax(quelle,0.);

			quelle=rQInt(2*i+2,k2)+mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i+1](vPosEner,fk);
#ifdef DEBUG
			if(isnan(quelle))
			{
				Log::mD<<"quelle : "<<quelle<<" rqint : "<<rQInt(2*i+2,k2)<<"  flux :"<<mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i+1](vPosEner,fk)<<endl;
				Log::mD<<" k2 "<<k2<<" i+1 "<<i+1<<" vPosEner :"<<vPosEner<<endl;
				Error err("nan quelle","","");
				throw(err);
			}
#endif
			src(fk,3*i+2)=nmax(quelle,0.);
		}
	}

//	Log::mL<<src<<endl;


	int usrtau=1;
	int nusrang=0;
	int vmu=0;
	int vnphi=0;
	double vfbeam=0.,vumu0=1,vphi0=0;
	int vlamber=1,vldeltam=1;
	int vexorc=1,vonlyfl=1,vlinear=1;
	double vbsrc=0,vtsrc=0;
	double vaccur=0;
	int vmaxphi=1;
	ublas::vector<double> vdfdt(ndbl);
	ublas::vector<double> vrfldir(ndbl);
	vdfdt.clear();
	vrfldir.clear();
	ublas::vector<double> vnumu(nang);
	vnumu.clear();
	ublas::vector<double> vhl(nang+1);
	vhl.clear();
	ublas::matrix<double> vsrcu(nang,3*(nblayer));
	vsrcu.clear();
	///ublas::matrix<double> vsrcu(zerovsrcu);
	ublas::vector<double> vphi(vmaxphi);
	vphi.clear();


	int vprtn[7]={0,0,0,0,0,0,0};
//	Log::mL<<"Definition of precipitation fluxes"<<endl;
//	ublas::zero_vector<double> zeropf(nang/2);
	ublas::vector<double> precip_flux(nang/2);
	precip_flux.clear();
	if(mpElecFlux->mPrecipitationDefined)
	{
	//	vector<double> prec=mpElecFlux->mPrecipFluxDowncm_2s_1sr_1[vPosEner];
	//	std::reverse(prec.begin(),prec.end());
	//	precip_flux=StdToUblas(prec);
		precip_flux=ublas::matrix_row< ublas::matrix<double> > (mpElecFlux->mPrecipFluxDowncm_2s_1sr_1,vPosEner);
		//std::reverse(precip_flux.begin(),precip_flux.end());
	//	cout<<"precip_flux "<<precip_flux<<endl;
	}
	else
	{
//		Log::mL<<"FLux redefined"<<endl;
		precip_flux.resize(nang/2);
		precip_flux.clear();
	}

	//ublas::zero_vector<double> zerovuu(nang*(ndbl)*vmaxphi);
	ublas::vector<double> vuu(nang*(ndbl)*vmaxphi);
	vuu.clear();
	// Output
	ublas::vector<double> fluxdown(ndbl), fluxup(ndbl),ruavg(ndbl);
	fluxdown.clear();
	fluxup.clear();
	ruavg.clear();
//	ublas::matrix<double> zerorUou(ndbl,nang);
	ublas::matrix<double> rUou(ndbl,nang); // ATTENTION there is a swich between fortran and c++ concerning array notation
	rUou.clear();
	int posener=static_cast<int>(vPosEner+1);// Does not have an influence on the result : just for error msg
//	Log::mL<<"Call to disort"<<endl;

	for(unsigned jt=0;jt<mPhaseFunction(vPosEner).size1();++jt)
	{

		for(unsigned kt=0;kt<mPhaseFunction(vPosEner).size2();++kt)
		{
			assert(mPhaseFunction(vPosEner)(jt,kt)<1.1);
		}
	}
/*
	Log::mL<<"mdifftau"<<mDiffTau<<endl;
	Log::mL<<"mssalb"<<mSsalb<<endl;
	Log::mL<<"mphasefunction"<<mPhaseFunction<<endl;
	Log::mL<<"mutau"<<mUTau<<endl;
	Log::mL<<"fluxdown"<<fluxdown<<endl;
	Log::mL<<(*mpElecC)[vPosEner]<<  precip_flux<<endl;
*/
//	Log::mL<<fluxup<<endl;

	disort_(&posener, &nblayer,&mDiffTau(vPosEner,0) ,&mSsalb(vPosEner,0), &mPhaseFunction(vPosEner)(0,0),
			&usrtau,&ndbl,&mUTau(vPosEner,0),&nang,&nusrang,
			&vmu,&vnumu(0),&vnphi,&vphi(0),&vfbeam,&vumu0,
			&src(0,0),&vsrcu(0,0),&vphi0,(&precip_flux(0)),&vlamber,
	//		&src(0,0),&vsrcu(0,0),&vphi0,(&precip_flux(0)+nang/2*sizeof(double)),&vlamber,
		//	&src(0,0),&vsrcu(0,0),&vphi0,(&precip_flux(0))+static_cast<unsigned>(nang/2+1),&vlamber,
			//&mAlbedo,&vhl(0)+1,&vbsrc,&vtsrc,&vldeltam,
			&mAlbedo,&vhl(0),&vbsrc,&vtsrc,&vldeltam,
			&vexorc,&vonlyfl,&vaccur,vprtn,&nblayer,&ndbl,
			&nang,&nang,&vmaxphi,&vrfldir(0),&fluxdown(0),&fluxup(0),
			&vdfdt(0),&ruavg(0),&vuu(0),&rUou(0,0),&nango2,&vlinear,&mIsDisortChecking);
//	Log::mL<<"Fin call to disort"<<endl;
/*
#ifdef DEBUG

	for(int na=0;na<nang;++na)
	{
		for(int m=0;m<2*nalt-1;++m)
		{
			if(rUou(m,na)<-1E-42)
			{
				Log::mL<<"Error, negative for (m,na)="<<m<<" "<<na<<endl;
				Log::mL<<"Error for pos energ = "<<vPosEner<<endl;
				Log::mL<<"Value equal to "<<rUou(m,na)<<endl;
				rUou(m,na)=0.;
			}
		}
	}

#endif
*/

	unsigned nang2=static_cast<unsigned>(nang/2);
	for(unsigned i=0;i< static_cast<unsigned>(nalt-1);++i)
	{
		unsigned mt=2*i; // top of the layer
		unsigned mc=2*i+1; //center of the layer
		rQIntensity(vPosEner)(mt,nang2)=PI*4*ruavg(mt);
		if(rQIntensity(vPosEner)(mt,nang2)<0)
			rQIntensity(vPosEner)(mt,nang2)=0.;

		rQIntensity(vPosEner)(mc,nang2)=PI*4*ruavg(mc);
		if(rQIntensity(vPosEner)(mc,nang2)<0)
			rQIntensity(vPosEner)(mc,nang2)=0.;

		mFluxHemisphericDown(vPosEner,i)=fluxdown(mt);
		mFluxHemisphericUp(vPosEner,i)=fluxup(mt);
		mFluxHemisphericTot(vPosEner,i)=fluxdown(mt)+fluxup(mt);
		mFluxHemisphericNet(vPosEner,i)=fluxdown(mt)-fluxup(mt);

		for(unsigned k=0;k<nang2;++k)
		{
			unsigned k2=k+nang/2+1;
			unsigned km=nang/2-k-1;

			unsigned fk=k+nang/2;
			unsigned fmk=nang/2-k-1;

			if(rUou(mt,fk)>0)
				rQIntensity(vPosEner)(mt,k2)=rUou(mt,fk);
			if(rUou(mt,fmk)>0)
				rQIntensity(vPosEner)(mt,km)=rUou(mt,fmk);
			if(rUou(mc,fk)>0)
				rQIntensity(vPosEner)(mc,k2)=rUou(mc,fk);
			if(rUou(mc,fmk)>0)
				rQIntensity(vPosEner)(mc,km)=rUou(mc,fmk);
		}
	}
	unsigned mb=2*(nalt-1);
	if(ruavg(mb)>0)
		rQIntensity(vPosEner)(mb,nang2)=PI*4*ruavg(mb);
	mFluxHemisphericDown(vPosEner,nalt-1)=fluxdown(mb);
	mFluxHemisphericUp(vPosEner,nalt-1)=fluxup(mb);
	mFluxHemisphericTot(vPosEner,nalt-1)=fluxdown(mb)+fluxup(mb);
	mFluxHemisphericNet(vPosEner,nalt-1)=fluxdown(mb)-fluxup(mb);
	for(unsigned k=0;k<nang2;++k)
	{
		unsigned k2=k+nang2+1;
		unsigned km=nang2-k-1;

		unsigned fk=k+nang2;
		unsigned fmk=nang2-k-1;
		if(rUou(mb,fk)>0)
			rQIntensity(vPosEner)(mb,k2)=rUou(mb,fk);
		if(rUou(mb,fmk)>0)
			rQIntensity(vPosEner)(mb,km)=rUou(mb,fmk);
	}

#ifdef DEBUG

	for(int na=0;na<nang+1;++na)
	{
		for(int m=0;m<2*nalt-1;++m)
		{
			if(rQIntensity(vPosEner)(m,na)<-1E-42)
			{
				Log::mD<<"Error, negative for (m,na)="<<m<<" "<<na<<endl;
				Log::mD<<"Error for pos energ = "<<vPosEner<<endl;
				Log::mD<<"rqintensity Value equal to "<<rQIntensity(vPosEner)(m,na)<<endl;
				Error err("machin","bidule","truc");
				throw err;
			}
		}
	}

#endif





}



void ElectronImpactIonization::QMstream(unsigned vPosEner,ublas::matrix<double>& rQInt
		, ublas::vector< ublas::matrix<double> >& rQIntensity )
{
	// In this function, we compute the summation term in equation 3.3b in the
	// Lummerzheim thesis
	rQInt.clear();

	unsigned nalt=mpElectronDensitycm_3->size();
	int nang=(mpElecFlux->NbAngles());
	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	unsigned nbsp=mpSp->size();

	//ublas::vector<double> sumrds(2*nalt-1),sumfsec(2*nalt-1);
	ublas::vector<double> sumint(2*nalt-1),sumintw(2*nalt-1);
	ublas::vector<double> sumsec(2*nalt-1),sumdeg(2*nalt-1),sumcoul(2*nalt-1);
//	sumrds.clear();
//	sumfsec.clear();
	sumint.clear();
	sumintw.clear();
	sumsec.clear(); // temp to compute the energy in secondaries
	sumdeg.clear(); // temp to compute the energy in degraded
	sumcoul.clear(); // temp to compute the energy in coulombian

	// GG: jan 2011, comments following the writing of the uncertainty paper
	// In this loop, we compute the Rin term (Lummerzheim, eq just after 3.3b)
	// and we multiply it by Delta E
	for(unsigned ene=0;ene<=vPosEner;++ene)
	{
#ifdef OLDCODE
		double dde=(*(mpElecFlux->mpElecDdengeV))[ene];
#endif
	//	ublas::zero_vector<double> machin(2*nalt-1);
		ublas::vector<double> rdp(2*nalt-1);
		ublas::vector<double> rds(2*nalt-1);
		ublas::vector<double> fsec(2*nalt-1);
		rdp.clear();
		rds.clear();
		fsec.clear();

		// GG: jan 2011, comments following the writing of the uncertainty paper
		// Here, computation of the sum of the cross sections x density of the species
		for(unsigned sp=0;sp<nbsp;++sp)
		{

			if( (*mpSp)[sp]->CheckElecCrs() && (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
			{
				for(unsigned i=0;i<nalt;++i)
				{
					double dens=(*mpSp)[sp]->mTotDensitycm_3[i];
					double omsec=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(vPosEner+1,ene);
					double omdeg=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(ene,vPosEner+1);
#ifdef DEBUG
					double t1,t2;
					t1=rdp(2*i);
					t2=rds(2*i);
#endif
					rdp(2*i)+=omdeg*dens;
					rds(2*i)+=omsec*dens;
#ifdef DEBUG
					if(isnan(rdp(2*i)) || isnan(rds(2*i)))
					{
						Log::mD<<"Altitude nb "<<i<<" / "<<nalt<<endl;
						Log::mD<<"Energy nb "<<ene<<" / "<<nen<<" et pos : "<<vPosEner<<endl;
						Log::mD<<omdeg<<" et dems "<<dens<<endl;
						Log::mD<<" omsec "<<omsec<<endl;
						Log::mD<<"Previous ="<<t1<<" "<<t2<<endl;
						Log::mD<<"now ="<<rdp(2*i)<<" "<<rds(2*i)<<endl;
						Log::mD<<" sp ="<<sp<<" i="<<i<<endl;
						Error err("rdp,rds erreur","tmpval","tmvpval");
						throw err;
					}
#endif
					if(i!=nalt-1)
					{
						double densi=(*mpSp)[sp]->mTotDensitycm_3[i+1];
						rdp(2*i+1)+=omdeg*(dens+densi)/2.;
						rds(2*i+1)+=omsec*(dens+densi)/2.;

					}
				}
			}
		}



		// GG: jan 2011, comments following the writing of the uncertainty paper
		// Here, division of that sum by the opacity tau, denominator in the Lummerzheim eq
		// and multiplication by the width
		for(unsigned m=0;m<2*nalt-1;++m)
		{
			double tmpval=mCTot(vPosEner+1,m);
#ifdef OLDCODE
			double tmpval=mCTot(vPosEner+1,m)/dde;
#endif

#ifdef DEBUG
			if(isnan(tmpval))
			{
				Log::mD<<tmpval<<endl;
				Error err("tmpval erreur","tmpval","tmvpval");
				throw err;
			}
#endif
			double tmp1=rdp(m);
			double tmp2=rds(m);
			rdp(m)/=tmpval;
			rds(m)/=tmpval;
#ifdef DEBUG
			if(isnan(rdp(m)) || isnan(rds(m)))
			{
				Log::mD<<tmpval<<" "<<rdp(m)<<" "<<rds(m)<<endl;
				Log::mD<<" rdp "<<tmp1<<"  "<<tmp2<<endl;
				Error err("rdp,rds erreur","tmpval","tmvpval");
				throw err;
			}
#endif
		}
//		Log::mL<<"rdp"<<rdp<<endl;
//		Log::mL<<"rds"<<rds<<endl;
//		Log::mL<<"fsec"<<fsec<<endl;
		// GG: jan 2011, comments following the writing of the uncertainty paper
		// Here we compute the total flux at the energy ene at the considered altitude
		// This flux is actually a mean flux divided by steradian
		for(unsigned k=0;k<static_cast<unsigned>(nang+1);++k)
		{
			if(k!=static_cast<unsigned>(nang/2))
			{
				double weight;
				if(k<static_cast<unsigned>(nang/2))
				{
					weight=mpElecFlux->ReturnAngleWeight(k);
				}else
				{
					weight=mpElecFlux->ReturnAngleWeight(k-1);
				}
	// GG: jan 2011, comments following the writing of the uncertainty paper
	// Using the weight, we compute fsec, which is the mean value of the intensity
	// we have to remember that this intensity is in unit of cm_2 s_1 eV_1 sr_1
	// and the same for fsec
				for(unsigned m=0;m<2*nalt-1;++m)
				{
					fsec(m)+=rQIntensity(ene)(m,k)*weight/2.;
					sumint(m)+=rQIntensity(ene)(m,k);
					sumintw(m)+=rQIntensity(ene)(m,k)*weight/2.;
				}
			}
		}

		// GG: jan 2011, comments following the writing of the uncertainty paper
		// Here, we compute the flux of secondary electrons, tmp0, which is
		// the same at every angle (isotropic creation)
		// The rdp(m)*rQintensity corresponds to the degradation in energy during
		// inelastic collisions, for which the angle is unchanged
		// This is the first term of rQint. 
		// The second one corresponds to the energy coming from the degradation
		// due to Coulombian collision (last term in Lummerzheim eq 3.3b)
		for(unsigned m=0;m<2*nalt-1;++m)
		{ // rds is dimensionless, so tmp0 is in the same units
		 // it is the same with rdp. So, rQint has the same units as rQIntensity
		 // which is what is expected!
			double tmp0=rds(m)*fsec(m);
			sumsec(m)+=tmp0;
			sumdeg(m)+=rdp(m)*fsec(m);
			sumcoul(m) = fsec(m);
		//	sumrds(m)+=rds(m);
		//	sumfsec(m)+=fsec(m);
			if(tmp0<1E33 and tmp0>1E-33)
			{
				rQInt(m,nang/2)+=tmp0/2;
			}
			for(unsigned k=0;k<static_cast<unsigned>(nang+1);++k)
			{
				if(k!=static_cast<unsigned>(nang/2))
				{
					double tmp1=rdp(m)*rQIntensity(ene)(m,k)+tmp0;
					if(tmp1<1E33 and tmp1>1E-33)
					{
						rQInt(m,k)+=tmp1;
					}
				}
			}
		}
	}
//	ublas::zero_vector<double> zfac(2*nalt-1);
	ublas::vector<double> fac(2*nalt-1);
	fac.clear();
	// GG: jan 2011, comments following the writing of the uncertainty paper
	// Comoputation of the the energy coming from the degradation
	// due to Coulombian collision (last term in Lummerzheim eq 3.3b)
	double dde=(*(mpElecFlux->mpElecDdengeV))[vPosEner+1];
	// Here ne * Lambda
	for(unsigned i=0;i<nalt;++i)
	{
		fac(2*i)=mEnergyLosseV(i,vPosEner)/mCTot(vPosEner+1,2*i)/dde;

		if(i!=nalt-1)
		{
			fac(2*i+1)=(mEnergyLosseV(i,vPosEner)+mEnergyLosseV(i+1,vPosEner))/mCTot(vPosEner+1,2*i+1)/dde/2.;
		}
	}
	
	// Here multiplied by the flux
	for(unsigned m=0;m<2*nalt-1;++m)
	{
		for(unsigned k=0;k<static_cast<unsigned>(nang+1);++k)
		{
			if(k!=static_cast<unsigned>(nang/2))
			{
				double tmp0=fac(m)*rQIntensity(vPosEner)(m,k);
				if(tmp0<1E33 and tmp0>1E-33)
				{
					rQInt(m,k)+=tmp0;
				}
			}
		}
		sumcoul(m)*= fac(m);

	}


	for(unsigned i = 0; i < nalt -1 ; ++i)
	{
		double tmpde = ((*mpAltGridKm)[i] - (*mpAltGridKm)[i+1]) / 4.;
		mSecondarySumeV  += (sumsec(2 * i) + sumsec(2*i + 1) + sumsec(2*i+2)) * tmpde * (*(mpElecFlux->mpElecDdengeV))[vPosEner] * (*(mpElecFlux->mpElecCentEeV))[vPosEner];
		mDegradedSumeV   += (sumdeg(2 * i) + sumdeg(2*i + 1) + sumdeg(2*i+2)) * tmpde * (*(mpElecFlux->mpElecDdengeV))[vPosEner] * (*(mpElecFlux->mpElecCentEeV))[vPosEner];
		mCoulombianSumeV += (sumcoul(2 * i) + sumcoul(2*i + 1) + sumcoul(2*i+2))  * tmpde  * (*(mpElecFlux->mpElecDdengeV))[vPosEner] * (*(mpElecFlux->mpElecCentEeV))[vPosEner];


	}

	// GG JAN 2011
	
	//	tmptot_e_intensity=MathFunction::TrapzInt(*mpAltGridKm,tmpintensityenergy)*1E5;





/*	Log::mL<<"Fac : "<<fac<<endl;
	Log::mL<<"Elosse"<<endl;
	for(unsigned i=0;i<nalt;++i)
	{
		Log::mL<<mEnergyLosseV(i,vPosEner)<<"\t";
	}
	Log::mL<<endl;
*
*/
#ifdef DEBUG
/*	Log::mL<<"rQint"<<endl;
	for(unsigned i=0;i<rQInt.size1();++i)
	{
		ublas::matrix_row< ublas::matrix<double> > ro(rQInt,i);
		Log::mL<<i<<"\t"<<ro(mNbang2)<<"\t"<<sumrds(i)<<"\t"<<sumfsec(i)<<"\t"<<sumint(i)<<"\t"<<sumintw(i)<<endl;//<<ro<<endl;
//		for(unsigned j=0;j<ro.size();++j)
//		{
//			Log::mL<<ro(j)<<"\t";
//		}
//		Log::mL<<endl;

	}
	*/
#endif
}


void ElectronImpactIonization::ComputeElectronImpact(std::deque<Specie*> vSp,
					   const EFlux& vFlux,
					   std::deque<Specie*>& rResult)
{
	//Log::SetPriority(Log::CONFIG);
	Log::mD<<"We init the process!"<<endl;
	mpElecFlux=const_cast<EFlux*>(&vFlux);
	//		double de= (*(mpElecFlux->mpElecDdengeV))[n]/2.;
	//		if(vEnergyeV>(*(mpElecFlux->mpElecCentEeV))[n]-de)
	//		{// If our energy is greater than the 
	//		// bottom energy of the grid -> Ok
	//			return n;
	mpElecC=mpElecFlux->mpElecCentEeV;
	mpElecB=mpElecFlux->mpElecBotEeV;
	mpElecD=mpElecFlux->mpElecDdengeV;
	mpGAngle=mpElecFlux->ReturnAngle();
	mNbang=mpGAngle->mNbAngles;
	mNbang2=mpGAngle->mNbAngles/2;
/*	mSecondarySumeV=0.;
	mDegradedSumeV=0.;
	mCoulombianSumeV=0.;
*/

	//!!! the init should be after the const cast!!!
	Init();
	// Energy loss is computed, porter and ruther also
	// (so mPhaseFunction is almost computed)

	// Init the species
	mpSp=&vSp;

	mNbsp=mpSp->size();
	mNben=mpElecC->size();
	mCenterGrideV.resize(mNben);
	for(unsigned i=0;i<mNben;++i)
	{
		mCenterGrideV[i]=(*(mpElecC))[i]-(*(mpElecD))[i]/2.;
	}

	ublas::vector<unsigned> ntherm(mNbalt);

	for(unsigned sp=0;sp<mNbsp;++sp)
	{
		//Log::SetPriority(Log::INFO);
		//Log::mL<<"Pre-initialization"<<endl;
		if( (*mpSp)[sp]->CheckElecCrs() && (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
		{
		//	Log::mL<<"INIT ELECTRON IMPACT FOR "<<(*mpSp)[sp]->mName<<endl;
			(*mpSp)[sp]->InitElectronImpact(mNbalt);
		}else
		{
			Log::mD<<"HOUSTON : "<<(*mpSp)[sp]->mName<<endl;
			if(! (*mpSp)[sp]->CheckElecCrs())
				Log::mW<<"Pb with CHECK elec CRS"<<endl;
			if(!  (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
				Log::mW<<"Pb with mIsDefinedCRS"<<endl;

		}
	}



	// And the phase function, with default parameters
	ComputePhaseFunction();// Here the mPhaseFunction is computed.
	ComputeCTot();

	//Log::mL<<"Number of alt"<<mNbalt<<endl;;
	//	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	//	int nang=(mpElecFlux->NbAngles());

	//	ublas::zero_matrix<double> zqint(2*nalt-1,nang+1);
	ublas::matrix<double> qint(2*mNbalt-1,mNbang+1);
	qint.clear();
	// Definition of the intensity matrix
	//a	ublas::matrix<double> ztmptny(2*mNbalt-1,mNbang+1);
	ublas::matrix<double> tmptny(2*mNbalt-1,mNbang+1);
	tmptny.clear();
	ublas::vector< ublas::matrix<double> > qintensity(mNben);
	//	ublas::vector< ublas::matrix<double> > qintensity2(nen);
	//

	//	ublas::zero_matrix<double> ztmpvecmat(nang,3*(nalt-1));
	//	ublas::zero_matrix<double> ztmpvecmat2(2*nalt-1,nang);
	//	ublas::zero_vector<double> znuavg(2*nalt-1);

	for(unsigned i=0;i<mNben;++i)
	{
		qintensity(i)=tmptny;
	}
	Log::mI<<"Begin of the energy loop"<<endl;


	mFluxHemisphericDown.resize(mNben,mNbalt);
	mFluxHemisphericUp.resize(mNben,mNbalt);
	mFluxHemisphericTot.resize(mNben,mNbalt);
	mFluxHemisphericNet.resize(mNben,mNbalt);

	mFluxHemisphericDown.clear();
	mFluxHemisphericUp.clear();
	mFluxHemisphericTot.clear();
	mFluxHemisphericNet.clear();
	/*
	   for(unsigned t=0;t<mNben;++t)
	   {
	   ublas::matrix_row<  ublas::matrix<double> > ro((*mpSp)[0]->mpElecCrs->mOmDegradcm2,t);
	   Log::mL<<ro<<endl;
	   }

	   for(unsigned t=0;t<mNben;++t)
	   {
	   for(unsigned u=0;u<=t;++u)
	   {
	   Log::mL<<(*mpSp)[0]->mpElecCrs->mOmDegradcm2(u,t)<<"\t";
	   }
	   Log::mL<<endl;
	   for(unsigned u=0;u<=t;++u)
	   {
	   Log::mL<<(*mpSp)[0]->mpElecCrs->mOmDegradcm2(t,u)<<"\t";
	   }
	   Log::mL<<endl;
	   }
	   ifstream wfs("qprim.dat");
	   for(unsigned i=0;i<mNbalt;++i)
	   {
	   for(int n=mNben-1;n>-1;--n)
	   {

	   for(unsigned k=0;k<mNbang;++k)
	   {
	   wfs>>mpElecFlux->mFluxAcm_2s_1eV_1sr_1[i](n,k);
	   }


	   }

	   }
	   wfs.close();
	   */
	/*
	   ifstream gfs("Utau.dat");
	   for(unsigned i=0;i<mNben-1;++i)
	   {
	   for(unsigned j=0;j<2*mNbalt-1;++j)
	   {
	   double tmp;
	   gfs>>tmp;
	   cout<<mUTau(i,j)<<" -- "<<tmp<<" || j="<<j<<" e "<<i<<endl;
	   mUTau(i,j)=tmp;
	   }
	   }
	   gfs.close();
	   ifstream hfs("Difftau.dat");
	   */
	/*
	   ifstream mfs("Ssalb.dat");

	//ifstream ffs("pmom.dat");
	ifstream lfs("ctot.dat");
	for(unsigned i=0;i<mNben-1;++i)
	{
	for(unsigned j=0;j<mNbalt-1;++j)
	{
	//		   double tmp;
	//		   hfs>>tmp;
	//		   cout<<mDiffTau(i,j)<<" __ "<<tmp<<" || j="<<j<<" e "<<i<<endl;
	//		   mDiffTau(i,j)=tmp;
	double tmp1,tmp0;
	lfs>>tmp0;
	lfs>>tmp1;
	double tmp2;
	mfs>>tmp2;
	cout.precision(9);
	cout<<mSsalb(i,j)<<" ** "<<tmp2<<"\t -> "<<(mSsalb(i,j)-tmp2)/tmp2*100<<"%\t E="<<i<<"\t alt="<<j<<endl;
	cout<<mCTot(i,j)<<"--"<<tmp1<<" ? "<<tmp0<<endl;
	//  if(nabs( (mSsalb(i,j)-tmp2)/tmp2*100)>0.01  )
	//	 mSsalb(i,j)=tmp2;
	// if(i>40&& i<46)
	//	 mSsalb(i,j)=tmp2;
	if(mSsalb(i,j)>0.9999)
	mSsalb(i,j)=tmp2;



	//		   for(unsigned k=0;k<mNbang;++k)
	//		   {
	//			   ffs>>mPhaseFunction(i)(j,k);
	//		   }
	}
	}
	// ffs.close();
	*/
	/*
	//Log::mL<<(*mpSp)[0]->mpElecCrs->mOmDegradcm2<<endl;
	ifstream tfs("omdeg.dat");
	for(unsigned n=1;n<mNben;++n)
	{

	for(unsigned sp=0;sp<mNbsp;++sp)
	{
	for(unsigned j=0;j<=n;++j)
	{
	//	double tmp=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(j,n);
	//	double tmp3=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(j,n);
	double tmp2=0;
	tfs>>tmp2; //(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(j,n);
	//double tmp2=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(j,n);
	//	if(tmp2!=0 && nabs((tmp-tmp2)/tmp2)>0.1)
	//	if(tmp2!=0 && nabs((tmp-tmp2)/tmp2)>1 && tmp2>1E-20  )
	//		Log::mL<<"|"<<tmp<<" "<<tmp2<<" ??"<<tmp3<<" ??sp="<<sp<<" n "<<n<<" j"<<j<<"|\t";
	//	double percent=(tmp-tmp2)/tmp2*100;
	//	if(abs(percent)>100)
	//	{
	//		Log::mL<<percent<<"\t";
	//	if(tmp2!=0 && nabs((tmp-tmp2)/tmp2)>0.1)
	//			if(tmp2!=0 && nabs((tmp-tmp2)/tmp2)>1 && tmp2>1E-20  )
	(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(j,n)=tmp2;
	//	}
	}
	}
	//	Log::mL<<endl;
	for(unsigned sp=0;sp<mNbsp;++sp)
	{
	for(unsigned j=0;j<=n;++j)
	{
	//	double tmp=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(n,j);
	double tmp2=0;
	tfs>>tmp2;
	//if(tmp2!=0 && nabs((tmp-tmp2)/tmp2)>0.1 && tmp2>1E-20  )
	//	Log::mL<<"|"<<tmp<<" "<<tmp2<<"|\t";
	(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(n,j)=tmp2;
	//gfs>>(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(n,j);
	//		Log::mL<<0.000000000<<"\t";
	//	}
	}
	}
	//Log::mL<<endl;


	}
	Log::mL<<"suite"<<endl;
	tfs.close();
	//double omdeg=(*mpSp)[sp]->mpElecCrs->mOmDegradcm2(ene,vPosEner+1);
	ifstream ifs("momdeg");
	for(unsigned i=0;i<nen;++i)
	{
	for(unsigned j=0;j<nen;++j)
	{
	ifs>>(*mpSp)[0]->mpElecCrs->mOmDegradcm2(i,j);

	}
	}
	ifstream lfs("ctot.dat");
	for(unsigned i=0;i<mNben;++i)
	{
	for(unsigned j=0;j<2*mNbalt-1;++j)
	{
	double tmp=0;
	}

	}
	*/

	// We start a timer here
	timer t0;

	// We start a progress display
	//Log::SetPriority(Log::INFO);
	progress_display loop_disp(mNben-1,cout,"\n Energy Loop ","             ","             ");

	//qint*=0.; qint is already clear!
	for(unsigned n=0;n<mNben-1;++n)
	{
		Mstream(n,qint,qintensity);
		// And we call qmstr
		QMstream(n,qint,qintensity);
		/*
		ublas::vector<double> tmpintensityenergy(mNbalt);
		double tmptot_e_intensity=0;
		for(unsigned i=0;i<mNbalt;++i)
		{
			for(unsigned j=0;j<mNben;++j)
			{
				tmpintensityenergy[i]+=qintensity(j)(2*i,mNbang2)*(*mpElecC)[j]*(*mpElecD)[j];
			}
		}

		tmptot_e_intensity=MathFunction::TrapzInt(*mpAltGridKm,tmpintensityenergy)*1E5;
//		Log::mL<<"Intensite energy "<<tmptot_e_intensity<<endl;		
		Log::mL<<n<<"\tIntensite energy "<<tmptot_e_intensity<<endl;		
		*/
		for(unsigned k=0;k<mNbang+1;++k)
		{

			qint(2*(mNbalt-1),k)=0.;
		}
#ifdef DEBUG
		//		qint*=0;
#endif
		++loop_disp;
	}
	/*
	   ifstream jfs("matrixintensity");

	   for(unsigned m=0;m<2*mNbalt-1;++m)
	   {
	   for(int n=(int)mNben-1;n>-1;--n)
	   {
	   for(unsigned k=0;k<mNbang+1;++k)
	   {
	   if(k==mNbang2)
	   {
	   double tmp;
	   jfs>>tmp;
	   cout<<" old intensity "<<qintensity(n)(m,k)<<" <-> "<<tmp<<endl;;
	   qintensity(n)(m,k)=tmp;
	   }else
	   {
	   double tmp;
	   jfs>>tmp;
	   }
	   }
	   }
	   }
	   jfs.close();
	   */	
	//Log::SetPriority(Log::INFO,"ComputeIonization");
	Log::mI<<"Loop execution time : "<<t0.elapsed()<<endl;
	/*
	   Log::mL<<"The intensity"<<endl;
	   for(unsigned n=0;n<nen-1;++n)
	   {
	   Log::mL<<"energy :"<<(*(mpElecFlux->mpElecCentEeV))[n]<<endl;
	   Log::mL<<qintensity[n]<<endl;
	   Log::mL<<"=========================================="<<endl;
	   Log::mL<<"=========================================="<<endl;
	   Log::mL<<"=========================================="<<endl;

	   }
	   */	
	ublas::vector<double> intensityenergy(mNbalt);
	intensityenergy.clear();
	
	
	mNumberDensity1.resize(mNbalt);
	mNumberDensity2.resize(mNbalt);
	mNumberFlux1.resize(mNbalt);
	mNumberFlux2.resize(mNbalt);

	
	double tot_e_intensity=0;
	for(unsigned i=0;i<mNbalt;++i)
	{
		mNumberDensity1(i) = 0;
		mNumberDensity2(i) = 0;
		mNumberFlux1(i) = 0;
		mNumberFlux2(i) = 0;


		for(unsigned j=0;j<mNben;++j)
		{
			intensityenergy[i]+=qintensity(j)(2*i,mNbang2)*(*mpElecC)[j]*(*mpElecD)[j];
			
			mNumberDensity1(i) += 1.6860e-8 * qintensity(j)(2*i,mNbang2) * (*mpElecD)[j];
			mNumberDensity2(i) += 2 * 1.6860e-8 * mFluxHemisphericTot(j,i)  / sqrt((*mpElecC)[j])  * (*mpElecD)[j];
			mNumberFlux1(i) += mFluxHemisphericNet(j, i)   * (*mpElecD)[j];
			
			for(unsigned k = 0 ; k < mNbang2; ++k)
			{
				mNumberFlux2(i) += 2 * M_PI *  qintensity(j)(2*i,k) * (mpGAngle->mWeight[k]) * mpGAngle->mXmu[k]  * (*mpElecD)[j];
				mNumberFlux2(i) += 2 * M_PI *  qintensity(j)(2*i,k + mNbang2 + 1) * (mpGAngle->mWeight[k+ mNbang2]) * mpGAngle->mXmu[k+ mNbang2]  * (*mpElecD)[j];
			}
		}
	}

	tot_e_intensity=MathFunction::TrapzInt(*mpAltGridKm,intensityenergy)*1E5;




	//	PrintQintensity("qintensity_flux",qintensity);
	// We compute the ntherm parameter
	//
	//Log::mL<<"we compute the ntherm parameter"<<endl;
	for(unsigned i=0;i<mNbalt;++i)
	{
		ntherm(i)=0;
		double etherm=(*mpElectronTemperatureK)[i]*BOLTZMANNJ_K/eV_IN_JOULE;
		double denelc=(*mpElectronDensitycm_3)[i];

		//		Log::mL<<"altitude "<<i<<" etherm "<<etherm<<" denelc "<<denelc<<endl;
		if(denelc<1E-30||etherm<0)
		{
			ntherm(i)=mNben-1;
		}else
		{
			if(false)
			{// Test avec Trans* initial
				double emass=5.69E-16;// boltz/eV 10E-4
				unsigned nmax=PosEner(10);
				unsigned nn=PosEner(etherm,nmax);

				double prefix=2*denelc/(sqrt(emass))*pow(2*PI*etherm,-1.5);
				ntherm(i)=mNben;
				for(unsigned j=nmax;j<=nn;++j)
				{// Inverted in comparison with Trans*
					double E=(*(mpElecFlux->mpElecCentEeV))[j];
					double flx=prefix*E*exp(-E/etherm);

					if(flx>=qintensity(j)(2*i,mNbang2))
					{
						ntherm(i)=j;
						break;// break of the for loop
					}
				}
				if(ntherm(i)==mNben)
				{
					Log::mD<<"Pb with ncross!!!"<<endl;
					ntherm(i)=mNben-1;
				}
			}

			if(true)
			{
				unsigned nmax=PosEner(10);
				unsigned nn=PosEner(etherm,nmax);
				//				Log::mL<<"thermal energy : "<<etherm<<" for altitude "<<i<<" Giving the energy position "<<nn<<" and nmax "<<nmax<<endl;
				double prefix=denelc/(2*PI);
				ntherm(i)=mNben;
				for(unsigned j=nmax;j<=nn;++j)
				{// Inverted in comparison with Trans*
					double E=(*(mpElecFlux->mpElecCentEeV))[j];
					double flx=prefix*E*exp(-E/etherm)*sqrt(E/(PI*pow(etherm,3)));

					if(flx>=qintensity(j)(2*i,mNbang2))
					{
						ntherm(i)=j;
						break;// break of the for loop
					}
				}
				if(ntherm(i)>=mNben-1)
				{
		//			if(i==0)// Write only one time
		//				Log::mL<<"Pb with ncross; but ok (check nlev with jean)!!!"<<endl;
					ntherm(i)=mNben-1;
				}



			}
		}
	}

	// We compute the creation of species
//	Log::mL<<"ntherm="<<ntherm<<endl;
	
	//ublas::vector<double> newaltgrid = //ublas::inner_prod(*mpAltGridKm, mIsmgdpa);
	//ublas::vector<double> newaltgrid = (*mpAltGridKm) * mIsmgdpa;
	mTotalEnergy=mpElecFlux->FluxEnergy(*mpAltGridKm);
	Log::mI<<".........................................."<<endl;
	Log::mI<<"Total energy in the flux "<<tot_e_intensity<<endl;
	Log::mI<<".........................................."<<endl;
	Log::mI<<"Total energy :"<<mTotalEnergy<<endl;


	mReflectedEnergy=ReflectedEnergy(ntherm);
	Log::mI<<"Reflected energy = "<<mReflectedEnergy<<endl;
#ifdef ANALYSE_ELEC_FLUX		
	Log::mI<<"Reflected energy absolute value = "<<mTestFlux<<endl;
#endif
	mTransmittedEnergy=TransmittedEnergy(ntherm);
	Log::mI<<"Transmitted energy = "<<mTransmittedEnergy<<endl;

	mAbsorbedEnergy=Heating(ntherm,qintensity);
	Log::mI<<"Absorbed energy = "<<mAbsorbedEnergy<<endl;

	//ublas::zero_vector<double> zdep(nbsp);
	ublas::vector<double> inelastic_energy_dep(mNbsp);
	inelastic_energy_dep.clear();
	ublas::matrix<double> enrate;

//	Log::mL<<"We compute the creation of the species"<<endl;
	IonizeSpecies(qintensity,ntherm,inelastic_energy_dep,enrate);

	mInelasticEnergy=ublas::sum(inelastic_energy_dep);


	//we print the deposition
	Log::mI<<"Inelastic energy deposition : "<<inelastic_energy_dep<<endl;
	Log::mI<<"Inelastic energy deposition total : "<<mInelasticEnergy<<endl;
	ComputeDeposition();

	/*
	   ublas::vector<double> intensalt(mNbalt);
	   intensalt.clear();
	   for(unsigned i=0;i<mNbalt-1;++i)
	   {
	   for(unsigned j=0;j<qintensity.size();++j)
	   {
	   double intens=0;
	   for(unsigned k=0;k<mNbang2;++k)
	   {
	   intens+=qintensity[j](2*i+1,k)*mpElecFlux->ReturnAngleWeight(k);
	   intens+=qintensity[j](2*i+1,k+1+mNbang2)*mpElecFlux->ReturnAngleWeight(k+mNbang2);

	   }
	   intensalt[i]=intens*(*mpElecC)[j]*(*mpElecD)[j];
	//intensalt[i]=qintensity[j](2*i+1,mNbang2)*(*mpElecC)[j]*(*mpElecD)[j];
	}
	cout<<intensalt[i]<<endl;
	}

	double totqintenergy=MathFunction::TrapzInt(*mpAltGridKm,intensalt)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV
	*/	

	//	Log::mL<<"Energy inside the qint : "<<totqintenergy<<endl;


	SpecieUtils::SpeciesToResu(*mpSp,rResult);
	//	Log::mL<<"Vector of species : size ="<<vSp.size()<<endl;
	//	Log::mL<<"Vector of resu : size ="<<rResult.size()<<endl;
}


void ElectronImpactIonization::ComputeDeposition()
{
	mResultTotalEnergy=mAbsorbedEnergy+mInelasticEnergy;
	mElos1=0.;
	mElos2=0.;
	mElos3=0.;
	mIlos1=0.;
	mIlos2=0.;
	mIlos3=0.;
	if(mProdTot>0)
	{
		mElos1=mResultTotalEnergy/mProdTot;
		mElos2=(mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdTot;
		mElos3=mTotalEnergy/mProdTot;
		mIlos1=mResultTotalEnergy/mProdIonTot;
		mIlos2=(mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdIonTot;
		mIlos3=mTotalEnergy/mProdIonTot;
	}
	Log::mI<<"============================================="<<endl;
	Log::mI<<" Total input energy: "<<mTotalEnergy<<endl;
	Log::mI<<"---------------------------------------------"<<endl;
	Log::mI<<" Total Heating:      "<<mAbsorbedEnergy<<endl;
	Log::mI<<" Total inelastic:    "<<mInelasticEnergy<<endl;
	Log::mI<<" Total reflected:    "<<mReflectedEnergy<<endl;
#ifdef ANALYSE_ELEC_FLUX		
	Log::mI<<" Total reflected abs:"<<mTestFlux<<endl;
#endif
	Log::mI<<" Total transmitted:  "<<mTransmittedEnergy<<endl;
	Log::mI<<"---------------------------------------------"<<endl;

	Log::mI<<" Total absorbed energy (inelastic + heating) "<<mResultTotalEnergy<<endl;
	Log::mI<<" Total computed energy (tot abs   + reflect) "<< mResultTotalEnergy+mReflectedEnergy+mTransmittedEnergy<<endl;
	mError= -(mTotalEnergy -(mResultTotalEnergy+mReflectedEnergy+mTransmittedEnergy))/mTotalEnergy;
	Log::mI<<"---------------------------------------------"<<endl;
	Log::mI<<"\t Energy conservation error :  "<<ntostr(mError*100.)<<"% (>0 : gain)"<<endl;
	Log::mL<<"\t elos (1,2,3)"<<mElos1<<" "<<mElos2<<" "<<mElos3<<endl;
	Log::mL<<"\t ilos (1,2,3)"<<mIlos1<<" "<<mIlos2<<" "<<mIlos3<<endl;
	Log::mL<<"\t Production total: "<<mProdTot<<endl;
	Log::mL<<"\t Production total, ions: "<<mProdIonTot<<endl;
	Log::mI<<"---------------------------------------------"<<endl;
	double tot = mSecondarySumeV+mDegradedSumeV+mCoulombianSumeV;
	Log::mL<<"Energy from the secondaries      : "<<mSecondarySumeV<<" = "<< mSecondarySumeV*100./tot<<"%"<<endl; 
	Log::mL<<"Energy from the degraded      : "<<mDegradedSumeV<<" = "<< mDegradedSumeV*100./tot<<"%"<<endl;  
	Log::mL<<"Energy from the coulombian      : "<<mCoulombianSumeV<<" = "<< mCoulombianSumeV*100./tot<<"%"<<endl;


	Log::mI<<"============================================="<<endl;
	Log::AddVersionInfo("# Electron impact ionization error: "+ntostr(mError*100.)+"% (>0 : gain)");
}


void ElectronImpactIonization::IonizeSpecies(ublas::vector< ublas::matrix<double> > vQIntensity,ublas::vector<unsigned> vPositionThermic,ublas::vector<double>& rInelasticDepeV,ublas::matrix<double>& rEnrate)
{

	int nang=(mpElecFlux->NbAngles());
//	Log::mL<<"Ionization of the different species"<<endl;
	unsigned nalt=mpElectronDensitycm_3->size();
	ublas::zero_matrix<double> zenrate(mpSp->size(),nalt);
	rEnrate.resize(mpSp->size(),nalt);
	rEnrate.clear();
	ublas::matrix<double> frate(zenrate);
	ublas::vector<double> vect_elec_production(nalt);// To store the number of electron produced
	ublas::vector<double> vect_ion_production(nalt);// To store the number of ions produced
	frate.clear();
	vect_elec_production.clear();
	vect_ion_production.clear();
	for(unsigned i=0;i<nalt;++i)
	{
		//for(unsigned ene=nen-1;ene>=vPositionThermic[i];--ene)
		unsigned posmax=vPositionThermic[i];
		//	if(vPositionThermic[i]==(nen-1))
		//		--posmax;
		for(unsigned ene=0;ene<=posmax;++ene)
		{
			double energy=(*(mpElecFlux->mpElecCentEeV))[ene];
			double engdd=(*(mpElecFlux->mpElecDdengeV))[ene];
			std::deque<Specie*>::iterator it;
			unsigned sp=0;
			for(it=mpSp->begin();it!=mpSp->end();++it)
			{
				double dens=(*it)->mTotDensitycm_3[i];
				if( (*mpSp)[sp]->CheckElecCrs() && (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
				{
					for(unsigned j=0;j<(*it)->mProcessNames.size();++j)
					{

						double threshold=(*it)->mpElecCrs->mThresholdseV[j];

						if(energy>threshold)
						{
							double cinex=(*it)->mpElecCrs->mCrscm2[j][ene];
							double produ=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex)*dens;
							(*it)->mSpeciesProductioncm_3s_1[j][i]+=produ;

						}
					}

					for(unsigned j=0;j<(*it)->mpElecCrs->mExcitationCrsPosition.size();++j)
					{
						unsigned position=(*it)->mpElecCrs->mExcitationCrsPosition[j];
						//		Log::mL<<"Position "<<position<<" of "<<j<<" / "<<(*it)->mpElecCrs->mExcitationCrsPosition.size()<<endl;
						double thresh=(*it)->mpElecCrs->mThresholdseV[position];
						//			double thresh=threshexc;
						double cinex=(*it)->mpElecCrs->mCrscm2[position][ene];
						double cross=thresh*cinex;
						rEnrate(sp,i)+=(vQIntensity(ene)(2*i,nang/2)*engdd*cross)*dens;
						frate(sp,i)+=(engdd*cross)*dens;

					}
					for(unsigned j=0;j<(*it)->mpElecCrs->mIonizationCrsPosition.size();++j)
					{
						unsigned position=(*it)->mpElecCrs->mIonizationCrsPosition[j];
						//		Log::mL<<"Position "<<position<<" of "<<j<<" / "<<(*it)->mpElecCrs->mIonizationCrsPosition.size()<<endl;
						double elnumber=((*it)->mpElecCrs->mNumberOfElectrons[position]);
						double thresh=(*it)->mpElecCrs->mThresholdseV[position];
						double cinex=(*it)->mpElecCrs->mCrscm2[position][ene];
						double cross=thresh*cinex;

						double produ=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex)*dens*elnumber;
						//mProdTot+=produ;
						vect_elec_production[i]+=produ;
						vect_ion_production[i]+=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex)*dens*((*it)->mpElecCrs->mNumberOfIons[position]);
						(*it)->mElecElecProductioncm_3s_1[i]+=produ;
						rEnrate(sp,i)+=(vQIntensity(ene)(2*i,nang/2)*engdd*cross)*dens;
						frate(sp,i)+=(engdd*cross)*dens;
						if((*it)->mpElecCrs->mIsAuger.at(position))
						{// If the auger process is defined
							for(unsigned z=0;z<(*it)->mpElecCrs->mAugerEnergy.at(position).size();++z)
							{
							//	unsigned energypos2=Search((*it)->mpPhotoCrs->mAugerEnergy.at(i).at(j));//  We search the position for the auger electron
								double augerenergy=(*it)->mpElecCrs->mAugerEnergy.at(position).at(z);
								rEnrate(sp,i)-=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex*augerenergy)*dens*(*it)->mpElecCrs->mAugerEfficiency.at(position).at(z);
								//rProelecE(iAlt,energypos2)+=produ/(*mpEDdengeV)[energypos2]*(*it)->mpPhotoCrs->mAugerEfficiency.at(i).at(j) ;
								(*it)->mElecElecProductioncm_3s_1[i]+=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex)*dens*(*it)->mpElecCrs->mAugerEfficiency.at(position).at(z);
								vect_elec_production[i]+=(vQIntensity(ene)(2*i,nang/2)*engdd*cinex)*dens*(*it)->mpElecCrs->mAugerEfficiency.at(position).at(z);
							}
						}
					}
				}
				++sp;
			}
		}
		//	Log::mL<<"altitude "<<i<<" enrate for sp 0 "<<enrate(0,i)<<endl;
		//	Log::mL<<"altitudekm "<<(*mpAltGridKm)[i]<<" frate for sp 0 "<<frate(0,i)<<endl;
	}

	for(unsigned i=1;i<nalt;++i)
	{
		double ddz=((*mpAltGridKm)[i-1]-(*mpAltGridKm)[i])*mIsmgdpa[i]*1E5;// 1E5 : km to cm
		for(unsigned sp=0;sp<mpSp->size();++sp)
		{
			rInelasticDepeV(sp)+=(rEnrate(sp,i)+rEnrate(sp,i-1))*ddz/2.;
		}
		mProdTot+=(vect_elec_production(i)+vect_elec_production(i-1))*ddz/2.;
	}
	// This other techniques works also very well
	//double prodtot2=MathFunction::TrapzInt(*mpAltGridKm,vect_elec_production)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
	mProdIonTot=MathFunction::TrapzInt(*mpAltGridKm,vect_ion_production)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
//	Log::mL<<"Old version mProdTot: "<<mProdTot<<endl;//<<"  new version : "<<prodtot2<<endl;
//	Log::mL<<"Ion production : "<<mProdIonTot<<endl;
	for(unsigned sp=0;sp<mpSp->size();++sp)
	{
	//	rInelasticDepeV(sp)+=(rEnrate(sp,i)+rEnrate(sp,i-1))*ddz/2.;
		ublas::matrix_row< ublas::matrix<double> > ro(rEnrate,sp);
//		ublas::matrix<double> machin(ro);
		Log::mD<<" Sp: "<<sp<<" energy dep : "<< MathFunction::TrapzInt(*mpAltGridKm,ro)*1E5<<endl;
	}
//	Log::mL<<"fin truc"<<endl;
}


double ElectronImpactIonization::Heating(ublas::vector<unsigned> vNtherm,const ublas::vector< ublas::matrix<double> >& vQIntensity )
{
	double elhsum=0.;
	unsigned nalt=mpElectronDensitycm_3->size();
	unsigned nen=vQIntensity.size();
	//vector<double>* egrid = mpElecFlux->mpElecCentEeV;
	//vector<double>* wgrid = mpElecFlux->mpElecDdengeV;

	unsigned nang2=mNbang2;//static_cast<unsigned>((mpElecFlux->NbAngles())/2);

	ublas::vector<double> Qetherm(nalt),Qeflux(nalt);//lossener(nalt);
	Qetherm.clear();
	mQe.resize(nalt);
	mQe.clear();
	Qeflux.clear();
	// lossener is created to check what is the absorbed energy lose in our computational process
//	lossener.clear();

	for(unsigned i=0;i<nalt;++i)
	{
		unsigned nn=vNtherm(i);
		if(nn>nen-1)
			nn=nen-1;
		double etherm=(*mpElectronTemperatureK)[i]*BOLTZMANNJ_K/eV_IN_JOULE;
		Qetherm(i)=mEnergyLosseV(i,nn)*( (*mpElecC)[nn]-1.5*etherm)*vQIntensity(nn)(2*i,nang2);
		if(Qetherm[i]<0)
			Qetherm[i]=0;
		//Log::mL<<"Qetherm"<<Qetherm(i)<<endl;
		//Log::mL<<( (*egrid)[nn]-1.5*etherm)<<endl;
		//assert(!(Qetherm(i)<0));
		for(unsigned n=0;n<=vNtherm[i];++n)
		{
			Qeflux(i)+=mEnergyLosseV(i,n)*vQIntensity(n)(2*i,nang2)*(*mpElecD)[n];
		}
		mQe(i)=Qeflux(i)+Qetherm(i);
	}
	for(unsigned i=1;i<nalt;++i)
	{
		double ddz=((*mpAltGridKm)[i-1]-(*mpAltGridKm)[i])*mIsmgdpa[i]*1E5;// 1E5 : km to cm
		elhsum+=(mQe(i-1)+mQe(i))*ddz/2.;
	}


	return elhsum;
}

double ElectronImpactIonization::ReflectedEnergy(ublas::vector<unsigned> vNtherm)
{
	double sum,sume;
	sum=0;
	sume=0;
	IntFluxP(vNtherm(0),0,mFluxHemisphericUp,0,sum,sume);
	
#ifdef ANALYSE_ELEC_FLUX		
	double suma;
	suma=0;
	IntFluxPa(vNtherm(0),0,mFluxHemisphericUp,0,sum,suma);
	mTestFlux = suma;
#endif
	return sume;
}

double ElectronImpactIonization::TransmittedEnergy(ublas::vector<unsigned> vNtherm)
{
	double sum,sume;
	sum=0;
	sume=0;
	IntFluxP(vNtherm(0),0,mFluxHemisphericDown,mNbalt-1,sum,sume);
	Log::mL<<"Transmitted (without taking albedo into account : "<<sume<<" and with albedo "<<sume*(1-mAlbedo)<<endl;
	return sume*(1-mAlbedo);
}

void ElectronImpactIonization::Init()
{
	// We call ComputeEnergyLoss
	//Log::SetPriority(Log::DEBUGG);
	//	Log::mL<<"We call the compute energy loss function"<<endl;

	mpParameter->ExistsOrDie("/aero_main/electron/albedo","You should define the albedo");
	mpParameter->GetValue("/aero_main/electron/albedo",mAlbedo);
	if(mpParameter->Exists("/aero_main/electron/no_disort_check"))
	{
		mIsDisortChecking=0; // in that case, disort does not check!
	}

	if(mpParameter->Exists("/aero_main/electron/LEfactor"))
	{
		mpParameter->GetValue("/aero_main/electron/LEfactor",mFactUncertLE);
	}
	ComputeEnergyLoss();
}


void ElectronImpactIonization::ComputeEnergyLoss()
{
	// GG: jan 2011, comments following the writing of the uncertainty paper
	// Simple computation of the L(E) parameter
	
	
	//Log::SetPriority(Log::DEBUGG);
	unsigned nalt=mpElectronDensitycm_3->size();
	assert(nalt==mNbalt);
	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	
	//double boltzev=BOLTZMANNJ_K/eV_IN_JOULE;
	double boltzev = BOLTZ_DIV_eV;
	mEnergyLosseV.resize(nalt,nen);
	mEnergyLosseV.clear();
	for(unsigned i=0;i<nalt;++i)
	{
		//		Log::mL<<"Electron temp "<<(*mpAltGridKm)[i]<<"  "<<(*mpElectronTemperatureK)[i]<<endl;;
		//		Log::mL<<"Electron Dens "<<i<<"  "<<(*mpElectronDensitycm_3)[i]<<endl;;
		double eth=boltzev* (*mpElectronTemperatureK)[i];// Thermal energy
		for(unsigned j=0;j<nen;++j)
		{
			double eng=(*(mpElecFlux->mpElecCentEeV))[j];//energy
			mEnergyLosseV(i,j)=0;
			if( (eng<150) && (eng>eth))
			{
				double fac= mFactUncertLE * pow((eng-eth)/(eng-0.53*eth),2.36);
				if(mFactUncertLE > 1.1)
					std::cout<<"MERDE MERDE MERDE "<<mFactUncertLE<<endl;
				if(fac>0)
				{
					mEnergyLosseV(i,j)=fac*3.37E-12*pow((*mpElectronDensitycm_3)[i],0.97)/pow(eng,0.94);
				}
#ifdef DEBUG
				if(isnan(mEnergyLosseV(i,j)))
				{
					Log::mW<<"We've got a pb, mEnergyLosseV negative"<<endl;
				}
#endif 
				assert(!isnan(mEnergyLosseV(i,j)));

			}
		}
	}
}


void ElectronImpactIonization::ComputePhaseFunction()
{
	// GG: jan 2011, comments following the writing of the uncertainty paper
	// Here, we compute the phase function for the elastic scattering.
	// It corresponds to the p in equation 3.2 of Lummerzheim phd thesis
	
	unsigned nlayer=mNbalt-1;// nalt-1
	unsigned nen=mNben;
	unsigned nangle=mNbang+1;
/*
	vector<double> tmp0(nlayer,0);
	vector< vector<double> > tmp1(nangle,tmp0);
	mPhaseFunction.resize(nen,tmp1);*/

	ublas::matrix<double> tmp(nlayer,nangle);
	tmp.clear();
	mPhaseFunction.resize(nen);
	for(unsigned i=0;i<nen;++i)
	{
		mPhaseFunction(i)=tmp;
	}


//	Log::mL<<"Init porter"<<endl;
	InitPorter();// We initialize the Porter function!
//	Log::mL<<"Init ruther"<<endl;
	InitRuther();// We initialize the Rutherford function!

//	Log::mL<<"Fill the phase function"<<endl;
	for(unsigned e=0;e<nen;++e)
	{
	//	Log::mL<<"energy size :
		double eng=((*(mpElecFlux->mpElecCentEeV))[e]);
	//	Log::mL<<" The energy :"<<eng<<endl;
		ublas::vector<double> gls;// To store the phase function
		if(eng<RUTHERFORD_PORTER_THRESHOLD)
		{
		//	Log::mL<<"energy < RUTHERFORD"<<endl;
			for(unsigned i=0;i<nlayer;++i)
			{
				gls=Porter(i,e);
				mPhaseFunction[e](i,0)=1;
				for(unsigned k=1;k<nangle;++k)
				{
				//	mPhaseFunction[e][k][i]=gls[k];
					if(gls[k]<0)
					{
					//	Log::mL<<"Error Porter  gls<0"<<endl; // ok gls<0 is possible
						gls[k]=0.;
					}
					assert(!(gls[k]>1.));
					mPhaseFunction(e)(i,k)=gls[k];
				}
			}
		}else
		{
			for(unsigned i=0;i<nlayer;++i)
			{
				mPhaseFunction[e](i,0)=1;
				for(unsigned k=1;k<nangle;++k)
				{
					//mPhaseFunction[e][k][i]=mRutherford[e][k];
					if(mRutherford(e,k)<0)
					{
						Log::mD<<"Error Rutherford gls<0"<<endl;
						mRutherford(e,k)=0.;
					}
					assert(!(mRutherford(e,k)>1.));
					mPhaseFunction(e)(i,k)=mRutherford(e,k);
				}
			}
		}
	}
//	Log::mL<<"Finish the phase function"<<endl;
}


void ElectronImpactIonization::InitPorter()
{
	// Fill the porter vector [specie][energy][angle]
//	Log::mL<<"Welcome to porter!"<<endl;
	unsigned nbsp=mNbsp;
	unsigned nen=mNben;
	unsigned nangle=mNbang+1;
	// Default : isotropic
	//vector<double> tmp0(nangle,0);
	//tmp0[0]=1;

	//vector< vector<double> > tmp1(nen,tmp0);
	mPorter.resize(nbsp);//,tmp1);
	
	ublas::matrix<double> tmp(nen,nangle);
	tmp.clear();
	for(unsigned i=0;i<nbsp;++i)
	{
		mPorter[i]=tmp;
	}

	Log::mI<<"We load the species in Porter"<<endl;
	for(unsigned sp=0;sp<nbsp;++sp)
	{
		(*mpSp)[sp]->LoadPorter();
//		Log::mL<<"Specie = "<<(*mpSp)[sp]->mName<<endl;
	//	Log::mL<<" We return the specie!"<<(*mpSp)[sp]->mName<<endl;
		ublas::vector<double>* eporter=(*mpSp)[sp]->ReturnEPorter();
	//	vector<double>* gporter=(*mpSp)[sp]->ReturnGPorter();
		ublas::vector<double>* gporter=(*mpSp)[sp]->ReturnGPorter();
	//	vector<double>* bporter=(*mpSp)[sp]->ReturnBPorter();
		ublas::vector<double>* bporter=(*mpSp)[sp]->ReturnBPorter();
	//	vector<double>* aporter=(*mpSp)[sp]->ReturnAPorter();
		ublas::vector<double>* aporter=(*mpSp)[sp]->ReturnAPorter();
		unsigned nbporter=eporter->size();
		assert(gporter->size()==nbporter);
		assert(bporter->size()==nbporter);
		assert(aporter->size()==nbporter);

		assert(eporter->size()>0);
		// eporter should be increasing!
		assert((*eporter)[0]<(*eporter)[nbporter-1]);

		// HEHEHE OPTIMIZATION!
		int posporter=nbporter-1;

		for(unsigned e=0;e<nen;++e)
		{// The energy is decreasing here!
			double ene=(*(mpElecC))[e];

			if(ene<(*eporter)[nbporter-1])
			{
				while( posporter>0 &&( ( ene > (*eporter)[posporter] ) || (ene< (*eporter)[posporter-1]) ) )
				{
					--posporter;
				}
				if(posporter==0&&( ene > (*eporter)[posporter] ))
				{
					break;
				}
			}
		//	if(posporter<0)
		//	{
//				Log::mL<<"Break for energy = "<<ene<<endl;
		//		break;
		//	}
			ublas::vector<double> phase_function;
		//	cout<<"123123123123131232"<<endl;
			if(nangle>3) // 2 + 1
			{
				phase_function=RecDown(nangle,(*aporter)[posporter],(*bporter)[posporter],(*gporter)[posporter]);
		//		Log::mL<<(*aporter)[posporter]<<" "<<(*bporter)[posporter]<<" "<<(*gporter)[posporter]<<endl;
//				Log::mL<<"recdown"<<endl;
			}else
			{
				phase_function=RecUp(nangle,(*aporter)[posporter],(*bporter)[posporter],(*gporter)[posporter]);
			}

			for(unsigned k=0;k<nangle;++k)
			{
//				if(phase_function[k]<0.)
//					Log::mL<<"Error phase function"<<endl;
				mPorter(sp)(e,k)=phase_function[k];
			}
		}
	}
	Log::mD<<"End of the  Init porter function"<<endl;
}

void ElectronImpactIonization::InitRuther()
{
	// Fill the ruther vector [energy][angle]
	unsigned nen=mNben;
	unsigned nangle=mNbang+1;
	//vector<double> tmp0(nangle,0);
	//tmp0[0]=1;
	mRutherford.resize(nen,nangle);
	mRutherford.clear();
	for(unsigned i=0;i<nen;++i)
	{
		mRutherford(i,0)=1;
	}

	for(unsigned e=0;e<nen;++e)
	{
		double ene=(*(mpElecFlux->mpElecCentEeV))[e];
		if(ene>1000)
		{
			double corr=0.6*pow((1000/ene),0.09);
			double eps=6.22E-5/(2.+ene/511000.)/(ene/511000.)*corr;
			ublas::vector<double> phase_function=RecUp(nangle,0.,0.,eps,false);
			for(unsigned k=0;k<nangle;++k)
			{
				mRutherford(e,k)=phase_function[k];
			}
		}else if(ene>=477.)
		{
			double corr=0.6;
			double eps=6.22E-5/(2.+ene/511000.)/(ene/511000.)*corr;
			ublas::vector<double> phase_function=RecDown(nangle,0.,0.,eps,false);
			for(unsigned k=0;k<nangle;++k)
			{
				mRutherford(e,k)=phase_function[k];
			}
		}else if(ene>12)
		{
			double eps=0.5/(pow((ene/12.),0.75) -1.);
			ublas::vector<double> phase_function=RecDown(nangle,0.,0.,eps,false);
			for(unsigned k=0;k<nangle;++k)
			{
				mRutherford(e,k)=phase_function[k];
			}
		}
	}
}

ublas::vector<double> ElectronImpactIonization::Porter(unsigned vAlt,unsigned vEner)
{// Returns the computed Porter function for the considered altitude

	unsigned nbsp=mNbsp;
	unsigned nangle=mNbang+1;
	ublas::vector<double> resu(nangle);
	resu.clear();
	double totdensity=0;
	for(unsigned sp=0;sp<nbsp;++sp)
	{
		double density=((*mpSp)[sp])->mTotDensitycm_3[vAlt];
		totdensity+=density;
		for(unsigned k=0;k<nangle;++k)
		{
			resu[k]+=density*(mPorter(sp)(vEner,k));
		}
	}

	assert(totdensity>0);//Theoretically not a problem
	
	for(unsigned k=0;k<nangle;++k)
	{
		resu[k]/=totdensity;//+=density*mPorter[sp][vEner][k];
	}


	return resu;

}

ublas::vector<double> ElectronImpactIonization::RecUp(unsigned vNbAngle,double vAlpha,double vBeta, double vGamma,bool vBG)
{
	// Isotropic init
	ublas::vector<double> resu(vNbAngle);
	ublas::vector<double> g1(vNbAngle);
	ublas::vector<double> g2(vNbAngle);
	resu.clear();
	g1.clear();
	g2.clear();
	resu[0]=1;
	if(!vBG && vGamma<1E-42 && vGamma>-1E-42)
	{// To ensure 0 for vGamma!!!
		return resu;
	}
	double fnorm=1/(4*vGamma*(1+vGamma));
	if(vBG)
	{// Default : porter
		fnorm= (fnorm + vBeta/(vAlpha*(1+vAlpha)))/4;
	}

	g1[0]=1./(4*vGamma*(1+vGamma))/fnorm;
	g1[1]=(1.+2.*vGamma)/(4.*vGamma*(1.+vGamma))-.5*log(1.+1./vGamma);
	g1[1]/=fnorm;

	if(vBG)
	{
                g2[0]=vBeta/(4.*vAlpha*(1.+vAlpha))/fnorm;
                g2[1]=(1.+2.*vAlpha)/(4.*vAlpha*(1.+vAlpha))-.5*log(1.+1./vAlpha);
	}
	for(unsigned k=2;k<vNbAngle;++k)
	{
		g1[k]=((2.*k-1.)*(1.+2.*vGamma)*g1[k-1]-k*g1[k-2])/(k-1.);
		if(vBG)
		{
			g2[k]=(-(2.*k-1.)*(1.+2.*vAlpha)*g2[k-1]-k*g2[k-2])/(k-1.);
		}
	}

	for(unsigned k=0;k<vNbAngle;++k)
	{
		resu[k]=g1[k]+g2[k];
	}
	return resu;



}


ublas::vector<double> ElectronImpactIonization::RecDown(unsigned vNbAngle,double vAlpha,double vBeta, double vGamma,bool vBG)
{
	// Isotropic init
	const unsigned madd=9;
	const unsigned maxstr=32;
	if(vNbAngle>maxstr+1)
	{
		Error err("RecDown","stream","too many streams");
		throw err;
	}
	unsigned infty=madd+maxstr;
	unsigned mdim=madd+maxstr+1;
//	unsigned amax=infty;//vNbAngle+madd;
	ublas::vector<double> resu(vNbAngle);
	ublas::vector<double> g1(mdim);
	ublas::vector<double> g2(mdim);
	resu.clear();
	g1.clear();
	g2.clear();
	resu[0]=1.;
	if(!vBG && vGamma<1E-42 && vGamma>-1E-42)
	{// To ensure 0 for vGamma!!!
		return resu;
	}
	double fnorm=1/(4*vGamma*(1+vGamma));
	double fnorm1=1.;
	double fnorm2=0;
	if(vBG)
	{// Default : porter
		fnorm+= (vBeta/(vAlpha*(1+vAlpha)))/4.;
		fnorm1=1./(4.*vGamma*(1.+vGamma))/fnorm;
		fnorm2=vBeta/(4.*vAlpha*(1.+vAlpha))/fnorm;
	}

	double del=1./pow(static_cast<double>(vNbAngle-1+madd),4); // vNbAngle is equal to the number of stream + 1!!!

	bool boucle=true;
	while(boucle)
	{
		boucle=false;// Theoretically 1 loop!
		g1[infty]=0;
		g2[infty]=0;
		g1[infty-1]=del;

		if(vBG)
		{
			g2[infty-1]=del;
		}

		for(unsigned k=infty-1;k>0;--k)
		{
			g1[k-1]= ((2.*k+1.)*(1.+2.*vGamma)*g1[k]-k*g1[k+1])/(k+1.);
			if(vBG)
			{
				g2[k-1]=(-(2.*k+1.)*(1.+2.*vAlpha)*g2[k]-k*g2[k+1])/(k+1.);
			}
			if(g1[k-1]>1E10 or g2[k-1]>1e10)
			{
				--infty;
				if(infty<1)
				{
					Error err("RecDown","Boucle","Error in recursion, too many iterations");
					throw err;
				}
				boucle=true;// One more loop
				break;
			}
		}
		
	}
	for(int i=vNbAngle-1;i>-1;--i)
	{
		g1[i]*=fnorm1/g1[0];
		if(vBG)
		{
			g2[i]*=fnorm2/g2[0];
		}else
		{
			g2[i]=0.;
		}
	}

	for(unsigned k=0;k<vNbAngle;++k)
	{
		resu[k]=g1[k]+g2[k];
		//,double vAlpha,double vBeta, double vGamma,bool vBG)
		assert(!(resu[k]>1.1));
/*		if(resu[k]<0) // Ok, possible
		{
			Log::mL<<"probleme grave sur RESU!"<<endl;
		}
		*/
	}
/*
	for(unsigned k=0;k<vNbAngle;++k)
	{
		Log::mL<<resu[k]<<"  ";
	}
	Log::mL<<endl;
	*/
	return resu;



}


void ElectronImpactIonization::ComputeCTot()
{
	Log::mD<<"Welcome to compute CTot"<<endl;
	// GG: jan 2011, comments following the writing of the uncertainty paper
	// mCTot corresponds to the opacity generalized (containing the Coulombian collisions)
	// it corresponds to the sum in equation 3.1 in the Lummerzheim phd thesis.
	
	
	
	unsigned nalt=mpElectronDensitycm_3->size();
	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	unsigned nbsp=mpSp->size();
	//	mCTot.resize(nalt*2-1,0.);
	mCTot.resize(nen,2*nalt-1);
	mCTot.clear();
	//inutile
	/*mCTot*=0.;
	for(unsigned i=0;i<nen;++i)
	{
		for(unsigned j=0;j<2*nalt-1;++j)
		{
			mCTot(i,j)=0;
		}
	}*/
	// Inverse of the sinus magnetic dip angle!
	ublas::vector<double> ismgdpa=ReturnInverseSinMagneticDipAngle();

	mIsmgdpa=ismgdpa;
	//ublas::zero_matrix<double> ztau(nen,nalt-1);// temporary : collision depth
	ublas::matrix<double> tau(nen,nalt-1);// temporary : collision depth
	tau.clear();
	// actually, this is the real collision depth, but
	// utau (differential in tau) and dtauc are used in the code
	//ublas::matrix<double> twork(nen,nalt);
	mDiffTau.resize(nen,nalt-1);
	mDiffTau.clear();
	mUTau.resize(nen,2*nalt-1);
	mUTau.clear();
//	mOmegaFunction.resize(nen,nalt);
//	mOmegaFunction*=0.;
	mSsalb.resize(nen,nalt-1);
	mSsalb.clear();
	/*
	for(unsigned en=0;en<nen;++en)
	{
		for(unsigned w=0;w<nalt-1;++w)
		{
			mSsalb(en,w)=0.;
			mDiffTau(en,w)=0;
			mUTau(en,2*w)=0;
			mUTau(en,2*w+1)=0;
		}
	}*/
	
	//ublas::zero_vector<double> ztwork(nalt);// temporary for tau
	ublas::vector<double> twork(nalt);
	twork.clear();
	for(unsigned ene=0;ene<nen;++ene) // I think that i can put nen-1 instead of nen... To be checked!!! Because I compute a lot of things here
	{
	//	Log::mL<<"Energy : "<<ene<<"/ "<<nen<<" = "<<(*(mpElecFlux->mpElecCentEeV))[ene]<<endl;
		twork.clear();
		double tau0=0;
		for(unsigned i=0;i<nalt;++i)
		{
			//		Log::mL<<"ialt ="<<i<<endl;
			unsigned m=i*2;
			mCTot(ene,m)=mEnergyLosseV(i,ene)/(*(mpElecFlux->mpElecDdengeV))[ene];

			assert(!(mEnergyLosseV(i,ene)<0));
			assert((*(mpElecFlux->mpElecDdengeV))[ene]>0.);
			//assert(mCTot(ene,m)>0);
			assert(!(mCTot(ene,m)<0));
			unsigned mm=i*2+1;
			//		Log::mL<<"mctot"<<endl;
			if(i!=nalt-1)
			{
				mCTot(ene,mm)=(mEnergyLosseV(i,ene)+mEnergyLosseV(i+1,ene))/(2.*(*(mpElecFlux->mpElecDdengeV))[ene]);
#ifdef DEBUG
				if((mCTot(ene,mm)<0.))
				{
					Log::mD<<" Merde : ctot ("<<ene<<","<<mm<<") = "<<mCTot(ene,mm)<<endl;
					Log::mD<<" energylossev i="<<i<<", ene; i+1, ene : "<<mEnergyLosseV(i,ene)<<"\t"<<mEnergyLosseV(i+1,ene)<<endl;
					Log::mD<<" elecddeng *2 : "<<(2.*(*(mpElecFlux->mpElecDdengeV))[ene])<<endl;
				}
				if(isnan(mEnergyLosseV(i+1,ene)))
				{
					Log::mD<<" Merde pour mEnergyLossseV("<<i+1<<" , "<<ene<<endl;
					Log::mD<<" Pour rappel imax = "<<nalt<<" et jmax ="<<nen<<endl;
				}
#endif
				
				assert(!isnan(mEnergyLosseV(i+1,ene)));
				assert(!isnan(mCTot(ene,mm)));
				assert(!(mCTot(ene,mm)<0.));

				// Tau
				twork(i+1)=twork(i)+mCTot(ene,mm)*((*mpAltGridKm)[i]-(*mpAltGridKm)[i+1])*1E5;
				// Be careful!!! we do a shift in Tau HERE!
				tau(ene,i)=twork(i+1);
			}


			for(unsigned sp=0;sp<nbsp;++sp)
			{
				if( (*mpSp)[sp]->CheckElecCrs() && (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
				{
					double cel= (*mpSp)[sp]->mpElecCrs->mElasticCrscm2[ene];
					double cin= (*mpSp)[sp]->mpElecCrs->mCrsCincm2[ene];
					double dens=(*mpSp)[sp]->mTotDensitycm_3[i];
					mCTot(ene,m)+=(cel+cin)*dens;
					assert(!(mCTot(ene,m)<0.));
					assert(cel>=0);
					assert(cin>=0);
					assert(!(dens<0.));
					assert(!((cel+cin)*dens<0.));
		//			mOmegaFunction(ene,i)+=cel*dens;


					if(i!=nalt-1)
					{
						double densi=(*mpSp)[sp]->mTotDensitycm_3[i+1];
						mCTot(ene,mm)+=(cel+cin)*(dens+densi)/2.;
						assert((cel+cin)*(dens+densi)/2.>0);
						assert(densi>=0);
						assert(!(mCTot(ene,mm)<0.));
					//	mSsalb(ene,i)+=cel*(dens+densi)*0.5/mCTot(ene,mm);
					}
					double colden=(*mpSp)[sp]->mColDenscm_2[i];

					if(i!=0)
					{
						tau(ene,i-1)+=(cel+cin)*colden;
					}else
					{
						tau0+=(cel+cin)*colden;
					}
				}

			}


		}
		// Recalculer ici ssalb


	// GG: jan 2011, comments following the writing of the uncertainty paper
	// mSsalb corresponds to the albedo of single scattering defined in
	// equation 3.3a of Lummerzheim phd thesis
	
		for(unsigned i=0;i<nalt-1;++i)
		{
			mSsalb(ene,i)=0.;

			unsigned mm=i*2+1;
			for(unsigned sp=0;sp<nbsp;++sp)
			{
				if( (*mpSp)[sp]->CheckElecCrs() && (*mpSp)[sp]->mpElecCrs->mIsDefinedCrs)
				{
					double cel= (*mpSp)[sp]->mpElecCrs->mElasticCrscm2[ene];
					double dens=(*mpSp)[sp]->mTotDensitycm_3[i];
					double densi=(*mpSp)[sp]->mTotDensitycm_3[i+1];
					mSsalb(ene,i)+=cel*(dens+densi)*0.5/mCTot(ene,mm);
				}
			}
		//	mSsalb(ene,i)/=mCTot(ene,mm);
		//	Log::mL<<mSsalb(ene,i)<<endl;
		//	assert(mSsalb(ene,i)<1.); We force it in the following...
			mSsalb(ene,i)=nmin(mSsalb(ene,i),0.99999);
		}
		// End of the work on tau
	//	Log::mL<<"Work on tau rows (taureau?)!"<<endl;
		ublas::matrix_row< ublas::matrix<double> > ro(tau,ene);
		for(unsigned m=0;m<nalt-1;++m)//ro.size();++m)
		{
			//		mOmegaFunction(ene,m)/=mCTot(ene,2*m);
			//		if(m!=nalt-1)
			//		{
			ro(m)-=tau0;// We suppress the shift!
			ro(m)*=ismgdpa[m];
			if(m==0)
			{
				mUTau(ene,2*m)=0;
				mUTau(ene,2*m+1)=ro(m)/2.;
				mDiffTau(ene,m)=ro(m);
			}else
			{
				mUTau(ene,2*m)=ro(m-1);
				mUTau(ene,2*m+1)=(ro(m)+ro(m-1))/2.;
				mDiffTau(ene,m)=ro(m)-ro(m-1);
#ifdef DEBUG
				//				if(!(mDiffTau(ene,m)>0))
				//						{
				//							Log::mL<<"Difftau<0"<<mDiffTau(ene,m)<<endl;
				//						}
#endif
					assert((mDiffTau(ene,m)>0));
					//assert(!(mDiffTau(ene,m)<0));
		/*		if(!(mDiffTau(ene,m)>0))
				{
					mDiffTau(ene,m)=1E-42;

				}*/
			}
			//		}
		}
		mUTau(ene,2*nalt-2)=ro(nalt-2);
	}

	Log::mD<<"End of compute CTot main loop"<<endl;
	Log::mI<<"Check if CTot > 1E-20"<<endl;

	vector<double>::iterator it;
	for(unsigned i=0;i<nen-1;++i)
	{// Because it is not valid for nen-1
		ublas::matrix_row< ublas::matrix<double> > ro(mCTot,i);
		for(ublas::matrix_row< ublas::matrix<double> >::iterator it=ro.begin();it!=ro.end();++it)
		{
			if(*it<1E-20)
			{
				//Log::SetPriority(Log::WARNING);
				Log::mE<<"Something is wrong with ctot!!!"<<endl;
				Log::mD<<ro<<endl;
				Error err("CTot","ComputeCTot","The Ctot value is < 1E-20, at position " + ntostr(i) + " Altitude: " + ntostr((*mpAltGridKm)[i])  +"km. \n This is not possible (except in really good vacuum, which is out of the capabilities of this code");
				throw err;
				break;
			}
		}
	}



}



ublas::vector<double> ElectronImpactIonization::ReturnInverseSinMagneticDipAngle()
{
	unsigned nalt=mNbalt;
	if(mpParameter->Exists("/aero_main/electron/magnetic_dip"))
	{
		mpParameter->ExistsOrDie("/aero_main/electron/magnetic_dip/alt","You want a magnetic dip, you must define the alt grid");
		mpParameter->ExistsOrDie("/aero_main/electron/magnetic_dip/inverse_sinus","You want a magnetic dip, you must define the inverse sinus");

		ublas::vector<double> altgrid,inverse_sinus;
		mpParameter->Get1DArray("/aero_main/electron/magnetic_dip/alt",altgrid);
		mpParameter->Get1DArray("/aero_main/electron/magnetic_dip/inverse_sinus",inverse_sinus);

		return MathFunction::IntLin(altgrid,inverse_sinus,*mpAltGridKm);


	}

	// Vertical case
	ublas::vector<double> resu(nalt);
	for(unsigned i=0;i<nalt;++i)
	{
		resu[i]=1.;// the precipitation is vertical : inverse sinuse of 90 = 1. (90°: vertical 0°: horizontal)
	}
	
	return resu;



}


unsigned ElectronImpactIonization::PosEner(double vEnergyeV,unsigned vPosmin)
{
	unsigned nen=(mpElecFlux->mpElecCentEeV)->size();
	double enemin=(*(mpElecFlux->mpElecCentEeV))[nen-1];
	if(vEnergyeV<enemin)
	{
		return nen-1;
	}
	for(unsigned n=vPosmin;n<nen;++n)
	{
		if(vEnergyeV>mCenterGrideV[n])
			return n;
//		double de= (*(mpElecFlux->mpElecDdengeV))[n]/2.;
//		if(vEnergyeV>(*(mpElecFlux->mpElecCentEeV))[n]-de)
//		{// If our energy is greater than the 
//		// bottom energy of the grid -> Ok
//			return n;
//		}
	}

	Error err("ElecCrossSection::PosEner"," Impossible to find the position","It was impossible to find the position of your energy in the position. Please try to put the electron grid in lower values for the minimum of energy");
        throw err;	
	return 0;

}

void ElectronImpactIonization::IntFluxP(unsigned vPos1,unsigned vPos2, ublas::matrix<double> vFlux,unsigned vPos3,double& rSum, double& rSumE)
{

	rSum=0;
	rSumE=0;

	int n1=static_cast<int>(vPos1);
	int n2=static_cast<int>(vPos2);
	//for(int n=n1;n>=n2;--n)
	for(int n=n2;n<=n1;++n)
	{
		rSum+=vFlux(n,vPos3)*(*(mpElecFlux->mpElecDdengeV))[n];
		rSumE+=vFlux(n,vPos3)*(*(mpElecFlux->mpElecCentEeV))[n]
					*(*(mpElecFlux->mpElecDdengeV))[n];
#ifdef ANALYSE_ELEC_FLUX		
		Log::mI<<" Energy : "<<(*(mpElecFlux->mpElecCentEeV))[n]<<"\t Flux: "<< vFlux(n,vPos3) <<"\t tot :"<<rSum<< " " <<rSumE <<endl;
#endif
	}
}

#ifdef ANALYSE_ELEC_FLUX		
void ElectronImpactIonization::IntFluxPa(unsigned vPos1,unsigned vPos2, ublas::matrix<double> vFlux,unsigned vPos3,double& rSum, double& rSumE)
{

	rSum=0;
	rSumE=0;

	int n1=static_cast<int>(vPos1);
	int n2=static_cast<int>(vPos2);
	for(int n=n1;n>=n2;--n)
	{
		rSum+=nabs(vFlux(n,vPos3)*(*(mpElecFlux->mpElecDdengeV))[n]);
		rSumE+=nabs(vFlux(n,vPos3)*(*(mpElecFlux->mpElecCentEeV))[n]
					*(*(mpElecFlux->mpElecDdengeV))[n]);
	//	Log::mI<<" Energy : "<<(*(mpElecFlux->mpElecCentEeV))[n]<<"\t Flux: "<< vFlux(n,vPos3) <<"\t tot :"<<rSum<< " " <<rSumE <<endl;
	}
}
#endif
void ElectronImpactIonization::IntFlux(double vEne1eV,double vEne2eV,ublas::matrix<double> vFluxeV,unsigned vPosAngle,double& rSum, double& rSumE)
{
	int n2=static_cast<int>(PosEner(vEne2eV));
	int n1=static_cast<int>(PosEner(vEne1eV,n2));
	rSum=0;
	rSumE=0;

	for(int n=n1;n>=n2;--n)
	{
		rSum+=vFluxeV(n,vPosAngle)*(*(mpElecFlux->mpElecDdengeV))[n];
		rSumE+=vFluxeV(n,vPosAngle)*(*(mpElecFlux->mpElecCentEeV))[n]
					*(*(mpElecFlux->mpElecDdengeV))[n];
	}



}


void ElectronImpactIonization::PrintPWOM(std::string vFilename)
{
	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# PWOM parameters function of altitude "<<endl;
	of<<"# Altitude in km "<<endl;
	of<<"# Density in cm-3/ flux eV/cm2/s"<<endl;
	of<<"# Ndens1"<<endl;
	of<<"# Ndens2"<<endl;
	of<<"# Nflux1"<<endl;
	of<<"# Nflux2"<<endl;

	unsigned nalt = mNumberFlux2.size();

	for(unsigned i=0;i<nalt;++i)
	{
		of<<(*mpAltGridKm)[i]<<"\t"<<mNumberDensity1[i]<<"\t"<<mNumberDensity2[i]<<"\t"<<mNumberFlux1[i]<<"\t"<<mNumberFlux2[i]<<endl;
	
	}
	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();


}

void ElectronImpactIonization::PrintQintensity(std::string vFilename,ublas::vector< ublas::matrix<double> > qintensity)
{
	unsigned nen=qintensity.size();
	if(nen==0)
		return;

	unsigned nalt2=qintensity(0).size1();
	unsigned nbang=qintensity(0).size2();
	
	Log::mD<<"Bonjour impresion qintensity"<<endl<<"nalt2 : "<<nalt2<<" nbang : "<<nbang<<" nen "<<nen<<endl;
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
	of<<"# Energy in eV"<<endl;
	of<<"# Flux in eV/cm2/s/km"<<endl;


	ublas::vector<double> newalt(nalt2),newene(nalt2),newene2(nalt2);
	newalt.clear();
	newene.clear();
	newene2.clear();
	for(unsigned i=0;i<nalt2;++i)
	{
		double energy=0.;
		double energy2=0;
		for(unsigned j=0;j<nen;++j)
		{
			for(unsigned k=0;k<nbang;++k)
			{
				if(k!=(nbang-1)/2)
				{
					unsigned kprim=k;
					if(k>(nbang-1)/2)
						kprim--;
				energy+=mpElecFlux->ReturnAngleWeight(kprim)*qintensity(j)(i,k)* (*(mpElecFlux->mpElecCentEeV))[j] *  (*mpElecD)[j];
				}else
				{
					energy2+=qintensity(j)(i,k)* (*(mpElecFlux->mpElecCentEeV))[j]  *  (*mpElecD)[j] ;
				}
			}
		}	

		double akm=0.;

		if(i%2)
		{
			akm= ((*mpAltGridKm)[i/2]+(*mpAltGridKm)[i/2+1])/2.;
		}else
		{
			akm=(*mpAltGridKm)[i/2];
		}
		newalt[i]=akm;
		newene[i]=energy;
		newene2[i]=energy2;


		of<<akm<<"\t"<<energy<<"\t"<<energy2<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;

	double enern=MathFunction::TrapzInt(newalt,newene)*1E5;
	double enern2=MathFunction::TrapzInt(newalt,newene2)*1E5;

	Log::mD<<"energy qint : "<<enern<<" et par nango2 : "<<enern2<<endl;
	of<<"energy qint : "<<enern<<" et par nango2 : "<<enern2<<endl;

	of.close();

}



void ElectronImpactIonization::PrintFluxes(std::string vFilename,ublas::vector<double> vAlts,unsigned option)
{
	unsigned nen=mNben;
	ublas::matrix<double> *flux=NULL;
	string precision;
	switch(option)
	{
		case 0:
			flux=&mFluxHemisphericDown;
			precision="down";
			break;
		case 1:
			flux=&mFluxHemisphericUp;
			precision="up";
			break;
		case 2:
			flux=&mFluxHemisphericTot;
			precision="total";
			break;

		case 3:
			flux=&mFluxHemisphericNet;
			precision="net";
			break;
		default:
			Error err("ElectronImpactIonization::PrintFluxes","Wrong option","The option for fluxes is between 0 and 3.");
			throw err;


	}


	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Energy flux "<<precision<<" in function of energy "<<endl;
	of<<"# Energy in eV"<<endl;
	of<<"# Flux in Particle/(eV.cm2.s)"<<endl;



	//	unsigned nen=mFluxHemisphericDown.size1();
	//	unsigned nalt=mFluxHemisphericDown.size2();

	vector<unsigned> positions;
	if(0!=vAlts.size())
	{
		for(unsigned i=0;i<vAlts.size();++i)
		{
			unsigned pos=CloseInVector(vAlts[i],*mpAltGridKm);
			positions.push_back(pos);
			of<<"# Altitude: "<<vAlts[i]<<" km ("<<(*mpAltGridKm)[pos]<<")"<<endl;
		}
	}else
	{
		for(unsigned i=0;i<(mpAltGridKm->size());++i)
		{
			positions.push_back(i);
			of<<"# Altitude: "<<(*mpAltGridKm)[i]<<" km "<<endl;
		}
	}



	if( ! ( flux==NULL ) )
	{
		for(unsigned i=0;i<nen;++i)
		{
			of<< (*mpElecFlux->mpElecCentEeV)[i]<<"\t";

			for(unsigned j=0;j<positions.size();++j)
			{
				of<<(*flux)(i,positions[j])<<"\t";
			}

			of<<endl;
		}
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();


}

void ElectronImpactIonization::PrintBentFluxes(std::string vFilename,const Path & vPath,ublas::vector<double> vAkm, ublas::vector<double> vLkm,unsigned option)
{
	unsigned nen=mNben;
	ublas::matrix<double> *flux=NULL;
	string precision;
	switch(option)
	{
		case 0:
			flux=&mFluxHemisphericDown;
			precision="down";
			break;
		case 1:
			flux=&mFluxHemisphericUp;
			precision="up";
			break;
		case 2:
			flux=&mFluxHemisphericTot;
			precision="total";
			break;

		case 3:
			flux=&mFluxHemisphericNet;
			precision="net";
			break;
		default:
			Error err("ElectronImpactIonization::PrintFluxes","Wrong option","The option for fluxes is between 0 and 3.");
			throw err;


	}


	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Energy flux "<<precision<<" in function of energy "<<endl;
	of<<"# Energy in eV"<<endl;
	of<<"# Flux in eV/cm2/s/sr"<<endl;



	//	unsigned nen=mFluxHemisphericDown.size1();
	//	unsigned nalt=mFluxHemisphericDown.size2();

	vector<unsigned> positions;

	if(0 == vAkm.size() && 0 == vLkm.size())
	{

		for(unsigned i=0;i<vPath.mLengthKm.size();++i)
		{

			positions.push_back(i);
			of<<"# Length: "<<(vPath.mLengthKm)[i]<<" km ; altitude ="<<vPath.mAltitudeKm[i]<<" km"<<endl;
		}
	}else
	{
		for(unsigned i=0;i<vAkm.size();++i)
		{
			unsigned pos=CloseInVector(vAkm[i],vPath.mAltitudeKm);
			positions.push_back(pos);
			of<<"# Altitude: "<<vAkm[i]<<" km ("<<(vPath.mAltitudeKm)[pos]<<"; length ="<<vPath.mLengthKm[pos]<<")"<<endl;
		}
		for(unsigned i=0;i<vLkm.size();++i)
		{
			unsigned pos=CloseInVector(vLkm[i],vPath.mLengthKm);
			positions.push_back(pos);
			of<<"# Length: "<<vLkm[i]<<" km ("<<(vPath.mLengthKm)[pos]<<"; altitude ="<<vPath.mAltitudeKm[pos]<<")"<<endl;
		}
	}


	if( ! ( flux==NULL ) )
	{
		for(unsigned i=0;i<nen;++i)
		{
			of<< (*mpElecFlux->mpElecCentEeV)[i]<<"\t";

			for(unsigned j=0;j<positions.size();++j)
			{
				of<<(*flux)(i,positions[j])<<"\t";
			}

			of<<endl;
		}
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();


}

/*
 *
	mProdIonTot=MathFunction::TrapzInt(*mpAltGridKm,vect_ion_production)*1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
	if(mProdTot>0)
	{
		mElos1=mResultTotalEnergy/mProdTot;
		mElos2=(mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdTot;
		mElos3=mTotalEnergy/mProdTot;
		mIlos1=mResultTotalEnergy/mProdIonTot;
		mIlos2=(mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdIonTot;
		mIlos3=mTotalEnergy/mProdIonTot;
	}
*/

void ElectronImpactIonization::PrintEnergyPerProduction(std::deque<Specie*> vSpecies, std::string vFilename)
{

	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	std::deque<Specie*>::iterator sp;
	of<<"<energyperproduction>"<<endl;
	of<<"<!-- mElos1 = mResultTotalEnergy/mProdTot "<<endl;
	of<<"# mElos2 = (mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdTot "<<endl;
	of<<"# mElos3 = mTotalEnergy/mProdTot  -->"<<endl;
	for(sp=vSpecies.begin();sp!=vSpecies.end();++sp)
	{
		of<<"<"<<(*sp)->mName<<">"<<endl;
		for(unsigned j=0;j<(*sp)->mStates.size();++j)
		{

			if( ((*sp)->mStateProductioncm_3s_1[j]).size() == mpAltGridKm->size())
			{
				double stateprod = MathFunction::TrapzInt(*mpAltGridKm, (*sp)->mStateProductioncm_3s_1[j]) * 1E5;// 1E5 : conversion km -> cm, so, we have a result in eV/cm2
				if(stateprod > 0)
				{
					of<<"\t<"<<(*sp)->mName<<"("<<(*sp)->mStates[j]<<")>"<<endl;

					of<<"\t\t<Elos1>"<< mResultTotalEnergy / stateprod<<"</Elos1>"<<endl;
					of<<"\t\t<Elos2>"<< (mTotalEnergy-mReflectedEnergy-mTransmittedEnergy) / stateprod<<"</Elos2>"<<endl;
					of<<"\t\t<Elos3>"<< mTotalEnergy / stateprod<<"</Elos3>"<<endl;


					of<<"\t</"<<(*sp)->mName<<"("<<(*sp)->mStates[j]<<")>"<<endl;
				}
			}
		}
		of<<"</"<<(*sp)->mName<<">"<<endl;
	}
	of<<"</energyperproduction>"<<endl;
}

void  ElectronImpactIonization::PrintElectronHeating(std::string vFilename)
{

	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Electron Heating ("<<MathFunction::TrapzInt(mQe, *mpAltGridKm) * 1E5<<" eV)"<<endl;
	of<<"# Heating eV/cm3"<<endl;
	of<<"# Altitude"<<endl;
	of<<"# Coulombian heating"<<endl;
	for(unsigned i = 0; i< mQe.size(); ++i)
	{
		of<<(*mpAltGridKm)[i]<<"\t"<<mQe[i]<<endl;
	}
	of<<Log::msMessageLog<<endl;
	of.close();
}



void  ElectronImpactIonization::PrintEnergyPairs(std::string vFilename)
{

	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Energy Loss, Elos1, Elos2, Elos3, Ilos1, Ilos2, Ilos3"<<endl;
	of<<mError*100<<"\t"<<mElos1<<"\t"<<mElos2<<"\t"<<mElos3<<"\t"<<mIlos1<<"\t"<<mIlos2<<"\t"<<mIlos3<<endl;
	of<<"# Error in %"<<endl;
	of<<"# mElos1 = mResultTotalEnergy/mProdTot "<<endl;
	of<<"# mElos2 = (mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdTot "<<endl;
	of<<"# mElos3 = mTotalEnergy/mProdTot "<<endl<<endl;
	of<<"# mIlos1 = mResultTotalEnergy/mProdIonTot "<<endl;
	of<<"# mIlos2 = (mTotalEnergy-mReflectedEnergy-mTransmittedEnergy)/mProdIonTot "<<endl;
	of<<"# mIlos3 = mTotalEnergy/mProdIonTot "<<endl;
	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();

}


void  ElectronImpactIonization::PrintEnergyConservation(std::string vFilename)
{

	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"<xml>"<<endl;
	of<<"<!-- "<<"============================================="<<" -->"<<endl;
	of<<"<!-- "<<" Total input energy: "<<"--><TotalEnergy>"<<mTotalEnergy<<"</TotalEnergy>"<<endl;
	of<<"<!-- "<<"============================================="<<" -->"<<endl;
	of<<"<!-- "<<" Total Heating:      "<<"--><AbsorbedEnergy>"<<mAbsorbedEnergy<<"</AbsorbedEnergy>"<<endl;
	of<<"<!-- "<<" Total inelastic:    "<<"--><InelasticEnergy>"<<mInelasticEnergy<<"</InelasticEnergy>"<<endl;
	of<<"<!-- "<<" Total reflected:    "<<"--><ReflectedEnergy>"<<mReflectedEnergy<<"</ReflectedEnergy>"<<endl;
	of<<"<!-- "<<" Total transmitted:  "<<"--><TransmittedEnergy>"<<mTransmittedEnergy<<"</TransmittedEnergy>"<<endl;
	of<<"<!-- "<<"=============================================" <<" -->"<<endl;
	of<<"<!-- "<<" Total absorbed energy (inelastic + heating) " <<"--><ResultTotalEnergy>"<<mResultTotalEnergy<<"</ResultTotalEnergy>"<<endl;
	of<<"<!-- "<<" Total computed energy (tot abs   + reflect) " <<"--><TotCener>"<< mResultTotalEnergy+mReflectedEnergy+mTransmittedEnergy<<"</TotCener>"<<endl;
	of<<"<!-- "<<"=============================================" <<" -->"<<endl;
	of<<"<!-- "<<"\t Energy conservation error :  "<<"--> <Error>"<<ntostr(mError*100.)<<"</Error><!-- % (>0 : gain)-->"<<endl;
	of<<"<!-- "<<"\t elos (1,2,3)"<<"--><Elos1>"<<mElos1<<"</Elos1> <Elos2> "<<mElos2<<"</Elos2> <Elos3>"<<mElos3<<"</Elos3>"<<endl;
	of<<"<!-- "<<"\t ilos (1,2,3)"<<"--><Ilos1>"<<mIlos1<<" </Ilos1> <Ilos2>"<<mIlos2<<" </Ilos2> <Ilos3>"<<mIlos3<<"</Ilos3>"<<endl;
	of<<"<!-- "<<"\t Production total: "<<"--><ProdTot>"<<mProdTot<<"</ProdTot>"<<endl;
	of<<"<!-- "<<"\t Production total, ions: "<<"--><ProdIonTot>"<<mProdIonTot<<"</ProdIonTot>"<<endl;
	of<<"<!-- "<<"============================================= -->"<<endl;
	double tot = mSecondarySumeV+mDegradedSumeV+mCoulombianSumeV;
	of<<"<!-- "<<"Energy from the secondaries      : "<<"--> <SecondarySumeV>"<<mSecondarySumeV<<" </SecondarySumeV> <perc>"<< mSecondarySumeV*100./tot<<"%</perc>"<<endl; 
	of<<"<!-- "<<"Energy from the degraded      : "<<"--> <DegradedSumeV>"<<mDegradedSumeV<<"</DegradedSumeV>  <perc>"<< mDegradedSumeV*100./tot<<"%</perc>"<<endl;  
	of<<"<!-- "<<"Energy from the coulombian      : "<<"--> <CoulombianSumeV>"<<mCoulombianSumeV<<" </CoulombianSumeV> <perc>"<< mCoulombianSumeV*100./tot<<"%</perc>"<<endl;
	of<<"<!-- "<<Log::msMessageLog<<"-->"<<endl;
	of<<"</xml>"<<endl;
	of.close();
}







