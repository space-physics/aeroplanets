#include "proton.hpp"
using namespace std;
typedef std::deque<Specie* >::iterator SpecIt;
using boost::timer;
using boost::progress_timer;
using boost::progress_display;

#define PGRID (*(mpProtFlux->mpGrideV))
#define SIGE( i , k , h ) ( h == 0 ? ((vSp[ k ])->mpProtCrs)->mTotExcitationCrscm2[ i ] : ((vSp[ k ])->mpHCrs)->mTotExcitationCrscm2[ i ] )
#define WEX( k ) ((vSp[ k ]->mpProtCrs)->mWexcitation)
#define PHASE( k , j , jj , st ) ( vSp[ k ]->mProtonPhase[st]( j , jj ) )
#define SIGI( i , k , h ) ( h == 0 ? ((vSp[ k ])->mpProtCrs)->mTotIonizationCrscm2[ i ] : ((vSp[ k ])->mpHCrs)->mTotIonizationCrscm2[ i ] )
#define WIO( i , k ) (((vSp[ k ])->mpProtCrs)->mElossIoni[ i ])
#define SIGEL( i , k , h ) ( h == 0 ? ((vSp[ k ])->mpProtCrs)->mElasticCrscm2[ i ] : ((vSp[ k ])->mpHCrs)->mElasticCrscm2[ i ] )
//#define WELAS( i , j , jj ) ( ( (vSp[ k ]->mpProtCrs)->mElossElas[ i ]( j , jj ) ) > 1 ? ( (vSp[ k ]->mpProtCrs)->mElossElas[ i ]( j , jj ) ) : 0)
#define WELAS( i , j , jj )  ( (vSp[ k ]->mpProtCrs)->mElossElas[ i ]( j , jj ) )
#define SIGT( i , k , h ) ( h == 0 ? ((vSp[ k ])->mpProtCrs)->mTotalCrscm2[ i ] : ((vSp[ k ])->mpHCrs)->mTotalCrscm2[ i ] )
// Nb this is really the opposite of mXmu (I worked days on that!)
#define MU( j ) ( - mpGa->mXmu[ j ] )
#define WEIGHT( j ) ( mpGa->mWeight[ j ] )
#define ENE( i ) ( PGRID[ i ] )
#define PHASERED( sp, k , j , jj , i , st ) (PhaseR( sp , k , j , jj , i , st ))
#define SIGC( i , k , h ) ( h == 0 ? ((vSp[ k ])->mpProtCrs)->mTotExchangeCrscm2[ i ] : ((vSp[ k ])->mpHCrs)->mTotExchangeCrscm2[ i ] )
#define W10( k )  ((vSp[ k ])->mpProtCrs->mW10 )
#define W01( k , i ) (((vSp[ k ])->mpHCrs)->mElossHioni[ i ])

#define DENE2( i )  ( (*(mpProtFlux->mpWidthGrideV))[ i ] )
#define DENE( i )  ( mProtonWidthGrideV[ i ] )


ProtonHydrogenTransport::ProtonHydrogenTransport(XmlParameters* pParam, ublas::vector<double>* vAltGridKm, boost::shared_ptr<MathFunction::GaussianAngle> vpGa, ublas::vector<double>* vpdB_B):mpParameter(pParam), mpAltGridKm(vAltGridKm), mpdB_B(vpdB_B), mpGa(vpGa)

{
	mbIsRedis = false; // Will be overloaded in the InitPhaseFunction; set here as a safeguard
	// Initialization of the Gaussian angle
	mNbAngle = mpGa -> mNbAngles;
	
	mError1 = 0;
	mError2 = 0;
	mReflectedEnergy = 0;
	mAbsorbedEnergy = 0;
	mInelasticEnergy = 0;
	
	// initialization of the mirror symmetry type
	mMirrorType = 0; // will be updated in InitPhaseFunction()
	mBeamSpreading = 1E-3; // The epsilon in the old code; or the beam spreading parameter (see Galand 1998 or Cyril's Phd Thesis p 132)
	RetrieveInverseSinus(); // We compute the inverse of the sinus of the geometry (non vertical magnetic field lines)
	//mpParameter->ExistsOrDie("/aero_main/proton/nb_angles","You have to define the number of angle for the multistream computation. Even if you do not compute this....");
	//mpParameter->GetValue("/aero_main/proton/nb_angles",mNbAngle);
	//boost::shared_ptr<MathFunction::GaussianAngle> tmpga(new MathFunction::GaussianAngle(mNbAngle));
	//mpGa = tmpga;
	// Initialization of the phase function
	mEinput = 0;
	InitPhaseFunction();

	// Initialization of the tolerance
	mInitialTolerance = 1.E-20;
	mLoopMax = 6; // nb of loops to reach equilibrium. 6 is a good start 
	mCoeffTol = 8.; // The tolerance coefficient
	mbIsReflected = true;
	mpParameter->GetValueOrDefault("/aero_main/proton/init_tol",mInitialTolerance);
	mpParameter->GetValueOrDefault("/aero_main/proton/coeff_tol",mCoeffTol);
	mpParameter->GetValueOrDefault("/aero_main/proton/loop_max",mLoopMax);
	
	if(mpParameter->Exists("/aero_main/proton/no_bottom_reflection"))
	{
		mbIsReflected = false;
	}




}



void ProtonHydrogenTransport::RetrieveInverseSinus()
{


	if(mpParameter->Exists("/aero_main/proton/magnetic_dip"))
	{
		mpParameter->ExistsOrDie("/aero_main/proton/magnetic_dip/alt","You want a magnetic dip, you must define the alt grid");
		mpParameter->ExistsOrDie("/aero_main/proton/magnetic_dip/inverse_sinus","You want a magnetic dip, you must define the inverse sinus");

		ublas::vector<double> altgrid,inverse_sinus;
		mpParameter->Get1DArray("/aero_main/proton/magnetic_dip/alt",altgrid);
		mpParameter->Get1DArray("/aero_main/proton/magnetic_dip/inverse_sinus",inverse_sinus);

		mInvSin = MathFunction::IntLin(altgrid,inverse_sinus,*mpAltGridKm);
		return;

	}
	unsigned nalt = mpAltGridKm->size();
	if( mpParameter->Exists("/aero_main/planet/dB_B/dip_angle"))
	{
		double dipa = 90;
		mpParameter->GetValue("/aero_main/planet/dB_B/dip_angle",dipa);
		dipa *= PI / 180;
		mInvSin.resize(nalt);
		for(unsigned i=0;i<nalt;++i)
		{
			mInvSin[i]=1./ sin(dipa);// the precipitation is vertical : inverse sinuse of 90 = 1. (90°: vertical 0°: horizontal)
		}
		return;
	}


	// Vertical case
	mInvSin.resize(nalt);
	for(unsigned i=0;i<nalt;++i)
	{
		mInvSin[i]=1.;// the precipitation is vertical : inverse sinuse of 90 = 1. (90°: vertical 0°: horizontal)
	}

}

void ProtonHydrogenTransport::InitPhaseFunction()
{

	if(mpParameter->Exists("/aero_main/proton/use_redistribution"))
	{
		mbIsRedis = true;
		Log::mL<<"We use the redistribution "<<endl;
	}


	mpParameter->GetValueOrDefault("/aero_main/proton/beam_spreading",mBeamSpreading);
	Log::mI<<" You use a proton beam spreading parameter (epsilon) of "<<mBeamSpreading<<endl;
	mpParameter->GetValueOrDefault("/aero_main/proton/mirror_type",mMirrorType);
	InitMirror(mpdB_B,mMirrorType);
}

void ProtonHydrogenTransport::ComputeProtonImpact(std::deque<Specie*> vSp,
				boost::shared_ptr<PHFlux> vPFlux, boost::shared_ptr<PHFlux> vHFlux,
				std::deque<Specie*>& rResult, EFlux& rResultFlux)
{
	Log::mI<<"We init the proton/hydrogen transport process!"<<endl;
	mpProtFlux=(vPFlux);
	mpHFlux=(vHFlux);
	Log::mI<<"Partie 1"<<endl;

	// We initialize the proton width
	mProtonWidthGrideV.resize(mpProtFlux->mpWidthGrideV->size());
	for(size_t e = 0; e < mProtonWidthGrideV.size() - 1; ++e)
	{
		mProtonWidthGrideV[e] = PGRID(e) - PGRID( (e + 1) );

	}
	mProtonWidthGrideV[mProtonWidthGrideV.size() - 1] = 0;
	Log::mI<<"Partie 2"<<endl;
	/*
	mProtonWidthGrideV[0] = 0;
	for(size_t e = 1; e < mProtonWidthGrideV.size() ; ++e)
	{
		mProtonWidthGrideV[e] = PGRID(e-1) - PGRID( (e) );
	}
	/
	for(size_t e = 0; e < mProtonWidthGrideV.size() ; ++e)
	{
		Log::mL<<mProtonWidthGrideV[e] <<" "<< DENE2(e)<<endl;
	}*/
	


	// Initialization of the redistribution (does not take time when it has already been computed (once))
	for(SpecIt i = vSp.begin(); i != vSp.end() ; ++ i)
	{
		try{

			if((*i)->CheckProtCrs())
			{
				Log::mI<<"Species: "<<(*i)->mName<<endl;
				(*i)->InitProtonImpact(mpAltGridKm->size(), mpGa, mbIsRedis, mBeamSpreading);
				Log::mI<<"Species ComputeEloss"<<endl;
				((*i)->mpProtCrs)->ComputeElasticLoss((*i)->mMass, *mpGa, mbIsRedis);
			}else
			{
				Log::mW<<"No proton/hydrogen cross section for "<<(*i)->mName << endl;
			}
		}
		catch(...)
		{
			Error err("Proton impact ionization", "Proton::init","There was an error in the initialization of the species " + (*i)->mName + "for the proton ionization");
			throw err;
		}
		//Log::mI<<"Species End"<<endl;
		//	Log::mI<<"PHASE "<<(*i)->mProtonPhase<<endl;
	}
	Log::mI<<"Partie 3"<<endl;

	mElecProduction.resize(mpAltGridKm->size());
	mIonProduction.resize(mpAltGridKm->size());
	mElecProduction.clear();
	mIonProduction.clear();
	Log::mI<<"Matrices initialization"<<endl;
	InitMatrices(vSp);
	// We compute the input energy:
	mEinput = mpProtFlux->FluxEnergyDown() + mpHFlux->FluxEnergyDown() + mpProtFlux->FluxEnergyUp() + mpHFlux->FluxEnergyUp();
	//Log::mL<<mEinput<<" "<< mpProtFlux->FluxEnergyUp() + mpHFlux->FluxEnergyUp()<<endl;


	Log::mI<<"Transport"<<endl;
	Transport();
	Log::mI<<"Ionization of the species"<<endl;
	IonizeSpecies(vSp, rResultFlux);

	Log::mI<<"Finish the Proton/hydrogen transport process!"<<endl;
	SpecieUtils::SpeciesToResu(vSp, rResult);
	Log::mI<<"Proton transport finished; estimation of the energy conservation"<<endl;
	EnergyConservation(vSp);
	Log::mI<<"Proton process finished"<<endl;

}



void ProtonHydrogenTransport::InitMirror(ublas::vector<double>* vpdB_B, int type)
{
	mMirror.resize(mpAltGridKm->size(),mNbAngle);
	mMirror.clear();
	
	if(type == -1)
	{
		ofstream f("mirror.dat");
		Log::mI<<"Magnetic mirror with dissymmetric derivative"<<endl;
		for(size_t i =0 ; i < mpAltGridKm->size(); ++i)
		{
			double muj = mpGa->mXmu[0];
			mMirror(i,0) = (1 -muj * muj) * (*vpdB_B)[i] * 1E-5 / (muj * (muj - mpGa->mXmu[1]));
			f<<mMirror(i,0)<<"\t"<<i<<"\t"<<0<<"\t"<<muj<<"\t"<< mpGa->mXmu[1]<<"\t"<<(*vpdB_B)[i]<<endl;
			for(int j = 1; j < mNbAngle ; ++j)
			{// Nb The mu should be the opposite of Xmu, but the - signs cancels
				double muj = mpGa->mXmu[j];
				double mujm = mpGa->mXmu[j-1];
				mMirror(i,j) = (1 -muj * muj) * (*vpdB_B)[i] * 1E-5 / (muj * (muj - mujm));
				f<<mMirror(i,j)<<"\t"<<i<<"\t"<<j<<"\t"<<muj<<"\t"<< mujm<<"\t"<<(*vpdB_B)[i]<<endl;
			}
		}
		f.close();
	}
}


double ProtonHydrogenTransport::PhaseR(std::deque<Specie*> vSp, unsigned vK, unsigned vJ, unsigned vJJ, unsigned vE, unsigned vSt)
{


	if( ENE(vE) > 1E4 && vSp[ vK ]->mProtonRedistributionFunction[vSt] == 3)
	{

		if(vJ == vJJ and WEIGHT( vJ ) > 1E-30 )
			return 1 / WEIGHT( vJ );
		return 0.;
	}
	return PHASE( vK , vJ , vJJ , vSt ); 

	/*	if( ENE(vE) < 1E4 || vSp[ vK ]->mProtonRedistributionFunction[vSt] != 3)
		{
		return PHASE( vK , vJ , vJJ , vSt ); 
		}
		if(vJ == vJJ and WEIGHT( vJ ) > 1E-30 )
		return 1 / WEIGHT( vJ );
		return 0.;*/
}

void ProtonHydrogenTransport::InitMatrices(std::deque<Specie*> vSp)
{
	// Initialization mLoss, mBdiag, mBinf
	// [species][energy][FinalParticle, InitialParticle][FinalAngle, InitialAngle]
	// (Particle: Proton or Hydrogen)
	ublas::matrix<double> tmp(mNbAngle, mNbAngle);
	tmp.clear();
	ublas::matrix< ublas::matrix<double> > ttmp(2,2);
	ttmp(0,0) = tmp;
	ttmp(0,1) = tmp;
	ttmp(1,0) = tmp;
	ttmp(1,1) = tmp;
	unsigned ne = PGRID.size();
	ublas::vector< ublas::matrix< ublas::matrix<double> > > tttmp(ne);
	for(unsigned i = 0; i < ne ; ++i)
	{
		tttmp[i] = ttmp;
	}
	unsigned nsp = vSp.size();
	mLoss.resize(nsp);
	mBdiag.resize(nsp);
	mBinf.resize(nsp);

	mMatInf.resize(mpAltGridKm->size());
	mMatDiag.resize(mpAltGridKm->size());
	assert(mpAltGridKm->size() > 2 );
	for(unsigned i = 0 ; i < nsp ; ++i)
	{
		mLoss[i] = tttmp;
		mBdiag[i] = tttmp;
		mBinf[i] = tttmp;
	}
	for(unsigned iz = 0 ; iz < mMatInf.size() ; ++iz)
	{
		mMatInf[iz] = tttmp;
		mMatDiag[iz] = tttmp;
	}
#ifdef DEBUG
	unsigned losscall = 0;
	unsigned lossnnul1 = 0;
	unsigned lossnnul2 = 0;
	unsigned lossnnul3 = 0;
#endif

	bool compute_diagmirror = true; // True if the first species
	// In the following, we will have mLoss[k][i][h,h'][j,jj]
	for(unsigned k = 0 ; k < nsp ; ++k)
	{
		if(not vSp[k]->CheckProtCrs())
		{
			Log::mW<<"Skip: "<<vSp[k]->mName<<endl;
			continue;
		}
		for(unsigned i = 0; i < ne ; ++i)
		{
			for(unsigned h = 0; h < 2; ++h)
			{// For that part, crossed terms, h == h'
				for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
				{
					for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
					{
						mLoss[k][i](h,h)(j,jj) = SIGE( i , k , h ) * WEX( k ) * PHASE( k , j , jj , 0) +  SIGI( i , k , h ) * WIO( i , k ) * PHASE( k , j , jj , 0 ) + SIGEL( i , k , h ) * WELAS( i , j , jj ) * PHASERED( vSp, k , j , jj , i , 1  );
#ifdef DEBUG
						losscall++;
						if(nabs(mLoss[k][i](h,h)(j,jj)) > 1E-30)
							lossnnul1++;
#endif
						if ( i == 0 and jj == j )
						{
							mBdiag[k][i](h,h)(j,jj) = - SIGT( i , k , h ) / MU( j );

						}else if(i > 0)
						{

							mBdiag[k][i](h,h)(j,jj) = WEIGHT( jj ) / MU( j ) * ( SIGE( i , k , h ) * PHASE( k , j , jj , 0) +  SIGI( i , k , h ) * PHASE( k , j , jj , 0 ) + SIGEL( i , k , h ) * PHASERED( vSp,  k , j , jj , i , 1 ) );

							double tmp1 = 1 / ( ENE( i ) * ( log(ENE( i - 1 )) - log(ENE( i )) ));
							double tmp2 = 1 / ( ( ENE( i - 1 ) - ENE( i ) ));
							if(mLoss[k][i](h,h)(j,jj) < 1E-37 or mLoss[k][i-1](h,h)(j,jj) < 1E-37)
							{
								mBdiag[k][i](h,h)(j,jj) += WEIGHT( jj ) / MU( j ) * ( (mLoss[k][i-1](h,h)(j,jj) - mLoss[k][i](h,h)(j,jj)) * tmp1 -  mLoss[k][i](h,h)(j,jj) * tmp2);
							}else
							{
								mBdiag[k][i](h,h)(j,jj) += WEIGHT( jj ) / MU( j ) *  mLoss[k][i](h,h)(j,jj) * ( ( log(mLoss[k][i-1](h,h)(j,jj)) - log(mLoss[k][i](h,h)(j,jj)) ) * tmp1 - tmp2);
							}

							mBinf[k][i](h,h)(j,jj) = WEIGHT( jj ) / MU( j ) * mLoss[k][i](h,h)(j,jj) / (ENE( i - 1) - ENE( i ) );
							if ((mBinf[k][i](h,h)(j,jj)) > 1E30)
							{
							//	Log::mI<<"Couille dans la partie 205"<<endl;
								Log::mI<< mLoss[k][i](h,h)(j,jj) <<"  "<<(ENE( i - 1) - ENE( i ) )<<endl;
								Log::mI<<SIGE( i , k , h ) << "  "<< WEX( k ) << " " << PHASE( k , j , jj , 0) <<" " <<  SIGI( i , k , h ) <<" "<< WIO( i , k ) <<" "<<PHASE( k , j , jj , 0 ) << " "<< SIGEL( i , k , h ) <<" "<< WELAS( i , j , jj ) <<" "<< PHASERED( vSp, k , j , jj , i , 1  )<<endl;
							}

							if( jj == j )
							{ // really inside the other, hence the -=
								mBdiag[k][i](h,h)(j,jj) -=  SIGT( i , k , h ) / MU( j );
							}
						} // I suppose i == 0 amd j != jj does not have any significance
					}
				}

			}

			if( i > 0 )
			{
				// Here, terms with h' != h
				for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
				{
					for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
					{
						//  H - P
						mLoss[k][i](0,1)(j,jj) =  SIGC( i , k , 1 ) * W01( k , i) * PHASERED( vSp, k , j , jj , i , 2 );
#ifdef DEBUG
						losscall++;
						if(nabs(mLoss[k][i](0,1)(j,jj)) > 1E-30)
							lossnnul2++;
#endif 

						mBdiag[k][i](0,1)(j,jj) = WEIGHT( jj ) / MU( j ) * ( SIGC( i , k , 1 ) * PHASERED( vSp, k , j , jj , i , 2 ));

						double tmp1 = 1 / ( ENE( i ) * ( log(ENE( i - 1 )) - log(ENE( i )) ));
						double tmp2 = 1 / ( ( ENE( i - 1 ) - ENE( i ) ));
						if(mLoss[k][i](0,1)(j,jj) < 1E-37 or mLoss[k][i-1](0,1)(j,jj) < 1E-37)
						{
							mBdiag[k][i](0,1)(j,jj) += WEIGHT( jj ) / MU( j ) * ( (mLoss[k][i-1](0,1)(j,jj) - mLoss[k][i](0,1)(j,jj)) * tmp1 -  mLoss[k][i](0,1)(j,jj) * tmp2);
						}else
						{
							mBdiag[k][i](0,1)(j,jj) += WEIGHT( jj ) / MU( j ) *  mLoss[k][i](0,1)(j,jj) * ( ( log(mLoss[k][i-1](0,1)(j,jj)) - log(mLoss[k][i](0,1)(j,jj)) ) * tmp1 - tmp2);
						}
						//	mBinf[k][i](0,1)(j,jj) = WEIGHT( jj ) / MU( j ) * mLoss[k][i](0,1)(j,jj) / (ENE( i - 1) - ENE( i ) );
						mBinf[k][i](0,1)(j,jj) = WEIGHT( jj ) / MU( j ) * mLoss[k][i](0,1)(j,jj) * tmp2 ;

						if ((mBinf[k][i](0,1)(j,jj)) > 1E30)
						{
							//Log::mI<<"Couille dans la partie 240"<<endl;
							Log::mI<< mLoss[k][i](0,1)(j,jj) <<"  "<<(ENE( i - 1) - ENE( i ) )<<endl;

							Log::mI<<mLoss[k][i](0,1)(j,jj)  <<"  "<<  SIGC( i , k , 1 )  <<"  "<< W01( k , i)  <<"  "<< PHASERED( vSp, k , j , jj , i , 2 )<<endl;
						}
						//  P - H

						mLoss[k][i](1,0)(j,jj) =  SIGC( i , k , 0 ) * W10( k ) * PHASERED( vSp, k , j , jj , i , 2 );
#ifdef DEBUG
						losscall++;
						if(nabs(mLoss[k][i](1,0)(j,jj)) > 1E-30)
							lossnnul3++;
#endif 

						mBdiag[k][i](1,0)(j,jj) = WEIGHT( jj ) / MU( j ) * ( SIGC( i , k , 0 ) * PHASERED( vSp, k , j , jj , i , 2 ));

						//	double tmp1 = 1 / ( ENE( i ) * ( log(ENE( i - 1 )) - log(ENE( i )) ));
						//	double tmp2 = 1 / ( ( ENE( i - 1 ) - ENE( i ) ));
						//	tmp1 and tmp2 are the same as in H-P
						if(mLoss[k][i](1,0)(j,jj) < 1E-37 or mLoss[k][i-1](1,0)(j,jj) < 1E-37)
						{
							mBdiag[k][i](1,0)(j,jj) += WEIGHT( jj ) / MU( j ) * ( (mLoss[k][i-1](1,0)(j,jj) - mLoss[k][i](1,0)(j,jj)) * tmp1 -  mLoss[k][i](1,0)(j,jj) * tmp2);
						}else
						{
							mBdiag[k][i](1,0)(j,jj) += WEIGHT( jj ) / MU( j ) *  mLoss[k][i](1,0)(j,jj) * ( ( log(mLoss[k][i-1](1,0)(j,jj)) - log(mLoss[k][i](1,0)(j,jj)) ) * tmp1 - tmp2);
						}
						mBinf[k][i](1,0)(j,jj) = WEIGHT( jj ) / MU( j ) * mLoss[k][i](1,0)(j,jj) / (ENE( i - 1) - ENE( i ) );

						if ((mBinf[k][i](1,0)(j,jj)) > 1E30)
						{
							//Log::mI<<"Couille dans la partie 260"<<endl;
							Log::mI<< mLoss[k][i](1,0)(j,jj) <<"  "<<(ENE( i - 1) - ENE( i ) )<<endl;

							Log::mI<< mLoss[k][i](1,0)(j,jj)  <<"  "<< SIGC( i , k , 0 )  <<"  "<<W10( k ) <<"  "<< PHASERED( vSp, k , j , jj , i , 2 )<<endl;
						}
					}
				}
			}
		}

		for(unsigned iz=0; iz < mpAltGridKm->size() ; ++iz)
		{
			for(unsigned i = 0; i < ne ; ++i)
			{
				for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
				{
					// Here, we find the flux up or down, and we compute the multiplicative factor

					double deltalt = 1E5; // unit conversion km -> cm
					if(iz == 0)
					{
						if( j < static_cast<unsigned>(mNbAngle/2))
						{
							deltalt *= ( (*mpAltGridKm)[1] * mInvSin[1] - (*mpAltGridKm)[0] * mInvSin[0] );
						}else
						{
							deltalt *= ( (*mpAltGridKm)[0] * mInvSin[0] - (*mpAltGridKm)[1] * mInvSin[1] );
						}
					}else if( iz == mpAltGridKm->size()-1)
					{
						if( j < static_cast<unsigned>(mNbAngle/2))
						{
							deltalt *= ( (*mpAltGridKm)[iz] * mInvSin[iz] - (*mpAltGridKm)[iz - 1] * mInvSin[iz - 1] );
						}else
						{
							deltalt *= ( (*mpAltGridKm)[iz - 1] * mInvSin[iz - 1] - (*mpAltGridKm)[iz] * mInvSin[iz] );
						}
					}else
					{
						if( j < static_cast<unsigned>(mNbAngle/2))
						{
							deltalt *= ( (*mpAltGridKm)[iz] * mInvSin[iz] - (*mpAltGridKm)[iz - 1] * mInvSin[iz - 1] );
						}else
						{
							deltalt *= ( (*mpAltGridKm)[iz] * mInvSin[iz] - (*mpAltGridKm)[iz + 1] * mInvSin[iz + 1] );
						}
					}
					// Now, we can compute the matrix which does not take into account the species.
	//				deltalt = 1. ;

					for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
					{
						for(unsigned h = 0 ; h < 2 ; ++h)
						{ 
							for(unsigned hp = 0 ; hp < 2 ; ++hp)
							{
								mMatInf[iz][i](h,hp)(j,jj) +=  vSp[k]->mTotDensitycm_3[iz] * mBinf[k][i](h,hp)(jj,j) * deltalt;
								mMatDiag[iz][i](h,hp)(j,jj) += vSp[k]->mTotDensitycm_3[iz] * mBdiag[k][i](h,hp)(jj,j) * deltalt;
		//						mMatInf[iz][i](h,hp)(j,jj) +=  mBinf[k][i](h,hp)(jj,j) * deltalt;
		//						mMatDiag[iz][i](h,hp)(j,jj) += mBdiag[k][i](h,hp)(jj,j) * deltalt;
							//	mMatInf[iz][i](h,hp)(j,jj) +=   mBinf[k][i](h,hp)(jj,j)* deltalt / 2.;
							//	mMatDiag[iz][i](h,hp)(j,jj) +=  mBdiag[k][i](h,hp)(jj,j) * deltalt  ;
								/*if(nabs(mBinf[k][i](h,hp)(j,jj)) > 1E-300)
									mMatInf[iz][i](h,hp)(j,jj) +=  mBinf[k][i](h,hp)(j,jj) * deltalt;
								if(nabs(mBdiag[k][i](h,hp)(j,jj)) > 1E-300)
									mMatDiag[iz][i](h,hp)(j,jj) +=  mBdiag[k][i](h,hp)(j,jj) * deltalt;*/
							}
						}
					}
					//if(0 == k) // We want to add the mirror term only once
					if(compute_diagmirror) // We want to add the mirror term only once
					{
						compute_diagmirror = false;
						if( 0 == j)
						{
							mMatDiag[iz][i](0,0)(0,0) += mMirror(iz,j) * deltalt;
							mMatDiag[iz][i](0,0)(1,0) -= mMirror(iz,j) * deltalt;
						}else
						{
							mMatDiag[iz][i](0,0)(j,j) += mMirror(iz,j) * deltalt;
							mMatDiag[iz][i](0,0)(j-1,j) -= mMirror(iz,j) * deltalt;
						}
					}
				}
			}
		}
	}
/*
	unsigned binfnb = 0, bdiagnb = 0;
	unsigned pbinfnb = 0, pbdiagnb = 0;
	unsigned lossnb = 0, plossnb = 0;
	unsigned total = 0;
		ofstream fli("testloss.dat");
		ofstream flii("testinf.dat");
		ofstream fliii("testdiag.dat"); 
	for(unsigned k = 0 ; k < nsp ; ++k)
	{
		for(unsigned i = 0; i < ne ; ++i)
		{
			for(unsigned h = 0; h < 2; ++h)
			{
				for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
				{
					for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
					{
						for(unsigned hh = 0; hh < 2; ++hh)
						{
							total++;

									fli<<mLoss[k][i](h,hh)(j,jj)<<"\t"<<k<<"\t"<<i<<"\t"<<h<<"\t"<<hh<<"\t"<<j<<"\t"<<jj<<"\t"<< SIGE( i , k , hh ) <<"\t"<<  WEX( k ) <<"\t"<<  PHASE( k , j , jj , 0) <<"\t"<<SIGI( i , k , hh )  <<"\t"<<  WIO( i , k )  <<"\t"<<  PHASE( k , j , jj , 0 )  <<"\t"<< SIGEL( i , k , hh )  <<"\t"<<  WELAS( i , j , jj ) <<"\t"<< PHASERED( vSp, k , j , jj , i , 1  )<<endl;
									flii <<mBinf[k][i](h,hh)(j,jj)<<endl;//"\t"
									if(i > 0)
									fliii <<mBdiag[k][i](h,hh)(j,jj)<<"\t" << mLoss[k][i](h,hh)(j,jj) <<"\t  "<<(ENE( i - 1) - ENE( i ) )<<"\t"<< SIGC( i , k , hh )  <<"\t  "<<W10( k ) <<"\t  "<< PHASERED( vSp, k , j , jj , i , 2 )<<"\t"<<WEIGHT( jj )<<"\t"<<MU( j )<<endl;
									else
									fliii <<mBdiag[k][i](h,hh)(j,jj)<<"\t" << mLoss[k][i](h,hh)(j,jj) <<"\t  "<<(ENE( i ) )<<"\t"<< SIGC( i , k , hh )  <<"\t  "<<W10( k ) <<"\t  "<< PHASERED( vSp, k , j , jj , i , 2 )<<"\t"<<WEIGHT( jj )<<"\t"<<MU( j )<<endl;
									

							if(nabs(mLoss[k][i](h,hh)(j,jj)) > 1E-30)
								lossnb++;
							if((mLoss[k][i](h,hh)(j,jj)) > 1E-30)
								plossnb++;
							if(nabs(mBinf[k][i](h,hh)(j,jj)) > 1E-30)
								binfnb++;
							if(nabs(mBdiag[k][i](h,hh)(j,jj)) > 1E-30)
								bdiagnb++;
							if((mBinf[k][i](h,hh)(j,jj)) > 1E-30)
								pbinfnb++;
							if((mBdiag[k][i](h,hh)(j,jj)) > 1E-30)
								pbdiagnb++;
						}
					}
				}
			}
		}
	}	
	ofstream flii("testmatinf.dat");
	ofstream fliii("testmatdiag.dat"); 
	*/
	for(unsigned iz=0; iz < mpAltGridKm->size() ; ++iz)
	{
		unsigned nn = 0, nn1 = 0, nn2 = 0, nn3 = 0, nn4 = 0, nnv2 = 0;
		for(unsigned i = 0; i < ne ; ++i)
		{
			for(unsigned h = 0; h < 2; ++h)
			{
				for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
				{
					for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
					{
						for(unsigned hp = 0; hp < 2; ++hp)
						{

							if(nabs(mMatInf[iz][i](h,hp)(j,jj)) > 1E-300)
							{
								nn +=1;
								if(i>0)
								{
									nnv2 +=1;
									if(j > static_cast<unsigned>(mNbAngle/2) - 1)
									{
										if(jj > static_cast<unsigned>(mNbAngle/2) - 1)
										{
											nn2 +=1;
										}else
										{
											nn4 +=1;
										}

									}else
									{
										if(jj > static_cast<unsigned>(mNbAngle/2) - 1)
										{
											nn3 +=1;
										}else
										{
											nn1 +=1;
										}
									}
								}
							}
							if(nabs(mMatDiag[iz][i](h,hp)(j,jj)) > 1E-300)
							{
								nn +=1;
								nnv2 +=1;

								if(j > static_cast<unsigned>(mNbAngle/2) - 1)
								{
									if(jj > static_cast<unsigned>(mNbAngle/2) - 1)
									{
										nn2 +=1;
									}else
									{
										nn4 +=1;
									}

								}else
								{
									if(jj > static_cast<unsigned>(mNbAngle/2) - 1)
									{
										nn3 +=1;
									}else
									{
										nn1 +=1;
									}
								}
							}
							/*
							   flii<<	mMatInf[iz][i](h,hp)(j,jj)<<"\t"<<iz<<"\t"<<mpAltGridKm[iz]<<"\t"<< vSp[0]->mTotDensitycm_3[iz]<<"\t"<< vSp[1]->mTotDensitycm_3[iz]<<"\t"<< vSp[2]->mTotDensitycm_3[iz]<<endl;
							   fliii<<	mMatDiag[iz][i](h,hp)(j,jj)<<"\t"<<iz<<"\t"<<mpAltGridKm[iz]<<"\t"<< vSp[0]->mTotDensitycm_3[iz]<<"\t"<< vSp[1]->mTotDensitycm_3[iz]<<"\t"<< vSp[2]->mTotDensitycm_3[iz]<<endl;*/
						}
					}
				}
			}
		}
		//Log::mL<<"Altitude "<<iz<<" "<<(*mpAltGridKm)[iz]<<" km "<<nn<<" (v2) "<<nnv2<<" "<<nn1<<" "<<nn2<<" "<<nn3<<" "<<nn4<<endl;
	}
	
	/*

	//fli.close();
	flii.close();
	fliii.close();
#ifdef DEBUG
	Log::mL<<"Call to Loss :"<<losscall<<" not null "<< lossnnul1<<" - "<< lossnnul2<<" - "<< lossnnul3<<endl;
#endif 
	Log::mL<<mpGa->mWeight<<endl;*/
	/*
	   Log::mL<<"PHASEREDIST"<<endl;
	   for(unsigned k = 0; k < 3; ++k)
	   {
	   unsigned k = 2;
	   for(unsigned j = 0; j < static_cast<unsigned>(mNbAngle); ++j)
	   {
	   for(unsigned jj = 0; jj < static_cast<unsigned>(mNbAngle); ++jj)
	   {
	   for(unsigned i = 0; i < ne ; ++i)
	   {
	   Log::mL<<"phasered( sp=0, j "<<jj<<" jj "<<j<<" i "<<i<<" k "<<k<<")"<<PHASERED(vSp, 0, jj, j, i, k)<<" ENE "<<ENE(i)<<endl;
	   }
	   }
	   }
	   }
	   */
	//	Log::mL<<"PHASE"<<vSp[ 0 ]->mProtonPhase[0]<<endl;
	//	Log::mL<<"PHASE"<<vSp[ 0 ]->mProtonPhase[1]<<endl;
	//	Log::mL<<"PHASE"<<vSp[ 0 ]->mProtonPhase[2]<<endl;

/*	Log::mL<<"Total number of data: "<<total<<endl;
	Log::mL<<"Number of inferior and diagonal non-0 data "<<binfnb<<" "<<bdiagnb<<" "<<pbinfnb<<" "<<pbdiagnb<<endl;
	Log::mL<<"Number of loss non-0 data "<<lossnb<<" "<<plossnb<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"Binf"<<endl<<mBinf<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"Bdiag"<<endl<<mBdiag<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"Loss"<<endl<<mLoss<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
	Log::mL<<"=========================================================="<<endl;
Log::mL<<"=========================================================="<<endl;*/

	//Log::mL<<"Altitudes"<<endl<<*mpAltGridKm<<endl<<"Angle"<<endl<<mInvSin<<endl;
}

void ProtonHydrogenTransport::Transport()
{

//n	Log::mI<<"begin proton transport"<<endl;
	// The proton transport function basically fills the matrices
	// mpProtFlux->mFluxAcm_2s_1eV_2sr_1
	// mpHFlux->mFluxAcm_2s_1eV_2sr_1
	
	// The fluxes are initialized at the top of the atmosphere
	unsigned nalt = mpAltGridKm->size();
	//unsigned nang = mNbAngle;


	assert(nalt > 1);
	assert( (*mpAltGridKm)[0] > (*mpAltGridKm)[1]);
	// Loop for the downward flux
	unsigned bottom_alt = nalt - 2; // may vary if flux too low when going down
	// the boundary conditions are forced at the end of the vector, i.e. nalt -1

	bool up = true;
	bool down = false;

	// Interface with fortran
	int order = PGRID.size() * mNbAngle;
	int (krylovsize) = 10;
	double (time) = 1.;

	int lwsp = static_cast<int>(PGRID.size() * mNbAngle * (krylovsize + 4) +  5 * (krylovsize + 3) * (krylovsize + 3) + krylovsize + 7);
	int liwsp = krylovsize + 4;
	ublas::vector<double> wsp(lwsp);
	ublas::vector<int> iwsp(liwsp);

	double tol_init = mInitialTolerance;//1.E-20; // interface in the xml (TODO Kill Jar Jar Binks)
	double tol = 1.E-20;



	int iflag = 1;
	int verbose = 0;
	// Little parameters
	double coeff_tol = mCoeffTol; //8.; //  interface in the xml
	unsigned loop_max = mLoopMax; // nb of loops to reach equilibrium. 6 is a good start 
	bool is_bottom_reflected = mbIsReflected; // TODO in the XML
	double val_min = 5E-3;
	double val_min_low = 1E-20;


//	double ecart_max = max(ublas::index_norm_inf(mpProtFlux->mPrecipFluxDowncm_2s_1sr_1), ublas::index_norm_inf(mpHFlux->mPrecipFluxDowncm_2s_1sr_1)) * 3. / 100.;

//	Log::mL<<ublas::index_norm_inf(mpProtFlux->mPrecipFluxDowncm_2s_1sr_1)<< "  "<< ublas::index_norm_inf(mpHFlux->mPrecipFluxDowncm_2s_1sr_1)<<endl;
//
	double pmax = 0;
	double hmax = 0;
//	Log::mL<<"Protonflux"<< mpProtFlux->mPrecipFluxDowncm_2s_1sr_1<<endl;
//	Log::mL<<"Hflux"<< mpHFlux->mPrecipFluxDowncm_2s_1sr_1<<endl;
	for(unsigned ene = 0; ene < mpProtFlux->mPrecipFluxDowncm_2s_1sr_1.size1(); ++ene)
	{
		for(unsigned ang = 0; ang < mpProtFlux->mPrecipFluxDowncm_2s_1sr_1.size2();++ang)
		{
			double t =  mpProtFlux->mPrecipFluxDowncm_2s_1sr_1(ene,ang);
			if(t > pmax)
				pmax = t;
		}
	}
	for(unsigned ene = 0; ene < mpHFlux->mPrecipFluxDowncm_2s_1sr_1.size1(); ++ene)
	{
		for(unsigned ang = 0; ang < mpHFlux->mPrecipFluxDowncm_2s_1sr_1.size2();++ang)
		{
			double t =  mpHFlux->mPrecipFluxDowncm_2s_1sr_1(ene,ang);
			if(t > pmax)
				pmax = t;
		}
	}


	double ecart_max = max(pmax, hmax) * 3. / 100.;


	double ecart = ecart_max * 2;

	Log::mL<<"Initialization of the transport; Maximum difference between each loops:" << ecart_max<<endl;
	//Log::mL<< mpProtFlux->mFluxAcm_2s_1eV_1sr_1[0]<<endl;
	FluxToVecAlt0();
	unsigned  loop = 0;
	//Log::mI<<"Maximum flux : " <<ecart_max<<endl;
	/*
#ifdef DEBUG
	ofstream of("inputs_outputs_flux");
	of.precision(9);
	of.setf(ios::scientific);
	of<<endl;
#endif
*/
	//
	//
//XmlParameters matrix("/home/gg/Travail/palantir/dir.venus/guillaume/cppblas/data/Terre/matrix_all.xml");
	while(loop < loop_max && ecart > ecart_max) 
	{
		++loop;
		ecart = 0;
		cout<<" Loop number "<< loop <<" (maximum: "<<loop_max<<")"<<" order "<<order<<endl;
		// We start a timer here
		timer t0;
		progress_display loop_disp(nalt * 2,cout,"\n Proton Flux Loop ","                   ","                   ");

		ublas::vector<double> oldout_flux(PGRID.size() * mNbAngle);
		oldout_flux.clear();
		for(unsigned i = 1; i < nalt ;  ++i)
		{
			// Down - Down Term
			ublas::vector<double> valdd, valdu, fluxup;
			ublas::vector<unsigned> rowdd, coldd, rowdu, coldu;
			double anorm = BlocMatrix(i, 1, valdd, rowdd, coldd );
//			//TiXmlNode* nod = matrix.GetNode("//v[@opt=" + ntostr(i+1) +"]");
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/dd/flux", valdd);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/dd/col", coldd);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/dd/row", rowdd);
			// Up - Down Term (non Homogeneous: down flux scattered from up flux)
			BlocMatrix(i, 4, valdu, rowdu, coldu);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/du/flux", valdu);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/du/col", coldu);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/du/row", rowdu);

//			Log::mL<<"altitude "<<i<<" valddsize "<<valdd.size()<<"  valdusize "<<valdu.size()<<endl;
			fluxup = FluxToVec(i, up);
			ublas::vector<double> flux_down_from_up = VectMultSparceMat(fluxup, valdu, rowdu, coldu);
			ublas::vector<double> flux_down = FluxToVec(i - 1, down);
			ublas::vector<double> out_flux(PGRID.size() * mNbAngle);
			out_flux.clear();

	/*		if( i == 1 and loop == 1)
			{
				ofstream fil("testblocmat.dat");
				for(size_t w = 0; w < valdd.size(); ++w)
				{
					fil<<valdd[w]<<"\t"<<coldd[w]<<"\t"<<w+1<<endl;
				}
				fil<<"#"<<rowdd<<endl;
				fil<<"#"<<flux_down<<endl;
				fil<<"#"<<oldout_flux<<endl;
				for(unsigned opt = 0; opt < 3 ; ++opt)
				{
					fil<<"#";
					for( unsigned aee = 0; aee < PGRID.size();++aee)
					{
						fil<<mpProtFlux->FluxInt(i-1, aee, opt)<<" ";
					}
					fil<<endl;
				}

				fil.close();
			}


			if( i == 10 and loop == 1)
			{
				ofstream fil("testblocmat10.dat");
				for(size_t w = 0; w < valdd.size(); ++w)
				{
					fil<<valdd[w]<<"\t"<<coldd[w]<<"\t"<<w+1<<endl;
				}
				fil<<"#"<<rowdd<<endl;
				fil<<"#"<<flux_down<<endl;
				fil<<"#"<<oldout_flux<<endl;
				for(unsigned opt = 0; opt < 3 ; ++opt)
				{
					fil<<"#";
					for( unsigned aee = 0; aee < PGRID.size();++aee)
					{
						fil<<mpProtFlux->FluxInt(i-1, aee, opt)<<" ";
					}
					fil<<endl;
				}
				fil.close();
			}


			if( i == 40 and loop == 1)
			{
				ofstream fil("testblocmat40.dat");
				for(size_t w = 0; w < valdd.size(); ++w)
				{
					fil<<valdd[w]<<"\t"<<coldd[w]<<"\t"<<w+1<<endl;
				}
				fil<<"#"<<rowdd<<endl;
				fil<<"#"<<flux_down<<endl;
				fil<<"#"<<oldout_flux<<endl;
				for(unsigned opt = 0; opt < 3 ; ++opt)
				{
					fil<<"#";
					for( unsigned aee = 0; aee < PGRID.size();++aee)
					{
						fil<<mpProtFlux->FluxInt(i-1, aee, opt)<<" ";
					}
					fil<<endl;
				}
				fil.close();
			}




			if( i == 60 and loop == 1)
			{
				ofstream fil("testblocmat60.dat");
				for(size_t w = 0; w < valdd.size(); ++w)
				{
					fil<<valdd[w]<<"\t"<<coldd[w]<<"\t"<<w+1<<endl;
				}
				fil<<"#"<<rowdd<<endl;
				fil<<"#"<<flux_down<<endl;
				fil<<"#"<<oldout_flux<<endl;
				for(unsigned opt = 0; opt < 3 ; ++opt)
				{
					fil<<"#";
					for( unsigned aee = 0; aee < PGRID.size();++aee)
					{
						fil<<mpProtFlux->FluxInt(i-1, aee, opt)<<" ";
					}
					fil<<endl;
				}
				fil.close();
			}*/
			tol = tol_init;
			iflag = 1;
			time = 1.;
	//		Log::mI<<"downward nb "<<i<<" Altitude"<<(*mpAltGridKm)[i]<<endl<<valdd<<endl<<rowdd<<endl<<coldd<<endl;
/*#ifdef DEBUG
			of<<"#input flux_down from loop "<<loop<<" altitude "<<i<<endl;
			of<<flux_down<<endl<<flux_down_from_up<<endl;
#endif*/
			while(1 == iflag)
			{	
				out_flux.clear();

		/*		Log::mI<<"Call dgphiv l500"<<endl<<"The norm of the matrix is "<<anorm<<endl<<" Tol is "<<tol<<endl;
				Log::mI<<"downward nb "<<i<<" Altitude"<<(*mpAltGridKm)[i]<<endl;
				Log::mI<<valdd.size()<<" "<<rowdd.size()<<" "<<coldd.size()<<endl;
				Log::mI<<order<<" "<<flux_down_from_up.size()<<" "<<flux_down.size()<<" "<<out_flux.size()<<endl;
			ublas::vector<double> testgg = VectMultSparceMat(flux_down, valdd, rowdd, coldd);
	//		Log::mL<<"test special: "<<testgg<<endl;
			Log::mL<< "max and min multiplication test"<<*(max_element(testgg.begin(), testgg.end()))<<"  "<<*(min_element(testgg.begin(), testgg.end()))<<endl;
			Log::mL<<"Max valdd: "<<*(max_element(valdd.begin(), valdd.end()))<<endl;
			Log::mL<<"Min valdd: "<<*(min_element(valdd.begin(), valdd.end()))<<endl;
			Log::mL<<"Max flux_down: "<<*(max_element(flux_down.begin(), flux_down.end()))<<endl;
			Log::mL<<"Min flux_down: "<<*(min_element(flux_down.begin(), flux_down.end()))<<endl; */
				wsp.clear();
				iwsp.clear();
				if(valdd.size() == 0)
				{
					Log::mW<<"The bloc matrix is empty!!!"<<endl;
					Log::mL<<valdu<<endl;
					Log::mL<<valdd<<endl;
					Log::mL<<"Fin de valud et valdd"<<endl;
					iflag = 0;
				}else
				{
					dgphiv_(&order, &krylovsize, &time, &flux_down_from_up(0), &flux_down(0), &out_flux(0), &tol, &anorm, &wsp(0), &lwsp, &iwsp(0), &liwsp, &verbose, &iflag, &valdd(0), &rowdd(0), &coldd(0));
				}
				if( 1 ==  iflag)
				{
					if( tol < 1)
					{
						(tol) *= coeff_tol;
					}else{
						Error err("PROTON","PROTON TRANSPORT", "DGPHIV failed to compute the flux: the tolerance is too high! The number of krylov bases should probably be increased; the tolerance or its coefficient may be modified...");
						throw err;
					}
				}

			}
//Log::mI<<"End dgphiv"<<endl<<out_flux<<endl;
//Log::mI<<"End dgphiv"<<endl;a
//			Log::mL<<"Max out_flux: "<<*(max_element(out_flux.begin(), out_flux.end()))<<endl;
//			Log::mL<<"Min out_flux: "<<*(min_element(out_flux.begin(), out_flux.end()))<<endl;
					//Error err("PROTON","PROTON TRANSPORT", "TEST");
					//throw err;
			double t_ecart = -1;
			oldout_flux = out_flux;
			double fmax = out_flux[0];
			ublas::vector<double> flux_down_prev = FluxToVec(i, down);
			assert(flux_down_prev.size() == out_flux.size());
			for(unsigned j = 0; j < out_flux.size(); ++j)
			{
				if( out_flux[j] < val_min_low)
					out_flux[j] = val_min_low;
				if( out_flux[j] > fmax)
					fmax = out_flux[j];

				double diff = nabs(out_flux[j] - flux_down_prev[j]);
				if (diff > t_ecart)
				{
					t_ecart = diff; 
				}
			}
/*#ifdef DEBUG
			of<<"#out _down from loop "<<loop<<" altitude "<<i<<endl;
			of<<out_flux<<endl;
#endif*/
			ecart = max(t_ecart, ecart);
//Log::mI<<"Ecart  down alt "<<i<<" "<<t_ecart<<endl;
//Log::mI<<"Max outflux "<<fmax<<" val_min"<<val_min<<endl;
			if(fmax < val_min)
			{
				bottom_alt = i - 2;
				break;
			}
			VecToFlux(out_flux, i, down);
			++loop_disp;
		}

//Log::mI<<"Init Bottom"<<endl;
		if( is_bottom_reflected )
		{
			ublas::vector<double> tmp_flux(PGRID.size() * mNbAngle );
			tmp_flux.clear();
			ublas::vector<double> dn_flux = FluxToVec(bottom_alt + 1, down);
			for(unsigned e = 0; e < 2 * PGRID.size(); ++e)
			{
				for(unsigned a = 0 ; a < static_cast<unsigned>(mNbAngle) / 2 ; ++a)
				{ // We reverse the angles (mirror reflection)
					tmp_flux[ a + e * mNbAngle / 2 ] = dn_flux[ (e + 1) * mNbAngle / 2 - a - 1];
				}
			}
			VecToFlux(tmp_flux, bottom_alt + 1, up); 
		}else
		{
			ublas::vector<double> tmp_flux(PGRID.size() * mNbAngle);
			tmp_flux.clear();
			VecToFlux(tmp_flux, bottom_alt + 1, up); 
		}
//Log::mI<<"End Init Bottom"<<endl;


		// work on flux up
		double ecart_t = -1;
		for(int iz = bottom_alt; iz > -1; --iz)
		{
			assert(iz < static_cast<int>(nalt));
			unsigned i = static_cast<unsigned>(iz);
			//double anorm = BlocMatrix(i, 2, valdd, rowdd, coldd );
			// Up- Down Term (non Homogeneous)
			//double nonh = BlocMatrix(i, 4, valdu, rowdu, coldu);
			// FluxBoundUp, FluxBoundDown
			//dgphiv_()
			//
			// fluxiz put to > 0
			// Up - Up Term
			ublas::vector<double> valuu, valud, fluxdown;
			ublas::vector<unsigned> rowuu, coluu, rowud, colud;
			double anorm = BlocMatrix(i, 2, valuu, rowuu, coluu );

//	matrix.Get1DArray("//v" + ntostr(i+1) + "/uu/flux", valuu);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/uu/col", coluu);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/uu/row", rowuu);
			// Down - Up Term (non Homogeneous: up flux scattered from down flux)
			BlocMatrix(i, 3, valud, rowud, colud);
	//		Log::mL<<"altitude "<<i<<" valuusize "<<valuu.size()<<"  valdusize "<<valud.size()<<endl;

//	matrix.Get1DArray("//v" + ntostr(i+1) + "/ud/flux", valud);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/ud/col", colud);
//	matrix.Get1DArray("//v" + ntostr(i+1) + "/ud/row", rowud);
			fluxdown = FluxToVec(i, down);
			// Former iflux_C
			//
			ublas::vector<double> flux_up_from_down = VectMultSparceMat(fluxdown, valud, rowud, colud);

			// The up flux from the altitude below
			ublas::vector<double> flux_up = FluxToVec(i + 1, up);
			ublas::vector<double> out_flux(PGRID.size() * mNbAngle);
			out_flux.clear();
			tol = tol_init;
			
			if(*max_element(flux_up.begin(), flux_up.end()) > val_min_low)
			{
				iflag = 1;
/*#ifdef DEBUG
			of<<"#input flux_up from loop "<<loop<<" altitude "<<i<<endl;
			of<<flux_up<<endl<<flux_up_from_down<<endl;
#endif*/
				while(1 == iflag)
				{
					out_flux.clear();
					wsp.clear();
					iwsp.clear();
					if(valuu.size() == 0)
					{
						Log::mW<<"The bloc matrix is empty!!!"<<endl;
						Log::mL<<valud<<endl;
						Log::mL<<valuu<<endl;
						Log::mL<<"Fin de valud et valdd"<<endl;
						iflag = 0;
					}else
					{
						dgphiv_(&order, &krylovsize, &time, &flux_up_from_down(0), &flux_up(0), &out_flux(0), &tol, &anorm, &wsp(0), &lwsp, &iwsp(0), &liwsp, &verbose, &iflag, &valuu(0), &rowuu(0), &coluu(0));
					}
					if( 1 == iflag)
					{
						if( tol < 1)
						{
							(tol) *= coeff_tol;
						}else{
							Error err("PROTON","PROTON TRANSPORT", "DGPHIV failed to compute the flux: the tolerance is too high! The number of krylov bases should probably be increased; the tolerance or its coefficient may be modified...");
							throw err;
						}
					}

				}

				ecart_t = -1;
//		a		Log::mL<<"Get Flux prev"<<i<<endl;
				ublas::vector<double> flux_up_prev = FluxToVec(i, up);
//				Log::mL<<"assert"<<i<<endl;
				assert(flux_up_prev.size() == out_flux.size());
				for(unsigned j = 0; j < out_flux.size(); ++j)
				{
					if( out_flux[j] < val_min_low)
						out_flux[j] = val_min_low;
					double diff = nabs(out_flux[j] - flux_up_prev[j]);
					if (diff > ecart_t)
					{
						ecart_t = diff; 
					}
				}
/*#ifdef DEBUG
			of<<"#out up from loop "<<loop<<" altitude "<<i<<endl;
			of<<out_flux<<endl;
#endif*/
			}


	//		Log::mL<<"Vec to flux in the up part... Forced"<<endl;
	//		for(unsigned test = 0 ; test < out_flux.size(); ++ test)
	//			out_flux[test]= 10;
			VecToFlux(out_flux, i, up);
//aLog::mI<<"Ecart  up alt "<<i<<" "<<ecart_t<<endl;
			++loop_disp;
		}
		ecart = max(ecart, ecart_t);
		Log::mI<<endl<<"End loop "<<loop<<" Maximum difference computed: "<<ecart<<endl;
	}

	Log::mI<<"End Proton transport"<<endl;

//ForceFluxTot();

	mpProtFlux->AnisotropicFluxToAverage();
	mpHFlux->AnisotropicFluxToAverage();

//	mpProtFlux->PrintFluxInt("out_proton_hemi.dat", *mpAltGridKm);
//	mpHFlux->PrintFluxInt("out_H_hemi.dat", *mpAltGridKm);
}

ublas::vector<double> ProtonHydrogenTransport::FluxToVec(unsigned vAlt, bool vbIsUp)
{ // Get the flux for Proton and Hydrogen in the vector form. For the altitude vAlt. vbIsUp: true if flux up
	ublas::vector<double> resu(mNbAngle * PGRID.size());
	resu.clear();

	unsigned amin = 0;
	unsigned amax = static_cast<unsigned>(mNbAngle) / 2;

	if(vbIsUp)
	{
		amin = static_cast<unsigned>(mNbAngle) / 2;
		amax = static_cast<unsigned>(mNbAngle);
	}
	if(vAlt != 0)
	{ // We put the proton part in the first half of the vector
		unsigned shift = PGRID.size() * mNbAngle / 2;
		for(unsigned e = 0 ; e < PGRID.size() ; ++e)
		{
			for(unsigned a = amin; a < amax; ++a)
			{
				resu[ e * mNbAngle / 2 + a - amin] = mpProtFlux->mFluxAcm_2s_1eV_1sr_1[vAlt](e, a);
				// And the H in the second half
				resu[ shift +  e * mNbAngle / 2 + a - amin] = mpHFlux->mFluxAcm_2s_1eV_1sr_1[vAlt](e, a);
			}
		}/*
		for(unsigned e = 0 ; e < PGRID.size() ; ++e)
		{
			for(unsigned a = amin; a < amax; ++a)
			{
			}
		}*/
	}else
	{ // Top of the atmosphere (we read the input precipitation flux; it works for both protons and hydrogens)
		if(mpProtFlux->mPrecipFluxUpcm_2s_1sr_1.size1() == 0)
		{
			Log::mL<<"resize P flux up"<<endl;
			mpProtFlux->mPrecipFluxUpcm_2s_1sr_1.resize(PGRID.size(),static_cast<unsigned>(mNbAngle) / 2);
			mpProtFlux->mPrecipFluxUpcm_2s_1sr_1.clear();
		}
		for(unsigned e = 0 ; e < PGRID.size() ; ++e)
		{
			for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
			{
				if(vbIsUp)
				{
					resu[ e * mNbAngle / 2 + a ] = mpProtFlux->mPrecipFluxUpcm_2s_1sr_1(e, a); //static_cast<unsigned>(mNbAngle) / 2 - a - 1);
				}else
				{
					resu[ e * mNbAngle / 2 + a ] = mpProtFlux->mPrecipFluxDowncm_2s_1sr_1(e, a);
				}
			}
		}
		unsigned shift = PGRID.size() * mNbAngle / 2;
		if(mpHFlux->mPrecipFluxUpcm_2s_1sr_1.size1() == 0)
		{
			Log::mL<<"resize H flux up"<<endl;
			mpHFlux->mPrecipFluxUpcm_2s_1sr_1.resize(PGRID.size(),static_cast<unsigned>(mNbAngle) / 2);
			mpHFlux->mPrecipFluxUpcm_2s_1sr_1.clear();
		}
		for(unsigned e = 0 ; e < PGRID.size() ; ++e)
		{
			for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
			{

				if(vbIsUp)
				{

					resu[shift + e * mNbAngle / 2 + a ] = mpHFlux->mPrecipFluxUpcm_2s_1sr_1(e, a); //static_cast<unsigned>(mNbAngle) / 2 - a - 1);
				}else
				{
					resu[shift + e * mNbAngle / 2 + a ] = mpHFlux->mPrecipFluxDowncm_2s_1sr_1(e, a);
				}
			}
		}
	}
	return resu;
}

void ProtonHydrogenTransport::FluxToVecAlt0()
{
	for(unsigned e = 0 ; e < PGRID.size() ; ++e)
	{
		for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
		{
			mpProtFlux->mFluxAcm_2s_1eV_1sr_1[0](e, a) = mpProtFlux->mPrecipFluxDowncm_2s_1sr_1(e, a);
		}
	}
	for(unsigned e = 0 ; e < PGRID.size() ; ++e)
	{
		for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
		{

			mpHFlux->mFluxAcm_2s_1eV_1sr_1[0](e, a) = mpHFlux->mPrecipFluxDowncm_2s_1sr_1(e, a);
		}
	}
}

void ProtonHydrogenTransport::VecToFlux(const ublas::vector<double>& vInputFlux, unsigned vAlt, bool vbIsUp)
{
	unsigned amin = 0;
	unsigned amax = mNbAngle / 2;

	if(vbIsUp)
	{
		amin = mNbAngle / 2;
		amax = mNbAngle;
	}
	unsigned shift = PGRID.size() * mNbAngle / 2;
	for(unsigned e = 0 ; e < PGRID.size() ; ++e)
	{
		for(unsigned a = amin; a < amax; ++a)
		{
			// We put the proton part which is in the first half of the vector
			mpProtFlux->mFluxAcm_2s_1eV_1sr_1[vAlt](e, a) = vInputFlux[ e * mNbAngle / 2 + a - amin]  ;
			// And the H which is in the second half
			mpHFlux->mFluxAcm_2s_1eV_1sr_1[vAlt](e, a) = vInputFlux[ shift +  e * mNbAngle / 2 + a - amin] ;
		}
	}
	if( 0 == vAlt)
	{
		if(vbIsUp)
		{
			for(unsigned e = 0 ; e < PGRID.size() ; ++e)
			{
				for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
				{
					mpProtFlux->mPrecipFluxUpcm_2s_1sr_1(e, a)= vInputFlux[ (e) * mNbAngle / 2 + a ]  ;
				}
			}
			for(unsigned e = 0 ; e < PGRID.size() ; ++e)
			{
				for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
				{

					mpHFlux->mPrecipFluxUpcm_2s_1sr_1(e, a) = vInputFlux[ shift +  (e) * mNbAngle / 2 + a ] ;
				}
			}
		}else
		{
			for(unsigned e = 0 ; e < PGRID.size() ; ++e)
			{
				for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
				{
					mpProtFlux->mPrecipFluxDowncm_2s_1sr_1(e, a)= vInputFlux[ e * mNbAngle / 2 + a ]  ;
				}
			}
			for(unsigned e = 0 ; e < PGRID.size() ; ++e)
			{
				for(unsigned a = 0; a < static_cast<unsigned>(mNbAngle) / 2; ++a)
				{

					mpHFlux->mPrecipFluxDowncm_2s_1sr_1(e, a) = vInputFlux[ shift +  e * mNbAngle / 2 + a ] ;
				}
			}
		}
	}
}



ublas::vector<double> ProtonHydrogenTransport::VectMultSparceMat(const ublas::vector<double>& vX, const ublas::vector<double> &vVal, const ublas::vector<unsigned> & vrRowP, const ublas::vector<unsigned> & vrCol)
{
	ublas::vector<double> resu(vX.size());
	resu.clear();
//	Log::mL<<" Multiplication matricielle " <<endl;
	assert(vVal.size() == vrCol.size() );
	//Log::mL<<vX.size()<<" "<<vrRowP.size()<<endl;
	if( 0 == vrRowP.size())
		return resu;
	assert(vX.size() + 1 == vrRowP.size());

//	for(unsigned i = 0 ; i < vX.size(); ++i)
	if( 0 == vVal.size())
		return resu;
	for(unsigned i = 0 ; i < vX.size(); ++i)
	{
		for(unsigned j = vrRowP[i] - 1 ; j < vrRowP[ i + 1 ] - 1  ; ++j)
		{ // nb : the indices in vrRowP are Fortran style
		//	Log::mL<<vX.size()<<" "<< i <<" "<< j <<" "<< vVal[ j ] <<" "<< vrCol[j]<<endl;
			resu[ i ] +=  vVal[ j ] * vX[ vrCol[ j ] - 1 ];
			//cout<<" no pb" <<endl;
		}
	}
//	Log::mL<<" Fin multiplication matricielle " <<endl;
	return resu;
}



double ProtonHydrogenTransport::BlocMatrix(unsigned vAltPos,unsigned vItype, ublas::vector<double>& vrVal, ublas::vector<unsigned>& vrRowP, ublas::vector<unsigned>& vrCol)
{
	assert(vItype >=1);
	assert(vItype < 5);
	assert(ENE( 0 ) > ENE( 1 ));
	std::deque<double> values;
	std::deque<unsigned> row, col, rowp; 
	unsigned imin, imax, jmin, jmax;

	switch(vItype)
	{
		case 1: {
				imin = 0; imax = mNbAngle / 2;
				jmin = 0; jmax = mNbAngle / 2;
			}; break;
		case 2: {
				imin = mNbAngle / 2; imax = mNbAngle;
				jmin = mNbAngle / 2; jmax = mNbAngle;
			}; break;
		case 3: {
				imin = 0; imax = mNbAngle / 2;
				jmin = mNbAngle / 2; jmax = mNbAngle;
			}; break;
		case 4: {
				imin = mNbAngle / 2; imax = mNbAngle;
				jmin = 0; jmax = mNbAngle / 2;
			}; break;
		default:
			{
				Error err("PROTON","BLOC MATRIX","BLOC TYPE OUT OF RANGE: CODING ERROR");
				throw err;
			}
	}
	double max = mMatInf[vAltPos][0](0,0)(imin, jmin);
	unsigned k = 0;
//vrRowP[0] = 1;
//
#ifdef DEBUG
	double maxmatinf = mMatInf[vAltPos][0](0,0)(imin, jmin);
	double maxmatdiag = mMatDiag[vAltPos][0](0,0)(imin, jmin);

#endif 	
//
	unsigned frow = 0;
	//rowp.push_back(1);
	for(unsigned pfinal = 0 ; pfinal < 2 ; ++pfinal)
	{ // Final particle
		for(unsigned e = 0 ; e < PGRID.size() ; ++ e)
		{ // Final energy
			for(unsigned j = jmin; j < jmax; ++j)
			{ // Final angle
				for(unsigned pinit = 0 ; pinit < 2 ; ++pinit)
				{ // initial particle
					if( e != 0)
					{ // inferior matrix (energy redistribution, to the energy -1 due to continuous loss approximation)
						for(unsigned i = imin; i < imax; ++i)
						{ // initial angle
							double t = mMatInf[vAltPos][e](pfinal,pinit)(i,j);
							if(nabs(t) > 1E-300)
							{
								if(t > max)
									max = t;
#ifdef DEBUG
								if(t > maxmatinf)
									maxmatinf = t;
#endif
								// Hummm complex here!
								//vrVal[k] = t;
								//vrCol[k] = (e - 1) * mNbAngle / 2 + pinit * PGRID.size() * mNbAngle / 2 + i + 1; // +1 for fortran indices in dgphiv
								//#endif	


								//

								values.push_back(t);
								row.push_back(frow + 1);
								// Be careful: we may have to suppress the thing here
								unsigned icorrect = 0;
								if(imin == static_cast<unsigned>(mNbAngle / 2))
								{
									icorrect = i - mNbAngle / 2;
								}else
								{
									icorrect = i;
								}

								col.push_back( (e - 1) * mNbAngle / 2 + pinit * PGRID.size() * mNbAngle / 2 + icorrect + 1);
								++k;
							}
						}
					}
					// Diagonal matrix
					for(unsigned i = imin; i < imax; ++i)
					{ // initial angle
						double t = mMatDiag[vAltPos][e](pfinal,pinit)(i,j);
						if(nabs(t) > 1E-300)
						{
							if(t > max)
								max = t;
#ifdef DEBUG
							if(t > maxmatdiag)
								maxmatdiag = t;
#endif
							values.push_back(t);
							row.push_back(frow + 1);
							unsigned icorrect = 0;
							if(imin == static_cast<unsigned>(mNbAngle / 2))
							{
								icorrect = i - mNbAngle / 2;
							}else
							{
								icorrect = i;
							}
							col.push_back( e * mNbAngle / 2 + pinit * PGRID.size() * mNbAngle / 2 + icorrect + 1);
							//#endif
							++k;
						}
					}
				}
				frow += 1; // e * mNbAngle + pinit * mNbAngle / 2 + i +1
				if(0 == rowp.size())
				{ // The initialization is not at 0 : it is at the value of the first non-zero data?
				//				rowp.push_back(k + 1);
					rowp.push_back( 1);
				}
				rowp.push_back(k + 1);
			}
		}
	}
	if( 0 == k)
		rowp[0] = 0; // To avoid problems whith empty matrices


#ifdef DEBUG
//	Log::mL<<"Blocmatrix, max: "<<max<<" Max Matrix inf: " << maxmatinf<<"  Max Matrix Diag: "<<maxmatdiag<<endl;
#endif
//	Log::mI<<"Initialization of the different matrices: "<< col.size()<<" "<<values.size()<<endl;
	if(0 == values.size())
		return 0;
	assert(col.size() == row.size());
	assert(values.size() == row.size());
//	Log::mL<<rowp.size()<<" "<< mNbAngle * PGRID.size() +1<<endl;
	vrRowP.resize(rowp.size());
	vrCol.resize(col.size());
	vrVal.resize(values.size()); 
/*	ublas::vector<double> inval(values.size());
	ublas::vector<unsigned> incol(col.size());
	std::copy(col.begin(), col.end(),incol.begin());
	std::copy(values.begin(), values.end(),inval.begin());

	ublas::vector<unsigned> roww(row.size());
	std::copy(row.begin(), row.end(),roww.begin());
	int t1 = (int)frow;
	int t2 = (int)k;
	coocsr_(&(t1), &(t2), &inval(0), (int*)(&incol(0)), (int*)(&roww(0)), (&vrVal(0)), (int*)(&vrCol(0)), (int*)(&vrRowP(0)));
*/
//	Log::mI<<"ROWP"<<endl;
	std::copy(rowp.begin(), rowp.end(),vrRowP.begin());
//	Log::mI<<"COL Initialization of the different matrices"<<endl;
	std::copy(col.begin(), col.end(),vrCol.begin());
//	Log::mI<<"VAL Initialization of the different matrices"<<endl;
	std::copy(values.begin(), values.end(),vrVal.begin());
//	Log::mI<<"ROW Initialization of the different matrices"<<endl;
#ifdef DEBUG
/*	March 1st 2012 : it is equal here */
/*	ublas::vector<unsigned> roww(row.size());
	std::copy(row.begin(), row.end(),roww.begin());
	Log::mL<<"Start comparizon"<<endl;
	Log::mL<<"Number of row :"<<frow<<" number of energy times angle "<<mNbAngle * PGRID.size()<<endl;
	int t1 = (int)frow;
	int t2 = (int)k;
	ublas::vector<double> outval(vrVal.size());
	ublas::vector<unsigned> outcol(vrCol.size());
	ublas::vector<unsigned> outrowp(vrRowP.size());
	//coocsr_(&(t1), &(t2),&vrVal(0), (int*)&vrCol(0), (int*)&roww(0), (&outval(0)), (int*)(&outcol(0)), (int*)(&outrowp(0)));
	coocsr_(&(t1), &(t2),&vrVal(0), (int*)&roww(0), (int*)&vrCol(0), (&outval(0)), (int*)(&outcol(0)), (int*)(&outrowp(0)));
	Log::mL<<"row "<<vrRowP.size()<<" "<<outrowp.size()<<endl;
	for(unsigned z = 0 ; z < frow; ++z)
	{
		Log::mL<< vrRowP[z]<<" "<< vrRowP[z+1]<<"  ====   "<< outrowp(z)<<" "<< outrowp(z+1)<<endl;
		Log::mL<<"\tROWP"<<endl;
		for(unsigned zi = vrRowP[z] -1; zi < vrRowP[z + 1] -1 ; ++zi)
		{
			Log::mL<<"\t\t\t"<<vrVal[zi]<<" "<<vrCol[zi]<<" "<<roww[zi]<<endl;
		}
		Log::mL<<"\tROWP"<<endl;
		for(unsigned zi = outrowp[z] - 1; zi < outrowp[z + 1] -1 ; ++zi)
		{
			Log::mL<<"\t\t\t"<<outval[zi]<<" "<<outcol[zi]<<" "<<z<<endl;
		}
	
	}
	if(MathFunction::VectorCompare(vrRowP, outrowp))
	{
		Log::mL<<"Cool, Rowp and outrowp are equal"<<endl;
	}
	if(MathFunction::VectorCompare(vrVal, outval))
	{
		Log::mL<<"Cool, Val and outval are equal"<<endl;
	}

	if(MathFunction::VectorCompare(vrCol, outcol))
	{
		Log::mL<<"Cool, Col and outcol are equal"<<endl;
	}
	Error err("PROTON","PROTON TRANSPORT", "TEST ROW");
	throw err;
*/
#endif

	return max;
}





void ProtonHydrogenTransport::IonizeSpecies(std::deque<Specie*> vSp,EFlux& rResultFlux)
{

	unsigned nalt = mpAltGridKm->size();
	unsigned nbelece = (rResultFlux.mpElecCentEeV)->size();
	ublas::matrix<double> elecprode(nalt,nbelece);
	elecprode.clear();
//	ublas::vector<double> elecprod(nalt);

	for(unsigned i=0;i<nalt;++i)
	{
	//	Log::mL<<"IOnization, altitude "<<i<<" < "<<nalt<<endl;

		for(unsigned ene=0;ene<PGRID.size();++ene)
		{
			std::deque<Specie*>::iterator it;
		//	unsigned sp=0;
			double energy = ENE( ene );
//DENE(ene) = 1;
			for(it=vSp.begin();it!=vSp.end();++it)
			{
				double dens=(*it)->mTotDensitycm_3[i];
				if( (*it)->CheckProtCrs() && (*it)->mpProtCrs->mIsDefinedCrs && (*it)->mpHCrs->mIsDefinedCrs)
				{
					unsigned protprocnumber  = (*it)->mpProtCrs->mProcessNames.size();
					for(unsigned j = 0; j < protprocnumber ; ++j)
					{

						double threshold=(*it)->mpProtCrs->mThresholdseV[j];

						if(energy > threshold)
						{
							double cinex=(*it)->mpProtCrs->mCrscm2[j][ene];
							double produ= mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) *cinex * dens;
							(*it)->mSpeciesProductioncm_3s_1[j][i] += produ ;
						}
					}
					for(unsigned j = 0 ; j < (*it)->mpHCrs->mProcessNames.size() ; ++j)
					{

						double threshold=(*it)->mpHCrs->mThresholdseV[j];

						if(energy > threshold)
						{
							double cinex=(*it)->mpHCrs->mCrscm2[j][ene];
							double produ= mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) *cinex * dens;
							(*it)->mSpeciesProductioncm_3s_1[j + protprocnumber][i] += produ;

						}
					}
					// The ionization by protons
					for(unsigned j = 0;j < (*it)->mpProtCrs->mIonizationCrsPosition.size(); ++j)
					{
						unsigned position = (*it)->mpProtCrs->mIonizationCrsPosition[j];
						double elnumber = ((*it)->mpProtCrs->mNumberOfElectrons[position]);
						double cinex = (*it)->mpProtCrs->mCrscm2[position][ene];
						mElecProduction[i] += (mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						mIonProduction[i] += (mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * ((*it)->mpProtCrs->mNumberOfIons[position]);
						(*it)->mProtElecProductioncm_3s_1[i] +=  (mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						if((*it)->mpProtCrs->mIsAuger.at(position))
						{// If the auger process is defined
							for(unsigned z=0;z<(*it)->mpProtCrs->mAugerEnergy.at(position).size();++z)
							{
								//double augerenergy=(*it)->mpProtCrs->mAugerEnergy.at(position).at(z);
								(*it)->mProtElecProductioncm_3s_1[i]+= mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex * dens * (*it)->mpProtCrs->mAugerEfficiency.at(position).at(z);
								mElecProduction[i]+= mpProtFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex * dens * (*it)->mpProtCrs->mAugerEfficiency.at(position).at(z);
							}
						}
					}

					// The ionization by H impact
					for(unsigned j = 0; j < (*it)->mpHCrs->mIonizationCrsPosition.size(); ++j)
					{
						unsigned position = (*it)->mpHCrs->mIonizationCrsPosition[j];
						double elnumber = ((*it)->mpHCrs->mNumberOfElectrons[position]);
						//double thresh = (*it)->mpHCrs->mThresholdseV[position];
						double cinex = (*it)->mpHCrs->mCrscm2[position][ene];
						mElecProduction[i] += (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						mIonProduction[i] += (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * ((*it)->mpHCrs->mNumberOfIons[position]);
						(*it)->mHElecProductioncm_3s_1[i] +=  (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						if((*it)->mpHCrs->mIsAuger.at(position))
						{// If the auger process is defined
							for(unsigned z=0;z<(*it)->mpHCrs->mAugerEnergy.at(position).size();++z)
							{
								//double augerenergy=(*it)->mpHCrs->mAugerEnergy.at(position).at(z);
								(*it)->mHElecProductioncm_3s_1[i]+= mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex * dens * (*it)->mpHCrs->mAugerEfficiency.at(position).at(z);
								mElecProduction[i]+= mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex * dens * (*it)->mpHCrs->mAugerEfficiency.at(position).at(z);
							}
						}
					}

					// The ionization of H
					for(unsigned j = 0; j < (*it)->mpHCrs->mExchangeCrsPosition.size(); ++j)
					{
						unsigned position = (*it)->mpHCrs->mExchangeCrsPosition[j];
						double elnumber = ((*it)->mpHCrs->mNumberOfElectrons[position]); // elnumber should be 1, but we let the possibility to tweak that process
						//double thresh = (*it)->mpHCrs->mThresholdseV[position];
						double cinex = (*it)->mpHCrs->mCrscm2[position][ene];
						mElecProduction[i] += (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						mIonProduction[i] += (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * ((*it)->mpHCrs->mNumberOfIons[position]);
						(*it)->mHElecProductioncm_3s_1[i] +=  (mpHFlux->mFluxcm_2s_1eV_1(i, ene)  * DENE( ene) * cinex) * dens * elnumber;
						// There is no Auger/multiple creation for the charge exchange 
					}

					//The electron production
					//
					for(unsigned e = 0 ; e < nbelece ; ++e )
					{
						if( ENE(ene) > (*rResultFlux.mpElecCentEeV)[e])
						{ // to simplify
							double Ea = sqrt( (*it)->mpProtCrs->mVion * ENE( ene ) /1836.  ) / 0.91;
							double f1 = 1. / Ea * exp( - ((*rResultFlux.mpElecCentEeV)[e]) / Ea);
							
							double Bg0 = sqrt( ENE( ene ) / 1836. /  ((*rResultFlux.mpElecCentEeV)[e]));
							double Ag0 = sqrt( ((*rResultFlux.mpElecCentEeV)[e]) / Bg0);
							double erfp = erf(Ag0 * (1 + Bg0));
							double erfm = erf(Ag0 * nabs(1 - Bg0));
							// beta = 1.
							double g0 = sqrt(PI * 1836. / ENE( ene ) ) * (erfp - erfm)/4. ;

							elecprode(i,e) += DENE( ene ) * dens * (
									mpProtFlux->mFluxcm_2s_1eV_1(i, ene) * f1 * (*it)->mpProtCrs->mTotIonizationCrscm2[ene] + 
									mpHFlux->mFluxcm_2s_1eV_1(i, ene) * f1 * (*it)->mpHCrs->mTotIonizationCrscm2[ene] +
									mpHFlux->mFluxcm_2s_1eV_1(i, ene) * g0 * (*it)->mpHCrs->mTotExchangeCrscm2[ene]
									);
						}
					}
				}
				//++sp;
			}
		}
	}// finish the work on the altitude grid
	// And we init the electron flux
	Log::mL<<"We init the electron flux"<<endl;
	//Log::mL<<mpProtFlux->mFluxcm_2s_1eV_1<<endl;
	rResultFlux.InitIsotropic(elecprode, mElecProduction ,vSp);
/*	Log::mL<<"Electron production: "<<mElecProduction<<endl;
	Log::mL<<"Ion production: "<<mIonProduction<<endl;*/
	Log::mL<<"Energy in the electron flux: "<<rResultFlux.FluxEnergy(*mpAltGridKm)<<" eV"<<endl;
}


void ProtonHydrogenTransport::EnergyConservation(std::deque<Specie*> vSp)
{

	// Computation of general heating
	// depending on the different species
	// and on the redistribution function
	//  - > Computes Qe
	//  - > Computes eoutput)min
	//   
	//
	//
	// Heating layer by layer
	unsigned nalt = mpAltGridKm->size();
//	unsigned nbelece = PGRID.size();

	ublas::vector<double> Eoutmin(nalt), Qe(nalt);
	Eoutmin.clear();
	Qe.clear();
	ublas::vector<double> Qee(nalt), Qei(nalt), Qec(nalt), Qeel(nalt), Qep(nalt), Qeh(nalt);
	Qee.clear();
	Qei.clear();
	Qec.clear();
	Qeel.clear();
	Qep.clear();
	Qeh.clear();
	unsigned nang = mNbAngle;
	Log::mL<<"Computation of Qe"<<endl;
	for(unsigned a=0;a<nalt;++a)
	{
		ublas::vector<double> li(PGRID.size());
		li.clear();
		ublas::vector<double> lie(PGRID.size()), lii(PGRID.size()), lic(PGRID.size()), liel(PGRID.size()), lip(PGRID.size()), lih(PGRID.size());
		lie.clear();
		lii.clear();
		lic.clear();
		liel.clear();
		lip.clear();
		lih.clear();
		for(unsigned i=0;i<PGRID.size();++i)
		{
			for(unsigned k = 0; k < vSp.size(); ++k)
			{

				if(not vSp[k]->CheckProtCrs())
				{
					Log::mW<<"Skip: "<<vSp[k]->mName<<endl;
					continue;
				}
				ublas::vector<double> we(nang);
				we.clear();
				for(unsigned j = 0; j < nang; ++j)
				{
					for(unsigned jj = 0; jj < nang; ++jj)
					{
						if( ENE(i) < 1E4 || vSp[ k ]->mProtonRedistributionFunction[1] != 3)
						{
							we(j) += WELAS( i, j, jj ) * PHASE( k , j , jj , 1 ) * WEIGHT( jj ); 
						}else{
							if( jj == j )
							{
								we(j) += WELAS( i, j, jj ); // We multiply and divide by weight here!
							}
						}
					}
					double IP = mpProtFlux->mFluxAcm_2s_1eV_1sr_1[a](i,j);
					double IH = mpHFlux->mFluxAcm_2s_1eV_1sr_1[a](i,j);
					li(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( 
							(SIGE( i , k , 0 ) * WEX( k ) +
							 SIGI( i , k , 0 ) * WIO( i , k ) +
							 SIGC( i , k , 0 ) * W10( k ) +
							 SIGEL( i , k , 0 ) * we(j) //WELAS( i , j , jj )
							)
							* IP+ 
							(SIGE( i , k , 1 ) * WEX( k ) +
							 SIGI( i , k , 1 ) * WIO( i , k ) +
							 SIGC( i , k , 1 ) * W01( k , i) +
							 SIGEL( i , k , 1 ) * we(j)//WELAS( i , j , jj )
							)
							* IH);
					lii(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( (  SIGI( i , k , 0 ) * WIO( i , k )  ) * IP+ ( SIGI( i , k , 1 ) * WIO( i , k ) ) * IH);
					lie(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( (SIGE( i , k , 0 ) * WEX( k ) ) * IP+ (SIGE( i , k , 1 ) * WEX( k )) * IH);
					lic(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( SIGC( i , k , 0 ) * W10( k )  * IP+ (SIGC( i , k , 1 )* W01( k , i)) * IH);
					liel(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( SIGEL( i , k , 0 ) * we(j) * IP+ (SIGEL( i , k , 1 ) * we(j)) * IH);
					lip(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ( (SIGE( i , k , 0 ) * WEX( k ) + SIGI( i , k , 0 ) * WIO( i , k ) + SIGC( i , k , 0 ) * W10( k ) + SIGEL( i , k , 0 ) * we(j) ) * IP);
					lih(i) += vSp[k]->mTotDensitycm_3[a] * WEIGHT( j ) * ((SIGE( i , k , 1 ) * WEX( k ) + SIGI( i , k , 1 ) * WIO( i , k ) + SIGC( i , k , 1 ) * W01( k , i) + SIGEL( i , k , 1 ) * we(j)) * IH);


					if( 0 == k && PGRID.size() - 1 == i)
					{
						Eoutmin[a] += WEIGHT( j ) * nabs(MU( j )) * (IP + IH);
					}
				}
			}
		}
	/*	if( 0 == a || 9 == a)
		{
			Log::mI<<"LI(0)"<<endl<<li<<endl;
		}*/
		for(unsigned i=1;i<PGRID.size();++i)
		{
			Qe[a] += 2 * PI * ( li(i-1) + li(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qei[a] += 2 * PI * ( lii(i-1) + lii(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qee[a] += 2 * PI * ( lie(i-1) + lie(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qec[a] += 2 * PI * ( lic(i-1) + lic(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qeel[a] += 2 * PI * ( liel(i-1) + liel(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qep[a] += 2 * PI * ( lip(i-1) + lip(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
			Qeh[a] += 2 * PI * ( lih(i-1) + lih(i) ) / 2. * ( ENE( i-1 ) - ENE( i ));
		}
	}
	// Energy in the flux : HEATING (eV cm-2 s-1)
	double QeI = 0.;
	// Energy deposited for E < Emin
	double Eoutput_min = 2 * PI * Eoutmin[0] * ENE( PGRID.size() - 1 ) * ( ENE( PGRID.size() - 2 ) - ENE( PGRID.size() - 1 ));


	ublas::vector<double> eta(nalt);
	eta.clear();

	double QeIi =0, QeIe =0, QeIc =0, QeIel =0, QeIp =0, QeIh =0;



	//Log::mL<<"Computation of QeI"<<endl;
	for(unsigned a = 1 ; a < nalt; ++a)
	{
		QeI += (Qe[a - 1] + Qe[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIi += (Qei[a - 1] + Qei[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIe += (Qee[a - 1] + Qee[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIc += (Qec[a - 1] + Qec[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIel += (Qeel[a - 1] + Qeel[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIp += (Qep[a - 1] + Qep[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
		QeIh += (Qeh[a - 1] + Qeh[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;


		Eoutput_min += 2 * PI * Eoutmin[a] * ENE( PGRID.size() - 1 ) * ( ENE( PGRID.size() - 2 ) - ENE( PGRID.size() - 1 ));

		ublas::vector<double> ddi(PGRID.size());
		ddi.clear();
		for(unsigned i=1;i<PGRID.size();++i)
		{

			for(unsigned j = 0; j < nang; ++j)
			{

					double pu = mpProtFlux->mFluxAcm_2s_1eV_1sr_1[a-1](i,j);
					double hu = mpHFlux->mFluxAcm_2s_1eV_1sr_1[a-1](i,j);
					double p = mpProtFlux->mFluxAcm_2s_1eV_1sr_1[a](i,j);
					double h = mpHFlux->mFluxAcm_2s_1eV_1sr_1[a](i,j);
				if(j > nang / 2)
				{
					ddi[i] += abs( MU( j ) * WEIGHT( j ) ) * (p + h - pu - hu) / 1E5 / ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]);
				}else
				{
					ddi[i] += abs( MU( j ) * WEIGHT( j ) ) * (pu + hu - p - h) / 1E5 / ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]);
				}
			}
			ddi[i] *= ENE( i );
		}
		for(unsigned i=1;i<PGRID.size();++i)
		{
			eta[a] += 2 * PI * (ddi[i - 1] + ddi[i]) / 2. * (ENE(i - 1) - ENE( i ));
			if(eta[a] < 1E-37)
				eta[a] = 1E-37;
		}
	}

//	Log::mL<<"Computation of Qei_90"<<endl;
	double Qei_90 = 0;
	for(unsigned a = 1 ; a < nalt; ++a)
	{
		Qei_90 +=  (eta[a - 1] + eta[a]) / 2. * ((*mpAltGridKm)[a - 1] * mInvSin[a - 1] - (*mpAltGridKm)[a] * mInvSin[a]) * 1E5;
	}



	// Energy in the precipitating flux (input)
//	Log::mL<<"Computation of Ein"<<endl;
//	double Ein = mpProtFlux->FluxEnergyDown() + mpHFlux->FluxEnergyDown();
	// Energy in the escaping flux (output)
//	Log::mL<<"Computation of Eout"<<endl;
	double Eou = mpProtFlux->FluxEnergyUp() + mpHFlux->FluxEnergyUp();

//	Log::mL<<"Nouvel input!"<<mEinput<<endl;
	double Ein = mEinput;

//	Log::mL<<"Computation of Ein2"<<endl;
	double Ein2 = mpProtFlux->FluxEnergyDown2() + mpHFlux->FluxEnergyDown2();

	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;


	mError1 = (Ein - Eou - Eoutput_min - QeI) / Ein;
	mError2 = (Ein - Eou - Eoutput_min - Qei_90) / Ein;
	mReflectedEnergy = Eou;
	mInelasticEnergy = QeI;
	mAbsorbedEnergy = Eoutput_min;



	Log::mL<<"Energy conservation (1st method) "<< ntostr(100. * (Ein - Eou - Eoutput_min - QeI) / Ein)<<"%"<<endl;
	Log::mL<<"Energy conservation (2nd method) "<< ntostr(100. * (Ein - Eou - Eoutput_min - Qei_90) / Ein)<<"%"<<endl;

	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"Energy in the precipitated flux: "<<Ein<<endl;
	Log::mL<<"Energy precipitating in the first layer: "<<Ein2<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"Energy in the escaping flux: "<<Eou<<endl;
	Log::mL<<"Energy deposited below Emin: "<<Eoutput_min<<endl;
	Log::mL<<"Energy transformed into Heating, ionization...: "<<QeI<<endl;
	Log::mL<<"Energy transformed into Heating, ionization... 2nd method: "<<Qei_90<<endl;

	Log::mL<<"Energy in Qe "<<QeI<<" : "<<endl;
	Log::mL<<"\t From excitation "<<QeIe<<endl;
	Log::mL<<"\t From ionization "<<QeIi<<endl;
	Log::mL<<"\t From exchange "<<QeIc<<endl;
	Log::mL<<"\t From elastic "<<QeIel<<endl;
	Log::mL<<"\t All of those from protons "<<QeIp<<endl;
	Log::mL<<"\t All of those from hydrogen "<<QeIh<<endl;


	Log::mL<<"==============================================="<<endl;
	Log::mL<< mpProtFlux->FluxEnergyDown() <<" "<< mpHFlux->FluxEnergyDown()<<" "<< mpProtFlux->FluxEnergyUp() <<" "<< mpHFlux->FluxEnergyUp()<<endl;
	Log::mL<< mpProtFlux->FluxEnergyDown2() <<" "<< mpHFlux->FluxEnergyDown2()<<" "<< mpProtFlux->FluxEnergyUp2() <<" "<< mpHFlux->FluxEnergyUp2()<<endl;
/*	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"Flux down alt 0 "<<endl<<FluxToVec(0, false)<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"Flux down alt 9 "<<endl<<FluxToVec(9, false)<<endl;
	Log::mL<<"==============================================="<<endl;
//	Log::mL<<"Flux up alt 0 "<<endl<<FluxToVec(0, true)<<endl;
	Log::mL<<"==============================================="<<endl;
//	Log::mL<<"Flux up alt 0 "<<endl<<FluxToVec(9, true)<<endl;
*/	
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;
	Log::mL<<"==============================================="<<endl;
}


void ProtonHydrogenTransport::PrintFluxes(std::string vFilename,ublas::vector<double> vAlts,unsigned vOption, bool vIsProton)
{
	unsigned nen=PGRID.size();
	//	ublas::matrix<double> *flux=NULL;
	string precision;
	switch(vOption)
	{
		case 0:
			//	flux=&mFluxHemisphericDown;
			precision="down";
			break;
		case 1:
			//	flux=&mFluxHemisphericUp;
			precision="up";
			break;
		case 2:
			//	flux=&mFluxHemisphericTot;
			precision="total";
			break;

		case 3:
			//	flux=&mFluxHemisphericNet;
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

	of<<"# Energy flux "<<precision<<" in function of energy, for the ";
	if(vIsProton)
	{
		of<<"protons"<<endl;
	}else
	{
		of<<"H"<<endl;
	}
	of<<"# Energy in eV"<<endl;
	of<<"# Flux in eV/cm2/s/sr"<<endl;



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



	for(unsigned i=0;i<nen;++i)
	{
		of<< ENE(i)<<"\t";

		for(unsigned j=0;j<positions.size();++j)
		{
			//of<<(*flux)(i,positions[j])<<"\t"
			double flx = 0;
			if(vIsProton)
			{
				flx = mpProtFlux->FluxInt(positions[j],i, vOption);

			}else
			{
				flx = mpHFlux->FluxInt(positions[j],i, vOption);
			}
			of<<flx<<"\t";
		}
		of<<endl;
	}

	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}



void ProtonHydrogenTransport::PrintEnergyFluxes(std::string vFilename,ublas::vector<double> vAlts, bool vIsProton, bool vIsLength)
{
	if(vIsProton)
	{
		mpProtFlux->PrintEnergyFlux(vFilename, vAlts, vIsLength);
	}else
	{
		mpHFlux->PrintEnergyFlux(vFilename, vAlts, vIsLength);
	}

	PrintFluxTot(vFilename + "_FLUXtot.dat");
}


void  ProtonHydrogenTransport::PrintFluxTot(std::string vFilename)
{
	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<endl;

	for(size_t alt = 0; alt < mpAltGridKm->size(); ++alt)
	{
		for(size_t a = 0; a < 2; ++a)
		{

			for(size_t sp = 0; sp < 2; ++sp)
			{

				for(size_t i = 0; i < PGRID.size(); ++i)
				{

					for(size_t ang = 0; ang < static_cast<size_t>(mNbAngle / 2); ++ang)
					{
						size_t shift = a * mNbAngle / 2;
						if(sp == 0)
						{
							of<<mpProtFlux->mFluxAcm_2s_1eV_1sr_1(alt)(i,ang + shift)<<"\t";
						}else
						{
							of<<mpHFlux->mFluxAcm_2s_1eV_1sr_1(alt)(i,ang + shift)<<"\t";
						}
					}	
				}	
			}	
		}	
	}
	of<<endl;
	of<<Log::msMessageLog<<endl;
	of.close();
}

void  ProtonHydrogenTransport::ForceFluxTot()
{
	std::string vFilename = "FLUXtot.txt";
	ifstream ifs(vFilename.c_str());

	for(size_t alt = 0; alt < mpAltGridKm->size(); ++alt)
	{
		for(size_t a = 0; a < 2; ++a)
		{

			for(size_t sp = 0; sp < 2; ++sp)
			{

				for(size_t i = 0; i < PGRID.size(); ++i)
				{

					for(size_t ang = 0; ang < static_cast<size_t>(mNbAngle / 2); ++ang)
					{
						size_t shift = a * mNbAngle / 2;
						if(sp == 0)
						{
							ifs>>mpProtFlux->mFluxAcm_2s_1eV_1sr_1(alt)(i,ang + shift);
						}else
						{
							ifs>>mpHFlux->mFluxAcm_2s_1eV_1sr_1(alt)(i,ang + shift);
						}
					}	
				}	
			}	
		}	
	}
	ifs.close();
	Log::mW<<"Flux FORCED"<<endl;
}
void  ProtonHydrogenTransport::PrintEnergyConservation(std::string vFilename)
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
	of<<"<!-- "<<" Total input energy: "<<"--><TotalEnergy>"<<mEinput<<"</TotalEnergy>"<<endl;
	of<<"<!-- "<<"============================================="<<" -->"<<endl;
	of<<"<!-- "<<" Total Heating:      "<<"--><AbsorbedEnergy>"<<mAbsorbedEnergy<<"</AbsorbedEnergy>"<<endl;
	of<<"<!-- "<<" Total inelastic:    "<<"--><InelasticEnergy>"<<mInelasticEnergy<<"</InelasticEnergy>"<<endl;
	of<<"<!-- "<<" Total reflected:    "<<"--><ReflectedEnergy>"<<mReflectedEnergy<<"</ReflectedEnergy>"<<endl;
//	of<<"<!-- "<<" Total transmitted:  "<<"--><TransmittedEnergy>"<<mTransmittedEnergy<<"</TransmittedEnergy>"<<endl;
//	of<<"<!-- "<<"=============================================" <<" -->"<<endl;
//	of<<"<!-- "<<" Total absorbed energy (inelastic + heating) " <<"--><ResultTotalEnergy>"<<mResultTotalEnergy<<"</ResultTotalEnergy>"<<endl;
//	of<<"<!-- "<<" Total computed energy (tot abs   + reflect) " <<"--><TotCener>"<< mResultTotalEnergy+mReflectedEnergy+mTransmittedEnergy<<"</TotCener>"<<endl;
	of<<"<!-- "<<"=============================================" <<" -->"<<endl;
	of<<"<!-- "<<"\t Energy conservation error :  "<<"--> <Error1>"<<ntostr(mError1*100.)<<"</Error1><!-- % (>0 : gain)-->"<<endl;
	of<<"<!-- "<<"\t Energy conservation error :  "<<"--> <Error2>"<<ntostr(mError2*100.)<<"</Error2><!-- % (>0 : gain)-->"<<endl;
/*	of<<"<!-- "<<"\t Production total: "<<"--><ProdTot>"<<mProdTot<<"</ProdTot>"<<endl;
	of<<"<!-- "<<"\t Production total, ions: "<<"--><ProdIonTot>"<<mProdIonTot<<"</ProdIonTot>"<<endl;
	of<<"<!-- "<<"============================================= -->"<<endl;
	double tot = mSecondarySumeV+mDegradedSumeV+mCoulombianSumeV;
	of<<"<!-- "<<"Energy from the secondaries      : "<<"--> <SecondarySumeV>"<<mSecondarySumeV<<" </SecondarySumeV> <perc>"<< mSecondarySumeV*100./tot<<"%</perc>"<<endl; 
	of<<"<!-- "<<"Energy from the degraded      : "<<"--> <DegradedSumeV>"<<mDegradedSumeV<<"</DegradedSumeV>  <perc>"<< mDegradedSumeV*100./tot<<"%</perc>"<<endl;  
	of<<"<!-- "<<"Energy from the coulombian      : "<<"--> <CoulombianSumeV>"<<mCoulombianSumeV<<" </CoulombianSumeV> <perc>"<< mCoulombianSumeV*100./tot<<"%</perc>"<<endl;
	*/
	of<<"<!-- "<<Log::msMessageLog<<"-->"<<endl;
	of<<"</xml>"<<endl;
	of.close();
}
