/**
 * \file marsatmotim.cpp
 * \brief Implements the martian atmosphere from MarTim lecture
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: marsatmotim.cpp 769 2009-12-02 23:27:04Z gronoff $
 */


#include "marsatmotim.hpp"
using namespace std;


MarsAtmoTim::MarsAtmoTim(std::string vCompoFilename, std::string vTempFilename):mCompoFilename(vCompoFilename),mTempFilename(vTempFilename)
{
	ReadCompo();
	ReadTemp();
}


void MarsAtmoTim::ReadTemp()
{

	ifstream ifs(mTempFilename.c_str());
	double dhelp,notempvar,tmplat,tmplon,tmpht;
	ifs>>tmplat>>tmplon>>tmpht>>dhelp>>notempvar;
	if( (tmplat!=mLatDim) || (tmplon!=mLonDim) || (tmpht!=mHDim)   )
	{
		Error err("MarsAtmoTim::ReadTemp","Check dimensions","Dimensions mismatch for composition and temperature files! The composition and temperature should come from a similar simulation (the same simulation if you want to be correct!), which is not the case here!");
		throw err;

	}


	for(unsigned m=0;m<mLatDim;++m)
	{
		for(unsigned n=0;n<mHDim;++n)
		{
			// temp
			for(unsigned i=0;i<mLonDim;++i)
			{
				ifs>>mTempNeutre[m][i][n];
			}
			// alt
			for(unsigned i=0;i<mLonDim;++i)
			{
				ifs>>mAltitude[m][i][n];
			}
			// 
			for(unsigned i=0;i<mLonDim;++i)
			{
				double tmp;
				ifs>>tmp;
			}
		}
	}
	ifs.close();

}


void MarsAtmoTim::ReadCompo()
{


	bfour::extent_gen bfextents;
	bcube::extent_gen bcextents;


	ifstream ifs(mCompoFilename.c_str());
	double dhelp,thingy;
	ifs>>mLatDim>>mLonDim>>mHDim>>dhelp>>thingy;
	
	// resize array and vectors
	mMassAtmo.resize(bcextents[mLatDim][mLonDim][mHDim]);
	mAltitude.resize(bcextents[mLatDim][mLonDim][mHDim]);
	mTempNeutre.resize(bcextents[mLatDim][mLonDim][mHDim]);
	mAtmosphere.resize(bfextents[MARTIM_NEUTRAL_SPECIES][mLatDim][mLonDim][mHDim]);
	mIonosphere.resize(bfextents[MARTIM_ION_SPECIES][mLatDim][mLonDim][mHDim]);


	for(unsigned m=0;m<mLatDim;++m)
	{
		for(unsigned n=0;n<mHDim;++n)
		{
			// mass
			for(unsigned i=0;i<mLonDim;++i)
			{
				ifs>>mMassAtmo[m][i][n];
			}

			// neutral atmo
			for(unsigned k=0;k<MARTIM_NEUTRAL_SPECIES;++k)
			{
				for(unsigned i=0;i<mLonDim;++i)
				{
					double tmp;
						ifs>>tmp;
					mAtmosphere[k][m][i][n]=pow(10,tmp);
				}
			}

			// den3d nden
			for(unsigned i=0;i<2*mLonDim;++i)
			{
				double tmp;
				ifs>>tmp;
			}
			// ions
			for(unsigned k=0;k<MARTIM_ION_SPECIES;++k)
			{
				for(unsigned i=0;i<mLonDim;++i)
				{
					double tmp;
						ifs>>tmp;
					mAtmosphere[k][m][i][n]=tmp;
				}
			}
			// mm and pp of neutral species
			for(unsigned i=0;i<14*mLonDim;++i)
			{
				double tmp;
				ifs>>tmp;
			}
		}
	}
	ifs.close();
}


double MarsAtmoTim::LsToSda(double vLsDeg)
{

	return asin(sin(25.3*PI/180)*sin(vLsDeg*PI/180.));
}


void MarsAtmoTim::ChgtCoord(double vLatDeg,double vLonDeg, double vPsiDeg, double& rLatDeg,double& rLongDeg)
{
	double lat=vLatDeg* DEG_TO_RAD;
	double lon=vLonDeg* DEG_TO_RAD;
	double psi=vPsiDeg* DEG_TO_RAD;
	double sinh=cos(psi)*sin(lat)+sin(psi)*cos(lat)*cos(lon);
	double cosh=sqrt(1-sinh*sinh);
	double rnlat=asin(sinh);
	rLatDeg=rnlat*RAD_TO_DEG;
	double sina=cos(lat)*sin(lon)/cosh;
	double cosa=(cos(psi)*cos(lat)*cos(lon)-sin(psi)*sin(lat))/cosh;
	double rnlong=atan2(sina,cosa);
	rLongDeg=rnlong*RAD_TO_DEG;

}

void MarsAtmoTim::InitRelPosition()
{
	//	mMartimLat=0;
	//	mMartimLong=0;
	double lattodegree=360./(static_cast<double>(mLatDim));
	double lotodegree=180./(static_cast<double>(mLonDim));

	int latmil=static_cast<int>(mLatDim/2);
	int lomil=static_cast<int>(mLonDim/2);

	int mlat=static_cast<int>(mLocalLatDeg/lattodegree)+latmil;
	int mlon=static_cast<int>(mLocalLongDeg/lotodegree)+lomil;

	if(mlat<0. || mlon<0.)
	{
		Error err("MarsTim atmosphere initialization","InitRelPosition","Problem with latitude and longitude conversion to martim grid");
		throw err;
	}
	mMartimLat=static_cast<unsigned>(mlat);
	mMartimLong=static_cast<unsigned>(mlon);

	
}


void MarsAtmoTim::InitLocalCoords(double vLsDeg,double vLocalLatDeg,double vLocalLongDeg)
{
	mLocalLatDeg=vLocalLatDeg;
	mLocalLongDeg=vLocalLongDeg;
	mSdaDeg=LsToSda(vLsDeg);
	mLsDeg=vLsDeg;
	mSZADeg=acos(-cos(mLocalLatDeg*PI/180)*cos(mLocalLongDeg*PI/180.)*cos(mSdaDeg*PI/180.)+sin(mLocalLatDeg*PI/180)*sin(mSdaDeg*PI/180.))*180./PI;

	mLatDeg=0.;
	mLongDeg=0.;
	
	ChgtCoord(mLocalLatDeg,mLocalLongDeg,-mSdaDeg,mLatDeg,mLongDeg);

	// we initialise the position in the array -> this is the closest approx
	// thanks to that, the other functions will be able to use the marstim function
	InitRelPosition();

}


void MarsAtmoTim::InitSubsolarCoords(double vLsDeg, double vLatDeg, double vLongDeg)
{
	mLatDeg=vLatDeg;
	mLongDeg=vLongDeg;
	mSdaDeg=LsToSda(vLsDeg);
	mLsDeg=vLsDeg;
	mLocalLatDeg=0.;
	mLocalLongDeg=0.;
	ChgtCoord(mLatDeg,mLongDeg,mSdaDeg,mLocalLatDeg,mLocalLongDeg);
	mSZADeg=acos(-cos(mLocalLatDeg*PI/180)*cos(mLocalLongDeg*PI/180.)*cos(mSdaDeg*PI/180.)+sin(mLocalLatDeg*PI/180)*sin(mSdaDeg*PI/180.))*180./PI;

	// we initialise the position in the array -> this is the closest approx
	// thanks to that, the other functions will be able to use the marstim function
	InitRelPosition();
}

std::deque< ublas::vector<double> > MarsAtmoTim::ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm)
{
	ublas::vector<double> tmpalt(vAltitudeGridKm.size());
	tmpalt.clear();
	deque< ublas::vector<double> > atmo_resu(MARTIM_NEUTRAL_SPECIES,tmpalt);
	//vector<double> altim(mAltitude[mMartimLat][mMartimLong],mAltitude[mMartimLat][mMartimLong]+mHDim);
	ublas::vector<double> altim(mHDim),tmpsp(mHDim);
	altim.clear();
	tmpsp.clear();

	// Smarter solution should exists... This one works
	for(unsigned k=0;k<mHDim;++k)
	{
		altim[k]=mAltitude[mMartimLat][mMartimLong][k];
	}


	for(unsigned i=0;i<MARTIM_NEUTRAL_SPECIES;++i)
	{
		//	vector<double> tmpsp(mAtmosphere[i][mMartimLat][mMartimLong],mAtmosphere[i][mMartimLat][mMartimLong]+mHDim);
		//
		for(unsigned k=0;k<mHDim;++k)
		{
			tmpsp[k]=mAtmosphere[i][mMartimLat][mMartimLong][k];
		}
		atmo_resu[i]=MathFunction::IntLog(altim,tmpsp,vAltitudeGridKm);
	}
	return atmo_resu;
	
}

std::deque< ublas::vector<double> > MarsAtmoTim::ReturnIono(const ublas::vector<double>& vAltitudeGridKm)
{
	ublas::vector<double> tmpalt(vAltitudeGridKm.size());
	tmpalt.clear();
	deque< ublas::vector<double> > atmo_resu(MARTIM_ION_SPECIES,tmpalt);
	//vector<double> altim(mAltitude[mMartimLat][mMartimLong],mAltitude[mMartimLat][mMartimLong]+mHDim);
	ublas::vector<double> altim(mHDim),tmpsp(mHDim);
	altim.clear();
	tmpsp.clear();

	// Smarter solution should exists... This one works
	for(unsigned k=0;k<mHDim;++k)
	{
		altim[k]=mAltitude[mMartimLat][mMartimLong][k];
	}


	for(unsigned i=0;i<MARTIM_ION_SPECIES;++i)
	{
		//	vector<double> tmpsp(mAtmosphere[i][mMartimLat][mMartimLong],mAtmosphere[i][mMartimLat][mMartimLong]+mHDim);
		//
		for(unsigned k=0;k<mHDim;++k)
		{
			tmpsp[k]=mIonosphere[i][mMartimLat][mMartimLong][k];
		}
		atmo_resu[i]=MathFunction::IntLog(altim,tmpsp,vAltitudeGridKm);
	}
	return atmo_resu;
	
}


ublas::vector<double> MarsAtmoTim::ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm)
{
	ublas::vector<double> altim(mHDim),tmpsp(mHDim);
	altim.clear();
	tmpsp.clear();
	for(unsigned k=0;k<mHDim;++k)
	{
		altim[k]=mAltitude[mMartimLat][mMartimLong][k];
	}
	for(unsigned k=0;k<mHDim;++k)
	{
		tmpsp[k]=mTempNeutre[mMartimLat][mMartimLong][k];
	}
	return MathFunction::IntLog(altim,tmpsp,vAltitudeGridKm);

}
