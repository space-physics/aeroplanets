/**
 * \file irimsis.cpp
 * \brief Implements the wrapper around Iri and Msis models
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: irimsis.cpp 1169 2010-10-22 21:28:32Z gronoff $ 
 */


#include "irimsis.hpp"
using namespace std;
void EarthIriMsis::ChgtCoord(double vLatDeg,double vLonDeg, double vPsiDeg, double& rLatDeg,double& rLongDeg)
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


void EarthIriMsis::DayToSda(double vDayOfYear)
{

           mSdaDeg= 23.5 * sin(2.0*PI*(vDayOfYear-80.0)/365.0);
}


EarthIriMsis::EarthIriMsis(int vModAtmo,double vF107,double vF107bar,int vYear,int vDayOfYear,double vUT,double vTinf,double vAp):mModAtmo(vModAtmo),mF107(vF107),mF107bar(vF107bar),mYear(vYear),mDayOfYear(vDayOfYear),mUT(vUT),mTinf(vTinf),mAp(vAp)
{
	DayToSda(mDayOfYear);

	mbIsMsisDone=false;
	mbIsIriDone=false;
}


void EarthIriMsis::InitLocalCoords(double vLocalLatDeg,double vLocalLongDeg)
{
	mLocalLatDeg=vLocalLatDeg;
	mLocalLongDeg=vLocalLongDeg;
	mSZADeg=acos(-cos(mLocalLatDeg*PI/180)*cos(mLocalLongDeg*PI/180.)*cos(mSdaDeg*PI/180.)+sin(mLocalLatDeg*PI/180)*sin(mSdaDeg*PI/180.))*180./PI;

	mLatDeg=0.;
	mLongDeg=0.;
	
	ChgtCoord(mLocalLatDeg,mLocalLongDeg,-mSdaDeg,mLatDeg,mLongDeg);
}


void EarthIriMsis::InitSubsolarCoords(double vLatDeg, double vLongDeg)
{
	mLatDeg=vLatDeg;
	mLongDeg=vLongDeg;
	mLocalLatDeg=0.;
	mLocalLongDeg=0.;
	ChgtCoord(mLatDeg,mLongDeg,mSdaDeg,mLocalLatDeg,mLocalLongDeg);
	mSZADeg=acos(-cos(mLocalLatDeg*PI/180)*cos(mLocalLongDeg*PI/180.)*cos(mSdaDeg*PI/180.)+sin(mLocalLatDeg*PI/180)*sin(mSdaDeg*PI/180.))*180./PI;
}

void EarthIriMsis::LaunchMsis(const ublas::vector<double>& vAltGrid)
{


	int neutspe=8;
	float dayno=static_cast<float>(mDayOfYear);
	int tamp=mYear/100;
	int nan=mYear-tamp*100;
	float ap=static_cast<float>(mAp);
	float f107=static_cast<float>(mF107);
	float f107bar=static_cast<float>(mF107bar);
	float ut=static_cast<float>(mUT);

	float alat=static_cast<float>(mLocalLatDeg);
	float along=static_cast<float>(mLocalLongDeg);


	int asize=vAltGrid.size();
	vector<float> z(asize);
	for(unsigned i=0;i<static_cast<unsigned>(asize);++i)
	{
		z[i]=static_cast<float>(vAltGrid[i]);
	}
	int sortie=0; // NO output
	int impress=6; // if some output : on the std output

	ublas::matrix<float> densout(asize,8);
	vector<float> fcdens(asize,1.);
	float tinf=static_cast<float>(mTinf);

	ublas::vector<float> tneutre(asize),xmasdens(asize);
	float fctemp=1.;
	atmntr_(&mModAtmo,&neutspe,&dayno,&nan,&ap,&f107bar,&f107,&alat,&along,&ut,&z[0],&asize,&sortie,&densout(0,0), &fcdens[0],&tinf,&tneutre(0),&fctemp,&xmasdens[0],&impress);


	mTempNeutre.resize(asize);
	for(unsigned i=0;i<static_cast<unsigned>(asize);++i)
	{
		mTempNeutre[i]=static_cast<double>(tneutre(i));
	}

	for(unsigned j=0;j<8;++j)
	{
		ublas::vector<double> tmp(asize);
		for(unsigned i=0;i<static_cast<unsigned>(asize);++i)
		{
			tmp[i]=static_cast<double>(densout(i,j));
		}
		mDensities.push_back(tmp);
	}
	mMsisAltGridKm=vAltGrid;

	mbIsMsisDone=true;
}

void EarthIriMsis::LaunchIri(const ublas::vector<double>& vAltGridKm)
{
	int nalt=vAltGridKm.size();
	mIriAltGridKm=vAltGridKm;
	float vmin=static_cast<float>(vAltGridKm[0]);
	float vmax=static_cast<float>(vAltGridKm[nalt-1]);
	if(vmin>vmax)
	{
		swap(vmin,vmax);
	}
	if(abs(vmax-vmin)<1.)
	{
		Error err("LaunchIri","altitudes","altitudes min and max should be different by several km");
		throw err;
	}

	float step=(vmax-vmin-0.1)/50.;

	int jf[12]={1,1,1,1,1,
		    1,1,1,1,1,
		    1,1};
	int jmag=0;

	float alat=static_cast<float>(mLocalLatDeg);
	float along=static_cast<float>(mLocalLongDeg);

	// The minus f107! Minus is important!
	float f107=-static_cast<float>(mF107);
	// Minus day of year! Minus is important!
	int dayno=-mDayOfYear;
	float ut=mUT+25; // Here, +25 is to account for UT, and not local time

	ublas::matrix<float> out(50,11);
	ublas::vector<float> outa(30);
	iris12_(jf,&jmag,&alat,&along,&f107,&dayno,&ut,&vmin,&vmax,&step,&out(0,0),&outa(0));

	ublas::vector<double> altiri(50);

	for(unsigned i=0;i<50;++i)
	{
		altiri[i]=static_cast<double>(vmin)+static_cast<double>(step)*static_cast<double>(i);
	}

	deque< ublas::vector<double> > tempm;
	for(unsigned j=0;j<11;++j)
	{

		ublas::vector<double> tmp(50);
		for(unsigned i=0;i<50;++i)
		{
			tmp[i]=static_cast<double>(out(i,j));
			if(tmp[i]<1E-42)
				tmp[i]=1E-42;
		}

		tempm.push_back(MathFunction::IntLog(altiri,tmp,vAltGridKm));
	}
	tempm[0]/=1E6; // m-3 -> cm-3

	mIriTemp.push_back(tempm[1]);
	mIriTemp.push_back(tempm[2]);
	mIriTemp.push_back(tempm[3]);
	ublas::vector<double> ntemp= ReturnNTemp(vAltGridKm);
	for(size_t i=0; i<vAltGridKm.size() ; ++i)
	{
		for(size_t k = 0 ; k<mIriTemp.size() ; ++ k)
		{
			if (mIriTemp[k][i] < ntemp[i])
			{
				mIriTemp[k][i] = ntemp[i];
			}
		}

	}
	mIonDensities.push_back(tempm[0]);
	mIonDensities.push_back(tempm[0]); // ! yes, it is really 6 times
	mIonDensities.push_back(tempm[0]); // The electron density !
	mIonDensities.push_back(tempm[0]);
	mIonDensities.push_back(tempm[0]);
	mIonDensities.push_back(tempm[0]);

	for(unsigned j=4;j<9;++j)
	{
		for(unsigned i=0;i<static_cast<unsigned>(nalt);++i)
		{
			mIonDensities[j-3][i]*=tempm[j][i]/100.; // percent
		}
	}
	mbIsIriDone=true;
}

std::deque< ublas::vector<double> > EarthIriMsis::ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm)
{

	if(!mbIsMsisDone)
	{
		LaunchMsis(vAltitudeGridKm);
	}

	if(MathFunction::VectorCompare(vAltitudeGridKm,mMsisAltGridKm))
	{
		return mDensities;
	}
	deque< ublas::vector<double> > Dens;

	for(unsigned i=0;i<mDensities.size();++i)
	{
		Dens.push_back(MathFunction::IntLog(mMsisAltGridKm, mDensities[i], vAltitudeGridKm));
	}
	return Dens;
}


ublas::vector<double> EarthIriMsis::ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm)
{

	if(!mbIsMsisDone)
	{
		LaunchMsis(vAltitudeGridKm);
	}


	if(MathFunction::VectorCompare(vAltitudeGridKm,mMsisAltGridKm))
	{
		return mTempNeutre;
	}
	return MathFunction::IntLog(mMsisAltGridKm,mTempNeutre,vAltitudeGridKm);
}

ublas::vector<double> EarthIriMsis::ReturnITemp(const ublas::vector<double>& vAltitudeGridKm)
{
	if(!mbIsIriDone)
	{
		LaunchIri(vAltitudeGridKm);
	}
	if(MathFunction::VectorCompare(vAltitudeGridKm,mIriAltGridKm))
	{
		return mIriTemp[1];
	}
	return MathFunction::IntLog(mIriAltGridKm,mIriTemp[1],vAltitudeGridKm);
}

ublas::vector<double> EarthIriMsis::ReturnETemp(const ublas::vector<double>& vAltitudeGridKm)
{
	if(!mbIsIriDone)
	{
		LaunchIri(vAltitudeGridKm);
	}
	if(MathFunction::VectorCompare(vAltitudeGridKm,mIriAltGridKm))
	{
		return mIriTemp[2];
	}
	ublas::vector<double> etemp= MathFunction::IntLog(mIriAltGridKm,mIriTemp[2],vAltitudeGridKm);
	ublas::vector<double> ntemp= ReturnNTemp(vAltitudeGridKm);
	for(unsigned i=0;i<vAltitudeGridKm.size();++i)
	{
		if(etemp[i]<ntemp[i])
			etemp[i]=ntemp[i];
	}
	return etemp;
}


std::deque< ublas::vector<double> > EarthIriMsis::ReturnIono(const ublas::vector<double>& vAltitudeGridKm)
{
	Log::mL<<"Return iono"<<endl;
	if(!mbIsIriDone)
	{
		Log::mL<<"launch iri"<<endl;
		LaunchIri(vAltitudeGridKm);
	}
	Log::mL<<"vector compare"<<endl;
	if(MathFunction::VectorCompare(vAltitudeGridKm,mIriAltGridKm))
	{
		Log::mL<<"return inside"<<endl;
		return mIonDensities;
	}
	deque< ublas::vector<double> > Dens;
	Log::mL<<"Int log"<<endl;

	for(unsigned i=0;i<mIonDensities.size();++i)
	{
		Log::mL<<"Ion density log for "<<i<<" / "<<mIonDensities.size()<<endl;
		Dens.push_back(MathFunction::IntLog(mIriAltGridKm,mIonDensities[i],vAltitudeGridKm));
	}
	return Dens;

}

