/**
 * \file chapman.hpp
 * \brief Implements the chapman function and subfunction based on Huestis function http://www-mpl.sri.com/software/chapman/chapman.html
 * Copyright G Gronoff August 2010
 * Last modification : $Id$
 */


#include "chapman.hpp"
using namespace std;



ublas::vector<double> ChapFunction::Chapman(ublas::vector<double> vX, double vChi)
{
	assert(!(vChi<0));
	assert(!(vChi>180));
	double chi=vChi;
	if(vChi>90)
		chi=180-vChi;

	ublas::vector<double> resu(vX.size());
	for(size_t i=0;i<vX.size();++i)
	{
		double atm8=0;
		assert(!(vX[i]<0));
		if(vX[i]<36)
		{
			atm8=ChapDeq(vX[i],chi);
		}else
		{
			atm8=ChapAsy(vX[i],chi);
		}

		if(vChi>90)
		{
			atm8=2*exp(vX[i]*2*pow(sin((90-chi)/(2*DEG_TO_RAD)),2))* CxK1(vX[i]*sin(chi/DEG_TO_RAD)) - atm8;
		}
		resu[i]=atm8;
	}
	return resu;
}



double ChapFunction::ChapAsy(double vX,double vChi)
{
	double sinchi=sin(vChi/DEG_TO_RAD);
	double s1=1+sinchi;
	double rx=sqrt(vX);
	double y0=rx*sqrt(2*pow(sin((90-vChi)/(2*DEG_TO_RAD)),2));
	double c[4]={0,0,0,0};
	double fact=1/sqrt(s1);
	c[0]=fact;
	fact/=s1;
	c[1]=fact*(0.5+sinchi);
	fact/=s1;
	c[2]=-fact*(0.125+.5*sinchi);
	fact/=s1;
	c[3]=fact*(0.0625+0.375*sinchi);
	double dn[4]={0,0,0,0};
	double xi[4]={0,0,0,0};
	Cgd3(y0,dn);
	fact=2*rx;

	for(size_t i=0;i<4;++i)
	{
		xi[i]=fact*dn[i];
		fact/=vX;
	}
	mFn[0]=c[0]*xi[0];
	for(size_t i=1;i<4;++i)
	{
		mFn[i]=mFn[i-1]+c[i]*xi[i];
	}
	return mFn[3];
}

double ChapFunction::ChapDeq(double vX,double vChi)
{
	double alpha=vX*sin(vChi/DEG_TO_RAD);
	mYI0=CyI0(vX,vChi);
	mYK0=CyK0(vX,vChi);
	mXI1=CxI1(alpha);
	mXK1=CxK1(alpha);
	return mXI1*mYK0+mXK1*mYI0;
	
}

double ChapFunction::CyI0(double vX,double vChi)
{

        double aI0[7]={ 1.0000000, 2.4999985E-01, 1.5625190E-02,
           4.3393973E-04, 6.8012343E-06, 6.5601736E-08,
           5.9239791E-10};
        double bI0[9]={ 3.9894228E-01, 4.9822200E-02, 3.1685484E-02,
          -8.3090918E-02, 1.8119815E+00,-1.5259477E+01,
           7.3292025E+01,-1.7182223E+02, 1.5344533E+02};
	double gg[7]={0,0,0,0,0,0,0};
	double qbeta[9]={0,0,0,0,0,0,0,0,0};
	double theta=(90-vChi)/(2*DEG_TO_RAD);
	double sint=sin(theta);
	double cost=cos(theta);
	double sinchit=sin(vChi/DEG_TO_RAD);
	double coschit=cos(vChi/DEG_TO_RAD);
	double sc1m=2*sint*sint;
	double alpha=vX*sinchit;

	double resu=0;
	if(alpha<3.75)
	{
		double x1=3.75/sinchit;
		Cgg06(vX,x1,gg);
		double rho=1;
		if(vX>1)
			rho=1/vX;
		double f=pow(sinchit/rho,2);
		double sum=aI0[6]*gg[6];
		for(int i=5;i>-1;--i)
		{
			sum= sum*f+aI0[i]*gg[i];
		}
		Cgq85(x1*sc1m,qbeta);
		double sum2=bI0[8]*qbeta[8];
		for(int i=7;i>-1;--i)
		{
			sum2= sum2/3.75+bI0[i]*qbeta[i];
		}
		resu=exp(-alpha)*coschit*sum+exp((vX-x1)*sc1m)*sum2*cost*sqrt(2/sinchit);
	}else
	{
		
		Cgq85(vX*sc1m,qbeta);
		double sum2=bI0[8]*qbeta[8];
		for(int i=7;i>-1;--i)
		{
			sum2= sum2/alpha+bI0[i]*qbeta[i];
		}
		resu=sum2*cost*sqrt(2/sinchit);
	}

	return resu;

}
double ChapFunction::CyK0(double vX,double vChi)
{
	return vX*vChi;
}
double ChapFunction::CxI1(double vX)
{
	return vX;
}
double ChapFunction::CxK1(double vX)
{
	return vX;
}

double ChapFunction::cXerfc(double vX)
{
	double T=1/(1+0.5*vX);
	return T*exp( -1.26551223 +T*(1.00002368 +T*( .37409196 +T*(  .09678418 +T*(-.18628806 +T*( .27886807 +T*(-1.13520398 +T*(1.48851587 +T*(-.82215223 +T*   .17087277) ))))))));
}
void ChapFunction::Cgd3(double vYo,double* vDn)
{
	vYo=vDn[0];
#ifdef DEBUG
	assert(vYo<1E42); // just dummy
#endif
}

double ChapFunction::Cgg06(double vX,double vX1,double* vGg)
{
	return vX+vX1+vGg[0];
}


void ChapFunction::Cgq85(double vX,double* vQbeta)
{
	vQbeta[0]=vX;
}



