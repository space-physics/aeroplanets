/**
 * \file mathflux.cpp
 * \brief Implements the mathematical function to initialize fluxes
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathflux.cpp 1681 2013-02-17 04:53:03Z gronoff $
 *
 */




#include "mathflux.hpp"
using namespace std;
#define PI M_PI

namespace MathFlux
{

	void InDirac(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,const ublas::vector<double>& vXmu,   ublas::matrix<double>  & rFxDown, ublas::matrix<double> &rFxUp)
	{
		unsigned nang2=vNang/2;
		unsigned nen= vCentEeV.size();
		unsigned ndirac=0;
		rFxUp.resize(nen,nang2); 
		rFxDown.resize(nen, nang2);
		rFxUp.clear();
		rFxDown.clear();


		double phi0=vEnTotErg*6.25E11;
		double epowlaw=0.;
		double powlaw=0.;
		if(vPowLaw)
		{
			if(vCentEeV[nen-1]<10)
			{// we have to get closer to 10. remember vCentEeV is decreasing
				for(unsigned ien=0;ien<nen-1;++ien)
				{
					if(vCentEeV[ien+1]<10)
					{
						epowlaw=vCentEeV[ien];
						powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
						break;
					}
				}
			}else{
				epowlaw=10.;
				powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
			}
		}
		// Search for the inferior energy closest to vEzeroeV
		if( (vCentEeV[nen-1]<vEzeroeV) and (vCentEeV[0]>vEzeroeV))
		{
			for(unsigned ien=0;ien<nen;++ien)
			{
					if(vCentEeV[ien]<vEzeroeV)
					{
						ndirac=ien;
						break;
					}
			}
		}else if(vCentEeV[nen-1]>vEzeroeV)
		{// GG add that part. Theoretically never happens
			ndirac=nen-1;
		}

		if(vNang==2)
		{
			for(unsigned ien=0;ien<nen;++ien)
			{
				double power=powlaw/(pow((vCentEeV[ien]/10.),2));
				if(ien==ndirac)
					power+=phi0;
				rFxDown(ien,0)=power;
				//fxup[ien][0]=0; par defaut
				if(vIsotro!=2)
					rFxUp(ien,0)=power;
			}
			return;
		}
		for(unsigned ien=0;ien<nen;++ien)
		{
			double power=powlaw/(pow((vCentEeV[ien]/10.),2));
			if(ien==ndirac)
				power+=phi0;

			if(vIsotro==2)
			{
				rFxDown(ien,0)=power;
			}else
			{
				for(unsigned iang=0;iang<nang2;++iang)
				{
					if(vIsotro==1)
					{
						rFxDown(ien,iang)=power;
						rFxUp(ien,iang)=power;
					}else
					{// We suppose that vIsotro=0: ie no supplementary options
						rFxDown(ien,iang)=power*vXmu[iang];
						rFxUp(ien,iang)=abs(power*vXmu[iang+nang2]);
					}
				}

			}
		}
	}





	void InGauss(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,const ublas::vector<double>& vXmu,            ublas::matrix<double> & rFxDown, ublas::matrix<double> & rFxUp)
	{
		unsigned nang2=vNang/2;
		unsigned nen= vCentEeV.size();
		rFxUp.resize(nen,nang2); 
		rFxDown.resize(nen, nang2);
		rFxUp.clear();
		rFxDown.clear();
		
		double ehalf=0.1*vEzeroeV;
		double phi0=vEnTotErg*6.25E11/(ehalf*vEzeroeV*sqrt(PI));
		double epowlaw=0.;
		double powlaw=0.;
		if(vPowLaw)
		{
			if(vCentEeV[nen-1]<10)
			{// we have to get closer to 10. remember vCentEeV is decreasing
				for(unsigned ien=0;ien<nen-1;++ien)
				{
					if(vCentEeV[ien+1]<10)
					{
						epowlaw=vCentEeV[ien];
						powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
						break;
					}
				}
			}else{
				epowlaw=10.;
				powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
			}
		}
		if(vNang==2)
		{
			for(unsigned ien=0;ien<nen;++ien)
			{
				double power=powlaw/(pow((vCentEeV[ien]/10.),2));
				power+=phi0*exp( - pow((vCentEeV[ien]-vEzeroeV)/ehalf ,2) );
				rFxDown(ien,0)=power*2;
				//fxup[ien][0]=0; par defaut
				if(vIsotro!=2)
					rFxUp(ien,0)=power*2;
			}
			return;
		}
		for(unsigned ien=0;ien<nen;++ien)
		{
			double power=powlaw/(pow((vCentEeV[ien]/10.),2));
			power+=phi0*exp( - pow((vCentEeV[ien]-vEzeroeV)/ehalf ,2) );

			if(vIsotro==2)
			{
				rFxDown(ien,0)=power;
			}else
			{
				for(unsigned iang=0;iang<nang2;++iang)
				{
					if(vIsotro==1)
					{
						rFxDown(ien,iang)=power;
						rFxUp(ien,iang)=power;
					}else
					{// We suppose that vIsotro=0: ie no supplementary options
						rFxDown(ien,iang)=power*vXmu[iang];
						rFxUp(ien,iang)=abs(power*vXmu[iang+nang2]);
					}
				}

			}
		}
	}


	
	void InMaxw(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,ublas::vector<double> vXmu,           ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp)
	{	
		unsigned nang2=vNang/2;
		unsigned nen= vCentEeV.size();
		
		rFxUp.resize(nen,nang2); 
		rFxDown.resize(nen, nang2);
		rFxUp.clear();
		rFxDown.clear();
		double phi0=vEnTotErg*6.25E11/(4*PI*pow(vEzeroeV,3));
		double epowlaw=0.;
		double powlaw=0.;
		if(vPowLaw)
		{
			if(vCentEeV[nen-1]<10)
			{// we have to get closer to 10. remember vCentEeV is decreasing
				for(unsigned ien=0;ien<nen-1;++ien)
				{
					if(vCentEeV[ien+1]<10)
					{
						epowlaw=vCentEeV[ien];
						powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
						break;
					}
				}
			}else{
				epowlaw=10.;
				powlaw=phi0*epowlaw*exp(-epowlaw/vEzeroeV);
			}
		}
		if(vNang==2)
		{
			for(unsigned ien=0;ien<nen;++ien)
			{
				double power=powlaw/(pow((vCentEeV[ien]/10.),2));
				power+=phi0*vCentEeV[ien]*exp(-vCentEeV[ien]/vEzeroeV);
				rFxDown(ien,0)=power*2;
				//fxup[ien][0]=0; par defaut
				if(vIsotro!=2)
					rFxUp(ien,0)=power*2;
			}
			return;
		}
		for(unsigned ien=0;ien<nen;++ien)
		{
			double power=powlaw/(pow((vCentEeV[ien]/10.),2));
			power+=phi0*vCentEeV[ien]*exp(-vCentEeV[ien]/vEzeroeV);
			if(power<0)
				power=0.;
			power/=3.;

			if(vIsotro==2)
			{
				rFxDown(ien,0)=power;
			}else
			{
				for(unsigned iang=0;iang<nang2;++iang)
				{
					if(vIsotro==1)
					{
						rFxDown(ien,iang)=power*2;
						rFxUp(ien,iang)=power*2;
					}else
					{// We suppose that vIsotro=0: ie no supplementary options
						rFxDown(ien,iang)=power*vXmu[iang]*3;
						rFxUp(ien,iang)=abs(power*vXmu[iang+nang2])*3;
					}
				}

			}
		}
	
	}


	void NormFlux(double vEnTotErg, ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp,ublas::vector<double> vWeight,ublas::vector<double> vXmu,ublas::vector<double> vEeV,ublas::vector<double> vDdengeV)
	{
		double qsum=vEnTotErg*6.25E11;
		double qtot=0;
		unsigned nang=vXmu.size();

		for(unsigned iang=0;iang<nang/2;++iang)
		{
			for(unsigned ien=0;ien<vEeV.size();++ien)
			{
				qtot+=rFxDown(ien,iang)*vWeight[iang]*vXmu[iang]*vEeV[ien]*vDdengeV[ien];
			}
		}
		qtot*=2*PI;
		double xnorm=qsum/qtot;

		// init small fluxes to 1E-10
		for(unsigned iang=0;iang<nang/2;++iang)
		{
			for(unsigned ien=0;ien<vEeV.size();++ien)
			{
				rFxDown(ien,iang)=max(rFxDown(ien,iang)*xnorm,0.);
				rFxUp(ien,iang)=max(rFxUp(ien,iang)*xnorm,0.);
			}
		}
		// verif
		qtot=0.;
		for(unsigned iang=0;iang<nang/2;++iang)
		{
			for(unsigned ien=0;ien<vEeV.size();++ien)
			{
				qtot+=rFxDown(ien,iang)*vWeight[iang]*vXmu[iang]*vEeV[ien]*vDdengeV[ien];
			}
		}
		qtot*=2*PI;


		cout.setf(ios::scientific,ios::floatfield);
		cout.precision(10);
		Log::mL<<"Energie totale integree :"<<qtot/6.25E11<<" desiree : "<<vEnTotErg<<" facteur de normalisation :"<<xnorm<<endl;
	}
	void FillFlux(ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp, ublas::vector<double> vFx)
	{
		unsigned nang = rFxDown.size2();
		unsigned nene = rFxDown.size1();
		Log::mL<<vFx<<endl;
		for(unsigned iang=0; iang < nang / 2; ++iang)
		{
			for(unsigned ien = 0; ien < nene; ++ien)
			{
				rFxDown(ien, iang) = vFx[ien]; 
				rFxUp(ien, iang) = vFx[ien];
			}
		}

	}
};

