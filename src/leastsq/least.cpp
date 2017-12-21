
/**
 * \file least.cpp
 * \brief Implementation for the least square system
 * Copyright G Gronoff July 2010
 * Last Modification : $Id$
 */
// We include cminpack here because we do not want the defined variables  to pollute the rest of the code
#include <minpack/cminpack.h>

#include "least.hpp"
using namespace std;

namespace LstSq
{
	int LMleast(LMlsq* pLm, double vTol, std::deque<double>& rParams, ublas::matrix<double> & rCovar)
	{
		int info=0; // initialized here but modified after
		int n=rParams.size();
		int m=static_cast<int>(pLm->Psize());
		if(m<n)
		{
			MathError err("LstSq::LMleast : Not enough points! You have more unknowns than measurements points, therefore, the Levenberg Marquardt algorithm cannot compute you points. Desole.");
			throw err;
		}
		double* fvec = new double[m];
		double* x = new double[n];
		for(unsigned i=0;i<static_cast<unsigned>(n);++i)
		{
			x[i]=rParams[i];
		}

		rCovar.resize(n,n);

		// Working arrays
		int lwa=m*n+5*n+m+1;

		double* wa=new double[lwa];
		int* iwa=new int[n+1];

		// Here, it is typically the lmdif1 definitions
		//
		int maxfev= (n+1)*200; // max number of calls
		int mp5n = m + n * 5;
		int mode=1;
		int nprint=0;
		double gtol=0.; 
		double epsfcn=0.;
    		const double factor = 100.;
		int nfev=0;  // initialized here but modified after

		// We set the xtol to the precision of the machine, vTol is the same if negative
  		double xtol = sqrt(dpmpar(1));
		if(vTol<=0)
			vTol=xtol;

		info=lmdif(Least,pLm,m,n,x,fvec,vTol,xtol,gtol,maxfev,epsfcn,&wa[1], mode, factor, nprint, &nfev, &wa[mp5n + 1], m, &iwa[1], &wa[n + 1], &wa[(n << 1) + 1], &wa[n * 3 + 1], &wa[(n << 2) + 1], &wa[n * 5 + 1]) ;
	// info=lmdif1((least),&chap,m,n,pretrieve,fvec,0.0001,iwa,wa,lwa);
	//
	// We compute the covariance
		double fnorm=enorm(m,fvec);
		double covfac=fnorm*fnorm/(m-n);
		covar(n,&wa[mp5n + 1],m,&iwa[1],vTol,&wa[(n << 1) + 1]);
		for(int i = 0; i<n ; ++i)
		{
			for(int j=0;j<n;++j)
			{
				rCovar(i,j)=covfac*(&wa[mp5n + 1])[(i)*m+j];
			}
		}
	
	
	// We retrieve the interesting parameters
		for(unsigned i=0;i<static_cast<unsigned>(n);++i)
		{
			rParams[i]=x[i];
		}

	// And we delete what we defined!
		delete[] fvec;
		delete[] iwa;
		delete[] wa;
		delete[] x;
		return info;
	}

	int Least(void* pP, int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag)
	{
		return ((LMlsq*)pP)->Exec(vDiffsize,vParamsize,pParam,pDiff,vFlag);
	}
};

