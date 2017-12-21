/**
 * \file mathfunction.cpp
 * \brief Implements basic functions, like integration...
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathfunction.cpp 1310 2011-09-14 18:38:47Z gronoff $
 */



#include "mathfunction.hpp"
using namespace std;
using namespace MathString;


namespace PhysTime{

	double CalToJul(int vYear, int vMonth, int vDay, int vHour, int vMin, double vSec)
	{
		if(vMonth < 2)
		{
			vYear -= 1;
			vMonth+=12;
		}
		double D = static_cast<double>(vDay) + static_cast<double>(vHour) / 24. + static_cast<double>(vMin) / 1440. + vSec / 86400.;
		int A = vYear / 100. ;
		int t = A / 4;
		int B = 2 - A + t;
		double ytmp = 365.25 * (static_cast<double>(vYear) + 4716);
		double mtmp = 30.6001 * (static_cast<double>(vMonth) + 1);
		int iytmp = static_cast<int>(ytmp);
		int imtmp = static_cast<int>(mtmp);
		return static_cast<double>(iytmp) + static_cast<double>(imtmp) + D + static_cast<double>(B) - 1524.5; 
	}

	double CalToJul(int vYear, int vMonth, int vDay, double vHour)
	{
		if(vMonth < 2)
		{
			vYear -= 1;
			vMonth+=12;
		}
		double D = static_cast<double>(vDay) + vHour / 24.;
		int A = vYear / 100. ;
		int t = A / 4;
		int B = 2 - A + t;
		double ytmp = 365.25 * (static_cast<double>(vYear) + 4716);
		double mtmp = 30.6001 * (static_cast<double>(vMonth) + 1);
		int iytmp = static_cast<int>(ytmp);
		int imtmp = static_cast<int>(mtmp);
		return static_cast<double>(iytmp) + static_cast<double>(imtmp) + D + static_cast<double>(B) - 1524.5; 
	}
};

namespace MathGrid
{

	ublas::vector<double> GridExp(double vMin, double vMax,int vNbPts)
	{
		return MathFunction::ExponentialDistrib(vMin,vMax,vNbPts);
	}


	ublas::vector<double> GridPow(double begin,double end,int nb_pts)
	{
		ublas::vector<double> tab(nb_pts);
		tab.clear();	
		double step=(end-begin)/((double)(nb_pts-1));
		for(unsigned i=0;i<tab.size();i++)
		{
			tab[i]=end - step*(double)(i);
			tab[i]=((tab[i]-2.*begin)*tab[i]+begin*end)/(end-begin);
		}
		return tab;
	}

	ublas::vector<double> GridCst(double begin,double end,int nb_pts)
	{
		ublas::vector<double> dtab(nb_pts);
		ublas::vector<double> tab(nb_pts);
		dtab.clear();
		tab.clear();

		double dzo=(end-begin)/((double)(nb_pts-1));
		tab[0]=end;
		dtab[0]=dzo;
		for(unsigned i=1;i<tab.size();i++)
		{
			tab[i]=tab[i-1]-dzo;
			dtab[i]=dzo;
		}
		tab[nb_pts-1]=begin;
		
		std::reverse(tab.begin(),tab.end());
		// if dtab is needed: uncomment and create the output
		//std::reverse(dtab.begin(),dtab.end());
		return tab;
	}



	// Subroutine for gridpolo
	double EnTab(int ntab,double spfac,double tmin)
	{
		double x=1;
		double sum=0;
		for(int i=0;i<ntab-1;i++)
		{
			sum+=x;
			//Avoid overflow that may appen at the beginning of the computation
			x*=spfac;
			if(x>1E25)
				break;
		}
		return tmin*(1+2*sum);
	}



	bool GridPolo(int ntab,double tmin,double tmax,ublas::vector<double>& tab,ublas::vector<double>&widthtab,double&spfac)
	{
		if(tmax<tmin)
			return false;
		if(tmin<0)
			return false;
		if(tmin==0||ntab==0)
		{
			return false;
		}
		if((tmax/tmin)<200)
		{
			cout<<"max should be at least 200x min"<<endl;
			cout<<"nb: it is an empirical statement,\n else, it takes too much time to compute" <<endl;
			return false;
		}
		// Initialisation
		tab.resize(ntab);
		tab.clear();
		widthtab.resize(ntab);
		widthtab.clear();
		tab[0]=tmin;
		tab[ntab-1]=tmax;
		double flntab= static_cast<double>(ntab);

		// First estimation
		double a,x;
		a = ( tab[ntab-1] - tab[0] ) / (2* tab[0] );
		x= 2. /  ((flntab-2.)*(flntab-1.));
		x*=(a-flntab+1);
		spfac=1+x;
		double oldspfac=-42;
		// Precision 0.01%
		double eps=0.0001;
		// Ajustement
		int iflag,niter;
		iflag=1;
		niter=0;

		while(iflag!=0)
		{
			niter++;

			if(nabs(oldspfac-spfac)<1E-15)
			{
				//Log::mL.SetPriority(Log::WARNING);
				Log::mW<<"WARNING: spfac not reached in power law grid"<<endl;
				iflag=0;
				//Log::mL.SetPriority(Log::DEBUGG);
			}
			double esup=EnTab(ntab,spfac,tmin);
			oldspfac=spfac;
			if(esup>tmax*(1+eps))
			{//spfac need reduction
				x/=1.2;
				spfac=1+x;
			}else if(esup<tmax*(1-eps))
			{//spfac needs augmentation
				x*=1.15;
				spfac=1+x;
			}else
			{
				iflag=0;
			}
			
		}
		//Set up the grid
		double de;
		de=2*tab[0];
		for(int i=1;i<ntab-1;i++) //tab[ntab-1] already defined
		{
			tab[i]=tab[i-1]+de;
			de*=spfac;
		}
		//Set up the width grid
		double dd;
		tab[ntab-1]=tmax;
		dd= tab[0];
		widthtab[0]=2*dd;

		for(int i=1;i<ntab;i++)
		{
			double ener=tab[i-1]+dd;
			dd=tab[i]-ener;
			widthtab[i]=2*dd;
		}

		return true;
	}



	ublas::vector<double> WidthExpGrid(ublas::vector<double> ttab)

	{
		bool croissant=true;
		if(*(ttab.begin())>*(ttab.end()-1))
		{
			croissant=false;
			std::reverse(ttab.begin(),ttab.end());
		}

		ublas::vector<double> widthtab(ttab.size());
		widthtab.clear();
		if(ttab.size()<3)
		{
			return widthtab;
		}
		double dsup=(ttab[1]+ttab[0])/2;
		double dinf=ttab[0]-dsup;
		if(dsup>ttab[0])
			dinf=0;
		widthtab[0]=dsup-dinf;
		dinf=dsup;
		for(unsigned i=1;i<ttab.size()-1;++i)
		{
			dsup=(ttab[i+1]+ttab[i])/2.;
			widthtab[i]=dsup-dinf;
			dinf=dsup;
		}
		unsigned ntab=ttab.size()-1;
		widthtab[ntab]=widthtab[ntab-1]+(ttab[ntab]-ttab[ntab-1])/(ttab[ntab-2]-ttab[ntab-1])*(widthtab[ntab-2]-widthtab[ntab-1]);
		if(croissant)
		{
			return widthtab;
		}
		std::reverse(widthtab.begin(),widthtab.end());
		return widthtab;

	}

	ublas::vector<double> WidthNonExpGrid(ublas::vector<double> ttab)

	{
		bool croissant=true;
		if(*(ttab.begin())>*(ttab.end()-1))
		{
			croissant=false;
			std::reverse(ttab.begin(),ttab.end());
		}

		ublas::vector<double> widthtab(ttab.size());
		widthtab.clear();
		if(ttab.size()<2)
		{
			return widthtab;
		}

		double dd=ttab[0];
		widthtab[0]=2*dd;
		bool iflag=false;
		double ener=0;
		for(unsigned i=1;i<ttab.size();++i)
		{
			ener=ttab[i-1]+dd;
			dd=ttab[i]-ener;
			if(dd<0)
			{
				cout<<"EXPGRID"<<endl;
				iflag=true;
			}
			widthtab[i]=2*dd;
		}
		if(iflag)
			return WidthExpGrid(ttab);
		if(croissant)
		{
			return widthtab;
		}
		std::reverse(widthtab.begin(),widthtab.end());
		return widthtab;

	}

	ublas::vector<double> WidthGrid(ublas::vector<double> tab,int itype)
	{
		if(itype==1)
		{
			return WidthExpGrid(tab);
		}
		return WidthNonExpGrid(tab);
	}


	ublas::vector<double> InterpolateNonChaoticFlux(const ublas::vector<double>& vFlux1, const ublas::vector<double>& vGridMin1, const ublas::vector<double>& vGridMax1, const ublas::vector<double>& vGridMin2, const ublas::vector<double>& vGridMax2, double &rSupplement)
	{
		// We check the input grids
		unsigned p_1_size=vGridMin1.size();
		assert(p_1_size==vFlux1.size());
		assert(p_1_size==vGridMax1.size());
		unsigned p_2_size=vGridMin2.size();
		assert(p_2_size==vGridMax2.size());
		rSupplement=0; // The supplement is initialized at 0



		// 1 - we create the MeanGrids, which will be used for the interpolation
		// 2 - we divide the flux 1 by the width of the box to create an interpolatable flux
		
		ublas::vector<double> nflux(p_1_size), mgrid1(p_1_size), mgrid2(p_2_size);
		nflux.clear();
		mgrid1.clear();
		mgrid2.clear();
		double numbertot1 = 0;
		double energytot1 = 0;
		for(size_t i = 0; i< p_1_size ; ++i)
		{
			nflux[i] = vFlux1[i] / (vGridMax1[i] - vGridMin1[i]);
			mgrid1[i] = (vGridMax1[i] + vGridMin1[i]) / 2.;
			numbertot1 += vFlux1[i];
			energytot1 += vFlux1[i] * mgrid1[i];
		}
		for(size_t i = 0; i< p_2_size ; ++i)
		{
			mgrid2[i] = (vGridMax2[i] + vGridMin2[i]) / 2.;
		}

		
		// 3 - we perform the interpolation
//		Log::mI<<mgrid1<<endl;
//		Log::mI<<nflux<<endl;
		
//#define TESTFLUXINTERPOLATION
#ifdef TESTFLUXINTERPOLATION
		ofstream of1("testf1");
		of1<<"# Flux initial"<<endl;
		of1<<"# Energy (eV)"<<endl;
		of1<<"# Flux (ph/eV)"<<endl;
		for(size_t i = 0; i< p_1_size ; ++i)
		{
			of1<<mgrid1[i]<<"\t"<<nflux[i]<<endl;
		}	
		of1.close();
#endif
		ublas::vector<double> resu = MathFunction::IntLog(mgrid1, nflux, mgrid2);
		
#ifdef TESTFLUXINTERPOLATION
		ofstream of2("testf2");
		of2<<"# Flux initial"<<endl;
		of2<<"# Energy (eV)"<<endl;
		of2<<"# Flux (ph/eV)"<<endl;
		for(size_t i = 0; i< p_2_size ; ++i)
		{
			of2<<mgrid2[i]<<"\t"<<resu[i]<<endl;
		}	
		of2.close();
#endif
		// 4 - we multiply the interpolated flux by the size of the new boxes to have the result flux
		double numbertot2 = 0;
		double energytot2 = 0;
		for(size_t i = 0; i< p_2_size ; ++i)
		{
			resu[i] *= (vGridMax2[i] - vGridMin2[i]);
			numbertot2 += resu[i];
			energytot2 += resu[i] * mgrid2[i];
		}
		// 5 - we integrate both grids to check the energy conservation
		// The integration in nb have been made in numbertot1 and numbertot2
		// The energy conservation is in energytot1 and energytot2
		rSupplement = energytot1 - energytot2;

		Log::mI<<" InterpolateNonChaoticFlux: The number is "<<numbertot2 - numbertot1<<" (To be compared with the total initial number :" << numbertot1 <<") for an energy of"<< energytot2<< " and the energy difference is : "<< rSupplement <<" Wich means a percentage difference of "<<(energytot1 - energytot2) / energytot1 * 100. <<"%"<<endl;
		
		

		return resu;
	}

	ublas::vector<double> RedistributeChaoticFlux(const ublas::vector<double>& vFlux1,const ublas::vector<double>& vGridMin1,const ublas::vector<double>& vGridMax1,const ublas::vector<double>& vGridMin2,const ublas::vector<double>& vGridMax2,double& rSupplement)
	{
		// Ok, this redistribution is very tricky
		// firstly, we check the entries
//		cout<<"We check the entries"<<endl;
		unsigned p_1_size=vGridMin1.size();
//		cout<<"We make asserts"<<endl;
		assert(p_1_size==vFlux1.size());
		assert(p_1_size==vGridMax1.size());
		unsigned p_2_size=vGridMin2.size();
		assert(p_2_size==vGridMax2.size());

//		cout<<"We initialize supplement"<<endl;

		//	unsigned j=0; // position on the second grid only for the non-chaotic case
		rSupplement=0; // The supplement is initialized at 0
//		cout<<"We initialize resu, size : "<<p_2_size<<endl;
		//vector<double> resu(p_2_size,0); // The output vector is initalized at 0
		ublas::vector<double> resu(p_2_size);
		resu.clear();
		
		//=vGridMin2;
		/*for(unsigned i=0;i<resu.size();++i)
		{
			resu[i]=0;
		}*/

//		cout<<"We do horribe initializaction of second parts"<<endl;
		// Sorry, but the grid must be increasing. 
		// So, I have to invert them if it is decreasing
		ublas::vector<double> tmp_g_min_1=vGridMin1;
		ublas::vector<double> tmp_g_max_1=vGridMax1;

		ublas::vector<double> tmp_g_min_2=vGridMin2;
		ublas::vector<double> tmp_g_max_2=vGridMax2;
		ublas::vector<double> flux1=vFlux1;
		// Even if we are in chaotic case, I do the job for the first line
		// it will be interesting to create the non chaotic case
		if( (*tmp_g_min_1.begin()) > (*(tmp_g_min_1.end()-1)))
		{
//			cout<<"We reverse the first vector"<<endl;
			std::reverse(tmp_g_min_1.begin(),tmp_g_min_1.end());
			std::reverse(tmp_g_max_1.begin(),tmp_g_max_1.end());
			// We also reverse the entry!!!
			std::reverse(flux1.begin(),flux1.end());
		}

		bool output_grid_reversed=false;

		if( (*tmp_g_min_2.begin()) > (*(tmp_g_min_2.end()-1)))
		{

			//cout<<"We reverse the second vector"<<endl;
			std::reverse(tmp_g_min_2.begin(),tmp_g_min_2.end());
			std::reverse(tmp_g_max_2.begin(),tmp_g_max_2.end());
			output_grid_reversed=true;
		}

//		cout<<"We start the loops"<<endl;
		// We start the loops
		//
		for(unsigned i=0;i<p_1_size;++i)
		{// Loop on the grid to redistribute

			// The extreme values for the input grid
			double ta=tmp_g_min_1[i];
			double tb=tmp_g_max_1[i];


			bool a_continue=true;
			unsigned j=0;
//			cout<<"Insertion "<<ta<<" -- "<<tb<<endl;
			if(tb< tmp_g_min_2[j])
			{
				rSupplement+=flux1[i]*(tb+ta)/2;
				a_continue=false;
			}
			while(j<p_2_size&&a_continue)
			{// If no matching was found, we continue

				while(j<p_2_size && tmp_g_max_2[j]<ta)
					++j;
//				cout<<j<<" et "<<p_2_size<<endl;

				if(j<p_2_size)
				{
					a_continue=false;
					double tap=ta;
					while(j<p_2_size && (tb>=tmp_g_max_2[j]))
					{// Our second grid is well formed
						// Therefore, the first time tb>tmp_g_max_2[j] we automatically have tmp_g_min_2[j]<=ta so this algorithm is correct
						resu[j]+=flux1[i]*(tmp_g_max_2[j]-tap)/(tb-ta);
						assert(tb-ta>0);
						tap=tmp_g_max_2[j];
						++j;
					}
					if(j<p_2_size)
					{// We finish the last part
						if(tb-ta>0)
							resu[j]+=flux1[i]*(tb-tap)/(tb-ta);
					}else
					{
						rSupplement+=flux1[i]*(tb-tap)/(tb-ta)*(tb+ta)/2;
						
					}
				}



			}

			if(j>=p_2_size&&a_continue)
			{// When really nothing was found
				// We put the energy in the grid
				rSupplement+=flux1[i]*(tb+ta)/2;
				a_continue=false;
			}
		}
		if(output_grid_reversed)
		{
			std::reverse(resu.begin(),resu.end());
		}
		return resu;
	}


	void AddGridLines(const ublas::vector<double>& vNewGridMin, const ublas::vector<double>& vNewGridMax, ublas::vector<double>& rGridMin, ublas::vector<double>& rGridMax)
	{
		assert(rGridMin.size()>0);
		assert(vNewGridMin.size()==vNewGridMax.size());
		assert(rGridMin.size()==rGridMax.size());
		// We work in the frame of increasing grids, 
		bool output_grid_reversed=false;
		if(*rGridMin.begin()> (*(rGridMin.end()-1)))
		{
			output_grid_reversed=true;
			std::reverse(rGridMin.begin(),rGridMin.end());
			std::reverse(rGridMax.begin(),rGridMax.end());
		}

		for(unsigned i=0;i<vNewGridMin.size();++i)
		{

			if(vNewGridMin[i]> (rGridMin[0]))
			{// The grid cut should be inside the grid

				if(vNewGridMax[i]< (*(rGridMax.end()-1)))
				{
					InsertGridCut(vNewGridMin[i],rGridMin,rGridMax);
					InsertGridCut(vNewGridMax[i],rGridMin,rGridMax);
				}
			}
		}



		if(output_grid_reversed)
		{// We put the grid in the previous order
			std::reverse(rGridMin.begin(),rGridMin.end());
			std::reverse(rGridMax.begin(),rGridMax.end());
		}


	}
	void InsertGridCut(const double & vCut, ublas::vector<double>& rGridMin,ublas::vector<double>&rGridMax)
	{ // This works only for increasing grids. But because it is a subfunction
	// We don't test that
		unsigned i=0;
		while(i<rGridMax.size() && vCut>rGridMax[i])
			++i;
	//	if(i==rGridMax.size())
		
		if(i<rGridMax.size())
		{
			// if rGridMax[i]==vCut, we do nothing
			if(vCut<rGridMax[i])
			{
			//	rGridMin.insert(rGridMin.begin()+i+1,vCut);
			//	rGridMax.insert(rGridMax.begin()+i,vCut);
				rGridMin=UblasInsert(rGridMin,i+1,vCut);
				rGridMax=UblasInsert(rGridMax,i,vCut);
			}
		}
	}


	ublas::vector<double> RedistributeEAAFlux(MathFunction::GaussianAngle* vpGangle,
			ublas::vector<double> vBotE1,
			ublas::vector<double> vTopE1,
			ublas::vector< ublas::matrix<double> > vOldElecDistro,
			ublas::vector<double> vBotE2,
			ublas::vector<double> vTopE2,
			ublas::vector< ublas::matrix<double> >& vNewElecDistro
			)
	{
		if(vOldElecDistro.size()==0)
		{
			ublas::vector<double> reste;
			return reste;
		}
		unsigned nbalt=vOldElecDistro[0].size1();
		unsigned nbang=vOldElecDistro[0].size2();
		unsigned nboe=vOldElecDistro.size();
		unsigned nbne=vBotE2.size();

		ublas::vector<double> reste(nbalt);
		reste.clear();
		ublas::matrix<double> restang(nbalt,nbang);


		vNewElecDistro.resize(nbne);
		for(unsigned i=0;i<nbne;++i)
		{
			vNewElecDistro[i].resize(nbalt,nbang);
			vNewElecDistro[i].clear();
		}

		for(unsigned i=0;i<nbalt;++i)
		{
			for(unsigned k=0;k<nbang;++k)
			{
				ublas::vector<double> odene(nboe);
				for(unsigned e=0;e<nboe;++e)
				{
					odene[e]=vOldElecDistro(e)(i,k);
				}
				double rest=0.;
				ublas::vector<double> nfene=RedistributeChaoticFlux(odene,vBotE1,vTopE1,vBotE2,vTopE2,rest);
				restang(i,k)=rest;
				for(unsigned e=0;e<nbne;++e)
				{
					vNewElecDistro(e)(i,k)=nfene(e);
				}
			}
		}
		//weight

		for(unsigned i=0;i<nbalt;++i)
		{
			for(unsigned k=0;k<nbang;++k)
			{
				reste(i)+=restang(i,k)*vpGangle->mWeight[k];
			}
		}


		return reste;


	}

}

