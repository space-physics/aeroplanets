/**
 * \file umathfunction.cpp
 * \brief Implements basic functions, like integration... with ublas vectors
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: umathfunction.cpp 1440 2012-03-05 22:24:56Z gronoff $
 */



#include "umathfunction.hpp"
using namespace std;
using namespace MathString;



namespace MathFunction
{


	void MoyGliss(const ublas::vector<double>& vData, int nMoy, ublas::vector<double>& rResult)
	{
		int ntab=vData.size();
		if(nMoy>ntab)
		{
			rResult.resize(ntab);
			rResult=vData;
			assert((&vData) != (&rResult));// check if this is really a copy!
			return;
		}
		rResult.resize(ntab);

		//		double parite= double(nmoy/2)-double(nmoy)/2;
		int ndeb,nfin;
		//		if(parite==0)
		if(!(nMoy%2))
		{
			ndeb=-(nMoy/2-1);
			nfin=nMoy/2;
		}else
		{
			ndeb=-nMoy/2;
			nfin=nMoy/2;
		}

		for(int i=-ndeb;i<ntab-nfin;i++)
		{
			rResult[i]=0;
			for(int j=ndeb;j<=nfin;j++)
			{
				rResult[i]=rResult[i]+vData[i+j];
			}
			rResult[i]=rResult[i]/static_cast<double>(nMoy);
		}
		// Premiers points:

		for(int i=0;i<-ndeb;i++)
		{
			rResult[i]=0;
			for(int j=i;j<=i+nMoy-1;j++)
			{
				rResult[i]=rResult[i]+vData[j];
			}
			rResult[i]=rResult[i]/static_cast<double>(nMoy);
		}
		//Derniers points
		for(int i=ntab-nfin;i<ntab;i++)
		{
			rResult[i]=0;
			for(int j=i-nMoy+1;j<=i;j++)
			{
				rResult[i]=rResult[i]+vData[j];
			}
			rResult[i]=rResult[i]/static_cast<double>(nMoy);
		}	



	}



	ublas::vector<double> ExponentialDistrib(double vMin, double vMax,int vNbPts)
	{
		double eps=1E-4;
		double dzo=(vMax-vMin)/(2*static_cast<double>(vNbPts));
		dzo=nmin(1.,dzo);
		double hsave=0;
		int itr=0;
		double resid=1;
		double ftab,hnu,scal;
		ftab=vNbPts-1;
		while(resid>eps)
		{
			itr++;
			hnu=ftab*dzo/log((vMax+ hsave)/(vMin+hsave)) -vMin;
			resid=nabs(hsave-hnu);
			scal=nmin( nabs(hsave),nabs(hnu));
			hsave=hnu;
			if (!scal==0)
			{
				resid=resid/scal;
			}
		}

		double ratio=(vMax+hsave)/(vMin+hsave);
		ublas::vector<double> distrib;
		distrib.resize(vNbPts);
		for(int i=0;i<vNbPts;i++)
		{
			distrib[vNbPts-i-1]= (vMin+hsave)*pow(ratio, (static_cast<double>(i))/ftab ) -hsave;
		}
		distrib[0]=vMax;
		distrib[vNbPts-1]=vMin;
		return distrib;
	}



	GaussianAngle::GaussianAngle(int vNbAngle)
	{
		if(! ((vNbAngle%2)==0))
		{
			//	Log::SetPriority(Log::ERROR);
			Log::mE<<"The number of gaussian angle should be even!!!!!"<<endl;
			MathError err("Gaussian Angle : error in the number of angles");
			throw err;
		}
		mNbAngles=vNbAngle;
		mNbmid=vNbAngle/2;
		Gaussa();

		ublas::vector<double> tmp=mGmu;
		std::reverse(tmp.begin(),tmp.end());

		mXmu.resize(mNbAngles);
		mWeight.resize(mNbAngles);
		mAngzb.resize(mNbAngles);
		//	Xmu.assign(tmp.begin(),tmp.end());
		std::copy(tmp.begin(),tmp.end(),mXmu.begin());
		for(unsigned i=0;i<mGmu.size();++i)
		{
			mXmu[i+mNbmid]=mGmu[i]*-1;
		}
		//cout<<"b"<<endl;
		ublas::vector<double> tmp2=mGwt;
		std::reverse(tmp2.begin(),tmp2.end());
		//	mWeight.assign(tmp2.begin(),tmp2.end());
		std::copy(tmp2.begin(),tmp2.end(),mWeight.begin());
		for(unsigned i=0;i<mGwt.size();++i)
		{
			mWeight[i+mNbmid]=mGwt[i];
		}

		mDThetaR.resize(mNbAngles,mNbAngles);
		for(unsigned i=0;i<mXmu.size();++i)
		{
			mAngzb[i]=acos(mXmu[i])*180./PI;

			for(unsigned j=0;j<mXmu.size();++j)
			{
				mDThetaR(i,j)=acos(mXmu[i]) - acos(mXmu[j]);
			}
		}
	}

	void GaussianAngle::Gaussa()
	{

		int m=mNbmid;
		mGmu.resize(m);
		mGwt.resize(m);
		double tol=1E-15;
		//if(m<=5)
		//	tol=1E-30;
		double en=m;
		double np1=m+1.;
		double nnp1=en*np1;
		double cona= (en-1.)/(8*en*en*en);
		int lim=m/2+1; // Attention, ici ca devrait faire comme en python

		if (m < 1)
		{
			cout<<"Bande de petits malins"<<endl;
			mGmu[0]=0.5;
			mGwt[0]=1;
			return;
		}

		double t,x,xi,tmp,pm2,pm1,p,ppr,p2pri;
		tmp=0;
		pm2=0;
		int ggtop=-1;

		for(int k=1; k<lim; ++k)
		{
			t=(4.*static_cast<double>(k)-1.)*PI/(4*en+2);
			x=cos(t+cona/tan(t));
			xi=x;
			ggtop=-1;

			while(abs(xi-x)>tol||ggtop==-1)
			{
				//	G4cout<<FunTrans::dabs(xi-x)<<" > "<<tol<<G4endl;
				x=xi;
				ggtop++;
				pm2=1.;
				pm1=x;
				p=0.;
				for(int nn = 2;nn < m+1; ++nn)
				{
					p=((2*nn-1)*x*pm1-(nn-1)*pm2)/nn;
					pm2=pm1;
					pm1=p;
				}
				tmp=1/(1-x*x);
				ppr=en*(pm2-x*p)*tmp;
				p2pri=(2.*x*ppr-nnp1*p)*tmp;
				xi=x-(p/ppr)*(1+(p/ppr)*p2pri/(2*ppr));
				//	cout<<" x , xi "<<x<<" - "<<xi<<endl;
				if (ggtop> 15)
				{
					exit(1);
				}
			}
			mGmu[k-1]=-x;
			mGwt[k-1]=2/(tmp*en*en*pm2*pm2);
			mGmu[m-k]=x;
			mGwt[m-k]=mGwt[k-1];
		}
		if(m%2!=0)
		{
			mGmu[lim - 1] = 0; // Minus 1 because we start at 0
			double prod = 1;
			for (int k = 3; k <= m; k += 2)
			{
				prod*=static_cast<double>(k)/static_cast<double>(k-1);
			}
			mGwt[lim - 1]=2/(prod*prod);
		}
		//		cout<<"Abscisses\tAngles\tWeight"<<endl;
		for(unsigned int i=0;i<mGmu.size();i++)
		{
			mGmu[i]=0.5*mGmu[i]+0.5;
			mGwt[i]=0.5*mGwt[i];
		}
	}


	ublas::vector<double> IntLin(ublas::vector<double> vOldx,ublas::vector<double> vOldy,ublas::vector<double>vNewx)
	{
		ublas::vector<double> gnewy(vNewx.size());
		gnewy.clear();
		bool reorder=false;
		if(vOldx.size()==0 or vNewx.size()==0)
		{
			return gnewy;
		}
		if(*(vOldx.begin())>*(vOldx.end()-1))
		{
			std::reverse(vOldx.begin(),vOldx.end());
			std::reverse(vOldy.begin(),vOldy.end());
		}
		if(*(vNewx.begin())>*(vNewx.end()-1))
		{
			std::reverse(vNewx.begin(),vNewx.end());
			reorder=true;
		}

		unsigned nin=vOldx.size();
		for(unsigned i=0;i<vNewx.size();++i)
		{
			unsigned j=0;
			bool pasfini=true;
			while( (j<nin) and pasfini)
			{
				/*aif (gnewx[i]<goldx[0])
				  {
				  gnewy[i]=goldy[1]-(goldx[1]-gnewx[i])*(goldy[1]-goldy[0])/(goldx[1]-goldx[0]);
				  cout<<"ici"<<
				  pasfini=false;
				  }else*/ 
				if(vNewx[i]>*(vOldx.end()-1))
				{
					gnewy[i]=*(vOldy.end()-1)-(*(vOldx.end()-1)-vNewx[i])*(*(vOldy.end()-1)-*(vOldy.end()-2))/(*(vOldx.end()-1)-*(vOldx.end()-2));
					pasfini=false;
				}else if(vNewx[i]==vOldx[j])
				{
					gnewy[i]=vOldy[j];
					pasfini=false;
				}else if ( vNewx[i]<vOldx[j])
				{
					if(j==0)
					{
						gnewy[i]=vOldy[j]-(vOldx[j]-vNewx[i])*(vOldy[1]-vOldy[0])/(vOldx[1]-vOldx[0]);
					}else
					{
						gnewy[i]=vOldy[j]-(vOldx[j]-vNewx[i])*(vOldy[j]-vOldy[j-1])/(vOldx[j]-vOldx[j-1]);
					}
					pasfini=false;
				}

				j++;
			}
		}

		if(reorder)
			std::reverse(gnewy.begin(),gnewy.end());
		return gnewy;

	}


	ublas::vector<double> IntLog(const ublas::vector<double>& vOldx,const ublas::vector<double>& vOldy,const ublas::vector<double>& vNewx)
	{
		ublas::vector<double> tmpy(vOldy.size());
		for(unsigned i=0;i<vOldy.size();++i)
		{
			// Attention: on ne fait pas de test pour savoir si c'est <=0...
			assert(vOldy[i]>0);
			tmpy[i]=log10(vOldy[i]);
		}
		ublas::vector<double> resu=IntLin(vOldx,tmpy,vNewx);

		for(unsigned i=0;i<resu.size();++i)
		{
			resu[i]=pow(10.,resu[i]);
		}
		return resu;

	}

	ublas::vector<double> IntLogLog(const ublas::vector<double>& vOldx,const ublas::vector<double>& vOldy,const ublas::vector<double>& vNewx)
	{
		ublas::vector<double> tmpy(vOldy.size()),tmpx(vOldx.size());
		for(unsigned i=0;i<vOldy.size();++i)
		{
			// Attention: on ne fait pas de test pour savoir si c'est <=0...
			assert(vOldy[i]>0);
			tmpy[i]=log10(vOldy[i]);
			assert(vOldx[i]>0);
			tmpx[i]=log10(vOldx[i]);
		}
		ublas::vector<double> tmpnx(vNewx.size());

		for(unsigned i=0;i<vNewx.size();++i)
		{
			assert(vNewx[i]>0);
			tmpnx[i]=log10(vNewx[i]);
		}
		ublas::vector<double> resu=IntLin(tmpx,tmpy,tmpnx);

		for(unsigned i=0;i<resu.size();++i)
		{
			resu[i]=pow(10.,resu[i]);
		}
		return resu;


	}

	double TrapzInt(const ublas::vector<double>& x,const ublas::vector<double>& y,int nb)
	{
		if(nb==-1)
		{
			nb=x.size()-1;

			if(x.size()!=y.size())
			{
				cout<<"Fatal error in hint : vector of different size"<<endl;
				exit(1);
			}


		}else
		{
			if(x.size()<(unsigned)nb)
			{
				cout<<"Fatal error in hint : x too short"<<endl;
				exit(1);
			}
			if(y.size()<(unsigned)nb)
			{
				cout<<"Fatal error in hint : y too short"<<endl;
				exit(1);
			}
		}

		double resu=0;
		for(int i=0;i<nb;++i)
		{
			resu+=0.5*(x[i+1]-x[i])*(y[i+1]+y[i]);
		}



		return nabs(resu);

	}


	bool VectorCompare(const ublas::vector<double>& vVec1,const ublas::vector<double> & vVec2, double vPercent,double vMin)
	{
		if(vVec1.size()!=vVec2.size())
		{
			cout<<"Size mismatch"<<endl;
			return false;
		}

		for(unsigned i=0;i<vVec1.size();++i)
		{
			if(vVec1[i]>vMin||vPercent>100)
			{
				if( ! (nabs(vVec1[i]-vVec2[i])/nmin(nabs(vVec1[i]),nabs(vVec2[i]))*100 < vPercent))
				{
					Log::mI<<" Value mismatch at position : "<<i<<endl;
					Log::mI<<" vec1 :"<<vVec1[i]<<" vec2 : "<<vVec2[i]<<" difference "<< nabs(vVec1[i]-vVec2[i])/nmin(nabs(vVec1[i]),nabs(vVec2[i]))*100 <<"%"<<endl;
					return false;
				}


			}else
			{
				if( (nabs(vVec1[i]-vVec2[i])/nmin(nabs(vVec1[i]),nabs(vVec2[i]))/vPercent*100.>1.))
				{
					Log::mI<<" vec1 :"<<vVec1[i]<<" vec2 : "<<vVec2[i]<<" difference "<< nabs(vVec1[i]-vVec2[i])/nmin(nabs(vVec1[i]),nabs(vVec2[i]))*100 <<"%"<<endl;
					//	return false;
					Log::mI<<"Position = "<<i<<endl;
					Log::mI<<"But, it is below the min value -> just a warning"<<endl;
				}
			}

		}


		return true;


	}


	bool VectorCompare2D(const ublas::matrix<double> & vVec1,const ublas::matrix<double> & vVec2, double vPercent,double vMin)
	{
		if( (vVec1.size1()!=vVec2.size1()) || (vVec1.size2()!=vVec2.size2()) )
		{
			cout<<"Size mismatch"<<endl;
			return false;
		}

		for(unsigned i=0;i<vVec1.size1();++i)
		{
			for(unsigned j=0;j<vVec1.size2();++j)
			{
				if(vVec1(i,j)>vMin||vPercent>100)
				{
					if( ! (nabs(vVec1(i,j)-vVec2(i,j))/nmin(nabs(vVec1(i,j)),nabs(vVec2(i,j)))*100 < vPercent))
					{
						Log::mI<<" Value mismatch at position : "<<i<<endl;
						Log::mI<<" vec1 :"<<vVec1(i,j)<<" vec2 : "<<vVec2(i,j)<<" difference "<< nabs(vVec1(i,j)-vVec2(i,j))/nmin(nabs(vVec1(i,j)),nabs(vVec2(i,j)))*100 <<"%"<<endl;
						return false;
					}


				}else
				{
					if( (nabs(vVec1(i,j)-vVec2(i,j))/nmin(nabs(vVec1(i,j)),nabs(vVec2(i,j)))/vPercent*100.>1.))
					{
						Log::mI<<" vec1 :"<<vVec1(i,j)<<" vec2 : "<<vVec2(i,j)<<" difference "<< nabs(vVec1(i,j)-vVec2(i,j))/nmin(nabs(vVec1(i,j)),nabs(vVec2(i,j)))*100 <<"%"<<endl;
						//	return false;
						Log::mI<<"Position = "<<i<<" , "<<j<<endl;
						Log::mI<<"But, it is below the min value -> just a warning"<<endl;
					}
				}
			}

		}
		return true;


	}





	ublas::vector<double> IntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector<double>* >& vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative)
	{
		assert(vInputSZA.size()==vInputVals.size());
		assert(vInputSZA.size()>0);
		assert(vOutputAlt.size()==vOutputSZA.size());

		size_t outsize=vOutputAlt.size();
		size_t szasize=vInputSZA.size();
		std::deque< ublas::vector<double> > interpolated_altitudes;
		ublas::vector<double> isza(szasize);

		for(unsigned i=0;i<vInputSZA.size();++i)
		{
			isza[i]=vInputSZA[i];
			assert(vInputAlt.size()==vInputVals[i]->size());
			ublas::vector<double> newval=MathFunction::IntLin(vInputAlt,*(vInputVals[i]),vOutputAlt);
			interpolated_altitudes.push_back(newval);
		}
		if(interpolated_altitudes.size()<1)
		{
			MathError err("You should give SZA to interpole for the IntLinPath function");
			throw err;
		}
		if(interpolated_altitudes.size()==1)
		{

			if(vNoNegative)
			{
				NoNegative(interpolated_altitudes[0]);
			}

			return interpolated_altitudes[0];
		}

		ublas::vector<double> resultat(outsize);

		for(unsigned  i=0;i<outsize;++i)
		{
			ublas::vector<double> outarray(1);
			outarray[0]=vOutputSZA[i];
			ublas::vector<double> valarray(szasize);
			for(unsigned j=0;j<szasize;++j)
			{
				valarray[j]=interpolated_altitudes[j][i];
			}
			ublas::vector<double> rez=MathFunction::IntLin(isza,valarray,outarray);
			resultat[i]=rez[0];
		}
		if(vNoNegative)
		{
			NoNegative(resultat);
		}
		return resultat;
	}

	ublas::vector<double> IntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector<double>* >& vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA)
	{
		assert(vInputSZA.size()==vInputVals.size());
		assert(vInputSZA.size()>0);
		assert(vOutputAlt.size()==vOutputSZA.size());

		size_t outsize=vOutputAlt.size();
		size_t szasize=vInputSZA.size();
		std::deque< ublas::vector<double> > interpolated_altitudes;
		ublas::vector<double> isza(szasize);

		for(unsigned i=0;i<vInputSZA.size();++i)
		{
			isza[i]=vInputSZA[i];
			assert(vInputAlt.size()==vInputVals[i]->size());
			ublas::vector<double> newval=MathFunction::IntLog(vInputAlt,*(vInputVals[i]),vOutputAlt);
			interpolated_altitudes.push_back(newval);
		}
		if(interpolated_altitudes.size()<1)
		{
			MathError err("You should give SZA to interpole for the IntLinPath function");
			throw err;
		}
		if(interpolated_altitudes.size()==1)
		{
			return interpolated_altitudes[0];
		}

		ublas::vector<double> resultat(outsize);

		for(unsigned  i=0;i<outsize;++i)
		{
			ublas::vector<double> outarray(1);
			outarray[0]=vOutputSZA[i];
			ublas::vector<double> valarray(szasize);
			for(unsigned j=0;j<szasize;++j)
			{
				valarray[j]=interpolated_altitudes[j][i];
			}
			ublas::vector<double> rez=MathFunction::IntLin(isza,valarray,outarray);
			resultat[i]=rez[0];
		}
		return resultat;
	}

	ublas::vector<double> NoPointIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, std::deque< ublas::vector<double> > vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative)
	{
		std::deque< ublas::vector<double>* > input_vals_nopointer;
		for(size_t i=0;i<vInputVals.size();++i)
		{
			input_vals_nopointer.push_back(&(vInputVals[i]));
		}
		return MathFunction::IntLinPath(vInputAlt,vInputSZA,input_vals_nopointer,vOutputAlt,vOutputSZA,vNoNegative);
	}

	ublas::vector<double> NoPointIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, std::deque< ublas::vector<double> > vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA)
	{
		std::deque< ublas::vector<double>* > input_vals_nopointer;
		for(size_t i=0;i<vInputVals.size();++i)
		{
			input_vals_nopointer.push_back(&(vInputVals[i]));
		}
		ublas::vector<double> resultat=MathFunction::IntLogPath(vInputAlt,vInputSZA,input_vals_nopointer,vOutputAlt,vOutputSZA);
		return resultat;
	}
	std::deque< ublas::vector<double> >  SetIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque<  std::deque< ublas::vector<double> > >& vSetInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative)
	{

		std::deque< ublas::vector<double> > resultat;

		if(vSetInputVals.size()<1)
			return resultat;

		for(size_t i=0;i<vSetInputVals[0].size();++i)
		{
			std::deque< ublas::vector<double>* > inputval;
			for(size_t j=0;j<vSetInputVals.size();++j)
			{
				assert(vSetInputVals[j].size()==vSetInputVals[0].size());
				inputval.push_back(const_cast< ublas::vector<double>* >(&(vSetInputVals[j][i])));

			}
			resultat.push_back(MathFunction::IntLinPath(vInputAlt,vInputSZA,inputval,vOutputAlt,vOutputSZA,vNoNegative));
		}

		return resultat;

	}

	std::deque< ublas::vector<double> >  SetIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque<  std::deque< ublas::vector<double> > >& vSetInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA)
	{

		std::deque< ublas::vector<double> > resultat;

		if(vSetInputVals.size()<1)
			return resultat;

		for(size_t i=0;i<vSetInputVals[0].size();++i)
		{
			//std::deque< ublas::vector<double>* > inputval;
			std::deque< ublas::vector<double> > inputval;
			Log::mD<<"starting number"<<i<<endl;
			for(size_t j=0;j<vSetInputVals.size();++j)
			{
				assert(vSetInputVals[i].size()==vSetInputVals[0].size());

				//MinValue(vSetInputVals[i][j],1E-42);
				//inputval.push_back(const_cast< ublas::vector<double>* >(&(vSetInputVals[i][j])));
				//		inputval.push_back((vSetInputVals[i][j]));
				//		Log::mL<<j<<" -< -1 . "<<inputval.size()<<endl;
				//		MinValue(inputval.at(j),1E-42);
				//		Log::mL<<Next
				ublas::vector<double> tmp=vSetInputVals[j][i];
				MinValue(tmp,1E-42);
				inputval.push_back(tmp);


			}
			resultat.push_back(MathFunction::NoPointIntLogPath(vInputAlt,vInputSZA,inputval,vOutputAlt,vOutputSZA));
		}

		return resultat;

	}
	// matrix(row number, column number)
	ublas::matrix<double> MatIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::matrix<double>* > & vMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative)
	{
		if(vMatInputVals.size()<1)
		{
			ublas::matrix<double> tmpresu;
			return tmpresu;
		}
		ublas::matrix<double> resultat(vOutputAlt.size(),(*(vMatInputVals[0])).size2());

		for(unsigned i=0;i<vMatInputVals[0]->size2();++i)
		{
			std::deque< ublas::vector<double> > tmpdeque;
			for(unsigned j=0;j<vMatInputVals.size();++j)
			{
				ublas::matrix_column< ublas::matrix<double> > mi(*(vMatInputVals[j]),i);
				ublas::vector<double> vec=mi;
				tmpdeque.push_back(vec);
			}

			ublas::matrix_column< ublas::matrix<double> > mo(resultat,i);
			mo=MathFunction::NoPointIntLinPath(vInputAlt,vInputSZA,tmpdeque,vOutputAlt,vOutputSZA,vNoNegative);
		}

		return resultat;
	}

	ublas::matrix<double> MatIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::matrix<double>* > & vMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA)
	{
		if(vMatInputVals.size()<1)
		{
			ublas::matrix<double> tmpresu;
			return tmpresu;
		}
		ublas::matrix<double> resultat(vOutputAlt.size(),(*(vMatInputVals[0])).size2());

		for(unsigned i=0;i<vMatInputVals[0]->size2();++i)
		{
			std::deque< ublas::vector<double> > tmpdeque;
			for(unsigned j=0;j<vMatInputVals.size();++j)
			{
				ublas::matrix_column< ublas::matrix<double> > mi(*(vMatInputVals[j]),i);
				ublas::vector<double> vec=mi;
				MinValue(vec,1E-42);
				tmpdeque.push_back(vec);
			}

			ublas::matrix_column< ublas::matrix<double> > mo(resultat,i);
			mo=MathFunction::NoPointIntLogPath(vInputAlt,vInputSZA,tmpdeque,vOutputAlt,vOutputSZA);
		}

		return resultat;
	}



	ublas::vector< ublas::matrix<double> > VecMatIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector< ublas::matrix<double> >* > & vVecMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative)
	{
		if(vVecMatInputVals.size()<1)
		{
			ublas::vector< ublas::matrix<double> > tmpresu;
			return tmpresu;
		}else if((vVecMatInputVals[0])->size()<1)
		{
			ublas::vector< ublas::matrix<double> > tmpresu;
			return tmpresu;
		}

		ublas::vector< ublas::matrix<double> > resultat(vOutputAlt.size());
		for(unsigned i=0;i<vOutputAlt.size();++i)
		{
			resultat[i].resize((*(vVecMatInputVals[0]))[0].size1(),(*(vVecMatInputVals[0]))[0].size2());

		}

		size_t nbalt=(vVecMatInputVals[0])->size();
		for(unsigned i=0;i<(*(vVecMatInputVals[0]))[0].size1();++i)
		{
			for(unsigned j=0;j<(*(vVecMatInputVals[0]))[0].size2();++j)
			{
				std::deque< ublas::vector<double> > tmpdeque;
				for(unsigned k=0;k<vVecMatInputVals.size();++k)
				{
					ublas::vector<double> tmpalt(nbalt);
					for(unsigned l=0;l<nbalt;++l)
					{
						tmpalt[l]=(*(vVecMatInputVals[k]))[l](i,j);
					}

					tmpdeque.push_back(tmpalt);
				}
				ublas::vector<double> res=MathFunction::NoPointIntLinPath(vInputAlt,vInputSZA,tmpdeque,vOutputAlt,vOutputSZA,vNoNegative);
				for(unsigned l=0;l<vOutputAlt.size();++l)
				{
					resultat[l](i,j)=res[l];
				}
			}
		}
		return resultat;
	}

	ublas::vector< ublas::matrix<double> > VecMatIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector< ublas::matrix<double> >* > & vVecMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA)
	{
		if(vVecMatInputVals.size()<1)
		{
			ublas::vector< ublas::matrix<double> > tmpresu;
			return tmpresu;
		}else if((vVecMatInputVals[0])->size()<1)
		{
			ublas::vector< ublas::matrix<double> > tmpresu;
			return tmpresu;
		}

		ublas::vector< ublas::matrix<double> > resultat(vOutputAlt.size());
		for(unsigned i=0;i<vOutputAlt.size();++i)
		{
			resultat[i].resize((*(vVecMatInputVals[0]))[0].size1(),(*(vVecMatInputVals[0]))[0].size2());

		}

		size_t nbalt=(vVecMatInputVals[0])->size();
		for(unsigned i=0;i<(*(vVecMatInputVals[0]))[0].size1();++i)
		{
			for(unsigned j=0;j<(*(vVecMatInputVals[0]))[0].size2();++j)
			{
				std::deque< ublas::vector<double> > tmpdeque;
				for(unsigned k=0;k<vVecMatInputVals.size();++k)
				{
					ublas::vector<double> tmpalt(nbalt);
					for(unsigned l=0;l<nbalt;++l)
					{
						tmpalt[l]=(*(vVecMatInputVals[k]))[l](i,j);
					}
					MinValue(tmpalt,1E-42);
					tmpdeque.push_back(tmpalt);
				}
				ublas::vector<double> res=MathFunction::NoPointIntLogPath(vInputAlt,vInputSZA,tmpdeque,vOutputAlt,vOutputSZA);
				for(unsigned l=0;l<vOutputAlt.size();++l)
				{
					resultat[l](i,j)=res[l];
				}
			}
		}
		return resultat;
	}

	/*       _\|/_
		 (o o)
		 +----oOO-{_}-OOo-----------------------------+
		 |                                            |
		 |Functions for the definition of new profiles|
		 |                                            |
		 +-------------------------------------------*/

	ublas::vector<double> BatesWalkerProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vT0K, double vTexoK, double vS, double vRkm, double vGoms_2, double vMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());

		double mu=vS+1/(vRkm+vAlt0Km);
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*mu*vTexoK)*1E5; // 1E5 for cm->km

		for(size_t i=0;i<vAltKm.size();++i)
		{
			double dzeta= ( vAltKm[i] - vAlt0Km)*(vRkm+vAlt0Km)/(vAltKm[i]+vRkm);
			double tmi= vT0K/(vTexoK- (vTexoK-vT0K)*exp(-gamma*dzeta));
			resu[i]=vN0cm_3*pow(tmi,1+gamma)*exp(-mu*gamma*dzeta);

			if(isnan(resu[i])||resu[i]<1E-42)
			{
				resu[i]=1E-42;
				//	Log::mL<<resu[i]<<"\t"<<vN0cm_3<<"\t"<<pow(tmi,1+gamma)<<"\t"<<exp(-mu*gamma*dzeta);
				//	Log::mL<<1+gamma<<"\t"<<(vTexoK- (vTexoK-vT0K)*exp(-gamma*dzeta))<<"\t"<<vTexoK<<"\t"<<vT0K;
				//	MathError err("Probleme in resu");
				//	throw err;
			}
		}
		return resu;

	}

	ublas::vector<double> ExpProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vRkm, double vGoms_2, double vMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());
		// gamma=1/H
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		for(size_t i=0;i<vAltKm.size();++i)
		{
			resu[i]=vN0cm_3*exp(-(vAltKm[i] - vAlt0Km)*gamma);
		}
		return resu;
	}

	ublas::vector<double> GaussianProfile(ublas::vector<double> vAltKm, double vAltPeakKm, double vN0cm_3, double vDevKm)
	{
		ublas::vector<double> resu(vAltKm.size());
		for(size_t i=0;i<vAltKm.size();++i)
		{
			resu[i]=vN0cm_3*exp(-0.5*pow((vAltKm[i]-vAltPeakKm)/vDevKm,2) );
		}
		return resu;
	}

	ublas::vector<double> ChapmanProfile(ublas::vector<double> vAltKm, double vAltPeakKm, double vN0cm_3, double vSHkm)
	{
		ublas::vector<double> resu(vAltKm.size());
		for(size_t i=0;i<vAltKm.size();++i)
		{
			resu[i]=vN0cm_3*exp(0.5*( 1.+ (vAltPeakKm-vAltKm[i])/vSHkm - exp((vAltPeakKm-vAltKm[i])/vSHkm) ));
		}
		return resu;
	}



	ublas::vector<double> SplineDeriv(ublas::vector<double> vX, ublas::vector<double> vY,double vD0, double vDn)
	{
		assert(vX.size()==vY.size());
		size_t n = vX.size();
		ublas::vector<double> dy(n), u(n);
		dy.clear();
		u.clear();
		if(vD0>1E33)
		{
			dy[0]=0;
			u[0]=0;
		}else
		{
			dy[0]=-0.5;
			u[0]=(3./(vX[1]-vX[0]))*( (vY[1]-vY[0])/(vX[1]-vX[0]) -vD0 );
		}
		for(size_t i=1; i<n-1;++i)
		{// it is really n-1
			double sig= (vX[i]-vX[i-1])/(vX[i+1]-vX[i-1]);
			double p= sig*dy[i-1]+2.;
			dy[i]= (sig-1.)/p;
			u[i]= (vY[i+1]-vY[i])/(vX[i+1]-vX[i]) - (vY[i]-vY[i-1])/(vX[i]-vX[i-1]);
			u[i]=(6.*u[i]/(vX[i+1]-vX[i-1])-sig*u[i-1])/p;
		}
		double qn=0.;
		if(vDn>1E33)
		{
			u[n-1]=0;
		}else
		{
			qn=0.5;
			u[n-1]=(3./(vX[n-1]-vX[n-2]))*(vDn- (vY[n-1]-vY[n-2])/(vX[n-1]-vX[n-2])  );
		}
		dy[n-1]=(u[n-1]-qn*u[n-2])/(qn*dy[n-2]+1.);
		for(int i=n-2;i>-1;--i)
		{
			dy[i]=dy[i]*dy[i+1]+u[i];
		}
		return dy;

	}

	double SplineInterpP(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vDy, double vPt)
	{
		assert(vX.size()==vY.size());
		unsigned n=vY.size();

		if(vX[0]>vX[n-1])
		{
			std::reverse(vX.begin(),vX.end());
			std::reverse(vY.begin(),vY.end());
		}


		unsigned klo=0;
		unsigned khi=n-1;
		while( khi-klo>1)
		{
			unsigned k=(khi+klo) >>1;
			if(vX[k]<vPt)
				klo=k;
			else
				khi=k;
		}
		double h=vX[khi]-vX[klo];
		if(h==0)
		{
			MathError err("Problem with your X input in SplineInterpP");
			throw err;
		}
		double a=(vX[khi]-vPt)/h;
		double b=(vPt-vX[klo])/h;
		return a*vY[klo]+b*vY[khi]+((a*a*a-a)*vDy[klo]+(b*b*b-b)*vDy[khi])*h*h/6.;
	}
	ublas::vector<double> SplineInterpForce(ublas::vector<double> vX, ublas::vector<double> vY,double vD0, double vDn,ublas::vector<double> vNewx)
	{
		size_t n = vY.size();
		if(vX[0]>vX[n-1])
		{
			std::reverse(vX.begin(),vX.end());
			std::reverse(vY.begin(),vY.end());
		}
		ublas::vector<double> dy=SplineDeriv(vX,vY,vD0,vDn);
		ublas::vector<double> resu(vNewx.size());
		resu.clear();
		for(size_t i=0;i<vNewx.size();++i)
		{
			resu[i]=SplineInterpP(vX,vY,dy,vNewx[i]);
		}
		return resu;
	}

	ublas::vector<double> SplineInterp(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vNewx)
	{
		assert(vX.size()==vY.size());
		size_t n = vY.size();
		if(vX[0]>vX[n-1])
		{
			std::reverse(vX.begin(),vX.end());
			std::reverse(vY.begin(),vY.end());
		}
		double d0=(vY[1]-vY[0])/(vX[1]-vX[0]);
		double dn=(vY[n-1]-vY[n-2])/(vX[n-1]-vX[n-2]);
		ublas::vector<double> dy=SplineDeriv(vX,vY,d0,dn);
		ublas::vector<double> resu(vNewx.size());
		resu.clear();
		for(size_t i=0;i<vNewx.size();++i)
		{
			resu[i]=SplineInterpP(vX,vY,dy,vNewx[i]);
		}
		return resu;
	}

	ublas::vector<double> SplineInterpExp(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vNewx)
	{
		assert(vX.size()==vY.size());
		size_t n = vY.size();
		if(vX[0]>vX[n-1])
		{
			std::reverse(vX.begin(),vX.end());
			std::reverse(vY.begin(),vY.end());
		}
		double d0=(vY[1]-vY[0])/(vX[1]-vX[0]);
		double dn=(vY[n-1]-vY[n-2])/(vX[n-1]-vX[n-2]);
		ublas::vector<double> dy=SplineDeriv(vX,vY,d0,dn);
		assert(dy.size()==vX.size());
		ublas::vector<double> resu(vNewx.size());
		resu.clear();
		for(size_t i=0;i<vNewx.size();++i)
		{

			resu[i]=exp(SplineInterpP(vX,vY,dy,vNewx[i]));
		}
		return resu;
	}

	ublas::vector<double> DoubleExpProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vAlt1Km, double vN0cm_3, double vN1cm_3, double vTexoK, double vTmesoK, double vRkm, double vGoms_2, double vMassamu, double vMixMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());
		// gamma=1/H
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		double gamma2=vMixMassamu*vGoms_2*pow(1+vAlt1Km/vRkm,-2)/(SCALEH_CONST*vTmesoK)*1E5; // 1E5 for cm->km
		for(size_t i=0;i<vAltKm.size();++i)
		{
			resu[i]=vN0cm_3*exp(-(vAltKm[i] - vAlt0Km)*gamma)
				+ vN1cm_3*exp(-(vAltKm[i] - vAlt1Km)*gamma2);
		}
		return resu;

	}



	ublas::vector<double> ChapmanCos(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vSZAdeg, double vRkm, double vGoms_2, double vMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		double multiplicator=1000;
		if(vSZAdeg<89.9 and vSZAdeg >0)
		{
			multiplicator = 1/cos(vSZAdeg*PI/180);

		}else
		{
			Log::mW<<"Your SZA for the chapman cos profile is out of bounds"<<endl;
		}
		for(size_t i=0;i<vAltKm.size();++i)
		{
			double z = gamma * (vAltKm[i] - vAlt0Km);
			resu[i]=vN0cm_3* exp(1-z-exp(-z)*multiplicator );
		}
		return resu;
	}	

	ublas::vector<double> ChapmanVar(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vC, double vRkm, double vGoms_2, double vMassamu)
	{
		double c=vC*2;
		ublas::vector<double> resu(vAltKm.size());
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		for(size_t i=0;i<vAltKm.size();++i)
		{
			double z = gamma * (vAltKm[i] - vAlt0Km);
			resu[i]=vN0cm_3* exp(c*(1-z-exp(-z)) );
		}
		return resu;
	}
	ublas::vector<double> Epstein(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vRkm, double vGoms_2, double vMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());
		double gamma=vMassamu*vGoms_2*pow(1+vAlt0Km/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		for(size_t i=0;i<vAltKm.size();++i)
		{
			double z = gamma/2. * (vAltKm[i] - vAlt0Km);
			double sech = 2/ (exp(z)+exp(-z));
			resu[i]=vN0cm_3* sech * sech;
		}
		return resu;
	}


	ublas::vector<double> ExpHyperbolaProfile(ublas::vector<double> vAltKm, double vAltmesoKm, double vAltexoKm, double vNmesocm_3, double vNexocm_3, double vTmesoK, double vTexoK, double vC, double vRkm, double vGoms_2, double vMassamu, double vMixMassamu)
	{
		ublas::vector<double> resu(vAltKm.size());

		double gammaexo=vMassamu*vGoms_2*pow(1+vAltexoKm/vRkm,-2)/(SCALEH_CONST*vTexoK)*1E5; // 1E5 for cm->km
		double gammameso=vMixMassamu*vGoms_2*pow(1+vAltmesoKm/vRkm,-2)/(SCALEH_CONST*vTmesoK)*1E5; // 1E5 for cm->km

		double Aexo = atan(gammaexo); // In the right quadrant!
		double Ameso = atan(gammameso);

		double lnexo = log(vNexocm_3);
		double lnmesoexo = log(vNmesocm_3 / vNexocm_3);
		double Zi = - lnmesoexo / (gammaexo * gammameso) + vAltexoKm / gammameso - vAltmesoKm / gammaexo;
		double denom = 1. / gammameso - 1. / gammaexo;
		Zi /= denom;

		double lni = exp(lnexo + lnmesoexo / (denom * gammameso) - (vAltexoKm - vAltmesoKm) / denom);



		double adenom =  cos(Ameso + Aexo) + cos(Ameso - Aexo);
		for(size_t i=0;i<vAltKm.size();++i)
		{
			double zmzi = vAltKm[i] - Zi;
			double val = sin(Ameso + Aexo) * zmzi;
			if (Ameso == Aexo)
			{
				resu[i]= lni * exp(-val / adenom);
			}else if(Ameso < Aexo)
			{
				double sqdelt = sqrt(abs(zmzi * zmzi / 2 * (1 - cos(2 * (Ameso - Aexo))) + vC));
				resu[i]= lni * exp(-(val + sqdelt) / adenom);
			}else
			{
				double sqdelt = sqrt(abs(zmzi * zmzi / 2 * (1 - cos(2 * (Ameso - Aexo))) + vC));
				resu[i]= lni * exp(-(val - sqdelt) / adenom);
			}
		}
		return resu;
	}



};

