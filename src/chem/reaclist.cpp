/**
 * \file reaclist.cpp
 * \brief Implements the different chemical reactions. Inherit from ChemReact.
 * These reactions will be used in the Chem class to compute the densities, and
 * in the emission class to compute the emissions
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: reaclist.cpp 1525 2012-07-03 17:53:32Z gronoff $
 *
 */
#include "reaclist.hpp"


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-------------+
 |                            |
 |O(1D) radiative deactivation|
 |                            |
 +---------------------------*/

Reac0::Reac0(XmlParameters* vpParam):ChemReact(vpParam,0)
{
	mbIsEmit=true;
	mEmitFreqnm=630.0304;
	mReactant.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=5.63E-3;
	mUncertainty="20%";
	mSupplInfo="NIST - Uncertainty estimated with galavis 97";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac0::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac1::Reac1(XmlParameters* vpParam):ChemReact(vpParam,1)
{
	mbIsEmit=true;
	mEmitFreqnm=636.3776;
	mReactant.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=1.82E-3;
	mUncertainty="50%";
	mSupplInfo="NIST - Uncertainty estimated";
	ReadParameters();
}

ublas::vector<double> Reac1::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}




Reac2::Reac2(XmlParameters* vpParam):ChemReact(vpParam,2)
{
	mbIsEmit=true;
	mEmitFreqnm=639.1733;
	mReactant.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=8.60E-7;
	mUncertainty="30%";
	mSupplInfo="NIST - Uncertainty estimated";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac2::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}
/*       _\|/_
         (o o)
 +----oOO-{_}-OOo------------+
 |                           |
 |O(1D) chemical deactivation|
 |                           |
 +--------------------------*/

Reac3::Reac3(XmlParameters* vpParam):ChemReact(vpParam,3)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("O2",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=3.3E-11;
	mUncertainty="30%";
	mSupplInfo="JPL - additional factor: exp(50/Tn). ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac3::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue*exp(50/vTn[i]);
	return resu;
}


Reac4::Reac4(XmlParameters* vpParam):ChemReact(vpParam,4)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=2.2E-11;
	mUncertainty="15%";
	mSupplInfo=" additional factor: (Tn/300)**0.14. Main constant from Kalogerakis 06, dependance from Fox 01 ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac4::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue*pow(vTn[i]/300.,0.14);
	return resu;
}

Reac5::Reac5(XmlParameters* vpParam):ChemReact(vpParam,5)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("N2",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=2.15E-11;
	mUncertainty="10%";
	mSupplInfo=" additional factor: exp(110/Tn)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac5::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue*exp(110/vTn[i]);
	return resu;
}


Reac6::Reac6(XmlParameters* vpParam):ChemReact(vpParam,6)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=2.87E-10;
	mUncertainty="15%";
	mSupplInfo=" additional factor: (Te/300)**0.91";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac6::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue*pow(vTe[i]/300.,0.91);
	return resu;
}




Reac7::Reac7(XmlParameters* vpParam):ChemReact(vpParam,7)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("CO2",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=7.50E-11;
	mUncertainty="15%";
	mSupplInfo=" additional factor: exp(115/Tn)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac7::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue*exp(115/vTn[i]);
	return resu;
}



Reac8::Reac8(XmlParameters* vpParam):ChemReact(vpParam,8)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O","1D"));
	mCatalys.push_back(SpecieId("CO",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=3.60E-10;
	mUncertainty="15%";
	mSupplInfo=" Old reaction rate ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac8::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}



Reac9::Reac9(XmlParameters* vpParam):ChemReact(vpParam,9)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2+",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O","1D"));
	mEfficiency[SpecieId("O","1D")]=1;
	mMainValue=1.;
	mUncertainty="14%";
	mSupplInfo="Kella 97 Sum of two recombination channels. Depends on Te. (Main value set to 1) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac9::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	double k1=0.; // Reaction rate for 1 O(1D) created
	double k2=0.; // Reaction rate for 2 O(1D) created
	for(unsigned i=0;i<resu.size();++i)
	{
		if(vTe[i]<1200.)
		{
			k1=(0.86E-7+9.75E-9)*pow(300./vTe[i],0.70);
			k2=(6.05E-8)*pow(300./vTe[i],0.70);
		}else
		{
			k1=(3.25E-8+3.69E-9)*pow(1200./vTe[i],0.56);
			k2=(2.29E-8)*pow(1200./vTe[i],0.56);
		}
		resu[i]=mMainValue*(k1+2*k2);
	}
	return resu;
}




Reac10::Reac10(XmlParameters* vpParam):ChemReact(vpParam,10)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O","1D"));
	mEfficiency[SpecieId("O","1D")]=1;
	mMainValue=1.;
	mUncertainty="20%";
	mSupplInfo=" Complicated reaction. (Mantas; Main value set to 1) 20\% estimated through seff uncertainty and \% in results";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac10::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	double kmantas=0.; // Reaction rate for 1 O(1D) created
	for(unsigned i=0;i<resu.size();++i)
	{
          	kmantas=0.596*pow((vTe(i)),(0.5))*(9329+vTe(i))*exp(-22756./vTe(i));
		resu[i]=mMainValue*kmantas*(pow((51183+vTe[i]),-3));
	}
	return resu;
}



Reac11::Reac11(XmlParameters* vpParam):ChemReact(vpParam,11)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2",""));
	mReactant.push_back(SpecieId("N+",""));
	mProducts.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("NO+",""));
	mEfficiency[SpecieId("O","1D")]=1;
	mEfficiency[SpecieId("NO+","")]=1;
	mMainValue=1.8E-10;
	mUncertainty="50%";
	mSupplInfo="Langford 86 ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac11::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac12::Reac12(XmlParameters* vpParam):ChemReact(vpParam,12)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("e",""));
	mReactant.push_back(SpecieId("CO+",""));
	mProducts.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("C",""));
	mEfficiency[SpecieId("O","1D")]=1;
	mEfficiency[SpecieId("C","")]=1;
	mMainValue=0.25E-7;
	mUncertainty="15%";
	mSupplInfo="Dependance pow(300/vTe,0.55) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac12::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*pow(300/vTe[i],0.55);
	}
	return resu;
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |O1S emissions   |
 |                |
 +---------------*/



Reac13::Reac13(XmlParameters* vpParam):ChemReact(vpParam,13)
{
	mbIsEmit=true;
	mEmitFreqnm=557.734;
	mReactant.push_back(SpecieId("O","1S"));
	mProducts.push_back(SpecieId("O","1D"));
	mEfficiency[SpecieId("O","1D")]=1;
	mMainValue=1.26;
	mUncertainty="10%";
	mSupplInfo=" NIST - Uncertainty estimated; Problem: 5577/2972 differs from Slanger!";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac13::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac14::Reac14(XmlParameters* vpParam):ChemReact(vpParam,14)
{
	mbIsEmit=true;
	mEmitFreqnm=297.229;
	mReactant.push_back(SpecieId("O","1S"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=0.0754;
	mUncertainty="10%";
	mSupplInfo=" NIST - Uncertainty estimated; Problem: 5577/2972 differs from Slanger!";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac14::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-----------+
 |                          |
 |O1S  chemical deactivation|
 |                          |
 +-------------------------*/


Reac15::Reac15(XmlParameters* vpParam):ChemReact(vpParam,15)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O","1S"));
	mCatalys.push_back(SpecieId("O2",""));
	mProducts.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","1D")]=1.36/4.4;
	mEfficiency[SpecieId("O","")]=1-1.36/4.4;
	mMainValue=4.4E-12;
	mUncertainty="40%";
	mSupplInfo="Dependance exp(-865/vTn) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac15::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*exp(-865/vTn[i]);
	}
	return resu;
}


Reac16::Reac16(XmlParameters* vpParam):ChemReact(vpParam,16)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("O","1S"));
	mCatalys.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=2E-14;
	mUncertainty="50%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac16::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac17::Reac17(XmlParameters* vpParam):ChemReact(vpParam,17)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("O","1S"));
	mCatalys.push_back(SpecieId("CO",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=7.4E-14;
	mUncertainty="6%";
	mSupplInfo=" exp(-961/vTn)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac17::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*exp(-961/vTn[i]) ;
	}
	return resu;
}




Reac18::Reac18(XmlParameters* vpParam):ChemReact(vpParam,18)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O","1S"));
	mCatalys.push_back(SpecieId("CO2",""));
	mProducts.push_back(SpecieId("O","1D"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","1D")]=2.02/3.11;
	mEfficiency[SpecieId("O","")]=1-2.02/3.11;
	mMainValue=3.21E-11;
	mUncertainty="8%";
	mSupplInfo="Dependance exp(-1327/vTn) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac18::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*exp(-1327/vTn[i]);
	}
	return resu;
}



Reac19::Reac19(XmlParameters* vpParam):ChemReact(vpParam,19)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("O","1S"));
	mCatalys.push_back(SpecieId("N2",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=5E-17;
	mUncertainty="50%";
	mSupplInfo="Negligible, estimated uncertainty ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac19::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



Reac20::Reac20(XmlParameters* vpParam):ChemReact(vpParam,20)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O","1S"));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=8.656E-9;
	mUncertainty="50%";
	mSupplInfo="Estimated uncertainty (theoretical work) dependance (pow((vTe)/300,0.94)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac20::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*(pow((vTe[i])/300,0.94));
	}
	return resu;
}




Reac21::Reac21(XmlParameters* vpParam):ChemReact(vpParam,21)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2+",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O","1S"));
	mEfficiency[SpecieId("O","1S")]=1;
	mMainValue=1.;
	mUncertainty="40%";
	mSupplInfo="Kella 97 Sum of two recombination channels. Depends on Te. (Main value set to 1) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac21::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	double k1=0.; // Reaction rate for O(1S) creation 
	for(unsigned i=0;i<resu.size();++i)
	{
		if(vTe[i]<1200.)
		{
			k1=(9.75E-9)*pow(300./vTe[i],0.70);
		}else
		{
			k1=(3.69E-9)*pow(1200./vTe[i],0.56);
		}
		resu[i]=mMainValue*(k1);
	}
	return resu;
}


Reac22::Reac22(XmlParameters* vpParam):ChemReact(vpParam,22)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2+",""));
	mReactant.push_back(SpecieId("N",""));
	mProducts.push_back(SpecieId("O","1S"));
	mProducts.push_back(SpecieId("NO+",""));
	mEfficiency[SpecieId("O","1D")]=1;
	mEfficiency[SpecieId("NO+","")]=1;
	mMainValue=0.9*2.5E-11;
	mUncertainty="50%";
	mSupplInfo="Est uncert. Controversial reaction. Frederick/Kopp. See Gronoff 2008 ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac22::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}




Reac23::Reac23(XmlParameters* vpParam):ChemReact(vpParam,23)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O","1S"));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O","1D"));
	mEfficiency[SpecieId("O","1D")]=1;
	mMainValue=8.50E-9;
	mUncertainty="50%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac23::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |Work with N(2D) |
 |                |
 +---------------*/


Reac24::Reac24(XmlParameters* vpParam):ChemReact(vpParam,24)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O",""));
	mReactant.push_back(SpecieId("N2+",""));
	mProducts.push_back(SpecieId("N","2D"));
	mEfficiency[SpecieId("N","2D")]=1;
	mMainValue=1.4E-10;
	mUncertainty="0%";
	mSupplInfo="Supplementary pow((vTn+vTi)/600,-0.44) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac24::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*pow((vTn[i]+vTi[i])/600,-0.44);
	}
	return resu;
}

Reac25::Reac25(XmlParameters* vpParam):ChemReact(vpParam,25)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("e",""));
	mReactant.push_back(SpecieId("N2+",""));
	mProducts.push_back(SpecieId("N","2D"));
	mEfficiency[SpecieId("N","2D")]=1;
	mMainValue=1.8E-7;
	mUncertainty="0%";
	mSupplInfo="Supplementary pow((vTe)/300,-0.39) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac25::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*pow((vTe[i])/300,-0.39);
	}
	return resu;
}


Reac26::Reac26(XmlParameters* vpParam):ChemReact(vpParam,26)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("e",""));
	mReactant.push_back(SpecieId("NO+",""));
	mProducts.push_back(SpecieId("N","2D"));
	mEfficiency[SpecieId("N","2D")]=1;
	mMainValue=0.78*4.2E-7;
	mUncertainty="0%";
	mSupplInfo="Supplementary pow((vTe)/300,-0.85) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac26::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*pow((vTe[i])/300,-0.85);
	}
	return resu;
}



Reac27::Reac27(XmlParameters* vpParam):ChemReact(vpParam,27)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2",""));
	mReactant.push_back(SpecieId("N+",""));
	mProducts.push_back(SpecieId("N","2D"));
	mEfficiency[SpecieId("N","2D")]=1;
	mMainValue=2E-10;
	mUncertainty="0%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac27::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac28::Reac28(XmlParameters* vpParam):ChemReact(vpParam,28)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2",""));
	mReactant.push_back(SpecieId("N","2D"));
	mProducts.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("O","1D"));
	mEfficiency[SpecieId("O","1D")]=5./5.3;
	mEfficiency[SpecieId("O","")]=0.3/5.3;
	mMainValue=5.3E-12;
	mUncertainty="30%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac28::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



Reac29::Reac29(XmlParameters* vpParam):ChemReact(vpParam,29)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O",""));
	mReactant.push_back(SpecieId("N","2D"));
	mMainValue=1.4E-12;
	mUncertainty="80%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac29::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac30::Reac30(XmlParameters* vpParam):ChemReact(vpParam,30)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("e",""));
	mReactant.push_back(SpecieId("N","2D"));
	mProducts.push_back(SpecieId("N",""));
	mMainValue=5.5E-10;
	mUncertainty="0%";
	mSupplInfo="Additional sqrt(vTe/300)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac30::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*sqrt(vTe[i]/300);
	}
	return resu;
}





Reac31::Reac31(XmlParameters* vpParam):ChemReact(vpParam,31)
{
	mbIsEmit=true;
	mEmitFreqnm=520;
	mReactant.push_back(SpecieId("N","2D"));
	mProducts.push_back(SpecieId("N",""));
	mEfficiency[SpecieId("N","")]=1;
	mMainValue=1.06E-5;
	mUncertainty="0%";
	mSupplInfo="";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac31::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac32::Reac32(XmlParameters* vpParam):ChemReact(vpParam,32)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O2+",""));
	mReactant.push_back(SpecieId("N","2D"));
	mMainValue=2.5E-10;
	mUncertainty="0%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac32::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac33::Reac33(XmlParameters* vpParam):ChemReact(vpParam,33)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O+",""));
	mReactant.push_back(SpecieId("N","2D"));
	mMainValue=1.3E-10;
	mUncertainty="0%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac33::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac34::Reac34(XmlParameters* vpParam):ChemReact(vpParam,34)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO2",""));
	mReactant.push_back(SpecieId("N","2D"));
	mMainValue=3.60E-13;
	mUncertainty="50%";
	mSupplInfo=" Fox 2001 - Herron 1999 ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac34::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |VEGARD KAPLAN   |
 |                |
 +---------------*/

Reac35::Reac35(XmlParameters* vpParam):ChemReact(vpParam,35)
{
	mbIsEmit=true;
	mEmitFreqnm=270;
	mReactant.push_back(SpecieId("N2","A3S"));
	mProducts.push_back(SpecieId("N2",""));
	mEfficiency[SpecieId("N","")]=1;
	mMainValue=1/2.37;
	mUncertainty="10%";
	mSupplInfo="Vegard -kaplan : BAND!";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac35::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac36::Reac36(XmlParameters* vpParam):ChemReact(vpParam,36)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("O",""));
	mReactant.push_back(SpecieId("N2","A3S"));
	mProducts.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("O","1S"));
	mProducts.push_back(SpecieId("N2",""));
	mEfficiency[SpecieId("O","1S")]=0.47;
	mEfficiency[SpecieId("O","")]=0.53;
	mEfficiency[SpecieId("N2","")]=1;
	mMainValue=3.2E-11;
	mUncertainty="25%";
	mSupplInfo="fact : sqrt(T/298) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac36::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*sqrt(vTn[i]/298.);
	}
	return resu;
}



Reac37::Reac37(XmlParameters* vpParam):ChemReact(vpParam,37)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("N2","A3S"));
	mReactant.push_back(SpecieId("CO",""));
	mEfficiency[SpecieId("CO","")]=1;
	mEfficiency[SpecieId("N2","")]=1;


	mMainValue=1.60E-12;
	mUncertainty="28%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac37::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac38::Reac38(XmlParameters* vpParam):ChemReact(vpParam,38)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("N2","A3S"));
	mReactant.push_back(SpecieId("CO2",""));
	mEfficiency[SpecieId("CO2","")]=1;
	mEfficiency[SpecieId("N2","")]=1;
	mMainValue=9.9E-15;
	mUncertainty="12%";
	mSupplInfo=" ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac38::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}




Reac39::Reac39(XmlParameters* vpParam):ChemReact(vpParam,39)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("N2","A3S"));
	mReactant.push_back(SpecieId("O2",""));
	mEfficiency[SpecieId("O2","")]=1;
	mEfficiency[SpecieId("N2","")]=1;
	mMainValue=2.50E-12;
	mUncertainty="20%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac39::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-+
 |                |
 |Cameron bands   |
 |                |
 +---------------*/


Reac40::Reac40(XmlParameters* vpParam):ChemReact(vpParam,40)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO2+",""));
	mReactant.push_back(SpecieId("e",""));
	mEfficiency[SpecieId("O","")]=1;
	mEfficiency[SpecieId("CO","")]=1;
	mEfficiency[SpecieId("CO","a3Pi")]=0.29;
	mEfficiency[SpecieId("CO","A1Pi")]=0.05;
	mMainValue=3.5E-7;
	mUncertainty="30%";
	mSupplInfo=" sqrt(300/Te) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac40::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*sqrt(300/vTe[i]);;
	}
	return resu;
}




Reac41::Reac41(XmlParameters* vpParam):ChemReact(vpParam,41)
{
	mbIsEmit=true;
	mEmitFreqnm=200;
	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=1/(3.7E-3);
	mUncertainty="10%";
	mSupplInfo="Cameron : BAND!";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac41::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac42::Reac42(XmlParameters* vpParam):ChemReact(vpParam,42)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mCatalys.push_back(SpecieId("CO2",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=1E-11;
	mUncertainty="20%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac42::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac43::Reac43(XmlParameters* vpParam):ChemReact(vpParam,43)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mCatalys.push_back(SpecieId("CO",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=5.7E-11;
	mUncertainty="20%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac43::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}
Reac44::Reac44(XmlParameters* vpParam):ChemReact(vpParam,44)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mCatalys.push_back(SpecieId("NO",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=17E-11;
	mUncertainty="35%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac44::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}
Reac45::Reac45(XmlParameters* vpParam):ChemReact(vpParam,45)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mCatalys.push_back(SpecieId("N2",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=1.6E-11;
	mUncertainty="40%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac45::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac46::Reac46(XmlParameters* vpParam):ChemReact(vpParam,46)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("CO","a3Pi"));
	mProducts.push_back(SpecieId("CO",""));
	mCatalys.push_back(SpecieId("O2",""));
	mEfficiency[SpecieId("CO","")]=1;
	mMainValue=6E-11;
	mUncertainty="30%";
	mSupplInfo="  ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac46::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

/*******
 * O+(2P) reactions
 */


Reac47::Reac47(XmlParameters* vpParam):ChemReact(vpParam,47)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O+","4S"));
	mCatalys.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O+","4S")]=1;
	mMainValue=5.2E-11;
	mUncertainty="50%";
	mSupplInfo="Stephan Meier et al 2003 (and ref therein) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac47::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



Reac48::Reac48(XmlParameters* vpParam):ChemReact(vpParam,48)
{
	mbIsEmit=false;
	mEmitFreqnm=0;

	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O+","4S"));
	mCatalys.push_back(SpecieId("N2",""));
	mEfficiency[SpecieId("O+","4S")]=1;
	mMainValue=1.8E-10;
	mUncertainty="15%";
	mSupplInfo="Stephan Meier et al 2003 (and ref therein) ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac48::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}




Reac49::Reac49(XmlParameters* vpParam):ChemReact(vpParam,49)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("e",""));
	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O",""));
	mEfficiency[SpecieId("O","")]=1;
	mMainValue=1.89E-7;
	mUncertainty="50%";
	mSupplInfo="Supplementary sqrt((vTe)/300), Krishna2009 et ref therein ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac49::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue*sqrt((vTe[i])/300);
	}
	return resu;
}


Reac50::Reac50(XmlParameters* vpParam):ChemReact(vpParam,50)
{
	mbIsEmit=true;
	mEmitFreqnm=732.0;
	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O","2D"));
	mEfficiency[SpecieId("O+","2D")]=1;
	mMainValue=0.160;
	mUncertainty="20%";
	mSupplInfo="NIST - Uncertainty estimated ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac50::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}

Reac51::Reac51(XmlParameters* vpParam):ChemReact(vpParam,51)
{
	mbIsEmit=true;
	mEmitFreqnm=733.0;
	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O","2D"));
	mEfficiency[SpecieId("O+","2D")]=1;
	mMainValue=0.130;
	mUncertainty="20%";
	mSupplInfo="NIST - Uncertainty estimated ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac51::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac52::Reac52(XmlParameters* vpParam):ChemReact(vpParam,52)
{
	mbIsEmit=true;
	mEmitFreqnm=247.0;
	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("O+","4S"));
	mEfficiency[SpecieId("O+","4S")]=1;
	mMainValue=0.073;
	mUncertainty="20%";
	mSupplInfo="NIST - Uncertainty estimated ";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac52::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
		resu[i]=mMainValue;
	return resu;
}


Reac53::Reac53(XmlParameters* vpParam):ChemReact(vpParam,53)
{
	mbIsEmit=false;
	mEmitFreqnm=0;


	mReactant.push_back(SpecieId("CO2",""));
	mReactant.push_back(SpecieId("O+","2P"));
	mProducts.push_back(SpecieId("CO+",""));
	mProducts.push_back(SpecieId("O","3P"));
	mEfficiency[SpecieId("O","3P")]=1;
	mEfficiency[SpecieId("CO+","")]=1;
	mMainValue=1.E-9;
	mUncertainty="80%";
	mSupplInfo="No report of this reaction in the litterature. Estimation of O+ + CO2 for that reaction.";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac53::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}









/***************************************
 *
 * WORK ON THE O++ reactions
 *
 **************************************/


Reac54::Reac54(XmlParameters* vpParam):ChemReact(vpParam,54)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("N2",""));
	mProducts.push_back(SpecieId("N2+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("N2+","")]=1;
	mEfficiency[SpecieId("O+","")]=1;
	mMainValue=1.3E-9;
	mUncertainty="25%";
	mSupplInfo="(Johnsen & Biondi 1978, Gronoff et al. 2007, assumed at 300K)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac54::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac55::Reac55(XmlParameters* vpParam):ChemReact(vpParam,55)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("CO2",""));
	mProducts.push_back(SpecieId("CO2+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("CO2+","")]=1;
	mEfficiency[SpecieId("O+","")]=1;
	mMainValue=2.E-9;
	mUncertainty="25%";
	mSupplInfo="(Fox & Victor 1981, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac55::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



Reac56::Reac56(XmlParameters* vpParam):ChemReact(vpParam,56)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("O+","")]=2;
	mMainValue=1.06E-10;
	mUncertainty="40%";
	mSupplInfo="(Simon et al. 2005, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac56::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac57::Reac57(XmlParameters* vpParam):ChemReact(vpParam,57)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("CO",""));
	mProducts.push_back(SpecieId("CO+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("O+","")]=1;
	mEfficiency[SpecieId("CO+","")]=1;
	mMainValue=1.6E-9;
	mUncertainty="25%";
	mSupplInfo="(Fox & Victor 1981, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac57::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac58::Reac58(XmlParameters* vpParam):ChemReact(vpParam,58)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("He",""));
	mProducts.push_back(SpecieId("He+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("O+","")]=1;
	mEfficiency[SpecieId("He+","")]=1;
	mMainValue=1.1E-10;
	mUncertainty="25%";
	mSupplInfo="(Kimura et al. 1996, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac58::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac59::Reac59(XmlParameters* vpParam):ChemReact(vpParam,59)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("O2",""));
	mProducts.push_back(SpecieId("O2+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("O+","")]=1;
	mEfficiency[SpecieId("O2+","")]=1;
	mMainValue=1.7E-9;
	mUncertainty="25%";
	mSupplInfo="(Howorka et al. (1979), Simon et al. (2005))";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac59::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac60::Reac60(XmlParameters* vpParam):ChemReact(vpParam,60)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("O+","")]=1;
	mMainValue=2.1E-11;
	mUncertainty="25%";
	mSupplInfo="(Nakada & Singer 1968, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac60::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue * 4. / sqrt((vTe[i] + vTi[i]) / 2.);
	}
	return resu;
}


Reac61::Reac61(XmlParameters* vpParam):ChemReact(vpParam,61)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("O++",""));
	mReactant.push_back(SpecieId("H",""));
	mProducts.push_back(SpecieId("O+",""));
	mProducts.push_back(SpecieId("H+",""));
	mEfficiency[SpecieId("O+","")]=1;
	mEfficiency[SpecieId("H+","")]=1;
	mMainValue=1.1E-9;
	mUncertainty="25%";
	mSupplInfo="(Honvault et al. 1995, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac61::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


/***************************************
 *
 * WORK ON THE N2++ reactions
 *
 **************************************/


Reac62::Reac62(XmlParameters* vpParam):ChemReact(vpParam,62)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("N2++",""));
	mReactant.push_back(SpecieId("N2",""));
	mProducts.push_back(SpecieId("N2+",""));
	mEfficiency[SpecieId("N2+","")]=2;
	mMainValue=2.7E-9;
	mUncertainty="25%";
	mSupplInfo="(Lilensten et al. 2005, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac62::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac63::Reac63(XmlParameters* vpParam):ChemReact(vpParam,63)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("N2++",""));
	mReactant.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("NO+",""));
	mProducts.push_back(SpecieId("N+",""));
	mEfficiency[SpecieId("NO+","")]=1;
	mEfficiency[SpecieId("N+","")]=1;
	mMainValue=1.8E-9;
	mUncertainty="50%";
	mSupplInfo="(Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac63::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}


Reac64::Reac64(XmlParameters* vpParam):ChemReact(vpParam,64)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("N2++",""));
	mReactant.push_back(SpecieId("CO2",""));
	mProducts.push_back(SpecieId("CO2+",""));
	mProducts.push_back(SpecieId("N2+",""));
	mEfficiency[SpecieId("N2+","")]=1;
	mEfficiency[SpecieId("CO2+","")]=1;
	mMainValue=3.E-9;
	mUncertainty="50%";
	mSupplInfo="(Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac64::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}

Reac65::Reac65(XmlParameters* vpParam):ChemReact(vpParam,65)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("N2++",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("N2+",""));
	mEfficiency[SpecieId("N2+","")]=1;
	mMainValue=5.8E-7;
	mUncertainty="25%";
	mSupplInfo="(Seiersen et al. 2003, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac65::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue * sqrt(300. / vTe[i]);
	}
	return resu;
}


Reac66::Reac66(XmlParameters* vpParam):ChemReact(vpParam,66)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("N2++",""));
	mProducts.push_back(SpecieId("N+",""));
	mEfficiency[SpecieId("N+","")]=2;
	mMainValue=1/3.;
	mUncertainty="50%";
	mSupplInfo="(Mathur et al. 1995, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac66::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



/***************************************
 *
 * WORK ON THE CO2++ reactions
 *
 **************************************/

Reac67::Reac67(XmlParameters* vpParam):ChemReact(vpParam,67)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("CO2++",""));
	mProducts.push_back(SpecieId("CO+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("CO+","")]=1;
	mEfficiency[SpecieId("O+","")]=1;
	mMainValue=1/4.2;
	mUncertainty="20%";
	mSupplInfo="(Mathur et al. 1995, Mathur et al. 2004, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac67::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}



Reac68::Reac68(XmlParameters* vpParam):ChemReact(vpParam,68)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("CO2++",""));
	mReactant.push_back(SpecieId("e",""));
	mProducts.push_back(SpecieId("CO2+",""));
	mEfficiency[SpecieId("CO2+","")]=1;
	mMainValue=6.2E-7;
	mUncertainty="25%";
	mSupplInfo="(Seiersen et al. 2003, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac68::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue * sqrt(300. / vTe[i]);
	}
	return resu;
}

Reac69::Reac69(XmlParameters* vpParam):ChemReact(vpParam,69)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("CO2++",""));
	mReactant.push_back(SpecieId("CO2",""));
	mProducts.push_back(SpecieId("CO2+",""));
	mEfficiency[SpecieId("CO2+","")]=2;
	mMainValue=2.13E-10;
	mUncertainty="25%";
	mSupplInfo="(Franceschi et al. 2003, Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac69::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue * sqrt(vTn[i] / 300.);
	}
	return resu;
}



Reac70::Reac70(XmlParameters* vpParam):ChemReact(vpParam,70)
{
	mbIsEmit=false;
	mEmitFreqnm=0;
	mReactant.push_back(SpecieId("CO2++",""));
	mReactant.push_back(SpecieId("O",""));
	mProducts.push_back(SpecieId("CO2+",""));
	mProducts.push_back(SpecieId("O+",""));
	mEfficiency[SpecieId("CO2+","")]=1;
	mEfficiency[SpecieId("O+","")]=1;
	mMainValue=2.E-9;
	mUncertainty="50%";
	mSupplInfo="(Gronoff et al. 2007)";
	ReadParameters();// should be here!
}

ublas::vector<double> Reac70::GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)
{
	assert(vTe.size()==vTi.size());
	assert(vTe.size()==vTn.size());
	ublas::vector<double> resu(vTn.size());
	for(unsigned i=0;i<resu.size();++i)
	{
		resu[i]=mMainValue;
	}
	return resu;
}






