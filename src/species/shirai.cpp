/**
 * \file shirai.cpp
 * \brief Implements the shirai class, used to translate the parameters of Shirai into actual cross sections 
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: shirai.cpp 996 2010-05-14 22:06:46Z gronoff $
 */


#include "shirai.hpp"
using namespace std;


Shirai::Shirai(XmlParameters* vpParams,TiXmlNode* vNode,double vThresholdeV)
{
//	Log::SetPriority(Log::DEBUGG);
//	Log::mL<<"Load Shirai cross section"<<endl;
	mEquationType=110; // To raise an error if not defined later
	mArticleId=""; // To raise an error if not defined later
	mpParams=vpParams;
	mNode=vNode;
	mThreshkeV=vThresholdeV/1E3;
	mSigmacm2=1E-16;
	mErkeV=1.361E-2;
//	Log::mL<<"Load Crs"<<endl;
	if(! LoadCrs() )
	{
		Error err("Shirai","Impossible to load the parameters","Your parameters in your shirai-kind cross section are wrong. Please check your cross section!");
		throw err;
	}
//	Log::mL<<"Load crs finished"<<endl;
}

ublas::vector<double> Shirai::GetCrs(ublas::vector<double> vEeV)
{
	//	Log::SetPriority(Log::DEBUGG);
	//	Log::mL<<"Get Crs"<<endl;
	assert((*(vEeV.begin()))> (*(vEeV.end()-1))); // vE decreasing
	ublas::vector<double> resultat;

	if(mArticleId=="CH4")
	{
		//	Log::mL<<"CH4 type"<<endl;
		resultat=ReturnCH4(vEeV);
	}else if(mArticleId=="CO2")
	{
		//	Log::mL<<"CO2 type"<<endl;
		resultat=ReturnCO2(vEeV);
	}else if(mArticleId=="N2")
	{
		//	Log::mL<<"N2 type"<<endl;
		resultat=ReturnN2(vEeV);
	}else
	{

		Error err("Get Crs Shirai","Error in the definition of the article id","The article Id is not defined for your cross section, it should be CO2, N2, or CH4 (corresponding to main species in the title of the article).");
		throw err;
	}

	size_t i=vEeV.size()-1;
	if(vEeV[0]/1E3 < mThreshkeV)
	{
		resultat.clear();
		return resultat;
	}
	while(vEeV[i]/1E3 < mThreshkeV && i > 0) // 1E3: eV to keV
	{
		resultat[i]=0;
		--i;
	}
	return resultat;
}

bool Shirai::LoadCrs()
{

	mpParams->ExistsOrDie(mNode,"//params","You should define params for Shirai");
	mpParams->ExistsOrDie(mNode,"//Emin","You should define Emin for Shirai");
	mpParams->ExistsOrDie(mNode,"//Emax","You should define Emax for Shirai");
	mpParams->ExistsOrDie(mNode,"//Equation","You should define Equation for Shirai");

	mpParams->Get1DArray(mNode,"//params",mAvalues);
	mpParams->GetValue(mNode,"//Emin",mMinEnerkeV);
	mMinEnerkeV/=1E3;//eV -> keV
	mpParams->GetValue(mNode,"//Emax",mMaxEnerkeV);
	mMaxEnerkeV/=1E3; //eV -> keV

	mpParams->GetNKey(mNode,"//Equation","type",mEquationType);
	mArticleId=trim(mpParams->GetKey(mNode,"//Equation","article_id"));

	return true;
}

ublas::vector<double> Shirai::ReturnCH4(ublas::vector<double> vE)
{
	ublas::vector<double> vE2(vE.size());
	for(size_t i=0;i<vE.size();++i)
		vE2[i]=vE[i]-mThreshkeV*1E3;// E1 in Shirai papers
	switch(mEquationType)
	{
		case 1: return S1CH4(vE2);
			break;
		case 2: return S2CH4(vE2);
			break;
		case 3: return S3CH4(vE2);
			break;
		case 4: return S4CH4(vE2);
			break;
		case 5: return S5CH4(vE2);
			break;
		case 6: return S6CH4(vE2);
			break;
		case 7: return S7CH4(vE2);
			break;
		case 8: return S8CH4(vE2);
			break;
		case 9: return S9CH4(vE2);
			break;
		case 10: return S10CH4(vE2);
			 break;
		case 11: return S11CH4(vE2);
			 break;
		case 12: return S12CH4(vE2);
			 break;
		case 13: return S13CH4(vE2);
			 break;
		case 14: return S14CH4(vE);
			 break;

		default:
			 Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
			 throw err;
	}
	// Never reaches this point. Raise an error!
	Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
	throw err;
	return vE;
}

ublas::vector<double> Shirai::ReturnCO2(ublas::vector<double> vE)
{
	ublas::vector<double> vE2(vE.size());
	for(size_t i=0;i<vE.size();++i)
		vE2[i]=vE[i]-mThreshkeV*1E3;// E1 in Shirai papers
	switch(mEquationType)
	{
		case 1: return S1CO2(vE2);
			break;
		case 2: return S2CO2(vE2);
			break;
		case 3: return S3CO2(vE2);
			break;
		case 4: return S4CO2(vE2);
			break;
		case 5: return S5CO2(vE2);
			break;
		case 6: return S6CO2(vE2);
			break;
		case 7: return S7CO2(vE2);
			break;
		case 8: return S8CO2(vE2);
			break;
		case 9: return S9CO2(vE2);
			break;
		case 10: return S10CO2(vE);
			 break;

		default:
			 Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
			 throw err;
	}
	// Never reaches this point. Raise an error!
	Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
	throw err;
	return vE;
}

ublas::vector<double> Shirai::ReturnN2(ublas::vector<double> vE)
{
	ublas::vector<double> vE2(vE.size());
	for(size_t i=0;i<vE.size();++i)
		vE2[i]=vE[i]-mThreshkeV*1E3;// E1 in Shirai papers
	switch(mEquationType)
	{
		case 1: return S1N2(vE2);
			break;
		case 2: return S2N2(vE2);
			break;
		case 3: return S3N2(vE2);
			break;
		case 4: return S4N2(vE2);
			break;
		case 5: return S5N2(vE2);
			break;
		case 6: return S6N2(vE2);
			break;
		case 7: return S7N2(vE);
			break;

		default:
			Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
			throw err;
	}
	// Never reaches this point. Raise an error!
	Error err("Return Shirai","Error in the definition of the equation type","The equation type is not defined for your cross section");
	throw err;
	return vE;//
}

ublas::vector<double> Shirai::F1(ublas::vector<double> vE,double vC1,double vC2)
{
	ublas::vector<double> resu(vE.size());
	resu.clear();

	for(unsigned i=0;i<vE.size();++i)
	{
		resu[i]=mSigmacm2*vC1*pow((vE[i]*1E-3/mErkeV),vC2);
	}
	return resu;
}


ublas::vector<double> Shirai::F2(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4)
{
	ublas::vector<double> resu(vE.size());
	resu.clear();
	ublas::vector<double> tmp=F1(vE,vC1,vC2);

	for(unsigned i=0;i<vE.size();++i)
	{
		resu[i]=tmp[i]/ (1+pow((vE[i]/(1E3*vC3)),vC2+vC4));
	}
	return resu;
}
ublas::vector<double> Shirai::F3(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4,double vC5,double vC6)
{
	ublas::vector<double> resu(vE.size());
	resu.clear();
	ublas::vector<double> tmp=F1(vE,vC1,vC2);

	for(unsigned i=0;i<vE.size();++i)
	{
		resu[i]=tmp[i]/ (1+pow((vE[i]/(1E3*vC3)),vC2+vC4)+pow((vE[i]/(1E3*vC5)),vC2+vC6));
	}
	return resu;
}

ublas::vector<double> Shirai::F4(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4)
{
	ublas::vector<double> resu(vE.size());
	resu.clear();

	for(unsigned i=0;i<vE.size();++i)
	{
		resu[i]=mSigmacm2*vC1*( (log(vE[i]*1E-3/mErkeV) + vC2)/(mThreshkeV*vE[i]*1E-3* (1+vC3/pow( (vE[i]*1E-3-mThreshkeV) ,vC4)) ) );
	}
	return resu;
}
ublas::vector<double> Shirai::S1CH4(ublas::vector<double> vE)
{
	return F1(vE,mAvalues[0],mAvalues[1]);
}
ublas::vector<double> Shirai::S2CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S3CH4(ublas::vector<double> vE)
{

	return F1(vE,mAvalues[0],mAvalues[1])+ F2(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S4CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+mAvalues[4]*F2(vE/mAvalues[5],mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S5CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[3]);
}
ublas::vector<double> Shirai::S6CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7]);
}
ublas::vector<double> Shirai::S7CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[3])+F2(vE,mAvalues[7],mAvalues[8],mAvalues[9],mAvalues[3]);
}
ublas::vector<double> Shirai::S8CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7])+F2(vE,mAvalues[8],mAvalues[9],mAvalues[10],mAvalues[11]);
}
ublas::vector<double> Shirai::S9CH4(ublas::vector<double> vE)
{

	return F3(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S10CH4(ublas::vector<double> vE)
{

	return F1(vE,mAvalues[0],mAvalues[1])+ F3(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7]);
}
ublas::vector<double> Shirai::S11CH4(ublas::vector<double> vE)
{

 	return F3(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5])+mAvalues[6]*F3(vE/mAvalues[7],mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S12CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F3(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[3]);
}
ublas::vector<double> Shirai::S13CH4(ublas::vector<double> vE)
{

	return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F3(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[9])+F2(vE,mAvalues[10],mAvalues[11],mAvalues[12],mAvalues[9]);
}
ublas::vector<double> Shirai::S14CH4(ublas::vector<double> vE)
{

	return F4(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S1CO2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S2CO2(ublas::vector<double> vE)
{

		return F1(vE,mAvalues[0],mAvalues[1])+ F2(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S3CO2(ublas::vector<double> vE)
{ // 14 may 2010 correction from http://www-jt60.naka.jaea.go.jp/english/JEAMDL/code/00001_errata.html
		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7]);

}
ublas::vector<double> Shirai::S4CO2(ublas::vector<double> vE)
{

		return F1(vE,mAvalues[0],mAvalues[1])+ F2(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5])+ F2(vE,mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[9]);
}
ublas::vector<double> Shirai::S5CO2(ublas::vector<double> vE)
{

		return F3(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S6CO2(ublas::vector<double> vE)
{

		return F1(vE,mAvalues[0],mAvalues[1])+ F3(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7]);
}
ublas::vector<double> Shirai::S7CO2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F3(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[3]);
}
ublas::vector<double> Shirai::S8CO2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F3(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[9]);
}
ublas::vector<double> Shirai::S9CO2(ublas::vector<double> vE)
{

		return F1(vE,mAvalues[0],mAvalues[1])+ F2(vE,mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5])+F3(vE,mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[9],mAvalues[10],mAvalues[11]);
}
ublas::vector<double> Shirai::S10CO2(ublas::vector<double> vE)
{

	return F4(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S1N2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}
ublas::vector<double> Shirai::S2N2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[3]);
}
ublas::vector<double> Shirai::S3N2(ublas::vector<double> vE)
{
		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7]);

}
ublas::vector<double> Shirai::S4N2(ublas::vector<double> vE)
{

		return F2(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3])+F2(vE,mAvalues[4],mAvalues[5],mAvalues[6],mAvalues[7])+F2(vE,mAvalues[8],mAvalues[9],mAvalues[10],mAvalues[11]);
}
ublas::vector<double> Shirai::S5N2(ublas::vector<double> vE)
{

		return F3(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5]);
}
ublas::vector<double> Shirai::S6N2(ublas::vector<double> vE)
{

		return  F3(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3],mAvalues[4],mAvalues[5])+F2(vE,mAvalues[6],mAvalues[7],mAvalues[8],mAvalues[9]);
}
ublas::vector<double> Shirai::S7N2(ublas::vector<double> vE)
{

	return F4(vE,mAvalues[0],mAvalues[1],mAvalues[2],mAvalues[3]);
}


