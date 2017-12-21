/**
 * \file chem.cpp
 * \brief Implements the chemistry class
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: chem.cpp 1270 2011-06-06 14:05:50Z gronoff $
 */


#include "chem.hpp"
using namespace std;
typedef std::deque< boost::shared_ptr<CalcDensPhEq> >::iterator CalcDensIt;
typedef std::deque<Specie* >::iterator SpecIt;


Chem::Chem(XmlParameters* vpParam,
		Planete* vpPlanet,
		ublas::vector<double> vAltGridKm,
		std::deque<Specie*> vNeutralDenscm_3,
		std::deque<Specie*> vTotProdcm_3s_1,
		std::deque<Specie*> vPhotProdcm_3s_1,
		std::deque<Specie*> vElecProdcm_3s_1,
		std::deque<Specie*> vProtProdcm_3s_1,
		std::deque<Specie*> vCosmoProdcm_3s_1,
		ublas::vector<double> vTempNeutreK,
		ublas::vector<double> vTempElecK,
		ublas::vector<double> vTempIonK,
		ublas::vector<double> vElecDenscm_3,
		double vElecDensMult):mpPlanet(vpPlanet),mpParam(vpParam),mAltGridKm(vAltGridKm),mTotProdcm_3s_1(vTotProdcm_3s_1),mPhotProdcm_3s_1(vPhotProdcm_3s_1),mElecProdcm_3s_1(vElecProdcm_3s_1),mProtProdcm_3s_1(vProtProdcm_3s_1),mCosmoProdcm_3s_1(vCosmoProdcm_3s_1),mTempNeutreK(vTempNeutreK),mTempElecK(vTempElecK),mTempIonK(vTempIonK),mElecDenscm_3(vElecDenscm_3),mElecDensMult(vElecDensMult)
{

	InitDensities(vNeutralDenscm_3);// appel readDensities
	InitChemReact();
	InitCalcDensPhEq();

}


double Chem::GetMultiplicator(std::string vName)
{
	double mult=1.;
	string name=StrReplace(vName,"+","_PLUS");
	if(mpParam->Exists("/aero_main/chem/specie_multiplicator/"+name))
	{
		mpParam->GetValue("/aero_main/chem/specie_multiplicator/"+name,mult);
		if(mpParam->NbParams("/aero_main/chem/specie_multiplicator/"+name,"mult_elec")>0)
		{
			mult*=mElecDensMult;
			Log::mD<<"Electron proportionality activated"<<endl;
		}
					
		Log::mD<<"Multiplicator activated for "<<name<<" : "<<mult<<endl;
	}
	return mult;

}


void Chem::InitDensities(std::deque<Specie*> vNeutralDenscm_3)
{
	assert(vNeutralDenscm_3.size()>0);

	if(!mpParam->Exists("/aero_main/chem/del_model_species"))
	{
		for(SpecIt it=vNeutralDenscm_3.begin();it!=vNeutralDenscm_3.end();++it)
		{
			mChemDensComputcm_3.push_back(*it);
		}
	}else
	{
		deque<string> del=mpParam->GetParamDeQue("/aero_main/chem/del_model_species");
		for(unsigned i=0;i<del.size();++i)
		{
			del[i]=StrReplace(trim(del[i]),"_PLUS","+");
		}

		for(SpecIt it=vNeutralDenscm_3.begin();it!=vNeutralDenscm_3.end();++it)
		{
			unsigned rpos=0;
			if(!PosInVector(del,(*it)->mName,rpos))
			{
				mChemDensComputcm_3.push_back(*it);
			}
		}

	}
	InitNewDensities();
	ReadDensities();
}

void Chem::ReadDensities()
{


	string spfile_name=mpParam->GetSubFileName("/aero_main/atmosphere/neutral/species_file");
	if(mpParam->Exists("/aero_main/chem/add_species"))
	{

		deque<string> add_names=mpParam->GetParamDeQue("/aero_main/chem/add_species/species");
		// check if the new species are not in the old ones 
		if(add_names.size()>0)
		{
			ublas::vector<double> alt;
			ublas::vector<double> tmpdens;
			mpParam->ExistsOrDie("/aero_main/chem/add_species/alt","If you want additionnal species, you should provide an atmosphere altitude model");
			mpParam->Get1DArray("/aero_main/chem/add_species/alt",alt);
			for(deque<string>::iterator st=add_names.begin();st!=add_names.end();++st)
			{
				string name=StrReplace(trim(*st),"_PLUS","+");

				double multi=GetMultiplicator(name);


				if(SpecieUtils::PosOfSp(name,mChemDensComputcm_3)>-1)
				{
					unsigned spposition=SpecieUtils::PosOfSp(name,mChemDensComputcm_3);
					//Log::SetPriority(Log::WARNING);
					//Log::mL<<"Species defined in Chem::InitNewDensities"<<" The species "+*st+" has already been defined for your atmosphere. Its density will be overloaded"<<endl;;
					//Log::SetPriority(Log::DEBUGG);
					if(mChemDensComputcm_3[spposition]->mTotDensitycm_3.size()>0)
					{
						Log::mD<<"Density for your species : "<<mChemDensComputcm_3[spposition]->mTotDensitycm_3<<endl;
						Error err("Error, species defined","Chem::ReadDensities"," The species "+*st+" has already been defined for your atmosphere.");
						throw err;
					}
					mpParam->Get1DArray("/aero_main/chem/add_species/"+*st,tmpdens);

					mChemDensComputcm_3[spposition]->mTotDensitycm_3=MathFunction::IntLog(alt,tmpdens,mAltGridKm)*multi;
				}else
				{
					Specie* tmp=new Specie(name,spfile_name);
					mpParam->Get1DArray("/aero_main/chem/add_species/"+*st,tmpdens);

					if(mpParam->KeyExists( "/aero_main/chem/add_species/"+*st,"DeleteZero"))
					{
						Log::mD<<"The zero values in your species density are set to 1E-40"<<endl;
						for(unsigned i=0; i< tmpdens.size(); ++i)
						{
							if(tmpdens[i] < 1E-40)
								tmpdens[i] = 1E-40;
						}
					}




					tmp->mTotDensitycm_3=MathFunction::IntLog(alt,tmpdens,mAltGridKm)*multi;
					mDefinedHere.push_back(tmp);
					mChemDensComputcm_3.push_back(tmp);
				}
			}
		}
	}

	if(mpParam->Exists("/aero_main/chem/overload_species"))
	{

		deque<string> overload_names=mpParam->GetParamDeQue("/aero_main/chem/overload_species/species");
		// check if the new species are not in the old ones 
		if(overload_names.size()>0)
		{
			ublas::vector<double> alt;
			mpParam->ExistsOrDie("/aero_main/chem/overload_species/alt","If you want overloaded species, you should provide an atmosphere altitude model");
			mpParam->Get1DArray("/aero_main/chem/overload_species/alt",alt);
			for(deque<string>::iterator st=overload_names.begin();st!=overload_names.end();++st)
			{
				string name=StrReplace(trim(*st),"_PLUS","+");
				int pos=SpecieUtils::PosOfSp(name,mChemDensComputcm_3);
				if(pos<0)
				{
					Error err("Error, species undefined","Chem::overload_sppecie"," The species "+name+" has never  been defined for your atmosphere.");
					throw err;
				}

				ublas::vector<double> tmpdens;
				mpParam->Get1DArray("/aero_main/chem/overload_species/"+*st,tmpdens);
				mpParam->ExistsOrDie("/aero_main/chem/overload_species/"+*st+"/alt_min","Impossible to get alt_min for your overloaded species");
				mpParam->ExistsOrDie("/aero_main/chem/overload_species/"+*st+"/alt_max","Impossible to get alt_max for your overloaded species");
				double altmin,altmax;
				mpParam->GetValue("/aero_main/chem/overload_species/"+*st+"/alt_min",altmin);
				mpParam->GetValue("/aero_main/chem/overload_species/"+*st+"/alt_max",altmax);
				tmpdens=MathFunction::IntLog(alt,tmpdens,mAltGridKm);
				assert(mAltGridKm.size()==mChemDensComputcm_3[pos]->mTotDensitycm_3.size());
				double multi=GetMultiplicator(name);
				for(unsigned i=0;i<mAltGridKm.size();++i)
				{
					if(mAltGridKm[i]>=altmin && mAltGridKm[i]<=altmax)
					{
						mChemDensComputcm_3[pos]->mTotDensitycm_3[i]=tmpdens[i]*multi;
					}

				}
			}
		}
	}
}

void Chem::InitNewDensities()
{
	if(mpParam->Exists("/aero_main/chem/add_model_species"))
	{
		string spfile_name=mpParam->GetSubFileName("/aero_main/atmosphere/neutral/species_file");
		deque<string> sp_names=mpParam->GetParamDeQue("/aero_main/chem/add_model_species");

		for(unsigned i=0;i<sp_names.size();++i)
		{
			sp_names.at(i)=trim(sp_names.at(i));
		}

		int type=-1;
		if(mpParam->NbParams("/aero_main/chem/add_model_species","type")!=0)
		{
			mpParam->GetNKey("/aero_main/chem/add_model_species","type",type);
		}else
		{
			mpParam->GetNKey("/aero_main/atmosphere/neutral/model","type",type);
		}
		map< string, ublas::vector<double> > my_atmo;
		my_atmo=mpPlanet->AtmoModel(mAltGridKm,sp_names,type);
		for(unsigned i=0;i<sp_names.size();++i)
		{
			double multi=GetMultiplicator(sp_names.at(i));
			if(SpecieUtils::PosOfSp(sp_names.at(i),mChemDensComputcm_3)>-1)
			{


				unsigned spposition=SpecieUtils::PosOfSp(sp_names.at(i),mChemDensComputcm_3);
				//Log::SetPriority(Log::WARNING);
				//Log::mL<<"Species defined in Chem::InitNewDensities"<<" The species "+*st+" has already been defined for your atmosphere. Its density will be overloaded"<<endl;;
				//Log::SetPriority(Log::DEBUGG);
				if(mChemDensComputcm_3[spposition]->mTotDensitycm_3.size()>0)
				{
					Log::mD<<"Density for your species : "<<mChemDensComputcm_3[spposition]->mName<<mChemDensComputcm_3[spposition]->mTotDensitycm_3<<endl;
					/*
					   ublas::vector<double> resu;
					   SpecieUtils::GetSpecieDensity("O2+","",mTotProdcm_3s_1,resu);
					   Log::mL<<" Densite O2+ dans tot : "<<resu<<endl;
					   resu.clear();
					   SpecieUtils::GetSpecieDensity("O2+","",mPhotProdcm_3s_1,resu);
					   Log::mL<<" Densite O2+ dans phot : "<<resu<<endl;
					   resu.clear();
					   SpecieUtils::GetSpecieDensity("O2+","",mElecProdcm_3s_1,resu);
					   Log::mL<<" Densite O2+ dans elec : "<<resu<<endl;
					   resu.clear();
					   SpecieUtils::GetSpecieDensity("O2+","",mProtProdcm_3s_1,resu);
					   Log::mL<<" Densite O2+ dans prot : "<<resu<<endl;
					   resu.clear();
					   SpecieUtils::GetSpecieDensity("O2+","",mCosmoProdcm_3s_1,resu);
					   Log::mL<<" Densite O2+ dans cosmo : "<<resu<<endl;
					   resu.clear();*/





					Error err("Add_model_species : species defined","Chem::InitNewDensities"," The species "+sp_names.at(i)+" has already been defined for your atmosphere.");
					throw err;
				}
				mChemDensComputcm_3[spposition]->mTotDensitycm_3=my_atmo[sp_names.at(i)]*multi;
			}else
			{
				Specie* tmp=new Specie(sp_names.at(i),spfile_name);
				tmp->mTotDensitycm_3=my_atmo[sp_names.at(i)]*multi;
				mDefinedHere.push_back(tmp);
				mChemDensComputcm_3.push_back(tmp);
			}
		}
	}
	if(mpParam->Exists("/aero_main/chem/add_model_ion"))
	{
		string spfile_name=mpParam->GetSubFileName("/aero_main/atmosphere/neutral/species_file");
		deque<string> sp_names=mpParam->GetParamDeQue("/aero_main/chem/add_model_ion");
		
		for(unsigned i=0;i<sp_names.size();++i)
		{
			sp_names.at(i)=StrReplace(trim(sp_names.at(i)),"_PLUS","+");
		}
		
		int type=0;
		if(mpParam->NbParams("/aero_main/chem/add_model_ion","type")!=0)
		{
			mpParam->GetNKey("/aero_main/chem/add_model_ion","type",type);
		}
		map< string, ublas::vector<double> > my_atmo;
		my_atmo=mpPlanet->AtmoModel(mAltGridKm,sp_names,type); // the atmo model stores both ions and neutral species
		for(unsigned i=0;i<sp_names.size();++i)
		{
			Specie* tmp=new Specie(sp_names.at(i),spfile_name);
			double multi=GetMultiplicator(sp_names.at(i));
			if(SpecieUtils::PosOfSp(sp_names.at(i),mChemDensComputcm_3)>-1)
			{
				Error err("Add_model_ion: species defined","Chem::InitNewDensities"," The species "+sp_names.at(i)+" has already been defined for your atmosphere.");
				throw err;
			}
			tmp->mTotDensitycm_3=my_atmo[sp_names.at(i)]*multi;
			mDefinedHere.push_back(tmp);
			mChemDensComputcm_3.push_back(tmp);
		}
	}



}



void Chem::InitChemReact()
{
	// here, the list of chemical reactions
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac0(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac1(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac2(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac3(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac4(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac5(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac6(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac7(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac8(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac9(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac10(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac11(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac12(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac13(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac14(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac15(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac16(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac17(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac18(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac19(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac20(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac21(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac22(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac23(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac24(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac25(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac26(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac27(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac28(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac29(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac30(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac31(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac32(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac33(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac34(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac35(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac36(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac37(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac38(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac39(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac40(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac41(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac42(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac43(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac44(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac45(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac46(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac47(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac48(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac49(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac50(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac51(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac52(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac53(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac54(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac55(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac56(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac57(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac58(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac59(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac60(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac61(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac62(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac63(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac64(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac65(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac66(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac67(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac68(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac69(mpParam)));
	mChemList.push_back(boost::shared_ptr<ChemReact>(new Reac70(mpParam)));


	// 
}

void Chem::InitCalcDensPhEq()
{
	// Here, the list of chemical eq
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensN2A3S()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensO1S()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensN2D()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensO1D()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensCOa3Pi()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensOplus2P()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensOpp()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensCO2pp()));
	mCalcPhEqList.push_back(boost::shared_ptr<CalcDensPhEq>(new CalcDensN2pp()));

	for(CalcDensIt it=mCalcPhEqList.begin();it!=mCalcPhEqList.end();++it)
	{
		(*it)->Init(mChemDensComputcm_3,
				mTotProdcm_3s_1,
				mPhotProdcm_3s_1,
				mElecProdcm_3s_1,
				mProtProdcm_3s_1,
				mCosmoProdcm_3s_1,
				mTempNeutreK,
				mTempElecK,
				mTempIonK,
				&mElecDenscm_3,
				mChemList
			   );

	}
}


Chem::~Chem()
{
	/// Delete the species defined inside the present class
	for(unsigned i = 0;i< mDefinedHere.size();++i)
	{
		delete mDefinedHere[i];
	}
	/* Chuis con, ce sont des smart pointers
	for(unsigned i=0;i<mCalcPhEqList.size();++i)
	{
		delete mCalcPhEqList[i];
	}

	for(unsigned i=0;i<mChemList.size();++i)
	{
		delete mChemList[i];
	}*/
}






void Chem::SelectedPrintDensities(std::deque<SpecieId>& rStates, ublas::vector<double>  vAltKmGrid, std::string rFilename)
{
	if(FileExists(rFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<rFilename<<endl;
	}
	ofstream of(rFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Species Densities "<<endl;
	of<<"# Density in cm-3"<<endl;
	of<<"# Altitude in km"<<endl;
	deque<Specie*>::iterator sp;


	deque< ublas::vector<double> > the_densities;
	ublas::vector<double> zero_densities(vAltKmGrid.size());
	zero_densities.clear();

	// Firstly, we print the different names.
	// While doing that, we find the production vectors
	// and we store their pointeur in the_densities
	// this vector will be used to print

	for(unsigned i = 0;i < rStates.size(); ++i)
	{
		of<<"# "<<rStates[i].StandardName()<<"    ";
		ublas::matrix<double> tmp;
		std::deque<std::string> rsubtmp;
		std::string warn;
		cout<<"The density for"<<rStates[i].StandardName()<<endl;
		the_densities.push_back(GetDensity(rStates[i],tmp, rsubtmp, warn));
	//	the_densities.push_back(GetDens(rStates[i].mName, rStates[i].mState)); 
		cout<<"Finished that density"<<endl;
		if(the_densities[i].size() == 0)
		{
			the_densities[i] = zero_densities;
		}
		of<<endl;
	}
	for(unsigned i = 0; i < vAltKmGrid.size(); ++i)
	{
		of<<vAltKmGrid[i]<<"\t";
		for(unsigned j=0;j<the_densities.size();++j)
		{
			of<<(the_densities[j])[i]<<"\t";
		}
		of<<endl;
	}
	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}

ublas::vector<double> Chem::GetDensity(SpecieId vId,ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string& rWarnings)
{
	ublas::vector<double> resu;
	unsigned i=0;
	bool found = false;
	while( (!found) && (i < mCalcPhEqList.size()))
	{
		if(vId == mCalcPhEqList[i]->GetId())
		{
			found = true;
			Log::mD<<"The density has been found"<<endl;
			resu=mCalcPhEqList[i]->GetDensity(rSubResults,rSubResultsNames,rWarnings);
			mDensityWarnings.insert(make_pair((mCalcPhEqList[i]->GetId()).StandardName(),rWarnings));
			PutDens(vId.mName, vId.mState, resu);
			cout<<"on continue"<<endl;	
		}else
		{
			//Log::mL<<"Not matching : "<<vId.StandardName()<<" _ "<<(mCalcPhEqList[i]->GetId()).StandardName()<<endl;
			++i;
		}
	}
	cout<<"on continue2"<<endl;	
	if(!found)
	{
			Log::mD<<"The density has NOT been found"<<endl;
			Log::mD<<"Size of the list : "<<mCalcPhEqList.size()<<endl;
			rWarnings+="\nThe density has NOT been found\n";
	}
	cout<<"On donne le resultat"<<endl;
	return resu;
}


void Chem::PrintChemList(std::string vFilename)
{
	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	for(unsigned i=0;i<mChemList.size();++i)
	{
		of<<mChemList[i]->GetInfo()<<endl;
	}
	of.close();

}


bool Chem::CheckReacId()
{
	bool resu=true;

	for(unsigned i=0;i<mChemList.size();++i)
	{
		bool bon= (i==mChemList[i]->GetId());
		if(!bon)
		{
			resu=false;
			//Log::SetPriority(Log::ERROR);
			Log::mE<<" Error in your chemical reaction list. Position "<<i<<endl;
		}
	}
	return resu;
}



void Chem::PrintDensWarnings(std::string vFilename)
{
	if(FileExists(vFilename))
	{

		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	for(std::map<string,string>::iterator it=mDensityWarnings.begin();it!=mDensityWarnings.end();++it)
	{
		of<<it->first<<"\t : \t"<<it->second<<endl;
	}
	of.close();
}


void Chem::PrintChemDensModels(std::string vFilename)
{
	
	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Atmosphere considered "<<endl;
	of<<"# Density in cm3"<<endl;
	of<<"# Altitude in km"<<endl;

	deque<Specie*>::iterator it;


	unsigned i=0;

	of<<"# e-"<<endl;
	assert(mElecDenscm_3.size()==mAltGridKm.size());
	for(it=mChemDensComputcm_3.begin();it!=mChemDensComputcm_3.end();++it)
	{
	//	cout<<"Specie "<<i<<" => "<<(*it)->mName<<endl;
		of<<"# "<<(*it)->mName<<endl;

	//	cout<<"Sizes : "<<(*it)->mTotDensitycm_3.size()<<"  alt_grid size : "<<mAltGridKm.size()<<endl;
		assert((*it)->mTotDensitycm_3.size()==mAltGridKm.size());
		++i;
	}


	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		of<<mElecDenscm_3[i]<<"\t";

		for(unsigned j=0;j<mChemDensComputcm_3.size();++j)
		{
			of<<mChemDensComputcm_3[j]->mTotDensitycm_3[i]<<"\t";
		}
		of<<endl;
	}
}


ublas::vector<double> Chem::GetDens(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieDensity(vSpecieName,vSpecieState,mChemDensComputcm_3,resu);
	return resu;
}

void Chem::PutDens(std::string vSpecieName,std::string vSpecieState,ublas::vector<double> vDensitycm_3)
{
	if(SpecieUtils::PutSpecieDensity(trim(vSpecieName),trim(vSpecieState),mChemDensComputcm_3,vDensitycm_3))
	{
		Log::mD<< trim(vSpecieName)<<"("<<trim(vSpecieState)<<") Density added correctly"<<endl;
	}
}


ublas::vector<double> Chem::GetPhotProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mPhotProdcm_3s_1,resu);
	return resu;
}

ublas::vector<double> Chem::GetElecProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mElecProdcm_3s_1,resu);
	return resu;
}

ublas::vector<double> Chem::GetProtProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mProtProdcm_3s_1,resu);
	return resu;
}

ublas::vector<double> Chem::GetCosmoProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mCosmoProdcm_3s_1,resu);
	return resu;
}


ublas::vector<double> Chem::GetTotProd(std::string vSpecieName,std::string vSpecieState)
{
	ublas::vector<double> resu;
	SpecieUtils::GetSpecieProduction(vSpecieName,vSpecieState,mTotProdcm_3s_1,resu);
	return resu;
}
