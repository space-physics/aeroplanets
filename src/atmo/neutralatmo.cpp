/**
 * \file neutralatmo.cpp
 * \brief implements the neutral atmo classes
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: neutralatmo.cpp 1550 2012-08-22 19:29:22Z gronoff $
 */


#include "neutralatmo.hpp"
using namespace std;


NeutralAtmo::NeutralAtmo(XmlParameters* pParam,ublas::vector<double> vAltGridKm,Planete* pPlanet,ublas::vector<double>* pPhGrideV,ublas::vector<double>* pElGrideV, ublas::vector<double>* pElGridEngddeV,ublas::vector<double>* pPrGrideV)
{
	mpProtonGrideV=pPrGrideV;
	mpPhotonGrideV=pPhGrideV;
	mpElectronGrideV=pElGrideV;
	mpElectronGridWidtheV=pElGridEngddeV;
	mpParameter=pParam;
	mAltGridKm=vAltGridKm;
	mpPlanet=pPlanet;
	ReadParameters();
}


NeutralAtmo::~NeutralAtmo()
{
	// We delete the elements of the vector
	for(deque<Specie*>::iterator it=mAtmoSpecies.begin();it!=mAtmoSpecies.end();++it)
	{
		delete *it;
	}

}


ublas::vector<double> NeutralAtmo::Temperature()
{
	mpParameter->ExistsOrDie("/aero_main/atmosphere/neutral/temperature","You have to define your temperature model");

	if(mpParameter->Exists("/aero_main/atmosphere/neutral/temperature/use_model"))
	{
		if(mpPlanet->mIsTemperatureDefined)
		{
			return mpPlanet->mTemperatureModelGridK;
		}else
		{
			Log::mE<<"Temperature is not defined in your model"<<endl;
			Log::mE<<"Please improve the model or add a data model by adding <alt> and <T>, and their respective values inside"<<endl;
			Error err("NeutralAtmo::Temperature()","Is temperature defined","Please improve the model or add a data model by adding <alt> and <T>, and their respective values inside");
			throw err;
		
		}

	}

	mpParameter->ExistsOrDie("/aero_main/atmosphere/neutral/temperature/alt","You have to define the neutral temperature altitudes");
	mpParameter->ExistsOrDie("/aero_main/atmosphere/neutral/temperature/T","You have to define the neutral temperature parameter");

	ublas::vector<double> lalt;
	mpParameter->Get1DArray("/aero_main/atmosphere/neutral/temperature/alt",lalt);
	ublas::vector<double> ltemp;
	mpParameter->Get1DArray("/aero_main/atmosphere/neutral/temperature/T",ltemp);
//	ublas::vector<double> temp=MathFunction::IntLin(lalt,ltemp,mAltGridKm);
	
//	return temp;
//
//
//	WARNING : THE MULTIPLICATOR IS SET UP IN THE MAIN ATMOSPHERE FUNCTION
//
	return MathFunction::IntLog(lalt,ltemp,mAltGridKm);
}



void NeutralAtmo::ReadParameters()
{
	//	/atmosphere/neutral/neutral_species <- string vector
	// 	/atmosphere/neutral/model ()type <- int
	// 	if type == 0
	// 	/atmosphere/neutral/data
	// 				-> /alt for the altitude grid
	// 				-> /specie_name for the specie
	// 			then, we can do a linear interpolation
	//

	if(!mpParameter->Exists("/aero_main/atmosphere/neutral/neutral_species"))
	{
		Log::mE<<"Error /atmosphere/neutral/neutral_species does not exists"<<endl;
		Error err("NeutralAtmo::ReadParameters","/atmosphere/neutral/neutral_species does not exists","Solution : implement that option in your main xml file");
		throw err;
	}	

	deque<string> sp_names=mpParameter->GetParamDeQue("/aero_main/atmosphere/neutral/neutral_species");

	// ...... Init neutral species
	//
	if(sp_names.size()<1)
	{
		Log::mE<<"Error /atmosphere/neutral/neutral_species is null. Please put your species here"<<endl;
		Error err("NeutralAtmo::ReadParameters","Neutral species null","Error /atmosphere/neutral/neutral_species is null. Please put your species here");
		throw err;
	}

	mpParameter->ExistsOrDie("/aero_main/atmosphere/neutral/species_file","The file defining the species is necessary");
	string spfile_name=mpParameter->GetSubFileName("/aero_main/atmosphere/neutral/species_file");
	//string spfile_name=mainParameter->elem("/atmosphere/neutral/species_file");

//	spfile_name=mainParameter->subfilename(spfile_name);

	for(unsigned i=0;i<sp_names.size();++i)
	{

		Specie* tmp=new Specie(sp_names.at(i),mpPhotonGrideV,mpElectronGrideV,mpElectronGridWidtheV,mpProtonGrideV,spfile_name,mpParameter->GetMonteCarlo());
		mAtmoSpecies.push_back(tmp);

	}


	//
	if(!mpParameter->Exists("/aero_main/atmosphere/neutral/model"))
	{

		Log::mE<<"Error /atmosphere/neutral/model does not exists"<<endl;
		Error err("NeutralAtmo::ReadParameters","Neutral model","Error /atmosphere/neutral/model does not exists");
		throw err;
	}	
	int type=-1;
	mpParameter->GetNKey("/aero_main/atmosphere/neutral/model","type",type);
	map< string, ublas::vector<double> > my_atmo;
	if(type==0)
	{// We have an atmosphere data model
		ublas::vector<double> alt;
		mpParameter->Get1DArray("/aero_main/atmosphere/neutral/data/alt",alt);

		for(deque<string>::iterator it=sp_names.begin();it!=sp_names.end();++it)
		{
			ublas::vector<double> tmpdens;
			mpParameter->Get1DArray("/aero_main/atmosphere/neutral/data/"+*it,tmpdens);
			Log::mD<<"Species to interpolate : "<<*it<<"; size of the data : "<<tmpdens.size()<<" size of the altitude grid : "<<alt.size()<<endl;
			if(alt.size()!=tmpdens.size())
			{
				Error err("Neutral atmosphere:: read parameters","Size mismatch"," The size of your altitude data is different than the size of your species "+*it+" data");
				throw err;
			}
			my_atmo[*it]=MathFunction::IntLog(alt,tmpdens,mAltGridKm);
			// Now we interpolate the tmpdens/alt grid on the alt_grid and we put that value inside the corresponding specie.
		}
	}else
	{
		my_atmo=mpPlanet->AtmoModel(mAltGridKm,sp_names,type);
	}

	//////// WE ADD OTHER SPECIES !!!! For example measured O+...
	// read the vector of add species
	//
	 if(mpParameter->Exists("/aero_main/atmosphere/neutral/add_species"))
	 {

		 deque<string> add_names=mpParameter->GetParamDeQue("/aero_main/atmosphere/neutral/add_species/species");
		 // check if the new species are not in the old ones 
		 if(add_names.size()>0)
		 {
			 ublas::vector<double> alt;
			 mpParameter->ExistsOrDie("/aero_main/atmosphere/neutral/add_species/alt","If you want additionnal species, you should provide an atmosphere altitude model");
			 mpParameter->Get1DArray("/aero_main/atmosphere/neutral/add_species/alt",alt);


			 for(deque<string>::iterator st=add_names.begin();st!=add_names.end();++st)
			 {
				 string name=StrReplace(trim(*st),"_PLUS","+");
				 unsigned pos;
				 if(PosInVector(sp_names,name,pos))
				 {
					 // if so, warning+ overload
					 Log::mW<<"Species "<<name<<" overloaded!!!";
				 }else
				 {
					 // init the new species (not overload)

					 Specie* tmp=new Specie(name,mpPhotonGrideV,mpElectronGrideV,mpElectronGridWidtheV,mpProtonGrideV,spfile_name,mpParameter->GetMonteCarlo());
					 mAtmoSpecies.push_back(tmp);
				 }

				 // read the value!
				 ublas::vector<double> tmpdens;
				 mpParameter->Get1DArray("/aero_main/atmosphere/neutral/add_species/"+*st,tmpdens);
				 if(mpParameter->KeyExists( "/aero_main/atmosphere/neutral/add_species/"+*st,"DeleteZero"))
				 {
					 Log::mD<<"The zero values in your species density are set to 1E-40"<<endl;
					 for(unsigned i=0; i< tmpdens.size(); ++i)
					 {
						 if(tmpdens[i] < 1E-40)
							 tmpdens[i] = 1E-40;
					 }
				 }

				 my_atmo[name]=MathFunction::IntLog(alt,tmpdens,mAltGridKm);
			 }
		 }
	 }


	 ////// We add parametrized species
	 if(mpParameter->Exists("/aero_main/atmosphere/neutral/add_parametrized_species"))
	 {

		 deque<string> add_names=mpParameter->GetParamDeQue("/aero_main/atmosphere/neutral/add_parametrized_species/species");
		 // check if the new species are not in the old ones 
		 if(add_names.size()>0)
		 {

			 for(deque<string>::iterator st=add_names.begin();st!=add_names.end();++st)
			 {
				 string name=StrReplace(trim(*st),"_PLUS","+");
				 unsigned pos;
				 Specie* ttmp;
				 if(PosInVector(sp_names,name,pos))
				 {
					 // if so, warning+ overload
					 Log::mW<<"Species "<<name<<" overloaded!!!";
					 ttmp=mAtmoSpecies[pos];
				 }else
				 {
					 // init the new species (not overload)

					 Specie* tmp=new Specie(name,mpPhotonGrideV,mpElectronGrideV,mpElectronGridWidtheV,mpProtonGrideV,spfile_name,mpParameter->GetMonteCarlo());
					 mAtmoSpecies.push_back(tmp);
					 ttmp=tmp;
				 }
				 TiXmlNode* node=mpParameter->GetNode("/aero_main/atmosphere/neutral/add_parametrized_species/"+*st);

				 int model=0;

				 mpParameter->GetNKey(node,"//Model","type",model);
				 switch(model)
				 {
					 case 0:
						 {
							 // Not really physical here, but interesting to check
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 double T0=0;
							 mpParameter->GetValue(node,"//Model/To",T0);
							 double shape=0;
							 mpParameter->GetValue(node,"//Model/Shape",shape);
							 //unsigned pos=CloseInVector(alt0,mAltGridKm);
							 my_atmo[name]=MathFunction::BatesWalkerProfile(mAltGridKm,alt0,dens0,T0,Texo,shape,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass);

						 }
						 break;

					 case 1:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double SH=0;
							 mpParameter->GetValue(node,"//Model/SH",SH);
							 my_atmo[name]=MathFunction::ChapmanProfile(mAltGridKm,alt0,dens0,SH);
						 }
						 break;
					 case 2:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double stddev=0;
							 mpParameter->GetValue(node,"//Model/Dev",stddev);
							 my_atmo[name]=MathFunction::GaussianProfile(mAltGridKm,alt0,dens0,stddev);
						 }
						 break;
					 case 3:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 my_atmo[name]=MathFunction::ExpProfile(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass);

						 }
						 break;
					 case 4:
						 {
							 ublas::vector<double> alts,logval;
							 mpParameter->Get1DArray(node,"//Model/altitudes",alts);
							 mpParameter->Get1DArray(node,"//Model/logvalues",logval);
							 if(alts.size()!=logval.size())
							 {
								 Log::mE<<"Size mismatch: The size of Model/altitudes and Model/logvalues are not compatible. please check"<<endl;
								 Error err("ReadParameters: add_parametrized_species case4","Size mismatch","The size of Model//altitudes and Model/logvalues are not compatible. please check");
								 throw err;
							 }
							 my_atmo[name]=MathFunction::SplineInterpExp(alts,logval,mAltGridKm);
						 }
						 break;
					 case 5:
						 {
							 // Not really physical here, but interesting to check
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/AltT",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/DensT",dens0);
							 double alt1=0;
							 mpParameter->GetValue(node,"//Model/AltM",alt1);
							 double dens1=0;
							 mpParameter->GetValue(node,"//Model/DensM",dens1);
							 if(dens0<0)
								 dens0=1E-50;
							 if(dens1<0)
								 dens1=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 double Tmeso=0;
							 mpParameter->GetValue(node,"//Model/Tmeso",Tmeso);
							 double mixmassamu=0;
							 mpParameter->GetValue(node,"//Model/MixMassAmu",mixmassamu);
							 //unsigned pos=CloseInVector(alt0,mAltGridKm);
							 my_atmo[name]=MathFunction::DoubleExpProfile(mAltGridKm,alt0,alt1,dens0,dens1,Texo,Tmeso,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass,mixmassamu);

						 }
						 break;
					 case 6:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 double SZA =0;
							 mpParameter->GetValue(node,"//Model/SZA",SZA);
							 my_atmo[name]=MathFunction::ChapmanCos(mAltGridKm,alt0,dens0,Texo,SZA,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass);

						 }
						 break;
					 case 7:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 double C =0;
							 mpParameter->GetValue(node,"//Model/C",C);
							 my_atmo[name]=MathFunction::ChapmanVar(mAltGridKm,alt0,dens0,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass);

						 }
						 break;
					 case 8:
						 {
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/alt",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/dens",dens0);
							 if(dens0<0)
								 dens0=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 my_atmo[name]=MathFunction::Epstein(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass);

						 }
						 break;
					 case 9:
						 {
							 // Not really physical here, but interesting to check
							 double alt0=0;
							 mpParameter->GetValue(node,"//Model/AltT",alt0);
							 double dens0=0;
							 mpParameter->GetValue(node,"//Model/DensT",dens0);
							 double alt1=0;
							 mpParameter->GetValue(node,"//Model/AltM",alt1);
							 double dens1=0;
							 mpParameter->GetValue(node,"//Model/DensM",dens1);
							 if(dens0<0)
								 dens0=1E-50;
							 if(dens1<0)
								 dens1=1E-50;
							 double Texo=0;
							 mpParameter->GetValue(node,"//Model/Texo",Texo);
							 double Tmeso=0;
							 mpParameter->GetValue(node,"//Model/Tmeso",Tmeso);
							 double mixmassamu=0;
							 mpParameter->GetValue(node,"//Model/MixMassAmu",mixmassamu);
							 double C = 0;
							 mpParameter->GetValue(node,"//Model/C",C);
							 //unsigned pos=CloseInVector(alt0,mAltGridKm);
							 my_atmo[name]=MathFunction::ExpHyperbolaProfile(mAltGridKm,alt1,alt0,dens1,dens0,Tmeso,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,(ttmp)->mMass,mixmassamu);

						 }
						 break;


					 default:
						 Error err("ReadParameters: add_parametrized_species","Model error"," The model "+ntostr(model)+" is not a valid model for the neutral density profile");
				 }







				 //			ResetSpecieDens(name,int vModel, std::deque<double> vParams, std::deque<double> vAddParams, const ublas::vector<double>& vNeutralT );

				 // read the value!
				 //			 ublas::vector<double> tmpdens;
				 //			 mpParameter->Get1DArray("/aero_main/atmosphere/neutral/add_species/"+*st,tmpdens);
				 //			 my_atmo[name]=MathFunction::IntLog(alt,tmpdens,mAltGridKm);
			 }
		 }
	 }



	//////// MULTIPLICATOR
	for(deque<Specie*>::iterator it=mAtmoSpecies.begin();it!=mAtmoSpecies.end();++it)
	{ 
		(*it)->mTotDensitycm_3=my_atmo[(*it)->mName];

		//we get the multiplicative factor
		string name=StrReplace((*it)->mName,"+","_PLUS");
		if(mpParameter->Exists("/aero_main/atmosphere/neutral/multiplicator/"+name))
		{
			double mult=1.;
			mpParameter->GetValue("/aero_main/atmosphere/neutral/multiplicator/"+name,mult);
			Log::mD<<"Multiplicator activated for "<<name<<" : "<<mult<<endl;
			(*it)->mTotDensitycm_3*=mult;
			//for(unsigned w=0;w<(*it)->mTotDensitycm_3.size();++w)
			//{
			//	(*it)->mTotDensitycm_3[w]*=mult;

			//}
		}
	}

}


void NeutralAtmo::PrintNeutralAtmo(std::string vFilename)
{
	
	if(FileExists(vFilename))
	{

		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Neutral atmosphere considered "<<endl;
	of<<"# Density in cm3"<<endl;
	of<<"# Altitude in km"<<endl;

	deque<Specie*>::iterator it;


	unsigned i=0;

	for(it=mAtmoSpecies.begin();it!=mAtmoSpecies.end();++it)
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

		for(unsigned j=0;j<mAtmoSpecies.size();++j)
		{
			of<<mAtmoSpecies[j]->mTotDensitycm_3[i]<<"\t";
		}
		of<<endl;
	}
}

void NeutralAtmo::PrintBentNeutralAtmo(std::deque<Specie*> vAtmo,ublas::vector<double> vKm,std::string vFilename)
{
	
	if(FileExists(vFilename))
	{

		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Neutral atmosphere considered "<<endl;
	of<<"# Density in cm3"<<endl;
	of<<"# Length in km (Bent case)"<<endl;

	deque<Specie*>::iterator it;


	unsigned i=0;

	for(it=vAtmo.begin();it!=vAtmo.end();++it)
	{
	//	cout<<"Specie "<<i<<" => "<<(*it)->mName<<endl;
		of<<"# "<<(*it)->mName<<endl;

	//	cout<<"Sizes : "<<(*it)->mTotDensitycm_3.size()<<"  alt_grid size : "<<mAltGridKm.size()<<endl;
		assert((*it)->mTotDensitycm_3.size()==vKm.size());
		++i;
	}


	for(unsigned i=0;i<vKm.size();++i)
	{
		of<<vKm[i]<<"\t";

		for(unsigned j=0;j<vAtmo.size();++j)
		{
			of<<vAtmo[j]->mTotDensitycm_3[i]<<"\t";
		}
		of<<endl;
	}
}



void NeutralAtmo::PrintNeutralColdens(std::string vFilename)
{
	
	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Neutral atmosphere considered "<<endl;
	of<<"# Column Density in cm2"<<endl;
	of<<"# Altitude in km"<<endl;

	deque<Specie*>::iterator it;


	unsigned i=0;

	for(it=mAtmoSpecies.begin();it!=mAtmoSpecies.end();++it)
	{
	//	cout<<"Specie "<<i<<" => "<<(*it)->mName<<endl;
		of<<"# "<<(*it)->mName<<endl;

	//	cout<<"Sizes : "<<(*it)->mColDenscm_2.size()<<"  alt_grid size : "<<mAltGridKm.size()<<endl;
		assert((*it)->mColDenscm_2.size()==mAltGridKm.size());
		++i;
	}


	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";

		for(unsigned j=0;j<mAtmoSpecies.size();++j)
		{
			of<<mAtmoSpecies[j]->mColDenscm_2[i]<<"\t";
		}
		of<<endl;
	}
}

void NeutralAtmo::PrintBentNeutralColdens(std::deque<Specie*> vAtmo,ublas::vector<double> vKm,std::string vFilename)
{
	
	if(FileExists(vFilename))
	{

		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Neutral atmosphere considered "<<endl;
	of<<"# Column Density in cm2"<<endl;
	of<<"# Length in km (Bent case)"<<endl;

	deque<Specie*>::iterator it;


	unsigned i=0;

	for(it=vAtmo.begin();it!=vAtmo.end();++it)
	{
	//	cout<<"Specie "<<i<<" => "<<(*it)->mName<<endl;
		of<<"# "<<(*it)->mName<<endl;

	//	cout<<"Sizes : "<<(*it)->mTotDensitycm_3.size()<<"  alt_grid size : "<<mAltGridKm.size()<<endl;
		assert((*it)->mColDenscm_2.size()==vKm.size());
		++i;
	}


	for(unsigned i=0;i<vKm.size();++i)
	{
		of<<vKm[i]<<"\t";

		for(unsigned j=0;j<vAtmo.size();++j)
		{
			of<<vAtmo[j]->mColDenscm_2[i]<<"\t";
		}
		of<<endl;
	}
}



void NeutralAtmo::ComputeVerticalColumnDensity(const ublas::vector<double>& vNeutralT)
{
	Log::mD<<"hello coldens"<<endl;
	// Computation of the scale height
	//         initial values of cold are set by assuming neutrals in
	//      diffusive equilibrium
	//        scale height definition: H = kT/mg
	//
	//          H      = {k/(amu*g)}  .    {(Re+z)/Re}**2   . T(z)/atomas
	//          with k boltzman
	//          	 amu atomic mass unit
	//          	 g gravity
	//          	 Re : radius of the planet
	//          	 z : altitude were the computation is made
	//          	 T(z) : neutral temperature at z
	//          	 atomas : mass of the considered specie in amu
	// here scaleh is in fact H*atomas. It is in cm, because atomas is without unit
	assert(vNeutralT.size()!=0);
	assert(vNeutralT.size()==mAltGridKm.size());
	double scaleh=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[0])/(mpPlanet->mRKm),2)*vNeutralT[0];

	// We use the scale height to compute the column density at the top
	// but after, we integrate. Thus, in case of error at the top

	deque<Specie*>::iterator sp;

	//	cout<<"Hello column density"<<endl;

	unsigned nbalt=mAltGridKm.size();
	for(sp=mAtmoSpecies.begin();sp!=mAtmoSpecies.end();++sp)
	{
		assert((*sp)->mTotDensitycm_3.size()!=0);

		double coldensinit=(*sp)->mTotDensitycm_3[0]*scaleh/(*sp)->mMass;

		double scalhspinit=scaleh/(*sp)->mMass;
		(*sp)->mColDenscm_2.resize(nbalt);
		(*sp)->mScaleHcm.resize(nbalt);
		(*sp)->mColDenscm_2[0]=(coldensinit);
		(*sp)->mScaleHcm[0]=(scalhspinit);

		for(unsigned i=1;i<nbalt;++i)
		{
			(*sp)->mColDenscm_2[i]=(MathFunction::TrapzInt(mAltGridKm,(*sp)->mTotDensitycm_3,i)*1E5+coldensinit);

			scalhspinit=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[i])/(mpPlanet->mRKm),2)*vNeutralT[i]/(*sp)->mMass;

			(*sp)->mScaleHcm[i]=(scalhspinit);
		}
	}

	Log::mD<<"bye bye column density"<<endl;


}



void NeutralAtmo::ClearProductions()
{
	deque<Specie*>::iterator it;
	for(it=mAtmoSpecies.begin();it!=mAtmoSpecies.end();++it)
	{
		(*it)->ClearProd();
	}
}




void NeutralAtmo::ResetSpecieDens(std::string vName,int vModel, std::deque<double> vParams, std::deque<double> vAddParams,const ublas::vector<double>& vNeutralT)
{
	int position=SpecieUtils::PosOfSp(vName,mAtmoSpecies);
	if(position<0)
	{
		Error err("ResetSpecieDens","Specie not found","The species "+vName+" has not been found in the list of the different species for reset");
	}
	Specie* sp=mAtmoSpecies[position];

	// We reset the density
	switch(vModel)
	{
		case 0:
			{
				assert(vParams.size()==4);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(1);
				double T0=vParams.at(2);
				double shape=vParams.at(3);
				//unsigned pos=CloseInVector(alt0,mAltGridKm);
				(sp)->mTotDensitycm_3=MathFunction::BatesWalkerProfile(mAltGridKm,alt0,dens0,T0,Texo,shape,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass);

			}
			break;

		case 1:
			{
				assert(vParams.size()==3);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				(sp)->mTotDensitycm_3=MathFunction::ChapmanProfile(mAltGridKm,vParams.at(0),dens0,vParams.at(2));
			}
			break;
		case 2:
			{
				assert(vParams.size()==3);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				(sp)->mTotDensitycm_3=MathFunction::GaussianProfile(mAltGridKm,vParams.at(0),dens0,vParams.at(2));
			}
			break;
		case 3:
			{
				assert(vParams.size()==2);
				double alt0=vAddParams.at(0);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(1);
				(sp)->mTotDensitycm_3=MathFunction::ExpProfile(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass);
				
			}
			break;
		case 4:
			{
				unsigned siz=static_cast<unsigned>(vAddParams.at(0));
				ublas::vector<double> alt(siz),logval(siz);
				alt.clear();
				logval.clear();
				for(unsigned k=0;k<siz;++k)
				{
					alt[k]=vAddParams.at(k+1);
					logval[k]=vParams.at(k);
				}
				(sp)->mTotDensitycm_3=MathFunction::SplineInterpExp(alt,logval,mAltGridKm);
//				Log::mI<<mAltGridKm<<endl;
//				Log::mI<<(sp)->mTotDensitycm_3<<endl;
//				Log::mI<<MathFunction::SplineInterp(alt,logval,mAltGridKm)<<endl;
			}
			break;
		case 5:
			{
				assert(vParams.size()==4);
				assert(vAddParams.size()==3);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double alt1=vAddParams.at(1);
				double mixmassamu=vAddParams.at(2);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double dens1=vParams.at(1);
				if(dens1<0)
					dens1=1E-50;

				double Texo=vParams.at(2);
				double Tmeso=vParams.at(3);
				(sp)->mTotDensitycm_3=MathFunction::DoubleExpProfile(mAltGridKm,alt0,alt1,dens0,dens1,Texo,Tmeso,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass,mixmassamu);
			}
			break;
		case 6:
			{
				assert(vParams.size()==4);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				double SZA=vParams.at(3);
				(sp)->mTotDensitycm_3=MathFunction::ChapmanCos(mAltGridKm,alt0,dens0,Texo,SZA,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass);
			}
			break;
		case 7:
			{
				assert(vParams.size()==3);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				double C=vAddParams.at(3);
				(sp)->mTotDensitycm_3=MathFunction::ChapmanVar(mAltGridKm,alt0,dens0,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass);
			}
			break;
		case 8:
			{
				assert(vParams.size()==3);
				double alt0=vParams.at(0);
				double dens0=vParams.at(1);
				if(dens0<0)
					dens0=1E-50;
				double Texo=vParams.at(2);
				(sp)->mTotDensitycm_3=MathFunction::Epstein(mAltGridKm,alt0,dens0,Texo,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass);
			}
			break;
		case 9:
			{
				assert(vParams.size()==5);
				assert(vAddParams.size()==3);
				// Not really physical here, but interesting to check
				double alt0=vAddParams.at(0);
				double alt1=vAddParams.at(1);
				double mixmassamu=vAddParams.at(2);
				double dens0=vParams.at(0);
				if(dens0<0)
					dens0=1E-50;
				double dens1=vParams.at(1);
				if(dens1<0)
					dens1=1E-50;

				double Texo=vParams.at(2);
				double Tmeso=vParams.at(3);
				double C = vParams.at(4);
				(sp)->mTotDensitycm_3=MathFunction::ExpHyperbolaProfile(mAltGridKm,alt1,alt0,dens1,dens0,Tmeso,Texo,C,mpPlanet->mRKm,mpPlanet->mGms_2,(sp)->mMass,mixmassamu);
			}
			break;
		default:
			Error err("Reset neutral density","Model error"," The model "+ntostr(vModel)+" is not a valid model for the neutral density profile");
	}


	// And we reset the scale height and the column density
	double scaleh=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[0])/(mpPlanet->mRKm),2)*vNeutralT[0];
	unsigned nbalt=mAltGridKm.size();
	assert((sp)->mTotDensitycm_3.size()!=0);

	double coldensinit=(sp)->mTotDensitycm_3[0]*scaleh/(sp)->mMass;

	double scalhspinit=scaleh/(sp)->mMass;
	(sp)->mColDenscm_2.resize(nbalt);
	(sp)->mScaleHcm.resize(nbalt);
	(sp)->mColDenscm_2[0]=(coldensinit);
	(sp)->mScaleHcm[0]=(scalhspinit);

	for(unsigned i=1;i<nbalt;++i)
	{
		(sp)->mColDenscm_2[i]=(MathFunction::TrapzInt(mAltGridKm,(sp)->mTotDensitycm_3,i)*1E5+coldensinit);

		scalhspinit=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[i])/(mpPlanet->mRKm),2)*vNeutralT[i]/(sp)->mMass;

		(sp)->mScaleHcm[i]=(scalhspinit);
	}



}
void NeutralAtmo::ResetSpecieDensInterp(std::string vName, ublas::vector<double> vDenscm_3,const ublas::vector<double>& vNeutralT)
{
	int position=SpecieUtils::PosOfSp(vName,mAtmoSpecies);
	if(position<0)
	{
		Error err("ResetSpecieInterp","Specie not found","The species "+vName+" has not been found in the list of the different species for reset");
	}
	Specie* sp=mAtmoSpecies[position];
	(sp)->mTotDensitycm_3=vDenscm_3;
	double scaleh=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[0])/(mpPlanet->mRKm),2)*vNeutralT[0];
	unsigned nbalt=mAltGridKm.size();

	double coldensinit=(sp)->mTotDensitycm_3[0]*scaleh/(sp)->mMass;

	double scalhspinit=scaleh/(sp)->mMass;
	(sp)->mColDenscm_2.resize(nbalt);
	(sp)->mScaleHcm.resize(nbalt);
	(sp)->mColDenscm_2[0]=(coldensinit);
	(sp)->mScaleHcm[0]=(scalhspinit);
	for(unsigned i=1;i<nbalt;++i)
	{
		(sp)->mColDenscm_2[i]=(MathFunction::TrapzInt(mAltGridKm,(sp)->mTotDensitycm_3,i)*1E5+coldensinit);

		scalhspinit=SCALEH_CONST/(mpPlanet->mGms_2)*pow((mpPlanet->mRKm+mAltGridKm[i])/(mpPlanet->mRKm),2)*vNeutralT[i]/(sp)->mMass;

		(sp)->mScaleHcm[i]=(scalhspinit);
	}
}
