 /** 
 * \file crs.cpp
 * \brief Implements the cross section class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: crs.cpp 1346 2011-11-11 05:05:49Z gronoff $
 *
 */

#include "crs.hpp"
using namespace std;

CrossSection::CrossSection(ublas::vector<double> vGreV)
{
	mGrideV=vGreV;
	mIsDefinedTotalCrs=false;
	mDefinedTotalCrsDirectly=false;
	mMaxDefinedEnergyTotalCrs=0.;
	mExtrapolationLogLog=false;
	mParamLoaded=false;
	mIsDefinedCrs=false;
	mbTotalSumIneel=false;
	mbIsDefinedElastic=false;
}

CrossSection::~CrossSection()
{
	if(mParamLoaded)
	{
	//	Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Detach cross section"<<endl;
		mpParam->DetachFileParameter();
	}
}

bool CrossSection::LoadCrs(std::string vFilename,std::string vName,bool bIsMonteCarlo,bool bLogInterp)
{
	//XmlParameters*param=new XmlParameters(vFilename);
	mpParam=XmlParameters::AttachFileParameter(vFilename);
	mParamLoaded=true;
	mSpecies=vName;
	bool resu=LoadCrs(mpParam,vName,bIsMonteCarlo,bLogInterp);
	//delete param;
	return resu;
}


ublas::vector<double> CrossSection::ExtractCrs(XmlParameters* pParams,TiXmlNode* vNode,double vThreshold,double& rMaxEnergy)
{

	if(pParams->Exists(vNode,"//Shirai"))
	{
		Shirai sh(pParams,vNode,vThreshold);
		if(!pParams->Exists(vNode,"//SetMCActive"))
		{	
			return sh.GetCrs(mGrideV);
		}else
		{
			ublas::vector<double> crs= sh.GetCrs(mGrideV);
			double percent=0;
			pParams->GetValue(vNode,"//uncertainty",percent);
			return pParams->ReturnMC(crs,percent);
		}
	}


	ublas::vector<double> tenergy;
	pParams->Get1DArray(vNode,"//Egrid",tenergy);
	if(tenergy[0]>(*(tenergy.end()-1)))
	{
		rMaxEnergy=tenergy[0];
	}else
	{
		rMaxEnergy=(*(tenergy.end()-1));
	}
	ublas::vector<double> tcrs;
	pParams->Get1DArray(vNode,"//Cross",tcrs);

	ublas::vector<double> resu=MathFunction::IntLin(tenergy,tcrs,mGrideV);
	for(unsigned i=0; i<resu.size();++i)
	{
		if(resu[i]<0)
			resu[i]=0;
		if(mGrideV[i]<vThreshold)
			resu[i]=0;
	}
	return resu;
}


ublas::vector<double> CrossSection::ExtractCrsLog(XmlParameters* pParams,TiXmlNode* vNode,double vThreshold,double& rMaxEnergy)
{

	if(pParams->Exists(vNode,"//Shirai"))
	{// There is no problem of interpolation with shirai like cross sections
		Shirai sh(pParams,vNode,vThreshold);
		if(!pParams->Exists(vNode,"//SetMCActive"))
		{	
			Log::mD<<"Shirai node with no uncertainty : "<<endl;
			return sh.GetCrs(mGrideV);
		}else
		{
			ublas::vector<double> crs= sh.GetCrs(mGrideV);
			double percent=0;
			pParams->GetValue(vNode,"//uncertainty",percent);
			Log::mD<<"Shirai node with uncertainty : "<<percent<<endl;
			return pParams->ReturnMC(crs,percent);
		}

	}
	ublas::vector<double> tenergy;
	pParams->Get1DArray(vNode,"//Egrid",tenergy);
	if(tenergy[0]>(*(tenergy.end()-1)))
	{
		rMaxEnergy=tenergy[0];
	}else
	{
		rMaxEnergy=(*(tenergy.end()-1));
	}
	ublas::vector<double> tcrs;
	pParams->Get1DArray(vNode,"//Cross",tcrs);

	for(unsigned i=0;i<tcrs.size();++i)
	{
		if(tcrs[i]<=1E-42)
			tcrs[i]=1E-42;
	}

	ublas::vector<double> resu=MathFunction::IntLog(tenergy,tcrs,mGrideV);
	for(unsigned i=0; i<resu.size();++i)
	{
		if(resu[i]<0)
			resu[i]=0;
		if(mGrideV[i]<vThreshold)
			resu[i]=0;
	}
	if(mExtrapolationLogLog)
	{
		ublas::vector<double> resu_nouveau=MathFunction::IntLogLog(tenergy,tcrs,mGrideV);
		for(unsigned i=0; i<resu_nouveau.size();++i)
		{
			if(resu_nouveau[i]<0)
				resu_nouveau[i]=0;
			// It is for the upper values, therefore, no need for that
			//	if(mGrideV[i]<vThreshold)
			//		resu[i]=0;
		}
		for(unsigned i=0;i<resu.size();++i)
		{
			if(mGrideV[i]>rMaxEnergy)
			{
				resu[i]=resu_nouveau[i];
			}
		}

	}

	return resu;
}
bool CrossSection::LoadCrs(XmlParameters*pParams,std::string vName,bool bIsMonteCarlo, bool bLogInterp)
{
	mIsDefinedCrs=true;

	if(bLogInterp)
		mExtrapolationLogLog=true;
	if(bIsMonteCarlo)
	{
		pParams->SetMonteCarloActive();
	}
	bool total_is_the_sum=false;// for photoionisation

//	mbTotalSumIneel=false; // electron impact
	
	
	// To be sure to have the real name
	string name=trim(vName);

	if(!pParams->Exists("/crs/"+name))
	{
		//Log::SetPriority(Log::WARNING);
		Log::mW<<"Impossible to find the cross section of your species : ["<<name<<"]"<<endl;
		return false;
	}
	Log::mL<<"Name of your cross section species :"<<name<<endl;

	if(pParams->Exists("/crs/"+name+"/ZeroCrs"))
	{
		mIsDefinedTotalCrs=true;
		ublas::vector<double> tmp(mGrideV.size());
		tmp.clear();
		mTotalCrscm2=tmp;


		mIsDefinedCrs=false;
		//Log::SetPriority(Log::WARNING);
		Log::mW<<"======================================="<<endl;
		Log::mW<<endl<<endl;
		Log::mW<<"YOUR SPECIE "<<name<<" HAS BEEN DEFINED WITHOUT CROSS SECTION"<<endl;
		Log::mW<<"IT IS VALID WHEN YOU WANT TO NEGLET THIS SPECIE, AND WHEN YOU NEED IT ONLY FOR CHEMICAL COMPUTATIONS"<<endl;
		Log::mW<<endl<<endl;
		Log::mW<<"======================================="<<endl;

		return true;

	}

	if(pParams->Exists("/crs/"+name+"/TotalCrs"))
	{
		mIsDefinedTotalCrs=true;
		mDefinedTotalCrsDirectly=true;

		TiXmlNode* totalnode=pParams->GetNode("/crs/"+name+"/TotalCrs");


		if(bLogInterp)
		{
			mTotalCrscm2=ExtractCrsLog(pParams,totalnode,0.,mMaxDefinedEnergyTotalCrs);
		}else
		{
			mTotalCrscm2=ExtractCrs(pParams,totalnode,0.,mMaxDefinedEnergyTotalCrs);
		}

	}

	if(pParams->Exists("/crs/"+name+"/DisableTotalCrsWarning"))
	{
		if(mIsDefinedTotalCrs)
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"Warning on total crs disabled, while totalcrs defined\n"<<
				"this is seen as an error (why disabling warning \n"<<
				"while no reasons ??? ). Please disable <DisableTotalCrsWarning/>\n"<<
				"or disable the TotalCrs data"<<endl;
			Error err("LoadCrs","DisableTotalCrsWarning","Warning on total crs disabled, while totalcrs defined\nthis is seen as an error (why disabling warning \n while no reasons ??? ). Please disable <DisableTotalCrsWarning/>\nor disable the TotalCrs data");
			throw err;
		}
		mIsDefinedTotalCrs=true;
	}
	if(pParams->Exists("/crs/"+name+"/TotalCrsIsTheSum"))
	{
		if(mIsDefinedTotalCrs)
		{

			//Log::SetPriority(Log::WARNING);
			Log::mE<<"Warning on total crs disabled, while totalcrsisthesum defined\n"<<
				"this is seen as an error (why adding \n"<<
				"while no reasons ??? ). Please disable <TotalCrsIsTheSum/>\n"<<
				"or disable the TotalCrs data"<<endl;
			Error err("LoadCrs","TotalCrsIsTheSum","Warning on total crs disabled, while totalcrsisthesum defined\nthis is seen as an error (why disabling warning \n while no reasons ??? ). Please disable <TotalCrsIsTheSum/>");
			throw err;
		}
		total_is_the_sum=true;
	}

	if(pParams->Exists("/crs/"+name+"/TotalCrsElectron"))
	{
		if(mIsDefinedTotalCrs)
		{

			//Log::SetPriority(Log::WARNING);
			Log::mE<<"Warning on total crs disabled, while totalcrselectron defined\n"<<
				"this is seen as an error (why adding \n"<<
				"while no reasons ??? ). Please disable <TotalCrsElectron/>\n"<<
				"or disable the TotalCrs data"<<endl;
			Error err("LoadCrs","TotalCrsElectron","Warning on total crs disabled, while totalcrsisthesum defined\nthis is seen as an error (why disabling warning \n while no reasons ??? ). Please disable <TotalCrsElectron/>");
			throw err;
		}
		mbTotalSumIneel=true;
	}

	if(pParams->Exists("/crs/"+name+"/TotalCrsProton") || pParams->Exists("/crs/"+name+"/TotalCrsHydrogen"))
	{
		if(mIsDefinedTotalCrs)
		{

			//Log::SetPriority(Log::WARNING);
			Log::mE<<"Warning on total crs disabled, while totalcrsproton or totalcrshydrogen defined\n"<<
				"this is seen as an error (why adding \n"<<
				"while no reasons ??? ). Please disable <TotalCrsHydrogen/> and/or <TotalCrsProton>\n"<<
				"or disable the TotalCrs data"<<endl;
			Error err("LoadCrs","TotalCrsProton/Hydrogen","Warning on total crs disabled, while totalcrsisthesum defined\nthis is seen as an error (why disabling warning \n while no reasons ??? ). Please disable <TotalCrsHydrogen/Proton/>");
			throw err;
		}
		mbTotalSumIneel=true;
	}
	// Loop on Process
	vector<TiXmlNode*> processes=pParams->GetNodes("/crs/"+name+"/Process");
	vector<TiXmlNode*>::iterator it;

	//	cout<<"Processes numbers : "<<processes.size()<<endl;
	for(it=processes.begin();it!=processes.end();++it)
	{
		//rem : node=*it
		//	TiXmlNode* node=*it;
		string procname=pParams->GetKey(*it,"/","name");
		//		cout<<"procname = "<<procname<<endl;

		// On remplit le vecteur:
		mProcessNames.push_back(procname);


		double electrons=0;
		pParams->GetNKey(*it,"/","electrons",electrons);

		mNumberOfElectrons.push_back(electrons);
		double threshold=0;
		pParams->GetNKey(*it,"/","threshold",threshold);
		mThresholdseV.push_back(threshold);

		if(pParams->KeyExists(*it,"/","ions"))
		{
			double ions=0;
			pParams->GetNKey(*it,"/","ions",ions);
			mNumberOfIons.push_back(ions);
		}else
		{
			mNumberOfIons.push_back(static_cast<double>(electrons));
		}


		if(pParams->Exists(*it,"//Auger"))
		{

			vector<TiXmlNode*> augers=pParams->GetNodes(*it,"//Auger");

			vector<TiXmlNode*>::iterator jt;
			std::deque<double> tmpAugerEnergy;
			std::deque<double> tmpAugerEfficiency;
			for(jt=augers.begin();jt!=augers.end();++jt)
			{
				double tmpaugerenergy=0.;
				pParams->GetNKey(*jt,"/","energy",tmpaugerenergy);
				mIsAuger.push_back(true);
				tmpAugerEnergy.push_back(tmpaugerenergy);
				Log::mD<<" Auger process : "<<procname<<" energy for the auger electron: "<<tmpaugerenergy<<endl;
				if(pParams->NbParams(*jt,"/","fact")>0)
				{
					double tmpaugereff;
					pParams->GetNKey(*jt,"/","fact",tmpaugereff);
					tmpAugerEfficiency.push_back(tmpaugereff);
					Log::mD<<"Factor (auger efficiency for this process: "<<tmpaugereff<<endl;

				}else
				{
					tmpAugerEfficiency.push_back(1.);
				}

			}
			mAugerEnergy.push_back(tmpAugerEnergy);
			mAugerEfficiency.push_back(tmpAugerEfficiency);
		}else
		{
			mIsAuger.push_back(false);
			std::deque<double> tmpAugerEnergy;
			std::deque<double> tmpAugerEfficiency;
			mAugerEnergy.push_back(tmpAugerEnergy);
			mAugerEfficiency.push_back(tmpAugerEfficiency);
		}	




		// we work on the different species!

		vector<TiXmlNode*> species=pParams->GetNodes(*it,"//Species/Specie");
		vector<TiXmlNode*>::iterator jt;

		std::deque<double> tmp_sp_multiplicator;
		deque<SpecieId> tmp_sp_name;

		//		cout<<"Taille species : "<<species.size()<<endl;
		for(jt=species.begin();jt!=species.end();++jt)
		{
			string spname=pParams->GetKey(*jt,"/","name");
			string spstate=pParams->GetKey(*jt,"/","state");
			//			cout<<"spstate : "<<spstate<<endl;
			SpecieId tmpsp(spname,spstate);
			tmp_sp_name.push_back(tmpsp);
			double bratio=1.;
			if(pParams->NbParams(*jt,"/","number")==1)
			{
				pParams->GetNKey(*jt,"/","number",bratio);
			}
			tmp_sp_multiplicator.push_back(bratio);
		}


		mCreatedSpecies.push_back(tmp_sp_name);
		mCreatedSpeciesMultiplicator.push_back(tmp_sp_multiplicator);

		double tmp_max_energy;
		ublas::vector<double> ProcessCrs;
		if(bLogInterp)
		{
			ProcessCrs=ExtractCrsLog(pParams,*it,threshold,tmp_max_energy);
		}else
		{
			ProcessCrs=ExtractCrs(pParams,*it,threshold,tmp_max_energy);
		}
		mCrscm2.push_back(ProcessCrs);
		mMaxDefinedEnergyProcess.push_back(tmp_max_energy);

		bool ionization_defined=false;
		if(pParams->Exists(*it,"//Ionization"))
		{
			ionization_defined=true;

			mIonizationThresholdeV.push_back(threshold);
			// The following commented technique was a great idea
			// unfortunately, the address of ProcessCrs
			// and the address of mCrscm_2[position]
			// can change when new cross sections are added
			//	mExcitationCrscm2.push_back(&ProcessCrs);
			//	unsigned sizpos=mCrscm_2.size()-1;
			//	mIonizationCrscm2.push_back(&(mCrscm_2[sizpos]));
			//	mIonizationCrscm2.push_back(&(*(mCrscm_2.end()-1)));//&ProcessCrs);
			// Therefore, we use position...
			mIonizationCrsPosition.push_back(mCrscm2.size()-1);
		}
		bool excitation_defined = false;
		if(pParams->Exists(*it,"//Excitation"))
		{
			excitation_defined = true;
			if(ionization_defined)
			{
				Error err("LoadCrs","A specie is defined Ionization and Excitation","Please correct: a specie cannot be both of them, see documentation");
				throw err;
			}
			mExcitationThresholdeV.push_back(threshold);
			// The following commented technique was a great idea
			// unfortunately, the address of ProcessCrs
			// and the address of mCrscm_2[position]
			// can change when new cross sections are added
			//	mExcitationCrscm2.push_back(&ProcessCrs);
			//	mIonizationCrsPosition.push_back(mCrscm2.size()-1);
			// we use position
			mExcitationCrsPosition.push_back(mCrscm2.size()-1);
		}
		if(pParams->Exists(*it,"//Exchange"))
		{
			if(ionization_defined)
			{
				Error err("LoadCrs","A specie is defined Ionization and Exchange","Please correct: a specie cannot be both of them, see documentation");
				throw err;
			}
			if(excitation_defined)
			{
				Error err("LoadCrs","A specie is defined Excitation and Exchange","Please correct: a specie cannot be both of them, see documentation");
				throw err;
			}
			mExchangeThresholdeV.push_back(threshold);
			mExchangeCrsPosition.push_back(mCrscm2.size()-1);
		}
	}
	// Finalisation : on remplit totalcrs là où il faut.


	if(total_is_the_sum)
	{
		Log::mI<<"Total cross section is the sum of the others"<<endl;
		mTotalCrscm2.resize(mGrideV.size(),0);
		for(unsigned i=0;i<mGrideV.size();++i)
		{
			for(unsigned j=0;j<mCrscm2.size();++j)
			{
				mTotalCrscm2[i]+=mCrscm2[j][i];
			}
		}
		mIsDefinedTotalCrs=true;

	}

	return true;
}


void CrossSection::PrintCrs(std::string vFilename)
{
	// assert de base: on verifie que les nombres sont identiques
	// nb : si on n'est pas en mode debug: les assert ne servent a rien)
	unsigned nb=mProcessNames.size();
	assert(nb==mCreatedSpecies.size());
	assert(nb==mNumberOfElectrons.size());
	assert(nb==mThresholdseV.size());
	assert(nb==mCrscm2.size());
	if(mIsDefinedTotalCrs)
	{
		assert(mGrideV.size()==mTotalCrscm2.size());
	}


	assert(CheckCreatedSpecies());



	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::INFO);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Cross section for "<<mSpecies<<endl;

	of<<"# Energy in eV"<<endl;
	of<<"# Cross section in cm2"<<endl;
	if(mIsDefinedTotalCrs)
	{
		of<<"# Total Cross Section"<<endl;
	}
	if(mbIsDefinedElastic)
	{
		of<<"# Elastic Cross Section"<<endl;
	}
	for(unsigned i=0;i<mProcessNames.size();++i)
	{
		of<<"# Process "<<i<<" : "<<mProcessNames[i]<<endl;
	}



	if(nb>0)
	{
		unsigned nb_e=mCrscm2[0].size();

		for(unsigned i=0;i<nb_e;++i)
		{
			of<<mGrideV[i]<<"\t";
			if(mIsDefinedTotalCrs)
			{
				of<<mTotalCrscm2[i]<<"\t";
			}
			if(mbIsDefinedElastic)
			{
				of<<mElasticCrscm2[i]<<"\t";
			}
			for(unsigned j=0;j<nb;++j)
			{
				of<<mCrscm2[j][i]<<"\t";
			}
			of<<endl;
		}

	}
	/// Very important: display the  credits!
	of<<Log::msMessageLog<<endl;
	of.close();
}


bool CrossSection::CheckCreatedSpecies()
{
	if(mCreatedSpecies.size()!=mCreatedSpeciesMultiplicator.size())
	{
		Log::mE<<"Created species != created_species_multiplicator"<<endl;
		return false;
	}
	for(unsigned i=0;i<mCreatedSpecies.size();++i)
	{
		if(mCreatedSpecies[i].size()!=mCreatedSpeciesMultiplicator[i].size())
		{
			Log::mE<<"Created species != created_species_multiplicator"<<endl;
			return false;
		}
	}
	return true;
}


