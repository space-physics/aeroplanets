/**
 * \file emission.cpp
 * \brief Implements the emission object : it stores the different emits, calls the geometry system, creates the synthetic spectrum 
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emission.cpp 1241 2011-03-25 20:15:03Z gronoff $
 */

#include "emission.hpp"
using namespace std;


void Emission::InitPath()
{
	mPath.reset((new GeoObs(mpParam,mpPlanet->mRKm,mAltGridKm[0])));
}

void Emission::InitEmission()
{
	// List of emissions
	//
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitN2A3S(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitO1S(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitND(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitO1D(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitCOa3Pi(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitCO2pB(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitCHA2D(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitOplus2P(mAltGridKm)));
	mEmitList.push_back(boost::shared_ptr<Emit>(new EmitO3p3P(mAltGridKm)));

	// End of the definition
	for(size_t lsize=0;lsize<mEmitList.size();++lsize)
	{
		mColFileName.push_back("");
		mProfFileName.push_back("");
		mExtraFileName.push_back("");
		mUseAbsorptionList.push_back(false);
		mbCompute.push_back(false);
		std::map< std::string, std::map<double,double> > tmp;
		mAbsorptionValscm_2.push_back(tmp);
	}
}

Emission::Emission(XmlParameters* vpParam,boost::shared_ptr<Chem> vChem,
				Planete* vpPlanet,
				ublas::vector<double> vAltGridKm
				):mpParam(vpParam),mChem(vChem),mpPlanet(vpPlanet),mAltGridKm(vAltGridKm)
{
	InitPath();
	InitEmission();
	mbComputeAll=false;
}



void Emission::InitMultipleEmit()
{
	Log::mS<<"Init multiple emit"<<endl;
	mbComputeAll=mpParam->Exists("/aero_main/emissions/use_all_emissions");
	vector<TiXmlNode*> femit=mpParam->GetNodes("/aero_main/emissions/emission");
	vector<TiXmlNode*>::iterator sit;
	/* We erase previous calculations */
	for(size_t lsize=0;lsize<mEmitList.size();++lsize)
	{
		mColFileName.at(lsize)="";
		mProfFileName.at(lsize)="";
		mExtraFileName.at(lsize)="";
		mUseAbsorptionList.at(lsize)=false;
		mbCompute.at(lsize)=false;
		std::map< std::string, std::map<double,double> > tmp;
		mAbsorptionValscm_2.at(lsize)=tmp;
	}
	
	for(sit=femit.begin();sit!=femit.end();++sit)
	{
		unsigned id=0;

		if(mpParam->KeyExists(*sit,"/","id"))
		{
			mpParam->GetNKey(*sit,"/","id",id);
			/*	if(mpParam->Exists(*sit,"//column_file"))
				{
				mColFileName.at(id)=mpParam->Elem(*sit,"//column_file")+vSuffix;
				}
				if(mpParam->Exists(*sit,"//profile_file"))
				{
				mProfFileName.at(id)=mpParam->Elem(*sit,"//profile_file")+vSuffix;
				}
				if(mpParam->Exists(*sit,"//extra_file"))
				{
				mExtraFileName.at(id)=mpParam->Elem(*sit,"//extra_file")+vSuffix;
				}
				*/

			if(mpParam->Exists(*sit,"//use_absorption"))
			{
				mUseAbsorptionList.at(id)=true;
				vector< TiXmlNode* > myabs = mpParam->GetNodes(*sit,"//absorption");
				vector< TiXmlNode*>::iterator ait;
				for(ait=myabs.begin();ait!=myabs.end();++ait)
				{
					if( mpParam->KeyExists(*ait,"/","species") )
					{
						double lambnm = -1;
						if( mpParam->KeyExists(*ait,"/","lamb") )
						{
							mpParam->GetNKey(*ait,"/","lamb",lambnm);
						}
						string aspecies = mpParam->GetKey(*ait,"/","species");
						double value;
						mpParam->GetValue(*ait,"/",value);

						Log::mD<<"Absorption for the species "<<aspecies<<" at the wavelength "<<lambnm<<" (-1 for the total absorption for the process): "<< value<<" cm-1"<<endl;
						mAbsorptionValscm_2.at(id)[aspecies][lambnm]=value;

					}
				}
			}
			mbCompute.at(id)=true;
		}else
		{
			std::string name,state;

			name = mpParam->GetKey(*sit,"/","name");
			state = mpParam->GetKey(*sit,"/","state");
			std::map<double,double> freqbratio;
			vector< TiXmlNode* > myfreq = mpParam->GetNodes(*sit,"//freq");
			vector< TiXmlNode*>::iterator frt;
			for(frt=myfreq.begin();frt!=myfreq.end();++frt)
			{
				double frequency,bratio;
				mpParam->GetNKey(*frt,"/","freq",frequency);
				mpParam->GetNKey(*frt,"/","branching",bratio);
				freqbratio[frequency]=bratio;
			}
			id = mAbsorptionValscm_2.size();
			mEmitList.push_back(boost::shared_ptr<Emit>(new EmitAllowed(mAltGridKm,id,name,state,freqbratio)));
			mColFileName.push_back("");
			mProfFileName.push_back("");
			mExtraFileName.push_back("");
			mUseAbsorptionList.push_back(false);
			mbCompute.push_back(false);
			std::map< std::string, std::map<double,double> > tmp;
			mAbsorptionValscm_2.push_back(tmp);

			if(mpParam->Exists(*sit,"//use_absorption"))
			{
				mUseAbsorptionList.at(id)=true;
				vector< TiXmlNode* > myabs = mpParam->GetNodes(*sit,"//absorption");
				vector< TiXmlNode*>::iterator ait;
				for(ait=myabs.begin();ait!=myabs.end();++ait)
				{
					if( mpParam->KeyExists(*ait,"/","species") )
					{
						double lambnm = -1;
						if( mpParam->KeyExists(*ait,"/","lamb") )
						{
							mpParam->GetNKey(*ait,"/","lamb",lambnm);
						}
						string aspecies = mpParam->GetKey(*ait,"/","species");
						double value;
						mpParam->GetValue(*ait,"/",value);

						Log::mD<<"Absorption for the species "<<aspecies<<" at the wavelength "<<lambnm<<" (-1 for the total absorption for the process): "<< value<<" cm-1"<<endl;
						mAbsorptionValscm_2.at(id)[aspecies][lambnm]=value;

					}
				}
			}


		}
	}

}
void Emission::CompMultipleEmit()
{
	Log::mS<<"Start multiple emit"<<endl;
	mEmitSpectrumIntegratedR.erase(mEmitSpectrumIntegratedR.begin(),mEmitSpectrumIntegratedR.end());
	for(size_t lsize=0;lsize<mEmitList.size();++lsize)
	{
		Log::mD<<"Emission number "<<lsize<<endl;
		if(mbComputeAll||mbCompute.at(lsize))
		{
			mEmitList[lsize]->ComputeReaction(mChem,mPath,mUseAbsorptionList.at(lsize),mAbsorptionValscm_2.at(lsize),mColFileName.at(lsize),mProfFileName.at(lsize),mExtraFileName.at(lsize));
			Log::mD<<"Spectrum integrated"<<endl;
			for(std::map<double,double>::iterator it=(mEmitList[lsize]->mEmitSpectrumIntegratedR).begin();it!=(mEmitList[lsize]->mEmitSpectrumIntegratedR).end();++it)
			{
				if(mEmitSpectrumIntegratedR.find(it->first)== mEmitSpectrumIntegratedR.end())
				{
					mEmitSpectrumIntegratedR[it->first]=it->second;
				}else
				{
					mEmitSpectrumIntegratedR[it->first]+=it->second;
				}
			}		
		}

	}
	Log::mD<<"End compmultipleemit"<<endl;
}

void Emission::ComputeEmissions(std::string vSuffix)
{
	if(!mpParam->Exists("/aero_main/emissions/use_emissions"))
		return;
	bool all=mpParam->Exists("/aero_main/emissions/use_all_emissions");

	// We search for the filenames

	vector<TiXmlNode*> femit=mpParam->GetNodes("/aero_main/emissions/emission");
	vector<TiXmlNode*>::iterator sit;
	for(sit=femit.begin();sit!=femit.end();++sit)
	{
		unsigned id=0;
		if(mpParam->KeyExists(*sit,"/","id"))
		{
			mpParam->GetNKey(*sit,"/","id",id);
			if(mpParam->Exists(*sit,"//column_file"))
			{
				mColFileName.at(id)=mpParam->Elem(*sit,"//column_file")+vSuffix;
			}
			if(mpParam->Exists(*sit,"//profile_file"))
			{
				mProfFileName.at(id)=mpParam->Elem(*sit,"//profile_file")+vSuffix;
			}
			if(mpParam->Exists(*sit,"//extra_file"))
			{
				mExtraFileName.at(id)=mpParam->Elem(*sit,"//extra_file")+vSuffix;
			}

			if(mpParam->Exists(*sit,"//use_absorption"))
			{
				mUseAbsorptionList.at(id)=true;
				vector< TiXmlNode* > myabs = mpParam->GetNodes(*sit,"//absorption");
				vector< TiXmlNode*>::iterator ait;
				for(ait=myabs.begin();ait!=myabs.end();++ait)
				{
					if( mpParam->KeyExists(*ait,"/","species") )
					{
						double lambnm = -1;
						if( mpParam->KeyExists(*ait,"/","lamb") )
						{
							mpParam->GetNKey(*ait,"/","lamb",lambnm);
						}
						string aspecies = mpParam->GetKey(*ait,"/","species");
						double value;
						mpParam->GetValue(*ait,"/",value);

						Log::mD<<"Absorption for the species "<<aspecies<<" at the wavelength "<<lambnm<<" (-1 for the total absorption for the process): "<< value<<" cm-1"<<endl;
						mAbsorptionValscm_2.at(id)[aspecies][lambnm]=value;

					}
				}
			}
			mbCompute.at(id)=true;
		}else
		{
			std::string name,state;

			name = mpParam->GetKey(*sit,"/","name");
			state = mpParam->GetKey(*sit,"/","state");
			std::map<double,double> freqbratio;
			vector< TiXmlNode* > myfreq = mpParam->GetNodes(*sit,"//freq");
			vector< TiXmlNode*>::iterator frt;
			for(frt=myfreq.begin();frt!=myfreq.end();++frt)
			{
				double frequency,bratio;
				mpParam->GetNKey(*frt,"/","freq",frequency);
				mpParam->GetNKey(*frt,"/","branching",bratio);
				freqbratio[frequency]=bratio;
			}
			id = mAbsorptionValscm_2.size();
			mEmitList.push_back(boost::shared_ptr<Emit>(new EmitAllowed(mAltGridKm,id,name,state,freqbratio)));
			mColFileName.push_back("");
			mProfFileName.push_back("");
			mExtraFileName.push_back("");
			if(mpParam->Exists(*sit,"//column_file"))
			{
				mColFileName.at(id)=mpParam->Elem(*sit,"//column_file")+vSuffix;
			}
			if(mpParam->Exists(*sit,"//profile_file"))
			{
				mProfFileName.at(id)=mpParam->Elem(*sit,"//profile_file")+vSuffix;
			}
			if(mpParam->Exists(*sit,"//extra_file"))
			{
				mExtraFileName.at(id)=mpParam->Elem(*sit,"//extra_file")+vSuffix;
			}		


			mUseAbsorptionList.push_back(false);
			mbCompute.push_back(false);
			std::map< std::string, std::map<double,double> > tmp;
			mAbsorptionValscm_2.push_back(tmp);

			if(mpParam->Exists(*sit,"//use_absorption"))
			{
				mUseAbsorptionList.at(id)=true;
				vector< TiXmlNode* > myabs = mpParam->GetNodes(*sit,"//absorption");
				vector< TiXmlNode*>::iterator ait;
				for(ait=myabs.begin();ait!=myabs.end();++ait)
				{
					if( mpParam->KeyExists(*ait,"/","species") )
					{
						double lambnm = -1;
						if( mpParam->KeyExists(*ait,"/","lamb") )
						{
							mpParam->GetNKey(*ait,"/","lamb",lambnm);
						}
						string aspecies = mpParam->GetKey(*ait,"/","species");
						double value;
						mpParam->GetValue(*ait,"/",value);

						Log::mD<<"Absorption for the species "<<aspecies<<" at the wavelength "<<lambnm<<" (-1 for the total absorption for the process): "<< value<<" cm-1"<<endl;
						mAbsorptionValscm_2.at(id)[aspecies][lambnm]=value;

					}
				}
			}


		}
	}
	// We compute the emissions


	for(size_t lsize=0;lsize<mEmitList.size();++lsize)
	{
		if(all||mbCompute.at(lsize))
		{
			mEmitList[lsize]->ComputeReaction(mChem,mPath,mUseAbsorptionList.at(lsize),mAbsorptionValscm_2.at(lsize),mColFileName.at(lsize),mProfFileName.at(lsize),mExtraFileName.at(lsize));
			for(std::map<double,double>::iterator it=(mEmitList[lsize]->mEmitSpectrumIntegratedR).begin();it!=(mEmitList[lsize]->mEmitSpectrumIntegratedR).end();++it)
			{
				if(mEmitSpectrumIntegratedR.find(it->first)== mEmitSpectrumIntegratedR.end())
				{
					mEmitSpectrumIntegratedR[it->first]=it->second;
				}else
				{
					mEmitSpectrumIntegratedR[it->first]+=it->second;
				}
			}		
		}

	}
	// When all the selected emissions are computed, we can call the entire spectrum function
	//

	if(mpParam->Exists("/aero_main/emissions/spectrum"))
	{
		string filename=mpParam->Elem("/aero_main/emissions/spectrum")+vSuffix;
		PrintEntireSpectrum(filename);
	}
	if(mpParam->Exists("/aero_main/emissions/emit_list"))
	{
		string filename=mpParam->Elem("/aero_main/emissions/emit_list")+vSuffix;
		PrintInfo(filename);
	}

}


bool Emission::CheckEmitList()
{
	bool resu=true;

	Log::mD<<"check emit list"<<endl;
	Log::mD<<"size of the list"<<mEmitList.size()<<endl;
	for(unsigned i=0;i<mEmitList.size();++i)
	{
		bool bon= (i==mEmitList[i]->GetId());
		if(!bon)
		{
			resu=false;
			//Log::SetPriority(Log::ERROR);
			Log::mE<<" Error in your emit reaction list. Position "<<i<<endl;
		}
	}
	return resu;
}



void Emission::PrintEntireSpectrum(std::string vFilename)
{
	
	if(FileExists(vFilename))
	{
		//Log::SetPriority(Log::DEBUGG);
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Syntetic vertical spectrum (w/o radiative transfer) "<<endl;
	of<<"# Wavelength in nm"<<endl;
	of<<"# Intensity in R"<<endl;

	for(std::map<double,double>::iterator it=mEmitSpectrumIntegratedR.begin();it!=mEmitSpectrumIntegratedR.end();++it)
	{
		of<<it->first<<" \t"<<it->second<<endl;
	}
	
	of<<Log::msMessageLog<<endl;
	of.close();

}


void Emission::PrintInfo(std::string vFilename)
{
	
	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Informations about the emissions"<<endl;
	for(unsigned i=0;i<mEmitList.size();++i)
	{
		of<<"Number: "<<i<<" : \t"<<mEmitList[i]->Info()<<endl;
	}
	
	of<<Log::msMessageLog<<endl;
	of.close();

}


ublas::vector<double> Emission::GetVER(SpecieId vSpec,int vPos)
{
	Log::mD<<"Number  of ver emissions: "<< mEmitList.size()<<endl;
	for(size_t i=0;i<mEmitList.size();++i)
	{
		Log::mD<<"Your species is : "<<vSpec.StandardName()<<endl;
		Log::mD<<"Check if your species is : "<<mEmitList[i]->GetSpecieId().StandardName()<<endl;
		if(mEmitList[i]->GetSpecieId()==vSpec)
		{
			return mEmitList[i]->GetVER(vPos);
		}
	}

	Error err("Emission::GetVER","Emitting species not found","The species "+vSpec.StandardName()+" was not found in the list of emitting species");
	throw err;

	ublas::vector<double> dummy;
	return dummy; // To avoid warning
}


ublas::vector<double> Emission::GetLimb(SpecieId vSpec,int vPos)
{
	for(size_t i=0;i<mEmitList.size();++i)
	{
		if(mEmitList[i]->GetSpecieId()==vSpec)
		{
			return mEmitList[i]->GetLimb(vPos);
		}
	}

	Error err("Emission::GetVER","Emitting species not found","The species "+vSpec.StandardName()+" was not found in the list of emitting species");
	throw err;

	ublas::vector<double> dummy;
	return dummy; // To avoid warning
}

