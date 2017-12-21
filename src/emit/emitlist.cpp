/**
 * \file emitlist.cpp
 * \brief Implements the different emissions 
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emitlist.cpp 1491 2012-05-07 21:35:29Z gronoff $
 */
#include "emitlist.hpp"
using namespace std;

typedef ublas::vector<double> Vec; 
typedef ublas::matrix<double> Mat;
typedef ublas::matrix_row< Mat > Row;
typedef std::deque< std::string > StrL;


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-------------------------+
 |                                        |
 |The Vegard Kaplan emissions ( N2 A3Su+ )|
 |                                        |
 +---------------------------------------*/


std::string EmitN2A3S::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ N2(A3Su+) VK\n";
	str+=" Computation of the emissions of the Vegard Kaplan bands\n";
	str+=" This reactions also set up the N2(A3Su+) density used in O(1S) computing \n";
	str+=" When this reaction is called, it fills the N2(A3S_ density in the chem list\n It allows for the other emissions to use that density (eg O(1S) emissions)\n\n";
	return str;
}
EmitN2A3S::EmitN2A3S(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,0,"N2","A3S")
{}

void EmitN2A3S::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the density of N2(A3Su+)"<<endl;
	Vec n2a3sdens=vChem->GetDensity(SpecieId("N2","A3S"),mExtraInfo,mExtraInfoName,mWarnings);
//	Log::mL<<"n2a3sdens : "<<n2a3sdens<<endl;

//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
	if(mColEmitphcm_3s_1.size()==0)
		return;

	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}

//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);


//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of N2(A3S) emission"<<endl;



}



void EmitN2A3S::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Extra informations about N2A3S"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;



	}
}





/*       _\|/_
         (o o)
 +----oOO-{_}-OOo----+
 |                   |
 |The O(1S) Emissions|
 |                   |
 +------------------*/





std::string EmitO1S::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ O(1S)\n";
	str+=" Computation of the emissions of O(1S) - 557.7 nm and 297.2 nm\n";
	str+=" The reaction uses the Frederick/Kopp O2+ + N -> O(1S) + NO+ \n";
	str+=" It is negligible for the dayside, but must be checked for the nightside\n\n";

	str+=" When this reaction is called, it fills the O1S density in the chem list\n It allows for the other emissions to use that density (eg O1D emissions)\n\n";
	return str;
}
EmitO1S::EmitO1S(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,1,"O","1S")
{}

void EmitO1S::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	mbUseAbsorption=vbUseAbsorption;

		 // This function fills basically mEmitSpectrumcm_3s_1 and mColEmitphcm_3s_1
		 // then, it calls EmitSpectrumIntegrate to fill mIntegratedColEmitR, mEmitSpectrumIntegratedR
		 // and finally, it calls EmitSpectrumToLimb


	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the density of O(1S)"<<endl;
	Vec o1sdens=vChem->GetDensity(SpecieId("O","1S"),mExtraInfo,mExtraInfoName,mWarnings);

//	Log::mL<<"o1sdens : "<<o1sdens<<endl;

// Done inside GetDensity Now!
//	Log::mL<<"Put dens "<<endl;
//	vChem->PutDens("O","1S",o1sdens);
	

//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,4); // tot prod already stored here
	if(mColEmitphcm_3s_1.size()==0)
		return;
	mEmitSpectrumcm_3s_1[557.734]=Row(mExtraInfo,5); // tot prod already stored here
	mEmitSpectrumcm_3s_1[297.229]=Row(mExtraInfo,6); // tot prod already stored here
	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}
	//mColumnAbsorption_cm[297.229]=vChem->GetDens("CO2","")*1E-25;



//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of O1S emission"<<endl;


}

void EmitO1S::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Extra informations about O1S"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;



	}
}


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-----------+
 |                          |
 |The N(2D) emissions 520 nm|
 |                          |
 +-------------------------*/





std::string EmitND::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ N(2D) 520nm \n";
	str+=" Computation of the emissions of the 520 nm\n";
	str+=" This reactions also set up the N(2D) photoequilibrium \n";
	str+=" The longlife of the N(2D) state can be a problem for the validity of the equilibrium.\n";
	str+="The validity of the reaction N(2D)+O2 -> O(1D) is questioned in Link 91: it may not exists.";
	return str;
}
EmitND::EmitND(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,2,"N","2D")
{}

void EmitND::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the density of N(2D)"<<endl;
	Vec nddens=vChem->GetDensity(SpecieId("N","2D"),mExtraInfo,mExtraInfoName,mWarnings);
//	Log::mL<<"n2dens : "<<nddens<<endl;

//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
	if(mColEmitphcm_3s_1.size()==0)
		return;

	mEmitSpectrumcm_3s_1[520]=Row(mExtraInfo,0); // tot prod already stored here
	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}

//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}
	Log::mD<<"End of N(2D) emission"<<endl;
}



void EmitND::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Extra informations about N(2D)"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;



	}
}

/*       _\|/_
         (o o)
 +----oOO-{_}-OOo------------------+
 |                                 |
 |O(1D) productions : 630nm triplet|
 |                                 |
 +--------------------------------*/

std::string EmitO1D::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ O(1D) 630nm \n";
	str+=" Computation of the emissions of the 630 nm triplet\n";
	str+=" This reactions also set up the O(1D) photoequilibrium \n";
	return str;
}


EmitO1D::EmitO1D(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,3,"O","1D")
{}

void EmitO1D::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{

	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the density of O(1D)"<<endl;
	Vec o1ddens=vChem->GetDensity(SpecieId("O","1D"),mExtraInfo,mExtraInfoName,mWarnings);

//	Log::mL<<"o1ddens : "<<o1ddens<<endl;
	

//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here

	if(mColEmitphcm_3s_1.size()==0)
		return;
	mEmitSpectrumcm_3s_1[630.0304]=Row(mExtraInfo,1); // tot prod already stored here
	mEmitSpectrumcm_3s_1[636.3776]=Row(mExtraInfo,2); // tot prod already stored here
	mEmitSpectrumcm_3s_1[639.1733]=Row(mExtraInfo,3); // tot prod already stored here
	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}


//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of O1D emission"<<endl;


}

void EmitO1D::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Extra informations about O1D"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;
	}
}




/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-------------------+
 |                                  |
 |CO(a3Pi) emissions : cameron bands|
 |                                  |
 +---------------------------------*/



std::string EmitCOa3Pi::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ CO(a3Pi) bands around 200nm \n";
	return str;
}


EmitCOa3Pi::EmitCOa3Pi(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,4,"CO","a3Pi")
{}

void EmitCOa3Pi::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{

	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the density of CO(a3Pi)"<<endl;
	Vec coa3pidens=vChem->GetDensity(SpecieId("CO","a3Pi"),mExtraInfo,mExtraInfoName,mWarnings);

//	Log::mL<<"[ CO(a3Pi) ] : "<<coa3pidens<<endl;
	

//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[630.0304]=Row(mExtraInfo,1); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[636.3776]=Row(mExtraInfo,2); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[639.1733]=Row(mExtraInfo,3); // tot prod already stored here

	if(mColEmitphcm_3s_1.size()==0)
		return;
	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}

//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of CO(a3Pi) emission"<<endl;


}

void EmitCOa3Pi::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Extra informations about CO(a3Pi)"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;
	}
}



/*       _\|/_
         (o o)
 +----oOO-{_}-OOo---------------+
 |                              |
 |CO2+(B-X) doublet around 289nm|
 |                              |
 +-----------------------------*/



std::string EmitCO2pB::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ CO2+(B) doublet at 289 nm \n";
	return str;
}


EmitCO2pB::EmitCO2pB(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,5,"CO2+","B")
{}

void EmitCO2pB::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{

	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the production of CO2+(B)"<<endl;
//	Vec coa3pidens=vChem->GetDensity(SpecieId("CO","a3Pi"),mExtraInfo,mExtraInfoName,mWarnings);
//	Log::mL<<"[ CO(a3Pi) ] : "<<coa3pidens<<endl;
//	Log::mL<<"Density info size "<<mExtraInfo.size1()<<" "<<mExtraInfo.size2()<<endl;
//	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[630.0304]=Row(mExtraInfo,1); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[636.3776]=Row(mExtraInfo,2); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[639.1733]=Row(mExtraInfo,3); // tot prod already stored here

	mColEmitphcm_3s_1=vChem->GetTotProd("CO2+","B");

	if(mColEmitphcm_3s_1.size()==0)
		return;
//	Log::mL<<" Emissions : "<<mColEmitphcm_3s_1<<endl;


	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}

//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of CO2+(B) emission"<<endl;


}

void EmitCO2pB::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
/*
	of<<"# Extra informations about CO(a3Pi)"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;
	}
	*/
}




/*       _\|/_
         (o o)
 +----oOO-{_}-OOo---------------+
 |                              |
 |CH(A2D) band around 420-440nm |
 |                              |
 +-----------------------------*/



std::string EmitCHA2D::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ CH(A2D) band around 420-440 nm \n";
	return str;
}


EmitCHA2D::EmitCHA2D(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,6,"CH","A2D")
{}

void EmitCHA2D::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{

	mbUseAbsorption=vbUseAbsorption;

	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
//	Log::mL<<"Compute the production of CH(A2D)"<<endl;

	mColEmitphcm_3s_1=vChem->GetTotProd("CH","A2D");

	if(mColEmitphcm_3s_1.size()==0)
		return;
//	Log::mL<<" Emissions : "<<mColEmitphcm_3s_1<<endl;

	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}


//	Log::mL<<"emit spectrum to integrate"<<endl;
	EmitSpectrumIntegrate();
//	Log::mL<<"emit spectrum to limb"<<endl;
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;

	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}


	Log::mD<<"End of CH(A2D) emission"<<endl;


}

void EmitCHA2D::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
}










/*        _\|/_
          (o o)
 +-----oOO-{_}-OOo-----+
 |                     |
 |The O+(2P) Emissions |
 |                     |
 +---------------------*/





std::string EmitOplus2P::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ O+(2P)\n";
	str+=" Computation of the emissions of O+(2P) - 732.0 733.0 and 247.0 nm\n";
	str+=" When this reaction is called, it fills the O+(2P) density in the chem list\n It allows for the other emissions to use that density (eg O+(2D) emissions)\n\n";
	return str;
}
EmitOplus2P::EmitOplus2P(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,7,"O+","2P")
{}

void EmitOplus2P::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	Log::mI<<"Work on O+(2P)"<<endl;
	mbUseAbsorption=vbUseAbsorption;

		 // This function fills basically mEmitSpectrumcm_3s_1 and mColEmitphcm_3s_1
		 // then, it calls EmitSpectrumIntegrate to fill mIntegratedColEmitR, mEmitSpectrumIntegratedR
		 // and finally, it calls EmitSpectrumToLimb


	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
	Log::mI<<"get density of O+(2P)"<<endl;
	Vec o1sdens=vChem->GetDensity(SpecieId("O+","2P"),mExtraInfo,mExtraInfoName,mWarnings);
	Log::mI<<"density retrieved"<<endl;
	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
	if(mColEmitphcm_3s_1.size()==0)
		return;
	mEmitSpectrumcm_3s_1[732.0]=Row(mExtraInfo,2); // tot prod already stored here
	mEmitSpectrumcm_3s_1[733.0]=Row(mExtraInfo,3); // tot prod already stored here
	mEmitSpectrumcm_3s_1[247.0]=Row(mExtraInfo,4); // tot prod already stored here

	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}
	EmitSpectrumIntegrate();
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;
	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}
	Log::mD<<"End of O+(2P) emission"<<endl;


}

void EmitOplus2P::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);

	of<<"# Extra informations about O+(2P)"<<endl;
	of<<"# Productions cm_3s_1"<<endl;
	of<<"# Altitudes km "<<endl;
	for(StrL::iterator it=mExtraInfoName.begin();it!=mExtraInfoName.end();++it)
	{
		of<<"# "<<*it<<endl;
	}

	for(unsigned i=0;i<mAltGridKm.size();++i)
	{
		of<<mAltGridKm[i]<<"\t";
		for(unsigned j=0;j<mExtraInfo.size1();++j)
		{
			of<<mExtraInfo(j,i)<<"\t";
		}
		of<<endl;



	}
}



/*        _\|/_
          (o o)
 +-----oOO-{_}-OOo--------+
 |                        |
 |The O((3p)3P) Emissions |
 |                        |
 +------------------------*/





std::string EmitO3p3P::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ O((3p)3P)\n";
	return str;
}
EmitO3p3P::EmitO3p3P(ublas::vector<double> vAltGridKm):Emit(vAltGridKm,8,"O","(3p)3P")
{}

void EmitO3p3P::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	Log::mI<<"Work on O((3p)3P)"<<endl;
	mbUseAbsorption=vbUseAbsorption;

		 // This function fills basically mEmitSpectrumcm_3s_1 and mColEmitphcm_3s_1
		 // then, it calls EmitSpectrumIntegrate to fill mIntegratedColEmitR, mEmitSpectrumIntegratedR
		 // and finally, it calls EmitSpectrumToLimb


	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
	Log::mI<<"get density of O((3p)3P)"<<endl;
	mColEmitphcm_3s_1 = vChem->GetTotProd("O","(3p)3P");
	mEmitSpectrumcm_3s_1[844.6]=mColEmitphcm_3s_1; // tot prod already stored here


	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}
	EmitSpectrumIntegrate();
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;
	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}
	Log::mD<<"End of O((3p)3P) emission"<<endl;


}

void EmitO3p3P::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
}

/*        _\|/_
          (o o)
 +-----oOO-{_}-OOo--------+
 |                        |
 |  The Allowed Emissions |
 |                        |
 +------------------------*/





std::string EmitAllowed::Info()
{
	string str="============= Emission nb "+ntostr(mId)+" =============\n";
	str+="================================================ Allowed\n";
	str+="This emission is user-defined: the user asks for a species, and computes directly its emission";
	return str;
}
EmitAllowed::EmitAllowed(ublas::vector<double> vAltGridKm, unsigned vId, std::string vName, std::string vState, std::map<double,double> vFreqBratio):Emit(vAltGridKm,vId,vName,vState)
{
	mFreqBratio = vFreqBratio;
}

void EmitAllowed::ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName,std::string vProfileFileName,std::string vExtraInfoFilename)
{
	Log::mI<<"Work on "<<mEmittingSpecie.StandardName()<<endl;
	mbUseAbsorption=vbUseAbsorption;

		 // This function fills basically mEmitSpectrumcm_3s_1 and mColEmitphcm_3s_1
		 // then, it calls EmitSpectrumIntegrate to fill mIntegratedColEmitR, mEmitSpectrumIntegratedR
		 // and finally, it calls EmitSpectrumToLimb


	mColumnFileName=trim(vColumnFileName);
	mProfileFileName=trim(vProfileFileName);
	mExtraInfoFilename=trim(vExtraInfoFilename);
	mColEmitphcm_3s_1 = vChem->GetTotProd(mEmittingSpecie.mName,mEmittingSpecie.mState);

	if(mColEmitphcm_3s_1.size()!=0)
	{
	//mEmitSpectrumcm_3s_1[844.6]=mColEmitphcm_3s_1; // tot prod already stored here


	for(std::map<double,double>::iterator it = mFreqBratio.begin();it!=mFreqBratio.end();++it)
	{
		mEmitSpectrumcm_3s_1[it->first] = mColEmitphcm_3s_1*it->second;
	}

	//	mColEmitphcm_3s_1=Row(mExtraInfo,0); // tot prod already stored here
//	if(mColEmitphcm_3s_1.size()==0)
//		return;
//	mEmitSpectrumcm_3s_1[732.0]=Row(mExtraInfo,2); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[733.0]=Row(mExtraInfo,3); // tot prod already stored here
//	mEmitSpectrumcm_3s_1[247.0]=Row(mExtraInfo,4); // tot prod already stored here

	if(mbUseAbsorption)
	{
		ComputeAbsorption(vAbsorptions,vChem);
	}
	EmitSpectrumIntegrate();
	EmitSpectrumToLimb(vPath);

//	Log::mL<<"extra info check"<<endl;
	if(trim(vExtraInfoFilename)!="")
	{
		mExtraInfoFilename=vExtraInfoFilename;
		WriteExtraInfo();
	}
	}else{
		Log::mW<<"Your species :"<<mEmittingSpecie.StandardName()<<"Does not seems to be produced, and therefore, its emission are not computed. Check if the species is created in the cross section files."<<endl;
	}
}

void EmitAllowed::WriteExtraInfo()
{
	if(FileExists(mExtraInfoFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<mExtraInfoFilename<<endl;
	}
	ofstream of(mExtraInfoFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
}



