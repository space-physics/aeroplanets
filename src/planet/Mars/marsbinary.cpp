/**
 * \file marsbinary.cpp
 * \brief Implements the martian atmosphere from binary lecture
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: marsbinary.cpp 1111 2010-08-12 19:43:33Z gronoff $
 */



#include "marsbinary.hpp"
using namespace std;


DumpData::~DumpData()
{
	Log::mL<<"bye bye dumpdata"<<endl;
}


void DumpData::CopyVect(ublas::matrix<double>& rMat,unsigned line,const ublas::vector<double>& vVec,unsigned start,unsigned size)
{
	for(unsigned i=0;i<size;++i)
	{
		rMat(line,i)=vVec[start+i];
	}
}


DumpData::DumpData(std::string vFilename,unsigned vSize):mFilename(vFilename),mSizeOf(vSize),mNbAlt(0),mNbCol(0),mIAnnee(0),mIMois(0),mIJour(0),mIHeure(0),mIMinute(0),mSeconde(0.),mIPas(0),mF107(0.),mChiDeg(0.)
{
}

void DumpData::GetVal(std::ifstream&  rIfs,double& rTmp)
{

	rIfs.read(reinterpret_cast<char*>(&rTmp),mSizeOf);
}
void DumpData::GetVect(std::ifstream& rIfs,unsigned vSize,ublas::vector<double>& vVect)
{
	vVect.resize(vSize);
	for(unsigned i=0;i<vSize;++i)
	{
		rIfs.read(reinterpret_cast<char*>(&(vVect[i])),mSizeOf);
	}
}

void DumpData::ReadTableau()
{
	ifstream ifs(mFilename.c_str(),ios::in|ios::binary);

	if(!ifs)
	{
	//	Log::SetPriority(Log::ERROR);
		Log::mE<<"Impossible to open your file :"<<mFilename<<endl;
		Error err("DumpData::ReadTableau","Impossible to open the file","It is not possible to open your binary mars file...");
		throw err;

	}
	double tmp=0;

	GetVal(ifs,tmp);
	mNbAlt=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mNbCol=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mIAnnee=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mIMois=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mIJour=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mIHeure=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mIMinute=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mSeconde=tmp;
	GetVal(ifs,tmp);
	mIPas=static_cast<unsigned>(tmp);
	GetVal(ifs,tmp);
	mF107=tmp;
	GetVal(ifs,tmp);
	mChiDeg=tmp;
	Log::mL<<"nombre de points "<<mNbAlt<<endl;
	Log::mL<<"nombre de colonnes "<<mNbCol<<endl;

	ublas::vector<double> vide;
	GetVect(ifs,mNbCol*2-11,vide);

	mAltitudeKm.resize(mNbAlt);
	mDensitecm_3.resize(mNbAlt,4);
	mVz.resize(mNbAlt,4);
	mVy.resize(mNbAlt,5);
	mTemperatureK.resize(mNbAlt,5);
	mFluxQ.resize(mNbAlt,5);
	mNeutreDensitycm_3.resize(mNbAlt,5);
	mNeutreTemperatureK.resize(mNbAlt);
	mProdioncm_3s_1.resize(mNbAlt,4);
	mLossioncm_3s_1.resize(mNbAlt,4);
	mHeat.resize(mNbAlt);
	mChampB.resize(mNbAlt);
	mSupra.resize(mNbAlt,4);
	mUy.resize(mNbAlt);
	mUz.resize(mNbAlt);

	for(unsigned i=0;i<mNbAlt;++i)
	{
		ublas::vector<double> data;
		GetVect(ifs,mNbCol,data);
		//cout<<data.size()<<endl;
		mAltitudeKm[i]=(data[0]);
//		Log::mL<<"Altitude : "<<data[0]<<endl;
		CopyVect(mDensitecm_3,i,data,1,4);
		CopyVect(mVz,i,data,5,4);
		CopyVect(mTemperatureK,i,data,14,5);
		CopyVect(mFluxQ,i,data,19,5);
		CopyVect(mNeutreDensitycm_3,i,data,24,5);
		mNeutreTemperatureK[i]=data[29];
		CopyVect(mProdioncm_3s_1,i,data,30,4);
		CopyVect(mLossioncm_3s_1,i,data,34,4);
		mHeat[i]=data[38];
		mChampB[i]=data[39];
		CopyVect(mSupra,i,data,40,4);
		mUy[i]=data[50];
		mUz[i]=data[51];
	}

	mDensitecm_3/=1E6;
	mNeutreDensitycm_3/=1E6;
	Log::mL<<"Fin lecture"<<endl;
}



std::deque< ublas::vector<double> > DumpData::ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm)
{
	Log::mL<<"neutral atmo"<<endl;
	//vector<double> alt=UblasToStd(mAltitudeKm);
	ublas::vector<double> tmpsp(mNbAlt);
	tmpsp.clear();
	deque< ublas::vector<double> > atmo_resu(5,mAltitudeKm);

	for(unsigned i=0;i<5;++i)
	{
		for(unsigned k=0;k<mNbAlt;++k)
		{
			tmpsp[k]=mNeutreDensitycm_3(k,i);
			if(tmpsp[k]<1E-42)
				tmpsp[k]=1E-42;
		}
		Log::mL<<"position = "<<i<<endl;
		atmo_resu[i]=MathFunction::IntLog(mAltitudeKm,tmpsp,vAltitudeGridKm);
	}
	return atmo_resu;
}
std::deque< ublas::vector<double> > DumpData::ReturnIono(const ublas::vector<double>& vAltitudeGridKm)
{
	Log::mL<<"return iono"<<endl;
//	vector<double> alt=UblasToStd(mAltitudeKm);
	ublas::vector<double> tmpsp(mNbAlt);
	tmpsp.clear();
	deque< ublas::vector<double> > atmo_resu(4,mAltitudeKm);
	for(unsigned i=0;i<4;++i)
	{
		for(unsigned k=0;k<mNbAlt;++k)
		{
			tmpsp[k]=mDensitecm_3(k,i);
//			Log::mL<<tmpsp[k]<<endl;
			if(tmpsp[k]<1E-42)
				tmpsp[k]=1E-42;
		}
		atmo_resu[i]=MathFunction::IntLog(mAltitudeKm,tmpsp,vAltitudeGridKm);
	}
	return atmo_resu;
}


ublas::vector<double> DumpData::ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm)
{

	Log::mL<<"ntemp"<<endl;
	//vector<double> alt=UblasToStd(mAltitudeKm);
/*	vector<double> tmpsp(mNbAlt,0.);
	for(unsigned k=0;k<mNbAlt;++k)
	{
		tmpsp[k]=mNeutreT(k);
	}
	*/
	//vector<double> tmpsp=UblasToStd(mNeutreTemperatureK);
	return MathFunction::IntLog(mAltitudeKm,mNeutreTemperatureK,vAltitudeGridKm);
}

ublas::vector<double> DumpData::ReturnITemp(const ublas::vector<double>& vAltitudeGridKm)
{

	//vector<double> alt=UblasToStd(mAltitudeKm);
	ublas::vector<double> tmpsp(mNbAlt);
	tmpsp.clear();
	for(unsigned k=0;k<mNbAlt;++k)
	{
		tmpsp[k]=mTemperatureK(k,2);
	}
	
	return MathFunction::IntLog(mAltitudeKm,tmpsp,vAltitudeGridKm);
}

ublas::vector<double> DumpData::ReturnETemp(const ublas::vector<double>& vAltitudeGridKm)
{

	//vector<double> alt=UblasToStd(mAltitudeKm);
	ublas::vector<double> tmpsp(mNbAlt);
	tmpsp.clear();
	for(unsigned k=0;k<mNbAlt;++k)
	{
		tmpsp[k]=mTemperatureK(k,4);
	}
	
	return MathFunction::IntLog(mAltitudeKm,tmpsp,vAltitudeGridKm);
}

