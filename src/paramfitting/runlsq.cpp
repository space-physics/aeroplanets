#include "runlsq.hpp"
using namespace std;
namespace bub = boost::numeric::ublas;
std::deque< ublas::vector<double> > RunAtmoLsq::RetrieveData(double vTmpMult)
{
	double multiplicator=1.;
	if(mbUseCalib)
	{
		multiplicator = vTmpMult;
	}
	mpAtmo->MCompute();
	//     RetrieveData
	//
	std::deque< ublas::vector<double> > resu; //SpecieId spec(mAdjustSpecie,"");
	for(size_t i = 0; i < mAdjustSpecieId.size(); ++i)
	{
		ublas::vector<double> data = mpAtmo->RetrieveObservable(mAdjustSpecieId[i], mAdjustParams[i], mXinit[i]);
		resu.push_back(1/multiplicator*data);
	}
	return resu;
}

void RunAtmoLsq::PrintConvergence(std::string vFilename)
{
	if(FileExists(vFilename))
	{
		Log::mD<<"Low level warning : We overwrite"<<vFilename<<endl;
	}
	ofstream of(vFilename.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Parameter convergence"<<endl;
	for(size_t i = 0; i<mStepsPars.size();++i)
	{

		of<<i<<"\t";

		for(size_t j=0;j<mStepsPars.at(i).size();++j)
		{
			for(size_t k=0;k<mStepsPars.at(i).at(j).size();++k)
			{
				of<<mStepsPars.at(i).at(j).at(k)<<"\t";

			}
		}
		of<<endl;
	}
	of.close();
}

int RunAtmoLsq::Exec(int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag)
{
	assert(vDiffsize>vParamsize);
	++mNbCalls;
	if(vFlag==0)
	{
		Log::mD<<"Oups, vFlag is equal to 0"<<endl;
		return 1;
	}
	// loop on species
	unsigned pos=0;// pos in the parameters
	std::deque< std::deque<double> > tmppars;
	for(size_t i=0;i<mSpecieName.size();++i)
	{
		std::deque<double> param;
		Log::mI<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		switch(mSpecieModel[i])
		{
			case 0:
				for(size_t k=0;k<4;++k)
				{
					param.push_back(exp(pParam[pos]));
					++pos;
				}
				Log::mI<<mSpecieName[i]<<endl<<"Dens 0 \t"<<param.at(0)<<endl<<"Texo  \t"<<param.at(1)<<endl<<"T0 \t"<<param.at(2)<<endl<<"Shape  \t"<<param.at(3)<<endl;
				break;
			case 1:
				for(size_t k=0;k<3;++k)
				{
					if(k==1)
						param.push_back(exp(pParam[pos]));
					else
						param.push_back(pParam[pos]);
					++pos;
				}
				Log::mI<<mSpecieName[i]<<endl<<"Alt 0 \t"<<param.at(0)<<endl<<"Dens  \t"<<param.at(1)<<endl<<"SH \t"<<param.at(2)<<endl;
				break;
			case 2:
				for(size_t k=0;k<3;++k)
				{
					if(k==1)
						param.push_back(exp(pParam[pos]));
					else
						param.push_back(pParam[pos]);
					++pos;
				}
				Log::mI<<mSpecieName[i]<<endl<<"Alt 0 \t"<<param.at(0)<<endl<<"Dens  \t"<<param.at(1)<<endl<<"Dev \t"<<param.at(2)<<endl;
				break;
			case 3:
				for(size_t k=0;k<2;++k)
				{
					param.push_back(exp(pParam[pos]));
					++pos;
				}
				Log::mI<<mSpecieName[i]<<endl<<"Dens 0 \t"<<param.at(0)<<endl<<"Texo  \t"<<param.at(1)<<endl;
				break;
			case 4:
				Log::mI<<mSpecieName[i]<<endl<<"Altitude \t Value"<<endl;
				for(size_t k=0;k<mAdditionalParameters[i][0];++k)
				{
					param.push_back(pParam[pos]);
					Log::mI<<mAdditionalParameters[i][k+1]<<"\t"<<param[k]<<endl;
					++pos;
				}

				break;
			case 5:
				{
					for(size_t k=0;k<4;++k)
					{
						param.push_back(exp(pParam[pos]));
						++pos;
					}
					Log::mI<<mSpecieName[i]<<endl<<"Dens Thermosphere \t"<<param.at(0)<<endl
						<<"Dens Mesosphere \t"<<param.at(1)<<endl
						<<"Temp Exosphere  \t"<<param.at(2)<<endl
						<<"Temp Mesosphere  \t"<<param.at(3)<<endl;


				}
				break;
			case 6:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(pParam[pos]));
						++pos;
					}
					param.push_back(acos(pParam[pos]));
					++pos;
					Log::mI<<mSpecieName[i]<<endl<<"Alt 0 \t"<<param.at(0)<<endl
						<<"Dens  \t"<<param.at(1)<<endl
						<<"Temp Exosphere  \t"<<param.at(2)<<endl
						<<"SZA  \t"<<param.at(3)<<endl;
				}
				break;
			case 7:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(pParam[pos]));
						++pos;
					}
					Log::mI<<mSpecieName[i]<<endl<<"Alt 0 \t"<<param.at(0)<<endl
						<<"Dens  \t"<<param.at(1)<<endl
						<<"Temp Exosphere  \t"<<param.at(2)<<endl;
				}
				break;
			case 8:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(pParam[pos]));
						++pos;
					}
					Log::mI<<mSpecieName[i]<<endl<<"Alt 0 \t"<<param.at(0)<<endl
						<<"Dens  \t"<<param.at(1)<<endl
						<<"Temp Exosphere  \t"<<param.at(2)<<endl;
				}
				break;
			case 9:
				{
					for(size_t k=0;k<4;++k)
					{
						param.push_back(exp(pParam[pos]));
						++pos;
					}
					param.push_back(pParam[pos]);
					++pos;
					Log::mI<<mSpecieName[i]<<endl<<"Dens Thermosphere \t"<<param.at(0)<<endl
						<<"Dens Mesosphere \t"<<param.at(1)<<endl
						<<"Temp Exosphere  \t"<<param.at(2)<<endl
						<<"Temp Mesosphere  \t"<<param.at(3)<<endl
						<<"Exccentric C     \t"<<param.at(4)<<endl;


				}
				break;
			default:
				Error err("Error in the adjustement model","Model does not exists","The model number "+ntostr(mSpecieModel[i])+" does not exists - ERROR");
				throw err;
		}
		Log::mI<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		//     reset atmo & elec dens
		if(mSpecieName[i]=="e")
		{
			mpAtmo->ResetElectronDens(mSpecieModel[i],param,mAdditionalParameters[i]);
		}else
		{
			mpAtmo->ResetSpecieDens(mSpecieName[i],mSpecieModel[i],param,mAdditionalParameters[i]);
		}
		tmppars.push_back(param);
	}


	if(mbIsElecPrecipAdjusted)
	{
		Log::mI<<"---------------------- Electron precipitation adjustement"<<endl;
		std::deque<double> param;
		switch(mElecPrecipModel)
		{
			case 1:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(pParam[pos]);
						++pos;
					}
					Log::mI<<"Electron precipitation: Maxwell \nE tot \t"<<param.at(0)<< endl <<"E0 \t"<<param.at(1)<<endl;
				}
				break;
			case 2:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(pParam[pos]);
						++pos;
					}
					Log::mI<<"Electron precipitation: Gaussian \nE tot \t"<<param.at(0)<< endl <<"E0 \t"<<param.at(1)<<endl;
				}
				break;
			case 3:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(pParam[pos]);
						++pos;
					}
					Log::mI<<"Electron precipitation: Dirac \nE tot \t"<<param.at(0)<< endl <<"E0 \t"<<param.at(1)<<endl;
				}
				break;
			case 4:
				{
					for(size_t k = 0; k < mElecPrecipAdditionalParameters.size(); ++k)
					{
						param.push_back(pParam[pos]);
						Log::mI<<mElecPrecipAdditionalParameters[k]<<"\t"<<param.at(k)<<endl;
						++pos;
					}
				}
				break;

			case 0: // passthru
			default:
				Log::mI<<"Null precipitation"<<endl;


		}
		mpAtmo->ResetElecPrecip(mElecPrecipModel, param, mElecPrecipAdditionalParameters);
		tmppars.push_back(param); // We add the electron precipitation parameters in the step saving
	}


	double multiplicator=1.;
	if(mbUseCalib)
	{
		std::deque<double> tmpm;
		multiplicator=pParam[pos];
		Log::mI<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl<<"\t\t Multiplicator :  "<<multiplicator<<endl<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		++pos;
		tmpm.push_back(multiplicator);
		tmppars.push_back(tmpm);
	}
	mStepsPars.push_back(tmppars);
	//
	assert(pos==static_cast<unsigned>(vParamsize));
	//     Launch MCompute

	mpAtmo->MCompute();
	//     RetrieveData
// SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
	
	std::deque< ublas::vector<double> > datas;
	size_t datasize = 0;
	for(size_t j =0; j < mAdjustSpecieId.size(); ++j)
	{
		ublas::vector<double> data=mpAtmo->RetrieveObservable(mAdjustSpecieId[j], mAdjustParams[j], mXinit[j]);
		datasize += data.size();
		datas.push_back(data);
	}
	assert(datasize == static_cast<size_t>(vDiffsize));
	double tmp=0;
	double freedeg = static_cast<double>((vDiffsize - vParamsize));
	Log::mL<<"Diff option:"<<mDiffOption<<endl;
	switch(mDiffOption)
	{
		case 4: 
			{
				unsigned k = 0;
				Log::mL<<"Diff option 4: Standard convergence with no knowledge of error"<<endl;
				for(size_t j =0; j < datas.size(); ++j)
				{
					for(size_t i = 0; i < datas[j].size(); ++i)
					{
						pDiff[k]= (multiplicator * mYmeasu[j][i] - datas[j][i]) / sqrt(mTol[j] * mTol[j] + multiplicator * mYmeasu[j][i]);
						double tm = (multiplicator*mYmeasu[j][i] - datas[j][i]) / (multiplicator * mErrors[j][i]);
						tmp += tm*tm;
						Log::mI<<"pDiff["<<k<<"] = "<<pDiff[k]<<"\t measu "<<mYmeasu[j][i]<<"\t simu "<<datas[j][i]<<"\t measu2 "<<multiplicator*mYmeasu[j][i]<<"\t Errors: "<<multiplicator * mErrors[j][i]<<"\t chi2v "<<tm*tm/freedeg << "(" <<  1./freedeg << ")"<<endl;
						k += 1;
					}
				}
			}
			break;

		case 3: 
			{
				unsigned k = 0;
				Log::mL<<"Diff option 3: Standard convergence function"<<endl;
				Log::mI<<" Size of datas: "<<datas.size()<<endl;
				for(size_t j =0; j < datas.size(); ++j)
				{
					for(size_t i = 0; i < datas[j].size(); ++i)
					{
						pDiff[k]= (multiplicator * mYmeasu[j][i] - datas[j][i]) / (mTol[j] * mErrors[j][i]);
						double tm = (multiplicator * mYmeasu[j][i] - datas[j][i]) / (multiplicator * mErrors[j][i]);
						tmp += tm*tm;
						Log::mI<<"pDiff["<<k<<"] = "<<pDiff[k]<<"\t measu "<<mYmeasu[j][i]<<"\t simu "<<datas[j][i]<<"\t measu2 "<<multiplicator*mYmeasu[j][i]<<"\t Errors: "<<multiplicator * mErrors[j][i]<<"\t chi2v "<<tm*tm/freedeg << "(" <<  1./freedeg << ")"<<endl;
						k += 1;
					}
				}
			}
			break;

		case 2: // A percentage minimization; taking into account the errors... 
			// The best in that case?
			{
				unsigned k = 0;
				Log::mL<<"Diff option 2"<<endl;
				for(size_t j =0; j < datas.size(); ++j)
				{
					for(size_t i = 0; i < datas[j].size(); ++i)
					{
						pDiff[k] = (multiplicator * mYmeasu[j][i] - datas[j][i]) / (mTol[j]*mErrors[j][i] * (multiplicator * mYmeasu[j][i] + datas[j][i]));
						double tm = (multiplicator*mYmeasu[j][i] - datas[j][i]) / (multiplicator * mErrors[j][i]);
						tmp += tm*tm;
						Log::mI<<"pDiff["<<k<<"] = "<<pDiff[k]<<"\t measu "<<mYmeasu[j][i]<<"\t simu "<<datas[j][i]<<"\t measu2 "<<multiplicator*mYmeasu[j][i]<<"\t Errors: "<<multiplicator * mErrors[j][i]<<"\t chi2v "<<tm*tm/freedeg << "(" <<  1./freedeg << ")"<<endl;
						k += 1;
					}
				}
			}
			break;

		case 1: // A chi**2 minimization...
			// Theoretically, you shouldn't do it!!
			{
				unsigned k = 0;
				Log::mL<<"Diff option 1: Chi2 (bad!)"<<endl;
				for(size_t j =0; j < datas.size(); ++j)
				{
					for(size_t i = 0; i < datas[j].size(); ++i)
					{
						pDiff[k]=(multiplicator * mYmeasu[j][i] - datas[j][i]) / (mTol[j] * mErrors[j][i]);
						pDiff[k] *= pDiff[i];
						double tm = (multiplicator * mYmeasu[j][i] - datas[j][i]) / (multiplicator * mErrors[j][i]);
						tmp += tm*tm;
						Log::mI<<"pDiff["<<k<<"] = "<<pDiff[k]<<"\t measu "<<mYmeasu[j][i]<<"\t simu "<<datas[j][i]<<"\t measu2 "<<multiplicator*mYmeasu[j][i]<<"\t Errors: "<<multiplicator * mErrors[j][i]<<"\t chi2v "<<tm*tm/freedeg << "(" <<  1./freedeg << ")"<<endl;
						k += 1;
					}
				}
			}
			break;
		case 0: //passthrough
		default:
			{
				Log::mL<<"Diff option 0: default"<<endl;
				unsigned k = 0;
				for(size_t j =0; j < datas.size(); ++j)
				{
					for(size_t i = 0; i < datas[j].size(); ++i)
					{
						pDiff[k]=(log(multiplicator * mYmeasu[j][i]) - log(datas[j][i])) / (mTol[j] * mErrors[j][i]);
						double tm = (multiplicator * mYmeasu[j][i] - datas[j][i]) / (multiplicator * mErrors[j][i]);
						tmp += tm*tm;
						Log::mI<<"pDiff["<<k<<"] = "<<pDiff[k]<<"\t measu "<<mYmeasu[j][i]<<"\t simu "<<datas[j][i]<<"\t measu2 "<<multiplicator*mYmeasu[j][i]<<"\t Errors: "<<multiplicator * mErrors[j][i]<<"\t chi2v "<<tm*tm/freedeg << "(" <<  1./freedeg << ")"<<endl;
						k += 1;
					}
				}
			}
	}
		
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

	// Computation of the Chi Square v for the Chi square test
	mChiSquarev = tmp / freedeg;
	Log::mI<<"#  CHI2 reduced = "<<mChiSquarev<<endl;
	

	string postfix="_call_"+ntostr(mNbCalls)+".dat";
	if(mbPrintIonosphere)
	{
		mpAtmo->PrintIonosphere(mPrintIonoFile+postfix);
	}
	if(mbPrintThermosphere)
	{
		mpAtmo->PrintNeutralAtmo(mPrintThermoFile+postfix);
	}
	if(mbPrintDiffplot)
	{
		string filename=mPrintDiffFile+postfix;
		ofstream oaf(filename.c_str());
		oaf.precision(9);
		oaf.setf(ios::scientific);
		oaf<<"# Data Model comparison"<<endl;
		oaf<<"# Alt in Km"<<endl<<"# Prod in cm-3s-1"<<endl;
		oaf<<"# Measurements"<<endl;
		oaf<<"# Model"<<endl;
////////////////////////////////////////////////////////////////
		for(size_t j = 0; j < mXinit.size(); ++j)
		{
			for(size_t i = 0; i < mXinit[j].size(); ++i)
			{
				oaf<<mXinit[j][i]<<"\t"<<multiplicator*mYmeasu[j][i]<<"\t"<<datas[j][i]<<endl;
			}
		}
////////////////////////////////////////////////////////////////
		oaf.close();
	}
	return 0;
}
