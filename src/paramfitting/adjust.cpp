#include "adjust.hpp"
using namespace std;
namespace bub = boost::numeric::ublas;





void Adjustements::StartAdjustements()
{
	vector<TiXmlNode*> adjustements=mpParams->GetNodes("/aerofit/adjustment");

	for(vector<TiXmlNode*>::iterator it=adjustements.begin();it!=adjustements.end();++it)
	{
		Adjust(*it);
	}
}

void Adjustements::Adjust(TiXmlNode* vpNode)
{
	std::deque<std::string> names;
	std::deque<double> initparams;
	std::deque< std::deque<double> > nonadjustparams;
	std::deque<int> models;
	Log::mL<<"Adjustement parameters for the atmosphere"<<endl;
	if(mpParams->Exists(vpNode,"//Parameter"))
	{
		vector<TiXmlNode*> sptoadjust = mpParams->GetNodes(vpNode,"//Parameter");

		for(vector<TiXmlNode*>::iterator it=sptoadjust.begin();it!=sptoadjust.end();++it)
		{
			std::deque<double> noadjustp;
			names.push_back(trim(mpParams->Elem(*it,"//Specie")));
			int type=0;
			mpParams->GetNKey(*it,"//Model","type",type);
			models.push_back(type);
			switch(type)
			{
				case 0:
					{
						double alt0=0;
						double dens0=0;
						double texo=0;
						double tzero=0;
						double shape=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						noadjustp.push_back(alt0);
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Texo",texo);
						initparams.push_back(log(texo));
						mpParams->GetValue(*it,"//Model/To",tzero);
						initparams.push_back(log(tzero));
						mpParams->GetValue(*it,"//Model/Shape",shape);
						initparams.push_back(log(shape));
					}
					break;
				case 1:
					{
						double alt0=0;
						double dens0=0;
						double sH=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						initparams.push_back(alt0);
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/SH",sH);
						initparams.push_back(sH);
					}
					break;
				case 2:
					{

						double alt0=0;
						double dens0=0;
						double stdev=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						initparams.push_back(alt0);
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Dev",stdev);
						initparams.push_back(log(stdev));
					}
					break;

				case 3:
					{
						double alt0=0;
						double dens0=0;
						double texo=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						noadjustp.push_back(alt0);
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Texo",texo);
						initparams.push_back(log(texo));
					}
					break;
				case 4:
					{
						ublas::vector<double> alts,logval;
						mpParams->Get1DArray(*it,"//Model/altitudes",alts);
						mpParams->Get1DArray(*it,"//Model/logvalues",logval);
						if(alts.size()!=logval.size())
						{
							Log::mE<<"Size mismatch: The size of Model/altitudes and Model/logvalues are not compatible. please check"<<endl;
							Error err("Ajustement::adjust","Size mismatch","The size of Model//altitudes and Model/logvalues are not compatible. please check");
							throw err;
						}
						unsigned siz=alts.size();
						noadjustp.push_back(static_cast<double>(siz));
						for(unsigned k=0;k<siz;++k)
						{
							noadjustp.push_back(alts[k]);
							initparams.push_back(logval[k]);
						}


					}
					break;
				case 5:
					{
						double alt0=0;
						double dens0=0;
						double alt1=0;
						double dens1=0;
						double Texo=0;
						double Tmeso=0;
						double mixmassamu=0;
						mpParams->GetValue(*it,"//Model/AltT",alt0);
						mpParams->GetValue(*it,"//Model/AltM",alt1);
						mpParams->GetValue(*it,"//Model/MixMassAmu",mixmassamu);
						noadjustp.push_back(alt0);
						noadjustp.push_back(alt1);
						noadjustp.push_back(mixmassamu);

						mpParams->GetValue(*it,"//Model/DensT",dens0);
						mpParams->GetValue(*it,"//Model/DensM",dens1);
						if(dens0<0)
							dens0=1E-50;
						if(dens1<0)
							dens1=1E-50;
						initparams.push_back(log(dens0));
						initparams.push_back(log(dens1));
						mpParams->GetValue(*it,"//Model/Texo",Texo);
						mpParams->GetValue(*it,"//Model/Tmeso",Tmeso);
						initparams.push_back(log(Texo));
						initparams.push_back(log(Tmeso));
					}
					break;
				case 6:
					{
						double alt0=0;
						double dens0=0;
						double texo=0;
						double sza=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						initparams.push_back(log(alt0));
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Texo",texo);
						initparams.push_back(log(texo));
						mpParams->GetValue(*it,"//Model/SZA",sza);
						initparams.push_back(cos(sza*PI/180.));
					}
					break;
				case 7:
					{
						double alt0=0;
						double dens0=0;
						double texo=0;
						double c=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						initparams.push_back(log(alt0));
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Texo",texo);
						initparams.push_back(log(texo));
						mpParams->GetValue(*it,"//Model/C",c);
						noadjustp.push_back(c);
					}
					break;
				case 8:
					{
						double alt0=0;
						double dens0=0;
						double texo=0;
						mpParams->GetValue(*it,"//Model/alt",alt0);
						initparams.push_back(log(alt0));
						mpParams->GetValue(*it,"//Model/dens",dens0);
						initparams.push_back(log(dens0));
						mpParams->GetValue(*it,"//Model/Texo",texo);
						initparams.push_back(log(texo));
					}
					break;
				case 9:
					{
						double alt0=0;
						double dens0=0;
						double alt1=0;
						double dens1=0;
						double Texo=0;
						double Tmeso=0;
						double mixmassamu=0;
						double C=0;
						mpParams->GetValue(*it,"//Model/AltT",alt0);
						mpParams->GetValue(*it,"//Model/AltM",alt1);
						mpParams->GetValue(*it,"//Model/MixMassAmu",mixmassamu);
						noadjustp.push_back(alt0);
						noadjustp.push_back(alt1);
						noadjustp.push_back(mixmassamu);

						mpParams->GetValue(*it,"//Model/DensT",dens0);
						mpParams->GetValue(*it,"//Model/DensM",dens1);
						if(dens0<0)
							dens0=1E-50;
						if(dens1<0)
							dens1=1E-50;
						initparams.push_back(log(dens0));
						initparams.push_back(log(dens1));
						mpParams->GetValue(*it,"//Model/Texo",Texo);
						mpParams->GetValue(*it,"//Model/Tmeso",Tmeso);
						mpParams->GetValue(*it,"//Model/C",C);
						initparams.push_back(log(Texo));
						initparams.push_back(log(Tmeso));
						initparams.push_back(C);
					}
					break;

				default:
					Error err("Error in the adjustement parameters","Model does not exists","The model number "+ntostr(type)+" does not exists - ERROR");
					throw err;
			}
			nonadjustparams.push_back(noadjustp);
		}
	}

	Log::mL<<"Adjustement parameters for the electron precipitation"<<endl;
	bool isElecPrecip = false;
	deque<double> eaddparams;
	//deque<double> einitparams, eaddparams;
	unsigned eElecPrecipType = 0;
	if(mpParams->Exists(vpNode,"//ElectronPrecipParameter"))
	{
		isElecPrecip = true;

		TiXmlNode* eprecipn = mpParams->GetNode(vpNode,"//ElectronPrecipParameter");
		unsigned type=0;
		mpParams->GetNKey(eprecipn, "//Model", "type", type);
		switch(type)
		{
			case 1:
				{
					unsigned isotro=0, powlaw=0;
					double entot=0.,E0=0.;
					mpParams->GetValue(eprecipn, "//Model/entot", entot);
					mpParams->GetValue(eprecipn, "//Model/E0", E0);
					mpParams->GetValue(eprecipn, "//Model/powlaw", powlaw);
					mpParams->GetValue(eprecipn, "//Model/isotro", isotro);
					//einitparams.push_back(entot);
					//einitparams.push_back(E0);
					initparams.push_back(entot);
					initparams.push_back(E0);
					eaddparams.push_back(static_cast<double>(powlaw));
					eaddparams.push_back(static_cast<double>(isotro));
				}
				break;

			case 2:
				{
					unsigned isotro=0, powlaw=0;
					double entot=0.,E0=0.;
					mpParams->GetValue(eprecipn, "//Model/entot", entot);
					mpParams->GetValue(eprecipn, "//Model/E0", E0);
					mpParams->GetValue(eprecipn, "//Model/powlaw", powlaw);
					mpParams->GetValue(eprecipn, "//Model/isotro", isotro);
					initparams.push_back(entot);
					initparams.push_back(E0);
					eaddparams.push_back(static_cast<double>(powlaw));
					eaddparams.push_back(static_cast<double>(isotro));
				}
				break;

			case 3:
				{
					unsigned isotro=0, powlaw=0;
					double entot=0.,E0=0.;
					mpParams->GetValue(eprecipn, "//Model/entot", entot);
					mpParams->GetValue(eprecipn, "//Model/E0", E0);
					mpParams->GetValue(eprecipn, "//Model/powlaw", powlaw);
					mpParams->GetValue(eprecipn, "//Model/isotro", isotro);
					initparams.push_back(entot);
					initparams.push_back(E0);
					eaddparams.push_back(static_cast<double>(powlaw));
					eaddparams.push_back(static_cast<double>(isotro));
				}
				break;

			case 4:
				{
					ublas::vector<double> enes,logval;
					mpParams->Get1DArray(eprecipn, "//Model/energies", enes);
					mpParams->Get1DArray(eprecipn, "//Model/logflux", logval);
					if(enes.size() != logval.size())
					{
						Error err("Error in the electron precipitation parameters","Size mismatch ","The number of value in energy "+ntostr(enes.size())+" does not match the logflux" +ntostr(logval.size())+ " - ERROR");
					}
					for(size_t l = 0; l < enes.size(); ++l)
					{
						initparams.push_back(logval[l]);
						eaddparams.push_back(enes[l]);
					}
				}
				break;


			case 0:
				// passthru
			default:
				break;
		};
		eElecPrecipType = type;
	}
	Log::mL<<"End reading adjustement parameters for the electron precipitation"<<endl;

	// get measu
	//	ublas::vector<double> mYalts;
	//	ublas::vector<double> mYmeasu;

/* // Previous version, with only one possibility to get some data
	mpParams->Get1DArray(vpNode,"//data/altitudes",mYalts);
	mpParams->Get1DArray(vpNode,"//data/measurements",mYmeasu);
	if(mpParams->Exists(vpNode,"//data/errors"))
		mpParams->Get1DArray(vpNode,"//data/errors",mYerrors);
	if(mpParams->Exists(vpNode,"//data/use_mc"))
	{
		for(size_t i=0;i<mYerrors.size();++i)
		{
			mpParams->ApplyMCP(mYmeasu[i], -1, mYerrors[i]);
		}
	}

	//string adjustsp=mpParams->Elem(vpNode,"//data/Specie");

	string spname=mpParams->GetKey(vpNode, "//data/Specie", "name");
	string spstate=mpParams->GetKey(vpNode, "//data/Specie", "state");
	double tol=0;
	mpParams->GetValue(vpNode,"//data/tolerance/",tol);
	if(tol<0)
		tol=1;
	int adjustype=0;
	mpParams->GetNKey(vpNode,"//data","type",adjustype);

	std::deque<int> adjustparams;
	adjustparams.push_back(adjustype);


	if(mpParams->KeyExists(vpNode,"//data","position"))
	{
		int poslamb=0;
		mpParams->GetNKey(vpNode,"//data","position",poslamb);
		adjustparams.push_back(poslamb);
	}
*/
	std::vector<TiXmlNode*> datanodes = mpParams->GetNodes(vpNode, "//data");


	std::deque<double> toll;
	std::deque< std::deque<int> > adjustparamsl;
	std::deque< string > spnamel, spstatel;



	for(vector<TiXmlNode*>::iterator it = datanodes.begin(); it != datanodes.end(); ++it)
	{
		ublas::vector<double> tmpalt, tmpmeas, tmperr;
		mpParams->Get1DArray(*it, "//altitudes", tmpalt);
		mpParams->Get1DArray(*it, "//measurements", tmpmeas);
		if(mpParams->Exists(*it,"//errors"))
			mpParams->Get1DArray(*it,"//errors",tmperr);
		if(mpParams->Exists(*it, "//use_mc"))
		{
			for(size_t i=0 ; i < tmperr.size(); ++i)
			{
				mpParams->ApplyMCP(tmpmeas[i], -1, tmperr[i]);
			}
		}
		mMyalts.push_back(tmpalt);
		mMymeasu.push_back(tmpmeas);
		mMyerrors.push_back(tmperr);
		string spname=mpParams->GetKey(*it, "//Specie", "name");
		string spstate=mpParams->GetKey(*it, "//Specie", "state");
		spnamel.push_back(spname);
		spstatel.push_back(spstate);
		double tol=0;
		mpParams->GetValue(*it, "//tolerance/",tol);
		if(tol<0)
			tol=1;

		toll.push_back(tol);

		int adjustype=0;
		mpParams->GetNKey(*it, "/", "type", adjustype);
		std::deque<int> adjustparams;
		adjustparams.push_back(adjustype);
		if(mpParams->KeyExists(*it,"/","position"))
		{
			int poslamb=0;
			mpParams->GetNKey(*it,"/","position",poslamb);
			adjustparams.push_back(poslamb);
		}
		adjustparamsl.push_back(adjustparams);
	}

//	RunAtmoLsq* leastobj=new RunAtmoLsq(mYalts, mYmeasu, mYerrors, spname, spstate, adjustparams, mpAtmo, names, models, nonadjustparams, isElecPrecip, eElecPrecipType, eaddparams, tol);
	RunAtmoLsq* leastobj=new RunAtmoLsq(mMyalts, mMymeasu, mMyerrors, spnamel, spstatel, adjustparamsl, mpAtmo, names, models, nonadjustparams, isElecPrecip, eElecPrecipType, eaddparams, toll);
	if(mpParams->Exists(vpNode, "//convergence_option"))
	{
		unsigned cnv_opt = 0;
		mpParams->GetValue(vpNode, "//convergence_option", cnv_opt);
		leastobj->SetDiffOption(cnv_opt);
	}


	if(mpParams->Exists(vpNode,"//use_calib_adjust"))
	{
		leastobj->SetUseCalib();

		if(!mpParams->KeyExists(vpNode,"//use_calib_adjust","value"))
		{
			initparams.push_back(1.);
		}else
		{
			double keycalib=1.;
			mpParams->GetNKey(vpNode,"//use_calib_adjust","value",keycalib);
			initparams.push_back(keycalib);

		}
		mbUseCalib=true;
	}

	if(mpParams->Exists(vpNode,"//print_ionosphere_list"))
	{
		string fi=mpParams->Elem(vpNode,"//print_ionosphere_list");
		leastobj->SetIono(fi);
	}
	if(mpParams->Exists(vpNode,"//print_thermosphere_list"))
	{
		string fi=mpParams->Elem(vpNode,"//print_thermosphere_list");
		leastobj->SetThermo(fi);
	}
	if(mpParams->Exists(vpNode,"//print_compar_list"))
	{
		string fi=mpParams->Elem(vpNode,"//print_compar_list");
		leastobj->SetDiffPlot(fi);
	}

	ublas::matrix<double> comat;
	std::deque<double> origparams=initparams;
	int ret=LstSq::LMleast(leastobj,-1,initparams,comat);


	Log::mI<<"Return Value : "<<ret<<endl;
	if(ret>0)
	{
		ApplyModifs(names,initparams,models,nonadjustparams, isElecPrecip, eElecPrecipType, eaddparams);
		mFinalData = leastobj->RetrieveData(initparams.at(initparams.size()-1));
	}

	if(mpParams->Exists(vpNode,"//print_convergence_list"))
	{
		string fi=mpParams->Elem(vpNode,"//print_convergence_list");
		leastobj->PrintConvergence(fi);
	}
	FillReport(ret, names, origparams, initparams, nonadjustparams, models, isElecPrecip, eElecPrecipType, eaddparams, leastobj->GetNbCalls(),leastobj->GetChiSquarev(),comat);
	delete leastobj;
}

void Adjustements::ApplyModifs(std::deque<std::string> vNames, std::deque<double> vParams, std::deque<int> vModels,std::deque< std::deque<double> > vAdditionalParameters, bool vbIsElecPrecipAdjusted, unsigned vElecPrecipModel, std::deque<double> vElecPrecipAdditionalParameters)
{
	unsigned pos=0;// pos in the parameters
	for(size_t i=0;i<vNames.size();++i)
	{
		std::deque<double> param;
		switch(vModels[i])
		{
			case 0:
				for(size_t k=0;k<4;++k)
				{
					param.push_back(exp(vParams[pos]));
					++pos;
				}
				break;
			case 1:
				for(size_t k=0;k<3;++k)
				{
					if(k==1)
						param.push_back(exp(vParams[pos]));
					else
						param.push_back(vParams[pos]);
					++pos;
				}
				break;
			case 2:
				for(size_t k=0;k<3;++k)
				{
					if(k==1)
						param.push_back(exp(vParams[pos]));
					else
						param.push_back(vParams[pos]);
					++pos;
				}
				break;

			case 3:
				for(size_t k=0;k<2;++k)
				{
					param.push_back(exp(vParams[pos]));
					++pos;
				}
				break;
			case 4:
				for(size_t k=0;k<vAdditionalParameters[i][0];++k)
				{
					param.push_back(vParams[pos]);
					++pos;
				}
				break;
			case 5:
				{
					for(size_t k=0;k<4;++k)
					{
						param.push_back(exp(vParams[pos]));
						++pos;
					}
				}
				break;
			case 6:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(vParams[pos]));
						++pos;
					}
					param.push_back(acos(vParams[pos]));
					++pos;
				}
				break;
			case 7:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(vParams[pos]));
						++pos;
					}
				}
				break;
			case 8:
				{
					for(size_t k=0;k<3;++k)
					{
						param.push_back(exp(vParams[pos]));
						++pos;
					}
				}
				break;
			case 9:
				{
					for(size_t k=0;k<4;++k)
					{
						param.push_back(exp(vParams[pos]));
						++pos;
					}
					param.push_back((vParams[pos]));
					++pos;
				}
				break;

			default:
				Error err("Error in the adjustement model","Model does not exists","The model number "+ntostr(vModels[i])+" does not exists - ERROR");
				throw err;
		}
		//     reset atmo & elec dens
		if(vNames[i]=="e")
		{
			mpAtmo->ResetElectronDens(vModels[i],param,vAdditionalParameters[i]);
		}else
		{
			mpAtmo->ResetSpecieDens(vNames[i],vModels[i],param,vAdditionalParameters[i]);
		}
	}
	if(vbIsElecPrecipAdjusted)
	{
		std::deque<double> param;
		switch(vElecPrecipModel)
		{
			case 1:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(vParams[pos]);
						++pos;
					}
				}
				break;
			case 2:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(vParams[pos]);
						++pos;
					}
				}
				break;
			case 3:
				{
					for(size_t k=0;k<2;++k)
					{
						param.push_back(vParams[pos]);
						++pos;
					}
				}
				break;
			case 4:
				{
					for(size_t k = 0; k < vElecPrecipAdditionalParameters.size(); ++k)
					{
						param.push_back(vParams[pos]);
						++pos;
					}
				}
				break;

			case 0: // passthru
			default:
				Log::mI<<"Null precipitation"<<endl;
		}
		mpAtmo->ResetElecPrecip(vElecPrecipModel, param, vElecPrecipAdditionalParameters);
	}
	if(mbUseCalib)
	{
		++pos;
	}
	//
	assert(pos==static_cast<unsigned>(vParams.size()));

}


void Adjustements::FillReport(int vRet, std::deque<std::string> vNames, std::deque<double> vOrigParams, std::deque<double> vParams, std::deque< std::deque<double> > vAdditionalParameters, std::deque<int> vModels, bool vbIsElecPrecipAdjusted, unsigned vElecPrecipModel, std::deque<double> vElecPrecipAdditionalParameters, unsigned vNbCalls, double vChiSquarev, ublas::matrix<double> vComatrix)
{
	mXmlReport+="<xml>\n";
	assert(vOrigParams.size()==vParams.size());

	mReport+="\n\n\n================================\n New fitting report\n----------------------------\n\n";
	if(vRet<1||vRet>3)
	{
		mReport+="=========== FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		mReport+="=========== FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		mReport+="=========== FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		mXmlReport+="\t<Failed/>\n";

		switch(vRet)
		{
			case 0: mReport+=" improper input parameters.\n";break;
			case 4: mReport+="the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.\n";break;
			case 5: mReport+="number of calls to fcn has reached or exceeded maxfev.\n";break;
			case 6: mReport+="ftol is too small. No further reduction in the sum of squares is possible.\n";break;
			case 7: mReport+="xtol is too small. No further improvement in the approximate solution x is possible.\n";break;
			case 8: mReport+="gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.\n";break;
			default: mReport+=" Unknown error number "+ntostr(vRet)+"\n";
		}
	}else
	{
		mXmlReport+="\t<Success/>\n";
		mReport+="------ SUCCESS!\n\n";
	}
	mReport+="++++ After "+ntostr(vNbCalls)+" calls to your function\n";
	mXmlReport+="\t<nb_call>"+ntostr(vNbCalls)+"</nb_call>\n";
	mReport+="++++ The Chi2 reduced value is "+ntostr(vChiSquarev)+", it must be close to 1 for a good fit.\n";
	mReport+="+++++++++++++++++ For that, measurement errors should be defined correctly!\n";
	mXmlReport+="\t<chi2v>"+ntostr(vChiSquarev)+"</chi2v> <!-- should be close to 1 -->\n";

	unsigned pos=0;
	for(size_t i=0;i<vNames.size();++i)
	{
		mReport+="*************\n"+vNames[i]+"\n*************\n";

		mXmlReport+="\t<Specie name=\""+vNames[i]+"\" type=\""+ntostr(vModels[i])+"\">\n";

		switch(vModels[i])
		{	
			case 0:
				{
					double t0=0;
					double dens0=0;
					double texo=0;
					double shape=0;
					double ot0=0;
					double odens0=0;
					double otexo=0;
					double oshape=0;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					otexo=vOrigParams[pos];
					texo=vParams[pos];
					++pos;
					ot0=vOrigParams[pos];
					t0=vParams[pos];
					++pos;
					oshape=vOrigParams[pos];
					shape=vParams[pos];
					++pos;
					mReport+="Exponential: Bates Walker profile \n";
					mXmlReport+="\t\t<BatesWalker/>\n";
					mXmlReport+="\t\t<alt>"+ntostr(vAdditionalParameters.at(i).at(0))+"</alt>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(otexo))+"</Texo>\n";
					mXmlReport+="\t\t\t<To>"+ntostr(exp(ot0))+"</To>\n";
					mXmlReport+="\t\t\t<Shape>"+ntostr(exp(oshape))+"</Shape>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(texo))+"</Texo>\n";
					mXmlReport+="\t\t<To>"+ntostr(exp(t0))+"</To>\n";
					mXmlReport+="\t\t<Shape>"+ntostr(exp(shape))+"</Shape>\n";

					mReport+="\n altitude (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(0));
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Texo  : "+ntostr(exp(texo))+"\t\t ("+ntostr(exp(otexo))+")";
					mReport+="\n T0    : "+ntostr(exp(t0))+"\t\t ("+ntostr(exp(ot0))+")";
					mReport+="\n Shape : "+ntostr(exp(shape))+"\t\t ("+ntostr(exp(oshape))+")";
					mReport+="\n";
				}
				break;
			case 1:
				{
					double alt0=0;
					double dens0=0;
					double sH=0;
					double oalt0=0;
					double odens0=0;
					double osH=0;
					oalt0=vOrigParams[pos];
					alt0=vParams[pos];
					++pos;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					osH=vOrigParams[pos];
					sH=vParams[pos];
					++pos;
					mReport+="Chapman profile \n Alt 0 : "+ntostr(alt0)+"\t\t ("+ntostr(oalt0)+")";
					mXmlReport+="\t\t<Chapman/>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<alt>"+ntostr((oalt0))+"</alt>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<SH>"+ntostr((osH))+"</SH>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<alt>"+ntostr((alt0))+"</alt>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<SH>"+ntostr((sH))+"</SH>\n";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n SH    : "+ntostr(sH)+"\t\t ("+ntostr(osH)+")";
					mReport+="\n";
				}
				break;
			case 2:
				{

					double alt0=0;
					double dens0=0;
					double stdev=0;
					double oalt0=0;
					double odens0=0;
					double ostdev=0;
					oalt0=vOrigParams[pos];
					alt0=vParams[pos];
					++pos;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					ostdev=vOrigParams[pos];
					stdev=vParams[pos];
					++pos;
					mXmlReport+="\t\t<Gaussian/>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<alt>"+ntostr((oalt0))+"</alt>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Dev>"+ntostr((ostdev))+"</Dev>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<alt>"+ntostr((alt0))+"</alt>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Dev>"+ntostr((stdev))+"</Dev>\n";
					mReport+="Gaussian profile \n Alt 0 : "+ntostr(alt0)+"\t\t ("+ntostr(oalt0)+")";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Dev   : "+ntostr(stdev)+"\t\t ("+ntostr(ostdev)+")";
					mReport+="\n";
				}
				break;

			case 3:
				{
					double dens0=0;
					double texo=0;
					double odens0=0;
					double otexo=0;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					otexo=vOrigParams[pos];
					texo=vParams[pos];
					++pos;
					mXmlReport+="\t\t<ExpH/>\n";
					mReport+="Exponential: H  profile \n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(otexo))+"</Texo>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(texo))+"</Texo>\n";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Texo  : "+ntostr(exp(texo))+"\t\t ("+ntostr(exp(otexo))+")";
					mReport+="\n";
				}
				break;
			case 4:
				{
					unsigned siz=static_cast<unsigned>(vAdditionalParameters.at(i).at(0));


					mXmlReport+="\t\t<Spline/>\n";
					mXmlReport+="\t\t<altitudes>\n\t\t\t";
					std::deque<double> altlist,ologval,logval;
					for(size_t k=0;k<siz;++k)
					{
						double talt=vAdditionalParameters.at(i).at(k+1);
						altlist.push_back(talt);
						ologval.push_back(vOrigParams[pos]);
						logval.push_back(vParams[pos]);
						++pos;
						mXmlReport+=ntostr(talt)+"  ";

					}
					mXmlReport+="\t\t</altitudes>\n";
					mXmlReport+="\t\t<first_guess><logvalues>\n\t\t\t";
					for(size_t k=0;k<siz;++k)
					{
						mXmlReport+=ntostr(ologval[k])+" ";
					}
					mXmlReport+="\t\t</logvalues></first_guess>\n";

					mXmlReport+="\t\t<logvalues>\n\t\t\t";
					for(size_t k=0;k<siz;++k)
					{
						mXmlReport+=ntostr(logval[k])+" ";
					}
					mXmlReport+="\t\t</logvalues>\n";

					mReport+="Spline profile (altitudes, log(val), (first guess) \n";
					for(size_t k=0;k<siz;++k)
					{
						mReport+=ntostr(altlist[k])+" \t "+ntostr(logval[k])+" \t ( "+ntostr(ologval[k])+" )\n";
					}
				}
				break;
			case 5:
				{
					//double alt0=0;
					double dens0=0;
					double odens0=0;
					//double alt1=0;
					double dens1=0;
					double odens1=0;
					double Texo=0;
					double oTexo=0;
					double Tmeso=0;
					double oTmeso=0;
					//double mixmassamu=0;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					odens1=vOrigParams[pos];
					dens1=vParams[pos];
					++pos;
					oTexo=vOrigParams[pos];
					Texo=vParams[pos];
					++pos;
					oTmeso=vOrigParams[pos];
					Tmeso=vParams[pos];
					++pos;
					mReport+="Double Exponential: GS profile \n";
					mXmlReport+="\t\t<DoubleExp/>\n";
					mXmlReport+="\t\t<AltT>"+ntostr(vAdditionalParameters.at(i).at(0))+"</AltT>\n";
					mXmlReport+="\t\t<AltM>"+ntostr(vAdditionalParameters.at(i).at(1))+"</AltM>\n";
					mXmlReport+="\t\t<MixMassAmu>"+ntostr(vAdditionalParameters.at(i).at(2))+"</MixMassAmu>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<DensT>"+ntostr(exp(odens0))+"</DensT>\n";
					mXmlReport+="\t\t\t<DensM>"+ntostr(exp(odens1))+"</DensM>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(oTexo))+"</Texo>\n";
					mXmlReport+="\t\t\t<Tmeso>"+ntostr(exp(oTmeso))+"</Tmeso>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<DensT>"+ntostr(exp(dens0))+"</DensT>\n";
					mXmlReport+="\t\t<DensM>"+ntostr(exp(dens1))+"</DensM>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(Texo))+"</Texo>\n";
					mXmlReport+="\t\t<Tmeso>"+ntostr(exp(Tmeso))+"</Tmeso>\n";

					mReport+="\n Altitude thermosphere (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(0));
					mReport+="\n Altitude mesosphere (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(1));
					mReport+="\n Mean mass in mesosphere (amu, not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(2));
					mReport+="\n Dens T: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Dens M: "+ntostr(exp(dens1))+"\t\t ("+ntostr(exp(odens1))+")";
					mReport+="\n Texo  : "+ntostr(exp(Texo))+"\t\t ("+ntostr(exp(oTexo))+")";
					mReport+="\n Tmeso  : "+ntostr(exp(Tmeso))+"\t\t ("+ntostr(exp(oTmeso))+")";
					mReport+="\n";
				}
				break;

			case 6:
				{
					double alt0=0;
					double oalt0=0;
					double dens0=0;
					double texo=0;
					double odens0=0;
					double otexo=0;
					double sza=0;
					double osza=0;
					oalt0=vOrigParams[pos];
					alt0=vParams[pos];
					++pos;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					otexo=vOrigParams[pos];
					texo=vParams[pos];
					++pos;
					osza=vOrigParams[pos];
					sza=vParams[pos];
					++pos;

					mXmlReport+="\t\t<ChapCos/>\n";
					mReport+="Exponential: Chapman profile with SZA correction\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<alt>"+ntostr((oalt0))+"</alt>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(otexo))+"</Texo>\n";
					mXmlReport+="\t\t\t<SZA>"+ntostr(exp(osza))+"</SZA>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<alt>"+ntostr((alt0))+"</alt>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(texo))+"</Texo>\n";
					mXmlReport+="\t\t<SZA>"+ntostr(exp(sza))+"</SZA>\n";
					mReport+="\n Alt 0 : "+ntostr(alt0)+"\t\t ("+ntostr(oalt0)+")";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Texo  : "+ntostr(exp(texo))+"\t\t ("+ntostr(exp(otexo))+")";
					mReport+="\n SZA  : "+ntostr(exp(sza))+"\t\t ("+ntostr(exp(osza))+")";
					mReport+="\n";
				}
				break;

			case 7:
				{
					double alt0=0;
					double oalt0=0;
					double dens0=0;
					double texo=0;
					double odens0=0;
					double otexo=0;
					oalt0=vOrigParams[pos];
					alt0=vParams[pos];
					++pos;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					otexo=vOrigParams[pos];
					texo=vParams[pos];
					++pos;

					mXmlReport+="\t\t<ChapC/>\n";
					mReport+="Exponential: Chapman profile with 2C correction\n";
					mXmlReport+="\t\t<C>"+ntostr(vAdditionalParameters.at(i).at(0))+"</C>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<alt>"+ntostr((oalt0))+"</alt>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(otexo))+"</Texo>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<alt>"+ntostr((alt0))+"</alt>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(texo))+"</Texo>\n";
					mReport+="\n C parameter (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(0));
					mReport+="\n Alt 0 : "+ntostr(alt0)+"\t\t ("+ntostr(oalt0)+")";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Texo  : "+ntostr(exp(texo))+"\t\t ("+ntostr(exp(otexo))+")";
					mReport+="\n";
				}
				break;

			case 8:
				{
					double alt0=0;
					double oalt0=0;
					double dens0=0;
					double texo=0;
					double odens0=0;
					double otexo=0;
					oalt0=vOrigParams[pos];
					alt0=vParams[pos];
					++pos;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					otexo=vOrigParams[pos];
					texo=vParams[pos];
					++pos;

					mXmlReport+="\t\t<Epstein/>\n";
					mReport+="Exponential: Epstein profile\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<alt>"+ntostr((oalt0))+"</alt>\n";
					mXmlReport+="\t\t\t<dens>"+ntostr(exp(odens0))+"</dens>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(otexo))+"</Texo>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<alt>"+ntostr((alt0))+"</alt>\n";
					mXmlReport+="\t\t<dens>"+ntostr(exp(dens0))+"</dens>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(texo))+"</Texo>\n";
					mReport+="\n Alt 0 : "+ntostr(alt0)+"\t\t ("+ntostr(oalt0)+")";
					mReport+="\n Dens 0: "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Texo  : "+ntostr(exp(texo))+"\t\t ("+ntostr(exp(otexo))+")";
					mReport+="\n";
				}
				break;
			case 9:
				{
					//double alt0=0;
					double dens0=0;
					double odens0=0;
					//double alt1=0;
					double dens1=0;
					double odens1=0;
					double Texo=0;
					double oTexo=0;
					double Tmeso=0;
					double oTmeso=0;
					double C = 0;
					double oC = 0;
					//double mixmassamu=0;
					odens0=vOrigParams[pos];
					dens0=vParams[pos];
					++pos;
					odens1=vOrigParams[pos];
					dens1=vParams[pos];
					++pos;
					oTexo=vOrigParams[pos];
					Texo=vParams[pos];
					++pos;
					oTmeso=vOrigParams[pos];
					Tmeso=vParams[pos];
					++pos;
					oC = vOrigParams[pos];
					C = vParams[pos];
					++pos;
					mReport+="Exponential Hyperbolic: GS profile \n";
					mXmlReport+="\t\t<DoubleExp/>\n";
					mXmlReport+="\t\t<AltT>"+ntostr(vAdditionalParameters.at(i).at(0))+"</AltT>\n";
					mXmlReport+="\t\t<AltM>"+ntostr(vAdditionalParameters.at(i).at(1))+"</AltM>\n";
					mXmlReport+="\t\t<MixMassAmu>"+ntostr(vAdditionalParameters.at(i).at(2))+"</MixMassAmu>\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<DensT>"+ntostr(exp(odens0))+"</DensT>\n";
					mXmlReport+="\t\t\t<DensM>"+ntostr(exp(odens1))+"</DensM>\n";
					mXmlReport+="\t\t\t<Texo>"+ntostr(exp(oTexo))+"</Texo>\n";
					mXmlReport+="\t\t\t<Tmeso>"+ntostr(exp(oTmeso))+"</Tmeso>\n";
					mXmlReport+="\t\t\t<C>"+ntostr((oC))+"</C>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<DensT>"+ntostr(exp(dens0))+"</DensT>\n";
					mXmlReport+="\t\t<DensM>"+ntostr(exp(dens1))+"</DensM>\n";
					mXmlReport+="\t\t<Texo>"+ntostr(exp(Texo))+"</Texo>\n";
					mXmlReport+="\t\t<Tmeso>"+ntostr(exp(Tmeso))+"</Tmeso>\n";
					mXmlReport+="\t\t<C>"+ntostr(C)+"</C>\n";

					mReport+="\n Altitude thermosphere (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(0));
					mReport+="\n Altitude mesosphere (not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(1));
					mReport+="\n Mean mass in mesosphere (amu, not adjusted) :"+ntostr(vAdditionalParameters.at(i).at(2));
					mReport+="\n Dens T : "+ntostr(exp(dens0))+"\t\t ("+ntostr(exp(odens0))+")";
					mReport+="\n Dens M : "+ntostr(exp(dens1))+"\t\t ("+ntostr(exp(odens1))+")";
					mReport+="\n Texo   : "+ntostr(exp(Texo))+"\t\t ("+ntostr(exp(oTexo))+")";
					mReport+="\n Tmeso  : "+ntostr(exp(Tmeso))+"\t\t ("+ntostr(exp(oTmeso))+")";
					mReport+="\n Exc C  : "+ntostr(C)+"\t\t ("+ntostr(C)+")";
					mReport+="\n";
				}
				break;


			default:
				Error err("REPORT!!! Error in the adjustement parameters","Model does not exists","The model number "+ntostr(vModels[i])+" does not exists - ERROR");
				throw err;
		}

	}

       	//bool vbIsElecPrecipAdjusted, unsigned vElecPrecipModel, std::deque<double> vElecPrecipAdditionalParameters
	if(vbIsElecPrecipAdjusted)
	{
		mXmlReport += "\n\t<ElectronPrecipitation>\n";
		mReport += " ================ Electron precipitation adjustement ==\n";
		switch(vElecPrecipModel)
		{
			case 1:
				{
					double oE0 = vOrigParams[pos];
					double E0 = vParams[pos];
					++pos;
					double oEtot = vOrigParams[pos];
					double Etot = vParams[pos];
					++pos;
					mXmlReport+="\t\t<Maxwell/>\n";
					mReport+="Maxwell preciptations profile\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<E0>"+ntostr((oE0))+"</E0>\n";
					mXmlReport+="\t\t\t<Etot>"+ntostr(oEtot)+"</Etot>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<E0>"+ntostr((E0))+"</E0>\n";
					mXmlReport+="\t\t<Etot>"+ntostr(Etot)+"</Etot>\n";
					mReport+="\n E 0 : "+ntostr(E0)+"\t\t ("+ntostr(oE0)+")";
					mReport+="\n E tot: "+ntostr(Etot)+"\t\t ("+ntostr(oEtot)+")";
					mReport+="\n";
				}
				break;	
			
			case 2:
				{
					double oE0 = vOrigParams[pos];
					double E0 = vParams[pos];
					++pos;
					double oEtot = vOrigParams[pos];
					double Etot = vParams[pos];
					++pos;
					mXmlReport+="\t\t<Gaussian/>\n";
					mReport+="Maxwell preciptations profile\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<E0>"+ntostr((oE0))+"</E0>\n";
					mXmlReport+="\t\t\t<Etot>"+ntostr(oEtot)+"</Etot>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<E0>"+ntostr((E0))+"</E0>\n";
					mXmlReport+="\t\t<Etot>"+ntostr(Etot)+"</Etot>\n";
					mReport+="\n E 0 : "+ntostr(E0)+"\t\t ("+ntostr(oE0)+")";
					mReport+="\n E tot: "+ntostr(Etot)+"\t\t ("+ntostr(oEtot)+")";
					mReport+="\n";
				}
				break;	
			
			case 3:
				{
					double oE0 = vOrigParams[pos];
					double E0 = vParams[pos];
					++pos;
					double oEtot = vOrigParams[pos];
					double Etot = vParams[pos];
					++pos;
					mXmlReport+="\t\t<Diran/>\n";
					mReport+="Maxwell preciptations profile\n";
					mXmlReport+="\t\t<first_guess>\n";
					mXmlReport+="\t\t\t<E0>"+ntostr((oE0))+"</E0>\n";
					mXmlReport+="\t\t\t<Etot>"+ntostr(oEtot)+"</Etot>\n";
					mXmlReport+="\t\t</first_guess>\n";
					mXmlReport+="\t\t<E0>"+ntostr((E0))+"</E0>\n";
					mXmlReport+="\t\t<Etot>"+ntostr(Etot)+"</Etot>\n";
					mReport+="\n E 0 : "+ntostr(E0)+"\t\t ("+ntostr(oE0)+")";
					mReport+="\n E tot: "+ntostr(Etot)+"\t\t ("+ntostr(oEtot)+")";
					mReport+="\n";
				}
				break;	
			
			case 4:
				{
					unsigned siz= static_cast<unsigned>(vElecPrecipAdditionalParameters.size());
					mXmlReport+="\t\t<Spline/>\n";
					mXmlReport+="\t\t<Energies>\n\t\t\t";
					std::deque<double> enelist,ologval,logval;
					for(size_t k=0;k<siz;++k)
					{
						double tene=vElecPrecipAdditionalParameters.at(k);
						enelist.push_back(tene);
						ologval.push_back(vOrigParams[pos]);
						logval.push_back(vParams[pos]);
						++pos;
						mXmlReport+=ntostr(tene)+"  ";

					}
					mXmlReport+="\t\t</Energies>\n";
					mXmlReport+="\t\t<first_guess><logvalues>\n\t\t\t";
					for(size_t k=0;k<siz;++k)
					{
						mXmlReport+=ntostr(ologval[k])+" ";
					}
					mXmlReport+="\t\t</logvalues></first_guess>\n";

					mXmlReport+="\t\t<logvalues>\n\t\t\t";
					for(size_t k=0;k<siz;++k)
					{
						mXmlReport+=ntostr(logval[k])+" ";
					}
					mXmlReport+="\t\t</logvalues>\n";

					mReport+="Spline profile (Energies, log(val), (first guess) \n";
					for(size_t k=0;k<siz;++k)
					{
						mReport+=ntostr(enelist[k])+" \t "+ntostr(logval[k])+" \t ( "+ntostr(ologval[k])+" )\n";
					}

				}
				break;
			case 0: //passthru
			default:
				mXmlReport += "\t\t<Null/>\n";
				mReport += " Null precipitation \n";
		};
		mXmlReport += "\t</ElectronPrecipitation>\n";
	}


	if(mbUseCalib)
	{
		double multiplicator=vParams[pos];

		mReport+=" ------------------------\n Calibration modification : "+ntostr(multiplicator)+"\t"+ntostr(vOrigParams.at(pos));
		mReport+="\n it means that you have to multiply the data by this value to get the model\n--------------\n\n";

		mXmlReport+="\t<calibration>"+ntostr(multiplicator)+"</calibration>\n";
		++pos;
	}

	// Now, we work with the covariant matrix
	unsigned n = vComatrix.size1(); // No pb, it is a squared matrix

	mReport += " -----------------------------\n------Covariant matrix -----\n";
	mXmlReport += "<comat>\n";

	for(unsigned i = 0 ; i < n ; ++i)
	{
		mXmlReport += " <i"+ntostr(i)+">";
		for(unsigned j = 0 ; j < n ; ++j)
		{
			mReport += ntostr(vComatrix(i,j))+"\t"; 
			mXmlReport +=  " <j"+ntostr(j)+">" + ntostr(vComatrix(i,j)) + " </j"+ntostr(j)+">"; 
		}
		mXmlReport += " </i"+ntostr(i)+">\n";
		mReport += "\n";
	}
	mXmlReport += "</comat>\n";




	mXmlReport+="\t</Specie>\n";
	mXmlReport+="\n</xml>\n";
	assert(pos==static_cast<unsigned>(vParams.size()));


}



void Adjustements::PrintReport(string vSuffix)
{
	if(!mpParams->Exists("/aerofit/xmlreportfile"))
	{
		mpParams->ExistsOrDie("/aerofit/reportfile","We need reportfile for the output of the fitting");
	}
	if(mpParams->Exists("/aerofit/reportfile"))
	{
		string filename=mpParams->Elem("/aerofit/reportfile");
		if(vSuffix!="")
		{
			filename+=vSuffix;
		}
		if(FileExists(filename))
		{
			Log::mD<<"Low level warning : We overwrite"<<filename<<endl;
		}
		ofstream of(filename.c_str());
		of.precision(9);
		of.setf(ios::scientific);
		of<<mReport<<endl;
		of.close();
		Log::mI<<mReport<<endl;
	}
	if(mpParams->Exists("/aerofit/xmlreportfile"))
	{
		string fname=mpParams->Elem("/aerofit/xmlreportfile");
		if(vSuffix!="")
		{
			fname+=vSuffix;
		}
		if(FileExists(fname))
		{
			Log::mD<<"Low level warning : We overwrite"<<fname<<endl;
		}
		ofstream ofa(fname.c_str());
		ofa.precision(9);
		ofa.setf(ios::scientific);
		ofa<<mXmlReport<<endl;
		ofa.close();
	}



	if(mFinalData.size()>0 and mpParams->Exists("/aerofit/comparfile"))
	{	
		string filename2=mpParams->Elem("/aerofit/comparfile");
		if(vSuffix!="")
		{
			filename2+=vSuffix;
		}
		if(FileExists(filename2))
		{
			Log::mD<<"Low level warning : We overwrite"<<filename2<<endl;
		}
		ofstream oaf(filename2.c_str());
		oaf.precision(9);
		oaf.setf(ios::scientific);
		oaf<<"# Data Model comparison"<<endl;
		oaf<<"# Alt in Km"<<endl<<"# Prod in cm-3s-1"<<endl;
		oaf<<"# Measurements"<<endl;
		oaf<<"# Model"<<endl;
		assert(mFinalData.size()==mMyalts.size());
		for(size_t j = 0; j < mFinalData.size(); ++j)
		{
			assert(mFinalData[j].size()==mMyalts[j].size());
			for(size_t i = 0;i< mMyalts[j].size();++i)
			{
				oaf<<mMyalts[j][i]<<"\t"<<mMymeasu[j][i]<<"\t"<<mFinalData[j][i]<<endl;
			}
		}
		oaf.close();
	}
}



