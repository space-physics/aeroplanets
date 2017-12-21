/**
 * \file flux_model.cpp
 * \brief The implementation of the solar flux
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: flux_model.cpp 1550 2012-08-22 19:29:22Z gronoff $
 *
 */



#include "flux_model.hpp"
using namespace std;
using namespace boost::assign;


SolarFlux::SolarFlux(XmlParameters* pParam) : mpParameter(pParam)
{
}

SolarFlux::~SolarFlux()
{
	//Log::SetPriority(Log::DEBUGG);
	Log::mD<<"Bye bye solar flux"<<endl;
}


void SolarFlux::BoxMultiplication(ublas::vector<double>&  vFluxBoxescm_2s_1)
{
	if(!mpParameter->Exists("/aero_main/sun/model/BoxFactor"))
			return; // Ok, no problem
	Log::mI<<"You multiply your photoionization boxes by a user-defined amount"<<endl;
	ublas::matrix<double> posi;

	mpParameter->Get2DArray("/aero_main/sun/model/BoxFactor",posi);
	//Log::mL<<posi<<endl;
	if(posi.size1()!=2)
	{

			Error err("Photoionization::Solarflux::BoxMultiplication","Error in the data inside BoxFactor","You should define an array int double of position - factor (position beginning at 0). The other solution is not to use BoxFactor");
			throw err;
	}
	unsigned sizeboxes=vFluxBoxescm_2s_1.size();
	for(size_t i=0;i<posi.size2();++i)
	{
		unsigned position=static_cast<unsigned>(posi(0,i));
		double fact=posi(1,i);
		if(position<sizeboxes)
		{
			vFluxBoxescm_2s_1[position]*=fact;
		}else
		{
			Error err("Photoionization::Solarflux::BoxMultiplication","Error in the data inside BoxFactor","You want to change the flux at the  position "+ntostr(posi)+" but the maximum position is "+ntostr(sizeboxes));
			throw err;
		}

	}



}

void SolarFlux::RangeMultiplication(const ublas::vector<double>& vMainGrideV, ublas::vector<double>&  vFluxBoxescm_2s_1)
{
	vector<TiXmlNode*> ranges = mpParameter->GetNodes("/aero_main/sun/model/RangeFactor");
	for(vector<TiXmlNode*>::iterator it = ranges.begin(); it != ranges.end() ; ++it)
	{
		string unit = mpParameter->GetKey(*it,"/","unit");
		double minv = 0, maxv = 0;
		mpParameter->GetNKey(*it,"/","min",minv);
		mpParameter->GetNKey(*it,"/","max",maxv);
		if(unit!="nm" and unit!="ev" and unit!="eV")
		{
			Error err("Photoionization::Solarflux::RangeMultiplication","Error in the data inside RangeFactor","The unit should be nm, ev, or eV (without spaces), and you put --"+unit+"--");
			throw err;
		}
		if(unit == "nm")
		{
			double tmp = maxv;
			maxv = NM_TO_EV / minv;
			minv = NM_TO_EV / tmp;
		}
		double fact = 1;
		mpParameter->GetValue(*it,"/",fact); 

		for(unsigned i = 0 ; i < vMainGrideV.size(); ++i)
		{
			if(vMainGrideV[i] < maxv and vMainGrideV[i] > minv)
			{
				vFluxBoxescm_2s_1[i]*=fact;
			}
		}
	}
}


//void Solar39Boxes::f107_to_flux(double f107, double UA, double fudge,const std::vector<double>& mainGrid,std::vector<double>& resu)
	
void Solar39Boxes::F107ToFlux(double vF107, double vUA,double vFudge,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1)
{


	double f107_min=68;
	double f107_max=243;

	//	vector<double> modelGrid;
	//	modelGrid+=532.915,330.625,185.976,103.32,72.3241,55.7929,48.3746,43.6334,45.4609,40.8771,40.8138,38.3761,33.685,33.2101,29.274,26.6507,26.1744,23.6697,22.3649,21.2182,21.6033,20.3333,19.6885,19.8693,18.3933,17.6287,17.1216,16.2039,16.0933,15.7069,16.0146,15.0422,14.1812,13.4135,12.69,12.7247,12.0875,12.015,12.1032;
	vector<double> photon_Grid_eV_max;
	photon_Grid_eV_max+= 652.54842105,413.28066667,247.9684,123.9842,82.65613333,61.9921,48.38463909,43.64336266,49.59368,40.88705648,40.82381263,41.32806667,33.69495123,35.42405714,30.99605,26.6606599,27.55204444,24.79684,22.37488266,21.22818151,22.54258182,20.34327867,19.69846966,20.66403333,19.07449231,17.63867015,17.71202857,16.21390773,16.10327501,15.71692713,16.53122667,15.498025,14.58637647,13.77602222,12.70003705,13.05096842,12.09752876,12.02502069,12.39842;
	vector<double> photon_Grid_eV_min;
	photon_Grid_eV_min+= 413.28066667,247.9684,123.9842,82.65613333,61.9921,49.59368,48.37463909,43.63336266,41.32806667,40.87705648,40.81381263,35.42405714,33.68495123,30.99605,27.55204444,26.6506599,24.79684,22.54258182,22.36488266,21.21818151,20.66403333,20.33327867,19.68846966,19.07449231,17.71202857,17.62867015,16.53122667,16.20390773,16.09327501,15.70692713,15.498025,14.58637647,13.77602222,13.05096842,12.69003705,12.39842,12.08752876,12.01502069,11.80801905;

	vector<double> model_Fmax;
	model_Fmax+= 
		1.2330e+08,  2.4690e+08,
		1.1487E+09, 3.4330E+08, 4.8498E+09, 3.7013E+09, 5.9470E+08,
		3.1675E+09, 4.1358E+09, 2.4995E+09, 1.1280E+10, 5.6326E+09,
		1.3949E+09, 2.1965E+09, 9.9320E+08, 3.6210E+08, 1.6716E+09,
		1.5468E+09, 1.5904E+09, 4.8664E+09, 1.0213E+09, 1.4621E+09,
		3.0180E+09, 4.8200E+08, 4.5540E+08, 7.1650E+08, 4.2560E+08,
		4.3180E+08, 6.7090E+08, 1.5869E+09, 2.1809E+09, 5.0135E+09,
		1.3298E+10, 1.2035E+10, 1.3177E+10, 4.4204E+09,
		1.3125E+10, 9.0426E+09, 8.6669E+09	;
	vector<double> model_Fmin;
	model_Fmin+= 1.119e+07, 3.119e+07,
		3.834e+08, 1.346e+08, 18.418e+08, 9.235e+08 , 2.713e+08,
		1.000e+08, 8.405e+08, 2.350e+08 , 6.000e+09 , 8.661e+08,
		7.394e+08, 2.121e+08, 3.926e+08 , 1.800e+08 , 3.063e+08,
		5.085e+08, 7.992e+08, 1.580e+09 , 4.843e+08 ,4.500e+08,
		1.500e+09, 1.746e+08, 2.223e+08 , 3.915e+08 , 1.667e+08,
		1.997e+08, 2.425e+08, 7.931e+08 , 8.728e+08 , 19.311e+08,
		44.325e+08, 42.170e+08,59.570e+08, 17.850e+08, 
		43.750e+08, 31.840e+08,36.401e+08;

	ublas::vector<double> tmp_flux(39);
	tmp_flux.clear();
	assert(vMinGrideV.size()==vMaxGrideV.size());
	assert(photon_Grid_eV_max.size()==model_Fmax.size());
	assert(photon_Grid_eV_min.size()==model_Fmin.size());
	assert(tmp_flux.size()==model_Fmin.size());
	for(unsigned i=0;i<model_Fmin.size();++i)
	{
		tmp_flux[i]=((vF107-f107_max)*(model_Fmin[i]-model_Fmax[i])/(f107_min-f107_max)+model_Fmax[i])/(vUA*vUA);
		if(tmp_flux[i]<0)
			tmp_flux[i]=0;
	}

	BoxMultiplication(tmp_flux);
	double prev_energy=0;
	double tot_number=0;
	/*	ofstream of("init_energy.dat");
		of.precision(9);
		of.setf(ios::scientific);*/
	for(unsigned i=0;i<tmp_flux.size();++i)
	{
		double tmp_mean=(photon_Grid_eV_max[i]+photon_Grid_eV_min[i])/2.;
//		double width_mean=nabs(photon_Grid_eV_max[i]-photon_Grid_eV_min[i]);
		prev_energy+=tmp_mean*tmp_flux[i];
		tot_number+=tmp_flux[i];
		//		of<<tmp_mean<<"\t"<<tmp_flux[i]/width_mean<<endl;
	}
	//	of.close();
	//MathString::print1d(tmp_flux);
	// We compute the final result
	//rResu=MathFunction::intlin(modelGrid,tmpFlux,mainGrid);
	// To do that, an interpolation is not valid, we have to redistribute
	// unfortunately, the model is chaotic-like, because boxes have
	// overlapping! 
	double supplement=0;


	ublas::vector<double> photon_Grid_eV_min2(photon_Grid_eV_min.size());
	std::copy(photon_Grid_eV_min.begin(),photon_Grid_eV_min.end(),photon_Grid_eV_min2.begin());
	ublas::vector<double> photon_Grid_eV_max2(photon_Grid_eV_max.size());
	std::copy(photon_Grid_eV_max.begin(),photon_Grid_eV_max.end(),photon_Grid_eV_max2.begin());

	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		rResuPhcm_2s_1 = MathGrid::InterpolateNonChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
	}else{
		rResuPhcm_2s_1 = MathGrid::RedistributeChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
	}
	if(mpParameter->Exists("/aero_main/sun/model/flux_file"))
	{
		string file=mpParameter->Elem("/aero_main/sun/model/flux_file");
		ofstream off(file.c_str());
		off<<"#================= FLUX DATA =================="<<endl;
		off<<"# max grid ; min grid ; mean grid ; flux ; flux/eV ; energy ; energy/eV"<<endl;
		
		for(size_t i = 0 ; i< tmp_flux.size() ; ++ i)
		{
			double mean_energy=(photon_Grid_eV_max2[i]+photon_Grid_eV_min2[i])/2.;
			off<<photon_Grid_eV_max2[i]<<"\t"<<photon_Grid_eV_min2[i]<<"\t"<<mean_energy<<"\t";
			off<<tmp_flux[i]<<"\t"<<tmp_flux[i]/(photon_Grid_eV_max2[i]-photon_Grid_eV_min2[i])<<"\t"<<tmp_flux[i]*mean_energy<<"\t"<<tmp_flux[i]*mean_energy/(photon_Grid_eV_max2[i]-photon_Grid_eV_min2[i])<<endl;
		}
		off<<"#================= OUTPUT =================="<<endl;
		for(size_t i = 0 ; i< rResuPhcm_2s_1.size() ; ++ i)
		{
			double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
			off<<vMaxGrideV[i]<<"\t"<<vMinGrideV[i]<<"\t"<<mean_energy<<"\t";
			off<<rResuPhcm_2s_1[i]<<"\t"<<rResuPhcm_2s_1[i]/(vMaxGrideV[i]-vMinGrideV[i])<<"\t"<<rResuPhcm_2s_1[i]*mean_energy<<"\t"<<rResuPhcm_2s_1[i]*mean_energy/(vMaxGrideV[i]-vMinGrideV[i])<<endl;
		}

		off<<Log::msMessageLog<<endl;
		off.close();
	}

	double tot_energy1=0.;
	double tot_number2=0;
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
		double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
		tot_energy1+=rResuPhcm_2s_1[i]*mean_energy;
		tot_number2+=rResuPhcm_2s_1[i];
	}
	//	cout<<"Total energy through this process : "<<tot_energy1<<" eV *** to be compared with previous energy :"<<prev_energy<<endl;
	//	cout<<"Total photons through this process : "<<tot_number2<<" eV *** to be compared with previous number :"<<tot_number<<endl;
	if(supplement>10)
	{
		//Log::SetPriority(Log::WARNING);
		Log::mW<<" Warning, your grid losses solar flux! Typically, you can extend your grid border energies to avoid that problem, another solution could be to increase the number of energy bin."<<endl;
	}
	// We do not want to interpolate where the model is not valid
	//
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
	/*	if(vMinGrideV[i]>533)
		{// Theoretically, redistro avoid this
			rResuPhcm_2s_1[i]=0;
		}
		if(vMaxGrideV[i]<11.8)
		{// Theoretically, redistro avoid this
			rResuPhcm_2s_1[i]=0;
		}
		*/
		if(vMinGrideV[i]>49)
		{// < 25.7 nm: we take the Torr & Torr fudge factor into account
			rResuPhcm_2s_1[i]*=vFudge;
		}
	}

	// We perform the integration to know the total energy
	double tot_energy=0.;
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
		double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
		tot_energy+=rResuPhcm_2s_1[i]*mean_energy;
	}
	//Log::SetPriority(Log::INFO,"Solar39Boxes::F107ToFlux");
	Log::mI<<"Total energy through this process : "<<tot_energy<<" eV === to be compared with previous energy :"<<prev_energy<<endl;
	Log::mI<<"This mean a modification of the energy through grid interpolation of "<<(tot_energy-prev_energy)/prev_energy*100.<<"%"<<endl;
	Log::mI<<"Number of photons in your new grid : "<<tot_number2<<" in your old one : "<<tot_number<<" meaning a difference of "<<(tot_number2-tot_number)/tot_number*100.<<"%"<<endl;
}

Solar39Boxes::Solar39Boxes(XmlParameters* pParam) :  SolarFlux(pParam)
{
	Log::mI<<"You are using the standard 39 boxes flux"<<endl;
}

//(double UA,const std::vector<double>& mainGrid,const std::vector<double>& minGrid,const std::vector<double>& maxGrid,std::vector<double>& resuGrid)

void Solar39Boxes::RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)
{
	assert(vMaxGrideV.size()==vMinGrideV.size());
	assert(vMaxGrideV.size()==vMainGrideV.size());

	// we retrieve the f107 parameter
	
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You need to set up the f107 flux for that model");

	double f107;
	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	
	// we retrieve the fudge factor
	mpParameter->ExistsOrDie("/aero_main/sun/model/fudge","You need to set up the fudge factor for that model. Typically, this is 1. but I am lazy: I stop");

	double fudge;
	mpParameter->GetValue("/aero_main/sun/model/fudge",fudge);

// We effectively compute the result
	F107ToFlux(f107,vUA,fudge,vMinGrideV,vMaxGrideV,rResuFluxPhcm_2s_1);
	RangeMultiplication(vMainGrideV, rResuFluxPhcm_2s_1);
}




SolarUserDefined::SolarUserDefined(XmlParameters* pParam) :  SolarFlux(pParam)
{
	//Log::SetPriority(Log::INFO);
	Log::mI<<"You are using the user defined photon flux"<<endl;
}


void SolarUserDefined::RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)
{
	assert(vMaxGrideV.size()==vMinGrideV.size());
	assert(vMaxGrideV.size()==vMainGrideV.size());

	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You need to set up the f107 flux for that model");


	mpParameter->ExistsOrDie("/aero_main/sun/model/continuum/Egrid","You should define your continuum grid");
	mpParameter->ExistsOrDie("/aero_main/sun/model/continuum/PhotFlux","You should define your continuum flux");

// We effectively compute the result
	MeasurementToFlux(vUA,vMinGrideV,vMaxGrideV,rResuFluxPhcm_2s_1);
	RangeMultiplication(vMainGrideV, rResuFluxPhcm_2s_1);
}

void SolarUserDefined::MeasurementToFlux(double vUA,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1)
{
	ublas::vector<double> cont_egrid,cont_flux,resuflux, error;
	mpParameter->Get1DArray("/aero_main/sun/model/continuum/Egrid",cont_egrid);
	mpParameter->Get1DArray("/aero_main/sun/model/continuum/PhotFlux",cont_flux);

	double energie_firstgrid=0.;
	double energie_secondgrid=0.;

	if(mpParameter->Exists("/aero_main/sun/model/continuum/Error"))
	{
		mpParameter->Get1DArray("/aero_main/sun/model/continuum/Error", error);
		if(mpParameter->Exists("/aero_main/sun/model/continuum/use_mc"))
		{
			mpParameter->VectorApplyMCP(cont_flux,error,true,false);
		}else if(mpParameter->Exists("/aero_main/sun/model/continuum/use_mc_fact"))
		{
			mpParameter->VectorApplyMCP(cont_flux,error,true,true);
		}

		MinValue(cont_flux, 1E-42); // We check that the value of the flux is > 1E-42

	}



	for(unsigned i=0;i<cont_egrid.size();++i)
	{
		energie_firstgrid+=cont_flux[i]*cont_egrid[i];
	}

	if(*cont_egrid.begin()<*(cont_egrid.end()-1))
	{
		std::reverse(cont_egrid.begin(),cont_egrid.end());
		std::reverse(cont_flux.begin(),cont_flux.end());
	}


	double supplement=0;
	double supplement2=0;
	unsigned nben=cont_egrid.size();
	ublas::vector<double> obote(nben),otope(nben);
	assert(*cont_egrid.begin()>*(cont_egrid.end()-1));

	double diffo=(cont_egrid(0)-cont_egrid(1))/2.;
	obote(0)=cont_egrid(0)-diffo;
	if(obote(0) < 0 )
	{
		obote(0) = 0;
	}
	otope(0)=cont_egrid(0)+diffo;

	for(unsigned e=1;e<nben-1;++e)
	{
		double diff=(cont_egrid(e)-cont_egrid(e+1))/2.;
		obote(e)=cont_egrid(e)-diff;
		otope(e)=obote(e-1);
	}
	otope(nben-1)=obote(nben-2);

	double diffin=obote(nben-2)-cont_egrid(nben-1);
	obote(nben-1)=cont_egrid(nben-1)-diffin;
	if(obote(nben-1)<0)
		obote(nben-1)=0;

	//	cout<<obote<<endl;
	//	cout<<otope<<endl;
	//	cout<<vMinGrideV<<endl;
	//	cout<<vMaxGrideV<<endl;
	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		//rResuPhcm_2s_1 = MathGrid::InterpolateNonChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
		resuflux=MathGrid::InterpolateNonChaoticFlux(cont_flux,obote,otope,vMinGrideV,vMaxGrideV,supplement);
		Log::mI<<"Use of the interpolating system for the continuum; for the lines, it is the same"<<endl;
	}else{
		resuflux=MathGrid::RedistributeChaoticFlux(cont_flux,obote,otope,vMinGrideV,vMaxGrideV,supplement);
	}

	if(mpParameter->Exists("/aero_main/sun/model/lines/PhotFlux"))
	{

		mpParameter->ExistsOrDie("/aero_main/sun/model/lines/Egrid","You should define your lines");
		mpParameter->ExistsOrDie("/aero_main/sun/model/lines/width","You should define your lines width (one width for all the lines)");
		double linewidth=0;
		mpParameter->GetValue("/aero_main/sun/model/lines/width",linewidth);
		ublas::vector<double> meanlines,MinLineseV,MaxLineseV;
		mpParameter->Get1DArray("/aero_main/sun/model/lines/Egrid",meanlines);

		unsigned esize=meanlines.size();
		MinLineseV.resize(esize);
		MaxLineseV.resize(esize);

		ublas::vector<double> resuflux2,lineflux,lerror;
		mpParameter->Get1DArray("/aero_main/sun/model/lines/PhotFlux",lineflux);

		if(mpParameter->Exists("/aero_main/sun/model/lines/Error"))
		{
			mpParameter->Get1DArray("/aero_main/sun/model/lines/Error", lerror);
			if(mpParameter->Exists("/aero_main/sun/model/lines/use_mc"))
			{
				mpParameter->VectorApplyMCP(lineflux,lerror,true,false);
			}else if(mpParameter->Exists("/aero_main/sun/model/lines/use_mc_fact"))
			{
				mpParameter->VectorApplyMCP(lineflux,lerror,true,true);
			}

		}


		for(unsigned i=0;i<esize;++i)
		{
			MaxLineseV[i]=meanlines[i]+linewidth/2.;
			MinLineseV[i]=meanlines[i]-linewidth/2.;
			energie_firstgrid+=meanlines[i]*lineflux[i];
		}


		resuflux2=MathGrid::RedistributeChaoticFlux(lineflux,MinLineseV,MaxLineseV,vMinGrideV,vMaxGrideV,supplement2);

		resuflux+=resuflux2;
	}
	rResuPhcm_2s_1=resuflux;
	for(unsigned i=0;i<rResuPhcm_2s_1.size();++i)
	{
		energie_secondgrid+=rResuPhcm_2s_1[i]*(vMinGrideV[i]+vMaxGrideV[i])/2.;
	}
	double totaloss=supplement+supplement2;
	if(mpParameter->Exists("/aero_main/sun/model/flux_AU"))
	{
		double measAU=0;
		mpParameter->GetValue("/aero_main/sun/model/flux_AU",measAU);
		for(unsigned i=0;i<rResuPhcm_2s_1.size();++i)
		{
			rResuPhcm_2s_1[i]*=measAU*measAU/(vUA*vUA);
		}
		totaloss*=measAU*measAU/(vUA*vUA);
		energie_firstgrid*=measAU*measAU/(vUA*vUA);
		energie_secondgrid*=measAU*measAU/(vUA*vUA);
	}
	Log::mI<<"Energy outside the grid borders: "<<totaloss<<endl;
	Log::mI<<"Energy of the initial grid (AU corrected): "<<energie_firstgrid<<endl;
	Log::mI<<"Energy of the final grid (AU corrected): "<<energie_secondgrid<<endl;
	Log::mI<<"Error in the modification: "<<((totaloss+energie_secondgrid-energie_firstgrid)/(energie_firstgrid)) *100.<<"%"<<endl;
	//	cout.setf(ios::scientific,ios::floatfield);
	//	cout<<"loss (scientific format) for the photons"<<totaloss<<endl;
}










void SolarTransBoxes::F107ToFlux(double vF107, double vUA,double vFudge,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1)
{


	double f107_min=68;
	double f107_max=243;

	vector<double> photon_Grid_eV_max;

//	photon_Grid_eV_max+= 652.54842105,413.28066667,247.9684,123.9842,82.65613333,61.9921,48.37463909,43.63336266,49.59368,40.87705648,40.81381263,41.32806667,33.68495123,35.42405714,30.99605,26.6506599,27.55204444,24.79684,22.36488266,21.21818151,22.54258182,20.33327867,19.68846966,20.66403333,19.07449231,17.62867015,17.71202857,16.20390773,16.09327501,15.70692713,16.53122667,15.498025,14.58637647,13.77602222,12.69003705,13.05096842,12.08752876,12.01502069,12.39842,11.80801905,10.33201667,10.19883686,9.04994161,8.85601429,8.55063448,8.26561333,7.99898065,7.7490125,7.51419394,7.29318824,8.89566353,8.83852663,8.00833231,7.99500893,7.94261371,7.48154719;

	photon_Grid_eV_max+=652.54842105,413.28066667,247.9684,123.9842,82.65613333,61.9921,48.37563909,43.63436266,49.59368,40.87805648,40.81481263,41.32806667,33.68595123,35.42405714,30.99605,26.6516599,27.55204444,24.79684,22.36588266,21.21918151,22.54258182,20.33427867,19.68946966,20.66403333,19.07449231,17.62967015,17.71202857,16.20490773,16.09427501,15.70792713,16.53122667,15.498025,14.58637647,13.77602222,12.69103705,13.05096842,12.08852876,12.01602069,12.39842,11.80801905,10.33201667,10.19983686,9.04994161,8.85601429,8.55063448,8.26561333,7.99898065,7.7490125,7.51419394,7.29318824,8.89666353,8.83952663,8.00933231,7.99600893,7.94361371,7.48254719;

	vector<double> photon_Grid_eV_min;
//	photon_Grid_eV_min+= 413.28066667,247.9684,123.9842,82.65613333,61.9921 ,49.59368,48.37463909,43.63336266,41.32806667,40.87705648 ,40.81381263,35.42405714,33.68495123,30.99605,27.55204444 ,26.6506599,24.79684,22.54258182,22.36488266,21.21818151 ,20.66403333,20.33327867,19.68846966,19.07449231,17.71202857 ,17.62867015,16.53122667,16.20390773,16.09327501,15.70692713 ,15.498025,14.58637647,13.77602222,13.05096842,12.69003705 ,12.39842,12.08752876,12.01502069,11.80801905,10.33201667 ,9.04994161,10.19883686,8.85601429,8.55063448,8.26561333 ,7.99898065,7.7490125,7.51419394,7.29318824,7.08481143 ,8.89566353,8.83852663,8.00833231,7.99500893,7.94261371 ,7.48154719;
 	photon_Grid_eV_min+=413.28066667,247.9684,123.9842,82.65613333,61.9921,49.59368,48.37363909,43.63236266,41.32806667,40.87605648,40.81281263,35.42405714,33.68395123,30.99605,27.55204444,26.6496599,24.79684,22.54258182,22.36388266,21.21718151,20.66403333,20.33227867,19.68746966,19.07449231,17.71202857,17.62767015,16.53122667,16.20290773,16.09227501,15.70592713,15.498025,14.58637647,13.77602222,13.05096842,12.68903705,12.39842,12.08652876,12.01402069,11.80801905,10.33201667,9.04994161,10.19783686,8.85601429,8.55063448,8.26561333,7.99898065,7.7490125,7.51419394,7.29318824,7.08481143,8.89466353,8.83752663,8.00733231,7.99400893,7.94161371,7.48054719;

	vector<double> model_Fmin;
	model_Fmin+=1.11900000e+07,3.11900000e+07,3.83400000e+08,1.34600000e+08
		,1.84180000e+09,9.23500000e+08,2.71300000e+08,1.00000000e+08
		,8.40500000e+08,2.35000000e+08,6.00000000e+09,8.66100000e+08
		,7.39400000e+08,2.12100000e+08,3.92600000e+08,1.80000000e+08
		,3.06300000e+08,5.08500000e+08,7.99200000e+08,1.58000000e+09
		,4.84300000e+08,4.50000000e+08,1.50000000e+09,1.74600000e+08
		,2.22300000e+08,3.91500000e+08,1.66700000e+08,1.99700000e+08
		,2.42500000e+08,7.93100000e+08,8.72800000e+08,1.93110000e+09
		,4.43250000e+09,4.21700000e+09,5.95700000e+09,1.78500000e+09
		,4.37500000e+09,3.18400000e+09,3.64010000e+09,7.37300000e+09
		,1.60000000e+10,2.92700000e+11,1.42000000e+10,4.57000000e+10
		,-8.34000000e+09,1.06092000e+11,1.57000000e+11,2.62000000e+11
		,3.20150000e+11,4.48900000e+11,1.06000000e+09,8.26000000e+08
		,3.32000000e+09,1.66000000e+09,2.35000000e+09,7.80000000e+09;
	vector<double> model_Fmax;
	model_Fmax+=1.23300000e+08,2.46900000e+08,1.14870000e+09,3.43300000e+08
		,4.84980000e+09,3.70130000e+09,5.94700000e+08,3.16750000e+09
		,4.13580000e+09,2.49950000e+09,1.12800000e+10,5.63260000e+09
		,1.39490000e+09,2.19650000e+09,9.93200000e+08,3.62100000e+08
		,1.67160000e+09,1.54680000e+09,1.59040000e+09,4.86640000e+09
		,1.02130000e+09,1.46210000e+09,3.01800000e+09,4.82000000e+08
		,4.55400000e+08,7.16500000e+08,4.25600000e+08,4.31800000e+08
		,6.70900000e+08,1.58690000e+09,2.18090000e+09,5.01350000e+09
		,1.32980000e+10,1.20350000e+10,1.31770000e+10,4.42040000e+09
		,1.31250000e+10,9.04260000e+09,8.66690000e+09,2.39300000e+10
		,1.60000000e+10,3.94600000e+11,1.00670000e+11,2.03630000e+11
		,2.04280000e+11,5.04000000e+11,6.24200000e+11,9.48000000e+11
		,1.49000000e+12,2.41800000e+12,7.88000000e+09,3.53000000e+09
		,1.72900000e+10,7.32000000e+09,7.15000000e+09,2.86000000e+10;


	ublas::vector<double> tmp_flux(56);
	tmp_flux.clear();
	assert(vMinGrideV.size()==vMaxGrideV.size());
	assert(photon_Grid_eV_max.size()==model_Fmax.size());
	assert(photon_Grid_eV_min.size()==model_Fmin.size());
	assert(tmp_flux.size()==model_Fmin.size());
	for(unsigned i=0;i<model_Fmin.size();++i)
	{
		tmp_flux[i]=((vF107-f107_max)*(model_Fmin[i]-model_Fmax[i])/(f107_min-f107_max)+model_Fmax[i])/(vUA*vUA);
		if(tmp_flux[i]<0)
			tmp_flux[i]=0;
	}


	BoxMultiplication(tmp_flux);
	double prev_energy=0;
	double tot_number=0;
	for(unsigned i=0;i<tmp_flux.size();++i)
	{
		double tmp_mean=(photon_Grid_eV_max[i]+photon_Grid_eV_min[i])/2.;
		prev_energy+=tmp_mean*tmp_flux[i];
	//	if(photon_Grid_eV_min[i]<10.18)
			tot_number+=tmp_flux[i];
	}
	double supplement=0;


	ublas::vector<double> photon_Grid_eV_min2(photon_Grid_eV_min.size());
	std::copy(photon_Grid_eV_min.begin(),photon_Grid_eV_min.end(),photon_Grid_eV_min2.begin());
	ublas::vector<double> photon_Grid_eV_max2(photon_Grid_eV_max.size());
	std::copy(photon_Grid_eV_max.begin(),photon_Grid_eV_max.end(),photon_Grid_eV_max2.begin());

//	Log::mL<<"vMinGrideV"<<vMinGrideV<<endl;
//	Log::mL<<"vMaxGrideV"<<vMaxGrideV<<endl;

	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		rResuPhcm_2s_1 = MathGrid::InterpolateNonChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
	}else{
		rResuPhcm_2s_1 =MathGrid::RedistributeChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
	}

	if(mpParameter->Exists("/aero_main/sun/model/flux_file"))
	{
		string file=mpParameter->Elem("/aero_main/sun/model/flux_file");
		ofstream off(file.c_str());
		off<<"#================= Solar Trans-boxes model =================="<<endl;
		off<<"#================= FLUX DATA =================="<<endl;
		off<<"# max grid ; min grid ; mean grid ; flux ; flux/eV ; energy ; energy/eV"<<endl;
		for(size_t i = 0 ; i< tmp_flux.size() ; ++ i)
		{
			double mean_energy=(photon_Grid_eV_max2[i]+photon_Grid_eV_min2[i])/2.;
			off<<photon_Grid_eV_max2[i]<<"\t"<<photon_Grid_eV_min2[i]<<"\t"<<mean_energy<<"\t";
			off<<tmp_flux[i]<<"\t"<<tmp_flux[i]/(photon_Grid_eV_max2[i]-photon_Grid_eV_min2[i])<<"\t"<<tmp_flux[i]*mean_energy<<"\t"<<tmp_flux[i]*mean_energy/(photon_Grid_eV_max2[i]-photon_Grid_eV_min2[i])<<endl;
		}
		off<<"#================= OUTPUT =================="<<endl;
		for(size_t i = 0 ; i< rResuPhcm_2s_1.size() ; ++ i)
		{
			double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
			off<<vMaxGrideV[i]<<"\t"<<vMinGrideV[i]<<"\t"<<mean_energy<<"\t";
			off<<rResuPhcm_2s_1[i]<<"\t"<<rResuPhcm_2s_1[i]/(vMaxGrideV[i]-vMinGrideV[i])<<"\t"<<rResuPhcm_2s_1[i]*mean_energy<<"\t"<<rResuPhcm_2s_1[i]*mean_energy/(vMaxGrideV[i]-vMinGrideV[i])<<endl;
		}
		off<<Log::msMessageLog<<endl;
		off.close();
	}
	double tot_energy1=0.;
	double tot_number2=0;
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
		if(rResuPhcm_2s_1[i]<0)
			rResuPhcm_2s_1[i]=0;
		double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
		tot_energy1+=rResuPhcm_2s_1[i]*mean_energy;
		tot_number2+=rResuPhcm_2s_1[i];
	}
	// We do not want to interpolate where the model is not valid
	//
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
		if(vMinGrideV[i]>49)
		{// < 25.7 nm: we take the Torr & Torr fudge factor into account
			rResuPhcm_2s_1[i]*=vFudge;
		}
	}

	// We perform the integration to know the total energy
	double tot_energy=0.;
	for(unsigned i=0;i<vMinGrideV.size();++i)
	{
		double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
//		if(vMinGrideV[i]<10.18)
			tot_energy+=rResuPhcm_2s_1[i]*mean_energy;
	}
	//Log::SetPriority(Log::INFO,"SolarTransBoxes::F107ToFlux");
	//Log::mL<<rResuPhcm_2s_1<<endl;
	Log::mI<<"Total energy through this process : "<<tot_energy<<" eV === to be compared with previous energy :"<<prev_energy<<endl;
	Log::mI<<"This mean a modification of the energy through grid interpolation of "<<(tot_energy-prev_energy)/prev_energy*100.<<"%"<<endl;
	Log::mI<<"Number of photons in your new grid : "<<tot_number2<<" in your old one : "<<tot_number<<" meaning a difference of "<<(tot_number2-tot_number)/tot_number*100.<<"%"<<endl;
	if(supplement>1000 && (tot_energy-prev_energy)/prev_energy>0.01 )
	{
		//Log::SetPriority(Log::WARNING);
		Log::mL<<" Warning, your grid losses solar flux! Try to extends the extreme values of your energy grid."<<endl;
	}
}

SolarTransBoxes::SolarTransBoxes(XmlParameters* pParam) :  SolarFlux(pParam)
{
	Log::mI<<"You are using the Trans* boxes flux    650eV - 7eV ; 1.9 - 175nm"<<endl;
}

//(double UA,const std::vector<double>& mainGrid,const std::vector<double>& minGrid,const std::vector<double>& maxGrid,std::vector<double>& resuGrid)

void SolarTransBoxes::RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)
{
	assert(vMaxGrideV.size()==vMinGrideV.size());
	assert(vMaxGrideV.size()==vMainGrideV.size());
	// we retrieve the f107 parameter
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You need to set up the f107 flux for that model");
	double f107;
	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	// we retrieve the fudge factor
	mpParameter->ExistsOrDie("/aero_main/sun/model/fudge","You need to set up the fudge factor for that model. Typically, this is 1. but I am lazy: I stop");
	double fudge;
	mpParameter->GetValue("/aero_main/sun/model/fudge",fudge);
// We effectively compute the result
	F107ToFlux(f107,vUA,fudge,vMinGrideV,vMaxGrideV,rResuFluxPhcm_2s_1);
	RangeMultiplication(vMainGrideV, rResuFluxPhcm_2s_1);
}








void HEUVAC::GetF74113(ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin)
{
	vector<double> tmpw,tmpf;
	
	
	tmpw += 14.25,14.40,15.01,15.26,16.01,16.77,17.05,17.11, 
	18.62,18.97,21.60,21.80,22.10,28.47,28.79,29.52, 
	30.02,30.43,33.74,40.95,43.76,44.02,44.16,45.66, 
	46.40,46.67,47.87,49.22,50.36,50.52,50.69,52.30, 
	52.91,54.15,54.42,54.70,55.06,55.34,56.08,56.92, 
	57.36,57.56,57.88,58.96,59.62,60.30,60.85,61.07,
	61.63,61.90,62.30,62.35,62.77,62.92,63.16,63.30,
	63.65,63.72,64.11,64.60,65.21,65.71,65.85,66.30,
	67.14,67.35,68.35,69.65,70.0,70.54,70.75,71.0, 
	71.94,72.31,72.63,72.80,72.95,73.55,74.21,74.44, 
	74.83,75.03,75.29,75.46,75.73,76.01,76.48,76.83, 
	76.94,77.30,77.74,78.56,78.70,79.08,79.48,79.76, 
	80.0,80.21,80.55,80.94,81.16,81.58,81.94,82.43, 
	82.67,83.25,83.42,83.67,84.0,84.26,84.50,84.72, 
	84.86,85.16,85.50,85.69,85.87,86.23,86.40,86.77, 
	86.98,87.30,87.61,88.10,88.11,88.14,88.42,88.64, 
	88.90,89.14,89.70,90.14,90.45,90.71,91.0,91.48, 
	91.69,91.81,92.09,92.55,92.81,93.61,94.07,94.25, 
	94.39,94.81,94.90,95.37,95.51,95.81,96.05,96.49, 
	96.83,97.12,97.51,97.87,98.12,98.23,98.50,98.88, 
	99.44,99.71,99.99,100.54,100.96,101.57,102.15,103.01, 
	103.15,103.17,103.58,103.94,104.23,104.76,105.23,106.25, 
	106.57,106.93,108.05,108.46,109.50,109.98,110.56,110.62, 
	110.76,111.16,111.25,113.80,114.09,114.24,115.39,115.82,
	116.75,117.20,120.40,121.15,121.79,122.70,123.50,127.65,
	129.87,130.30,131.02,131.21,136.21,136.28,136.34,136.45,
	136.48,141.20,144.27,145.04,148.40,150.10,152.15,154.18,
	157.73,158.37,159.98,160.37,164.15,167.50,168.17,168.55,
	168.92,169.70,171.08,172.17,173.08,174.58,175.26,177.24,
	178.05,179.27,179.75,180.41,181.14,182.17,183.45,184.53,
	184.80,185.21,186.60,186.87,187.95,188.23,188.31,190.02,
	191.04,191.34,192.40,192.82,193.52,195.13,196.52,196.65, 
	197.44,198.58,200.02,201.13,202.05,202.64,203.81,204.25, 
	204.94,206.26,206.38,207.46,208.33,209.63,209.78,211.32, 
	212.14,213.78,214.75,215.16,216.88,217.0,218.19,219.13, 
	220.08,221.44,221.82,224.74,225.12,227.01,227.19,227.47, 
	228.70,230.65,231.55,232.60,233.84,234.38,237.33,239.87, 
	240.71,241.74,243.03,243.86,244.92,245.94,246.24,246.91, 
	247.18,249.18,251.10,251.95,252.19,253.78,256.32,256.38, 
	256.64,257.16,257.39,258.36,259.52,261.05,262.99,264.24, 
	264.80,270.51,271.99,272.64,274.19,275.35,275.67,276.15, 
	276.84,277.0,277.27,278.40,281.41,284.15,285.70,289.17, 
	290.69,291.70,292.78,296.19,299.50,303.31,303.78,315.02, 
	316.20,319.01,319.83,320.56,335.41,345.13,345.74,347.39, 
	349.85,356.01,360.80,364.48,368.07,399.82,401.14,401.94, 
	403.26,405.0,406.0,407.0,408.0,409.0,410.0,411.0, 
	412.0,413.0,414.0,415.0,416.0,417.0,417.24,418.0, 
	419.0,420.0,421.0,422.0,423.0,424.0,425.0,426.0, 
	427.0,428.0,429.0,430.0,430.47,431.0,432.0,433.0, 
	434.0,435.0,436.0,436.70,437.0,438.0,439.0,440.0, 
	441.0,442.0,443.0,444.0,445.0,446.0,447.0,448.0, 
	449.0,450.0,451.0,452.0,453.0,454.0,455.0,456.0, 
	457.0,458.0,459.0,460.0,461.0,462.0,463.0,464.0, 
	465.0,465.22,466.0,467.0,468.0,469.0,470.0,471.0, 
	472.0,473.0,474.0,475.0,476.0,477.0,478.0,479.0, 
	480.0,481.0,482.0,483.0,484.0,485.0,486.0,487.0, 
	488.0,489.0,489.50,490.0,491.0,492.0,493.0,494.0, 
	495.0,496.0,497.0,498.0,499.0,499.37,500.0,501.0, 
	502.0,503.0,504.0,507.93,515.60,520.66,525.80,537.02, 
	554.37,558.60,562.80,584.33,599.60,608.0,609.0,609.76, 
	610.0,611.0,612.0,613.0,614.0,615.0,616.0,616.60, 
	617.0,618.0,619.0,620.0,621.0,622.0,623.0,624.0, 
	624.93,625.0,626.0,627.0,628.0,629.0,629.73,630.0, 
	631.0,632.0,633.0,634.0,635.0,636.0,637.0,638.0, 
	639.0,640.0,640.41,640.93,641.0,641.81,642.0,643.0, 
	644.0,645.0,646.0,647.0,648.0,649.0,650.0,651.0, 
	652.0,653.0,654.0,655.0,656.0,657.0,657.30,658.0, 
	659.0,660.0,661.0,661.40,662.0,663.0,664.0,665.0, 
	666.0,667.0,668.0,669.0,670.0,671.0,672.0,673.0, 
	674.0,675.0,676.0,677.0,678.0,679.0,680.0,681.0, 
	682.0,683.0,684.0,685.0,685.71,686.0,687.0,688.0, 
	689.0,690.0,691.0,692.0,693.0,694.0,695.0,696.0, 
	697.0,698.0,699.0,700.0,701.0,702.0,703.0,703.36, 
	704.0,705.0,706.0,707.0,708.0,709.0,710.0,711.0, 
	712.0,713.0,714.0,715.0,716.0,717.0,718.0,718.50, 
	719.0,720.0,721.0,722.0,723.0,724.0,725.0,726.0, 
	727.0,728.0,729.0,730.0,731.0,732.0,733.0,734.0, 
	735.0,736.0,737.0,738.0,739.0,740.0,741.0,742.0, 
	743.0,744.0,745.0,746.0,747.0,748.0,749.0,750.0, 
	751.0,752.0,753.0,754.0,755.0,756.0,757.0,758.0, 
	758.68,759.0,759.44,760.0,760.30,761.0,761.13,762.0, 
	762.0,763.0,764.0,765.0,765.15,766.0,767.0,768.0, 
	769.0,770.0,770.41,771.0,772.0,773.0,774.0,775.0, 
	776.0,776.0,777.0,778.0,779.0,780.0,780.32,781.0, 
	782.0,783.0,784.0,785.0,786.0,786.47,787.0,787.71, 
	788.0,789.0,790.0,790.15,791.0,792.0,793.0,794.0, 
	795.0,796.0,797.0,798.0,799.0,800.0,801.0,802.0, 
	803.0,804.0,805.0,806.0,807.0,808.0,809.0,810.0, 
	811.0,812.0,813.0,814.0,815.0,816.0,817.0,818.0, 
	819.0,820.0,821.0,822.0,823.0,824.0,825.0,826.0, 
	827.0,828.0,829.0,830.0,831.0,832.0,833.0,834.0, 
	834.20,835.0,836.0,837.0,838.0,839.0,840.0,841.0, 
	842.0,843.0,844.0,845.0,846.0,847.0,848.0,849.0, 
	850.0,851.0,852.0,853.0,854.0,855.0,856.0,857.0, 
	858.0,859.0,860.0,861.0,862.0,863.0,864.0,865.0, 
	866.0,867.0,868.0,869.0,870.0,871.0,872.0,873.0, 
	874.0,875.0,876.0,877.0,878.0,879.0,880.0,881.0, 
	882.0,883.0,884.0,885.0,886.0,887.0,888.0,889.0, 
	890.0,891.0,892.0,893.0,894.0,895.0,896.0,897.0, 
	898.0,899.0,900.0,901.0,902.0,903.0,904.0,904.10, 
	905.0,906.0,907.0,908.0,909.0,910.0,911.0,912.0, 
	913.0,914.0,915.0,916.0,917.0,918.0,919.0,920.0, 
	920.96,921.0,922.0,923.0,923.15,924.0,925.0,926.0, 
	926.20,927.0,928.0,929.0,930.0,930.75,931.0,932.0, 
	933.0,933.38,934.0,935.0,936.0,937.0,937.80,938.0, 
	939.0,940.0,941.0,942.0,943.0,944.0,944.52,945.0, 
	946.0,947.0,948.0,949.0,949.74,950.0,951.0,952.0, 
	953.0,954.0,955.0,956.0,957.0,958.0,959.0,960.0, 
	961.0,962.0,963.0,964.0,965.0,966.0,967.0,968.0, 
	969.0,970.0,971.0,972.0,972.54,973.0,974.0,975.0, 
	976.0,977.0,977.02,978.0,979.0,980.0,981.0,982.0, 
	983.0,984.0,985.0,986.0,987.0,988.0,989.0,989.79, 
	990.0,991.0,991.55,992.0,993.0,994.0,995.0,996.0, 
	997.0,998.0,999.0,1000.0,1001.0,1002.0,1003.0,1004.0, 
	1005.0,1006.0,1007.0,1008.0,1009.0,1010.0,1010.20,1011.0, 
	1012.0,1013.0,1014.0,1015.0,1016.0,1017.0,1018.0,1019.0, 
	1020.0,1021.0,1022.0,1023.0,1024.0,1025.0,1025.72,1026.0, 
	1027.0,1028.0,1029.0,1030.0,1031.0,1031.91,1032.0,1033.0, 
	1034.0,1035.0,1036.0,1036.34,1037.0,1037.02,1037.61,1038.0, 
	1039.0,1040.0,1041.0,1042.0,1043.0,1044.0,1045.0,1046.0, 
	1047.0,1048.0,1049.0,1050.0,1051.0,1052.0,1053.0,1054.0,1055.0 ;

	tmpf += 0.0,0.0,0.0E+00,0.0E+00,1.0E+05,0.0E+00,0.0E+00, 
	2.0E+05,1.0E+05,6.0E+05,5.0E+05,2.0E+05,5.0E+05,8.0E+05, 
	3.20E+06,2.80E+06,7.0E+05,1.0E+06,1.80E+06,8.0E+05,2.40E+06, 
	1.0E+06,1.20E+06,6.0E+05,3.0E+06,3.20E+06,3.60E+06,4.80E+06, 
	2.0E+05,6.0E+06,6.0E+06,3.20E+06,4.0E+05,5.30E+06,2.30E+06, 
	2.0E+05,2.80E+06,8.90E+06,2.0E+06,4.80E+06,4.0E+06,3.20E+06, 
	2.90E+06,1.50E+06,1.40E+06,2.0E+06,2.90E+06,5.20E+06,2.40E+06, 
	4.0E+06,1.0E+05,1.70E+06,3.0E+06,1.0E+05,2.80E+06,4.50E+06, 
	3.30E+06,2.0E+05,2.20E+06,2.0E+06,2.40E+06,3.30E+06,2.50E+06, 
	6.0E+06,2.60E+06,3.10E+06,1.80E+06,1.10E+07,4.0E+05,2.60E+06, 
	2.40E+06,3.50E+06,3.10E+06,5.10E+06,1.90E+06,1.70E+06,2.70E+06, 
	2.0E+06,2.0E+06,1.0E+06,3.20E+06,4.10E+06,2.0E+06,3.0E+06, 
	2.0E+06,5.60E+06,2.0E+06,3.30E+06,2.60E+06,2.30E+06,3.40E+06, 
	2.40E+06,2.90E+06,1.50E+06,2.80E+06,1.80E+06,2.30E+06,2.30E+06, 
	3.60E+06,1.60E+06,2.10E+06,2.40E+06,3.40E+06,4.80E+06,6.70E+06, 
	3.50E+06,4.10E+06,3.50E+06,4.70E+06,3.0E+06,4.70E+06,3.30E+06, 
	3.10E+06,2.40E+06,4.70E+06,2.20E+06,2.30E+06,1.70E+06,1.40E+06, 
	5.90E+06,3.60E+06,1.90E+06,1.60E+06,1.0E+05,4.90E+06,1.0E+05, 
	1.50E+06,1.90E+06,4.50E+06,3.20E+06,3.50E+06,3.70E+06,2.30E+06, 
	2.40E+06,3.10E+06,1.50E+06,4.60E+06,3.80E+06,3.0E+06,2.40E+06, 
	3.0E+06,4.50E+06,6.80E+06,1.0E+05,1.40E+06,1.10E+06,1.0E+05, 
	5.0E+06,2.30E+06,2.30E+06,9.20E+06,1.50E+06,2.60E+06,5.0E+06, 
	2.40E+06,2.20E+06,4.50E+06,5.20E+06,2.10E+06,1.20E+06,1.30E+06, 
	1.60E+06,2.20E+06,6.70E+06,2.10E+06,2.70E+06,3.60E+06,1.60E+06, 
	1.0E+05,3.40E+06,6.0E+06,4.90E+06,1.30E+06,1.60E+06,5.0E+06, 
	1.80E+06,9.0E+05,1.50E+06,1.50E+06,1.60E+06,1.60E+06,1.0E+05, 
	1.0E+05,1.0E+05,1.20E+06,1.0E+05,3.80E+06,2.40E+06,2.10E+06, 
	1.0E+05,1.0E+05,1.90E+06,3.10E+06,2.10E+06,1.0E+05,1.0E+05, 
	9.0E+05,3.70E+06,2.10E+06,3.60E+06,3.70E+06,1.0E+05,4.20E+06, 
	3.80E+06,1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.0E+05,9.0E+06, 
	9.0E+05,1.19E+07,4.0E+07,1.92E+07,1.92E+07,1.28E+07,7.50E+06, 
	1.20E+07,1.10E+07,6.0E+06,3.20E+06,1.93E+07,3.56E+07,2.02E+07, 
	1.24E+07,1.87E+07,3.0E+08,9.60E+06,1.71E+07,2.56E+08,2.65E+07, 
	1.36E+08,1.82E+07,4.0E+05,1.84E+07,2.20E+08,2.50E+07,6.60E+07, 
	1.69E+07,1.04E+08,6.90E+06,4.45E+07,2.58E+07,7.0E+07,1.90E+06, 
	1.54E+07,1.91E+08,4.80E+07,1.99E+07,2.25E+07,6.40E+07,7.68E+07, 
	1.0E+08,1.74E+08,4.53E+07,9.90E+06,1.82E+07,2.33E+07,3.79E+07, 
	6.60E+07,9.60E+07,2.07E+07,3.90E+07,1.50E+07,9.90E+06,9.0E+05, 
	9.0E+05,9.0E+05,1.90E+06,1.90E+06,9.0E+05,5.40E+07,1.18E+07, 
	9.90E+06,7.60E+06,2.17E+07,3.28E+07,5.10E+07,6.50E+07,1.20E+07, 
	1.70E+07,2.01E+07,1.90E+06,3.74E+07,5.10E+07,9.03E+07,1.90E+06, 
	2.98E+07,1.84E+07,1.05E+07,1.28E+07,1.87E+07,1.0E+07,8.26E+07, 
	5.50E+07,1.69E+07,5.40E+07,1.38E+08,9.80E+07,5.0E+07,6.75E+07, 
	1.49E+07,5.0E+07,1.87E+07,1.99E+07,2.02E+07,9.60E+06,7.30E+07, 
	3.50E+07,2.90E+07,4.0E+08,6.0E+07,3.0E+07,2.0E+08,2.0E+08, 
	1.40E+08,5.60E+07,6.20E+07,1.0E+07,6.80E+07,6.40E+07,3.0E+07, 
	6.0E+07,2.50E+07,9.0E+07,3.0E+07,2.50E+07,5.60E+06,1.0E+07, 
	4.0E+07,6.0E+07,5.0E+07,1.50E+07,2.10E+08,2.40E+07,1.30E+07, 
	4.0E+07,1.80E+07,6.34E+07,9.40E+07,9.80E+06,8.0E+08,6.90E+09, 
	1.20E+08,1.0E+08,2.50E+07,1.30E+08,2.0E+07,1.40E+08,1.0E+08, 
	8.0E+07,1.50E+08,1.0E+08,1.10E+08,7.0E+07,1.20E+08,6.50E+08, 
	1.40E+07,3.10E+07,8.20E+07,4.80E+07,0.0E+00,0.0E+00,0.0E+00, 
	0.0E+00,0.0E+00,0.0E+00,0.0E+00,0.0E+00,0.0E+00,0.0E+00, 
	0.0E+00,0.0E+00,0.0E+00,2.70E+07,0.0E+00,0.0E+00,0.0E+00, 
	0.0E+00,0.0E+00,3.0E+05,3.0E+05,3.0E+05,3.0E+05,3.0E+05, 
	3.0E+05,3.0E+05,3.0E+05,7.40E+07,3.0E+05,3.0E+05,3.0E+05, 
	3.0E+05,3.0E+05,3.0E+05,1.10E+08,3.0E+05,3.0E+05,3.0E+05, 
	3.0E+05,3.0E+05,3.0E+05,3.0E+05,7.0E+05,7.0E+05,7.0E+05, 
	7.0E+05,7.0E+05,7.0E+05,7.0E+05,7.0E+05,7.0E+05,1.0E+06, 
	1.0E+06,1.0E+06,1.0E+06,1.0E+06,1.0E+06,1.0E+06,1.30E+06, 
	1.30E+06,1.30E+06,1.30E+06,1.30E+06,1.60E+06,2.90E+08,1.60E+06, 
	1.60E+06,2.0E+06,2.0E+06,2.0E+06,2.0E+06,2.30E+06,2.30E+06, 
	2.60E+06,2.60E+06,2.60E+06,3.0E+06,3.0E+06,3.30E+06,3.30E+06, 
	3.60E+06,3.60E+06,4.0E+06,4.30E+06,4.30E+06,4.60E+06,4.90E+06, 
	5.30E+06,5.30E+06,1.0E+07,5.60E+06,5.90E+06,6.30E+06,6.60E+06, 
	6.90E+06,7.60E+06,7.90E+06,8.30E+06,8.60E+06,9.20E+06,1.0E+08, 
	9.60E+06,1.02E+07,1.06E+07,1.12E+07,1.19E+07,1.10E+08,2.50E+07, 
	4.60E+07,6.50E+07,1.20E+08,7.20E+08,4.50E+07,7.0E+07,1.27E+09, 
	1.40E+08,1.0E+05,1.0E+05,5.30E+08,1.0E+05,1.0E+05,1.0E+05, 
	1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.57E+07,1.0E+05,1.0E+05, 
	1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.0E+05,2.0E+05,2.40E+08, 
	2.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05,1.59E+09,2.0E+05, 
	1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.0E+05,1.0E+05, 
	1.0E+05,1.0E+05,1.0E+05,1.09E+07,1.29E+07,1.0E+05,1.69E+07, 
	1.0E+05,1.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05, 
	2.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05,2.0E+05, 
	2.0E+05,3.0E+05,1.10E+07,3.0E+05,3.0E+05,3.0E+05,3.0E+05, 
	1.10E+07,3.0E+05,3.0E+05,3.0E+05,3.0E+05,3.0E+05,4.0E+05, 
	4.0E+05,4.0E+05,4.0E+05,4.0E+05,4.0E+05,4.0E+05,4.0E+05, 
	5.0E+05,5.0E+05,5.0E+05,5.0E+05,5.0E+05,5.0E+05,5.0E+05, 
	5.0E+05,5.0E+05,5.0E+05,5.0E+05,9.0E+07,5.0E+05,6.0E+05, 
	6.0E+05,6.0E+05,6.0E+05,6.0E+05,7.0E+05,7.0E+05,7.0E+05, 
	7.0E+05,7.0E+05,7.0E+05,7.0E+05,7.0E+05,8.0E+05,8.0E+05, 
	8.0E+05,9.0E+05,3.60E+08,9.0E+05,9.0E+05,9.0E+05,9.0E+05, 
	9.0E+05,9.0E+05,1.0E+06,1.0E+06,1.0E+06,1.10E+06,1.10E+06, 
	1.20E+06,1.20E+06,1.20E+06,1.30E+06,5.0E+07,1.30E+06,1.30E+06, 
	1.30E+06,1.40E+06,1.40E+06,1.50E+06,1.50E+06,1.50E+06,1.50E+06, 
	1.60E+06,1.60E+06,1.70E+06,1.70E+06,1.70E+06,1.90E+06,1.90E+06, 
	1.90E+06,2.0E+06,2.0E+06,2.10E+06,2.10E+06,2.10E+06,2.20E+06, 
	2.30E+06,2.30E+06,2.40E+06,2.40E+06,2.50E+06,2.60E+06,2.70E+06, 
	2.70E+06,2.80E+06,2.90E+06,2.90E+06,3.0E+06,3.10E+06,3.10E+06, 
	3.20E+06,3.30E+06,3.40E+06,3.0E+07,3.60E+06,2.30E+07,3.60E+06, 
	8.0E+07,3.70E+06,2.0E+07,3.0E+07,3.80E+06,3.90E+06,4.0E+06, 
	4.10E+06,1.70E+08,4.20E+06,4.30E+06,4.40E+06,4.50E+06,4.60E+06, 
	2.60E+08,4.80E+06,4.90E+06,5.0E+06,5.20E+06,5.20E+06,5.40E+06, 
	1.10E+07,5.60E+06,5.70E+06,5.80E+06,6.0E+06,1.40E+08,6.10E+06, 
	6.30E+06,6.40E+06,6.60E+06,6.80E+06,6.90E+06,1.30E+08,7.20E+06, 
	2.50E+08,7.30E+06,7.50E+06,7.60E+06,4.30E+08,7.90E+06,8.10E+06, 
	8.20E+06,8.50E+06,8.70E+06,8.90E+06,9.10E+06,9.40E+06,9.60E+06, 
	9.80E+06,1.02E+07,1.04E+07,1.07E+07,1.10E+07,1.12E+07,1.15E+07, 
	1.18E+07,1.21E+07,1.24E+07,1.27E+07,1.31E+07,1.34E+07,1.37E+07, 
	1.41E+07,1.45E+07,1.48E+07,1.52E+07,1.55E+07,1.60E+07,1.63E+07, 
	1.68E+07,1.72E+07,1.77E+07,1.82E+07,1.86E+07,1.91E+07,1.96E+07, 
	2.0E+07,2.06E+07,2.11E+07,2.16E+07,2.22E+07,2.28E+07,2.33E+07, 
	6.20E+08,2.40E+07,2.45E+07,2.51E+07,2.58E+07,2.65E+07,2.71E+07, 
	2.79E+07,2.85E+07,2.93E+07,3.0E+07,3.08E+07,3.16E+07,3.24E+07, 
	3.32E+07,3.40E+07,3.49E+07,3.58E+07,3.67E+07,3.77E+07,3.86E+07, 
	3.96E+07,4.06E+07,4.17E+07,4.27E+07,4.38E+07,4.49E+07,4.61E+07, 
	4.72E+07,4.84E+07,4.96E+07,5.09E+07,5.22E+07,5.35E+07,5.49E+07, 
	5.63E+07,5.77E+07,5.92E+07,6.07E+07,6.23E+07,6.39E+07,6.55E+07, 
	6.71E+07,6.88E+07,7.06E+07,7.24E+07,7.42E+07,7.61E+07,7.81E+07, 
	8.01E+07,8.21E+07,8.42E+07,8.63E+07,8.85E+07,9.08E+07,9.31E+07, 
	9.55E+07,9.79E+07,1.0E+08,1.03E+08,1.06E+08,1.08E+08,1.11E+08, 
	1.14E+08,1.17E+08,1.20E+08,1.23E+08,1.26E+08,1.29E+08,1.32E+08, 
	1.36E+08,1.10E+08,1.39E+08,1.43E+08,1.46E+08,1.50E+08,1.54E+08, 
	1.58E+08,1.62E+08,1.66E+08,2.30E+06,2.30E+06,2.30E+06,2.40E+06, 
	2.40E+06,2.50E+06,2.50E+06,2.60E+06,5.60E+07,2.60E+06,2.60E+06, 
	2.80E+06,8.0E+07,2.80E+06,2.90E+06,2.90E+06,1.10E+08,3.0E+06, 
	3.0E+06,3.10E+06,3.10E+06,1.30E+08,3.20E+06,3.20E+06,3.30E+06, 
	1.03E+08,3.30E+06,3.40E+06,3.40E+06,3.50E+06,1.80E+08,3.60E+06, 
	3.60E+06,3.70E+06,3.70E+06,3.90E+06,4.0E+06,4.0E+06,6.80E+07, 
	4.10E+06,4.20E+06,4.20E+06,4.30E+06,4.40E+06,3.0E+08,4.40E+06, 
	4.50E+06,4.60E+06,4.70E+06,4.70E+06,4.80E+06,5.0E+06,5.10E+06, 
	5.20E+06,5.30E+06,5.30E+06,5.40E+06,5.50E+06,5.60E+06,5.70E+06, 
	5.80E+06,5.90E+06,6.10E+06,6.20E+06,6.30E+06,6.40E+06,6.50E+06, 
	6.60E+06,6.0E+08,6.70E+06,6.80E+06,7.0E+06,7.20E+06,7.30E+06, 
	4.40E+09,7.40E+06,7.50E+06,7.70E+06,7.80E+06,7.90E+06,8.0E+06, 
	8.30E+06,8.40E+06,8.50E+06,8.70E+06,8.80E+06,9.0E+06,1.70E+08, 
	9.10E+06,9.40E+06,3.40E+08,9.50E+06,9.70E+06,9.80E+06,1.0E+07, 
	1.02E+07,1.03E+07,1.06E+07,1.08E+07,1.10E+07,1.11E+07,1.13E+07, 
	1.16E+07,1.18E+07,1.20E+07,1.22E+07,1.24E+07,1.27E+07,1.29E+07, 
	1.31E+07,8.0E+07,1.33E+07,1.36E+07,1.39E+07,1.41E+07,1.44E+07, 
	1.46E+07,1.49E+07,1.52E+07,1.54E+07,1.57E+07,1.60E+07,1.63E+07, 
	1.66E+07,1.68E+07,1.72E+07,3.50E+09,1.75E+07,1.78E+07,1.82E+07, 
	1.85E+07,1.88E+07,1.91E+07,2.10E+09,1.95E+07,1.98E+07,2.02E+07, 
	2.06E+07,2.09E+07,2.0E+08,2.13E+07,2.50E+08,1.04E+09,2.17E+07, 
	2.21E+07,2.26E+07,2.29E+07,2.33E+07,2.38E+07,2.42E+07,2.46E+07, 
	2.51E+07,2.55E+07,2.60E+07,2.64E+07,2.70E+07,2.74E+07,2.79E+07, 
	2.84E+07,2.89E+07,2.95E+07 ;

	vWavelengthnm.resize(tmpw.size());
	std::copy(tmpw.begin(),tmpw.end(),vWavelengthnm.begin());
	vFluxphcm_1s_1_bin.resize(tmpf.size());
	std::copy(tmpf.begin(),tmpf.end(),vFluxphcm_1s_1_bin.begin());
	ModF74113(vWavelengthnm, vFluxphcm_1s_1_bin);


}

void HEUVAC::ModF74113(const ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin)
{
	for(unsigned i = 0; i <  vWavelengthnm.size(); ++i)
	{
		if(vWavelengthnm[i] > 500 and vWavelengthnm[i] < 550 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.18;
		if(vWavelengthnm[i] > 550 and vWavelengthnm[i] < 600 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.9;
		if(vWavelengthnm[i] > 600 and vWavelengthnm[i] < 650 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.65;
		if(vWavelengthnm[i] > 650 and vWavelengthnm[i] < 700 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.7;
		if(vWavelengthnm[i] > 700 and vWavelengthnm[i] < 750 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.1;
		if(vWavelengthnm[i] > 750 and vWavelengthnm[i] < 800 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.12;
		if(vWavelengthnm[i] > 800 and vWavelengthnm[i] < 850 and vFluxphcm_1s_1_bin[i] > 1E8)
			vFluxphcm_1s_1_bin[i] *= 1.01;

		if(vWavelengthnm[i] > 150 and vWavelengthnm[i] < 250)
			vFluxphcm_1s_1_bin[i] *= 2.0;
		if(vWavelengthnm[i] > 25 and vWavelengthnm[i] < 150)
			vFluxphcm_1s_1_bin[i] *= 3.0;
		if(vWavelengthnm[i] > 0 and vWavelengthnm[i] < 25)
			vFluxphcm_1s_1_bin[i] *= 9.0;
	}

}


unsigned HEUVAC::BiSplt(const std::vector<double>& vArr, double vAl)
{
	//""" find val index in an increasing array"""
	unsigned s = 0;
	unsigned e = vArr.size() - 1;

	if(vAl < vArr[0])
		return 0;
	
	while( s < e )
	{
		unsigned mid = (s + e) / 2;

		if(vAl == vArr[mid])
			return mid;
		if(vAl < vArr[mid])
			e = mid - 1;
		if(vAl > vArr[mid])
			s = mid + 1;
	}
	if(vAl < vArr[s])
		return s-1;
	return s;
}




void HEUVAC::HeuvacFlux(double vF107, double vF107a, ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin)
{
	GetF74113(vWavelengthnm,vFluxphcm_1s_1_bin); // WARNING : the wavelength are in angstrom here!
	vector<double> zcells, acells, zlines, alines;

	// Boundary wavelengths for the 50 A bins
	zcells += 0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050;
	// A coeffs for the 50 A bins. Note that the 0-50 A scaling factor
	// is from Hinteregger's model
	acells += .05,1.0017E-02,7.1250E-03,1.3375E-02,
	   1.9450E-02,2.6467E-02,2.5000E-03,3.6542E-02,7.4083E-03,
	   2.0225E-02,8.7583E-03,3.6583E-03,1.1800E-02,4.2667E-03,
	   4.7500E-03,4.7667E-03,4.8167E-03,5.6750E-03,4.9833E-03,
	   4.4167E-03,4.3750E-03,4.3750E-03;
	// Wavelengths of the discrete lines. 335.41 A line added with
	//  same scaling factor as 284.15 A line (both Fe lines). Note that
	//  the ACELL factor for the 300-350 bin had to be reduced to 2.5000E-03
	//  to preserve the total flux in that bin.
	zlines += 256.3,284.15,303.31,303.78,335.41,368.07,465.22,
	   554.37,584.33,609.76,629.73,703.36,765.15,770.41,787.71,
	   790.15,977.02,1025.72,1031.91;
	// A coefficients for the discrete lines
	alines += 2.78E-03,1.38E-01,2.50E-02,3.33E-03,1.38E-01,6.59E-03,
	   7.49E-03,3.27E-03,5.16E-03,1.62E-02,3.33E-03,3.04E-03,
	   3.85E-03,1.28E-02,3.28E-03,3.28E-03,3.94E-03,5.18E-03,
	   5.28E-03;
	double f = 0.5 * (vF107 + vF107a) - 80.;
	assert(vWavelengthnm.size() == vFluxphcm_1s_1_bin.size());
	for(size_t i = 0; i < vWavelengthnm.size(); ++i)
	{
		double fxfact = 1;
		unsigned j = BiSplt(zlines, vWavelengthnm[i]);
		if(nabs(vWavelengthnm[i] - zlines[j]) < 0.05)
		{
			fxfact = 1. + alines[j] * f;
		}
		else
		{
			j = BiSplt(zcells, vWavelengthnm[i]);
			fxfact = 1 + acells[j] * f;
		}
		if(fxfact < 0.8)
			fxfact = 0.8;
		vFluxphcm_1s_1_bin[i] *= fxfact;
	}
	vWavelengthnm /=10.; // the wavelength are actually in nm

}




HeuvacFlux::HeuvacFlux(XmlParameters* pParam) :  SolarFlux(pParam)
{
	Log::mI<<"You are using the HEUVAC flux; 1.4 - 105.5nm"<<endl;
}

//(double UA,const std::vector<double>& mainGrid,const std::vector<double>& minGrid,const std::vector<double>& maxGrid,std::vector<double>& resuGrid)

void HeuvacFlux::RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)
{
	assert(vMaxGrideV.size()==vMinGrideV.size());
	assert(vMaxGrideV.size()==vMainGrideV.size());
	// we retrieve the f107 parameter
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You need to set up the f107 flux for that model");
	double f107;
	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107av","You need to set up the f107av flux for that model");
	double f107av;
	mpParameter->GetValue("/aero_main/sun/model/f107av",f107av);
// We effectively compute the result
	HeuvacF107ToFlux(f107,f107av,vUA,vMinGrideV,vMaxGrideV,rResuFluxPhcm_2s_1);
	RangeMultiplication(vMainGrideV, rResuFluxPhcm_2s_1);
}


void HeuvacFlux::HeuvacF107ToFlux(double vF107, double vF107av, double vUA, const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1)
{
	ublas::vector<double> wavelengthnm, fluxphcm_1s_1_bin;
	HEUVAC::HeuvacFlux(vF107, vF107av, wavelengthnm, fluxphcm_1s_1_bin);

	ublas::vector<double> wmax(wavelengthnm.size()),wmin(wavelengthnm.size());

	ublas::vector<double> mine(wmin.size()), maxe(wmin.size());
	for(size_t i=0; i<(wavelengthnm.size()); ++i)
	{
		mine[i] = NM_TO_EV;
		maxe[i] = NM_TO_EV;
		if(wavelengthnm.size()-1 == i)
		{

			wmax[i] = 1.5 * wavelengthnm[i] - 0.5 * wavelengthnm[i-1] ;
		}else{
			wmax[i] =  0.5 * (wavelengthnm[i] +  wavelengthnm[i+1]) ;
		}
		if(0 == i)
		{
			wmin[i] = 1.5 * wavelengthnm[0] - 0.5 * wavelengthnm[1] ;
		}else{
			wmin[i] = wmax[i-1];
		}
	}

	BoxMultiplication(fluxphcm_1s_1_bin);
	mine = ublas::element_div(mine, wmax);
	maxe = ublas::element_div(maxe, wmin);

	double supplement = 0;
	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		//rResuPhcm_2s_1 = MathGrid::InterpolateNonChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
		rResuPhcm_2s_1 =MathGrid::InterpolateNonChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
	}else{
		//rResuPhcm_2s_1 =MathGrid::RedistributeChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement);
		rResuPhcm_2s_1 =MathGrid::RedistributeChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
	}
	rResuPhcm_2s_1 /= (vUA*vUA);
	

	if(mpParameter->Exists("/aero_main/sun/model/flux_file"))
	{
		string file=mpParameter->Elem("/aero_main/sun/model/flux_file");
		ofstream off(file.c_str());
		off<<"#================= HEUVAC model =================="<<endl;
		off<<"#================= FLUX DATA =================="<<endl;
		off<<"# max grid ; min grid ; mean grid ; flux ; flux/eV ; energy ; energy/eV"<<endl;
		for(size_t i = 0 ; i< fluxphcm_1s_1_bin.size() ; ++ i)
		{
			double mean_energy=(maxe[i]+mine[i])/2.;
			off<<maxe[i]<<"\t"<<mine[i]<<"\t"<<mean_energy<<"\t";
			off<<fluxphcm_1s_1_bin[i]<<"\t"<<fluxphcm_1s_1_bin[i]/(maxe[i]-mine[i])<<"\t"<<fluxphcm_1s_1_bin[i]*mean_energy<<"\t"<<fluxphcm_1s_1_bin[i]*mean_energy/(maxe[i]-mine[i])<<endl;
		}
		off<<"#================= OUTPUT =================="<<endl;
		for(size_t i = 0 ; i< rResuPhcm_2s_1.size() ; ++ i)
		{
			double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
			off<<vMaxGrideV[i]<<"\t"<<vMinGrideV[i]<<"\t"<<mean_energy<<"\t";
			off<<rResuPhcm_2s_1[i]<<"\t"<<rResuPhcm_2s_1[i]/(vMaxGrideV[i]-vMinGrideV[i])<<"\t"<<rResuPhcm_2s_1[i]*mean_energy<<"\t"<<rResuPhcm_2s_1[i]*mean_energy/(vMaxGrideV[i]-vMinGrideV[i])<<endl;
		}
		off<<Log::msMessageLog<<endl;
		off.close();
	}
}



HeuvacFluxLow::HeuvacFluxLow(XmlParameters* pParam) :  SolarFlux(pParam)
{
	Log::mI<<"You are using the HEUVAC flux Low ; 1.4 - 165.7nm"<<endl;
}

//(double UA,const std::vector<double>& mainGrid,const std::vector<double>& minGrid,const std::vector<double>& maxGrid,std::vector<double>& resuGrid)

void HeuvacFluxLow::RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)
{
	assert(vMaxGrideV.size()==vMinGrideV.size());
	assert(vMaxGrideV.size()==vMainGrideV.size());
	// we retrieve the f107 parameter
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107","You need to set up the f107 flux for that model");
	double f107;
	mpParameter->GetValue("/aero_main/sun/model/f107",f107);
	mpParameter->ExistsOrDie("/aero_main/sun/model/f107av","You need to set up the f107av flux for that model");
	double f107av;
	mpParameter->GetValue("/aero_main/sun/model/f107av",f107av);
// We effectively compute the result
	HeuvacLowF107ToFlux(f107,f107av,vUA,vMinGrideV,vMaxGrideV,rResuFluxPhcm_2s_1);
	RangeMultiplication(vMainGrideV, rResuFluxPhcm_2s_1);
}


void HeuvacFluxLow::HeuvacLowF107ToFlux(double vF107, double vF107av, double vUA, const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1)
{

	// We work on the standard HEUVAC model
	ublas::vector<double> wavelengthnm, fluxphcm_1s_1_bin;
	HEUVAC::HeuvacFlux(vF107, vF107av, wavelengthnm, fluxphcm_1s_1_bin);

	ublas::vector<double> wmax(wavelengthnm.size()),wmin(wavelengthnm.size());

	ublas::vector<double> mine(wmin.size()), maxe(wmin.size());
	for(size_t i=0; i<(wavelengthnm.size()); ++i)
	{
		mine[i] = NM_TO_EV;
		maxe[i] = NM_TO_EV;
		if(wavelengthnm.size()-1 == i)
		{

			wmax[i] = 1.5 * wavelengthnm[i] - 0.5 * wavelengthnm[i-1] ;
		}else{
			wmax[i] =  0.5 * (wavelengthnm[i] +  wavelengthnm[i+1]) ;
		}
		if(0 == i)
		{
			wmin[i] = 1.5 * wavelengthnm[0] - 0.5 * wavelengthnm[1] ;
		}else{
			wmin[i] = wmax[i-1];
		}
	}
	mine = ublas::element_div(mine, wmax);
	maxe = ublas::element_div(maxe, wmin);
	double supplement = 0;
	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		rResuPhcm_2s_1 =MathGrid::InterpolateNonChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
	}else{
		rResuPhcm_2s_1 =MathGrid::RedistributeChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
	}
//////////////////////////////
	double f107_min=68;
	double f107_max=243;
	vector<double> photon_Grid_eV_max;
	photon_Grid_eV_max+= 12.39842, 11.80801905,10.33201667,
		10.19983686,9.04994161,8.85601429,8.55063448,8.26561333,
		7.99898065,7.7490125,7.51419394,7.29318824,8.89666353,
		8.83952663,8.00933231,7.99600893,7.94361371,7.48254719;
	vector<double> photon_Grid_eV_min;
 	photon_Grid_eV_min+= 11.80801905, 10.33201667,9.04994161,
		10.19783686,8.85601429,8.55063448,8.26561333,7.99898065,
		7.7490125,7.51419394,7.29318824,7.08481143,8.89466353,
		8.83752663,8.00733231,7.99400893,7.94161371,7.48054719;
	vector<double> model_Fmin;
	model_Fmin+= 3.64010000e+09,7.37300000e+09 ,1.60000000e+10,
		2.92700000e+11,1.42000000e+10,4.57000000e+10,-8.34000000e+09,1.06092000e+11,
		1.57000000e+11,2.62000000e+11,3.20150000e+11,4.48900000e+11,1.06000000e+09,
		8.26000000e+08,3.32000000e+09,1.66000000e+09,2.35000000e+09,7.80000000e+09;
	vector<double> model_Fmax;
	model_Fmax+= 8.66690000e+09,2.39300000e+10 ,1.60000000e+10,
		3.94600000e+11,1.00670000e+11,2.03630000e+11,2.04280000e+11,5.04000000e+11,
		6.24200000e+11,9.48000000e+11,1.49000000e+12,2.41800000e+12,7.88000000e+09,
		3.53000000e+09,1.72900000e+10,7.32000000e+09,7.15000000e+09,2.86000000e+10;
	ublas::vector<double> tmp_flux(model_Fmax.size());
	tmp_flux.clear();
	assert(photon_Grid_eV_max.size()==model_Fmax.size());
	assert(photon_Grid_eV_min.size()==model_Fmin.size());
	assert(tmp_flux.size()==model_Fmin.size());
	for(unsigned i=0;i<model_Fmin.size();++i)
	{
		tmp_flux[i]=((vF107-f107_max)*(model_Fmin[i]-model_Fmax[i])/(f107_min-f107_max)+model_Fmax[i]);
		if(tmp_flux[i]<0)
			tmp_flux[i]=0;
	}
	/*double prev_energy=0;
	double tot_number=0;
	for(unsigned i=0;i<tmp_flux.size();++i)
	{
		double tmp_mean=(photon_Grid_eV_max[i]+photon_Grid_eV_min[i])/2.;
		prev_energy+=tmp_mean*tmp_flux[i];
			tot_number+=tmp_flux[i];
	}*/
	double supplement2=0;
	ublas::vector<double> photon_Grid_eV_min2(photon_Grid_eV_min.size());
	std::copy(photon_Grid_eV_min.begin(),photon_Grid_eV_min.end(),photon_Grid_eV_min2.begin());
	ublas::vector<double> photon_Grid_eV_max2(photon_Grid_eV_max.size());
	std::copy(photon_Grid_eV_max.begin(),photon_Grid_eV_max.end(),photon_Grid_eV_max2.begin());
	ublas::vector<double> newresu ;
	if(mpParameter->Exists("/aero_main/sun/model/InterpolateFlux"))
	{
		//	rResuPhcm_2s_1 =MathGrid::InterpolateNonChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
		newresu = MathGrid::InterpolateNonChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement2);
	}else{
		//	rResuPhcm_2s_1 =MathGrid::RedistributeChaoticFlux(fluxphcm_1s_1_bin,mine, maxe,vMinGrideV,vMaxGrideV,supplement);
		newresu = MathGrid::RedistributeChaoticFlux(tmp_flux,photon_Grid_eV_min2,photon_Grid_eV_max2,vMinGrideV,vMaxGrideV,supplement2);
	}


	assert(newresu.size() == rResuPhcm_2s_1.size());
	for(size_t i = 0; i< rResuPhcm_2s_1.size(); ++i)
	{
		if(rResuPhcm_2s_1[i] <=1E-42 && newresu[i] > 0)
		{
			rResuPhcm_2s_1[i] = newresu[i];
		}
	}

	
	
	rResuPhcm_2s_1 /= (vUA*vUA);
	

	if(mpParameter->Exists("/aero_main/sun/model/flux_file"))
	{
		string file=mpParameter->Elem("/aero_main/sun/model/flux_file");
		ofstream off(file.c_str());
		off<<"#================= HEUVAC Low model =================="<<endl;
		off<<"#================= FLUX DATA =================="<<endl;
		off<<"# max grid ; min grid ; mean grid ; flux ; flux/eV ; energy ; energy/eV"<<endl;
		for(size_t i = 0 ; i< fluxphcm_1s_1_bin.size() ; ++ i)
		{
			double mean_energy=(maxe[i]+mine[i])/2.;
			off<<maxe[i]<<"\t"<<mine[i]<<"\t"<<mean_energy<<"\t";
			off<<fluxphcm_1s_1_bin[i]<<"\t"<<fluxphcm_1s_1_bin[i]/(maxe[i]-mine[i])<<"\t"<<fluxphcm_1s_1_bin[i]*mean_energy<<"\t"<<fluxphcm_1s_1_bin[i]*mean_energy/(maxe[i]-mine[i])<<endl;
		}
		off<<"#================= OUTPUT =================="<<endl;
		for(size_t i = 0 ; i< rResuPhcm_2s_1.size() ; ++ i)
		{
			double mean_energy=(vMinGrideV[i]+vMaxGrideV[i])/2.;
			off<<vMaxGrideV[i]<<"\t"<<vMinGrideV[i]<<"\t"<<mean_energy<<"\t";
			off<<rResuPhcm_2s_1[i]<<"\t"<<rResuPhcm_2s_1[i]/(vMaxGrideV[i]-vMinGrideV[i])<<"\t"<<rResuPhcm_2s_1[i]*mean_energy<<"\t"<<rResuPhcm_2s_1[i]*mean_energy/(vMaxGrideV[i]-vMinGrideV[i])<<endl;
		}
		off<<Log::msMessageLog<<endl;
		off.close();
	}
}
