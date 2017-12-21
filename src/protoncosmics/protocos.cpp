/**
 * \file protocos.cpp
 * \brief Implementation for the proton-cosmics impact.
 * The effects of protons and cosmics rays is computed here with
 * planetocosmics. We use the output of the code averaged here.
 * Copyright G Gronoff Dec 2009
 * Last Modification : $Id: protocos.cpp 1529 2012-07-06 14:58:19Z gronoff $
 */

#include "protocos.hpp"
using namespace std;

ProtoCosIonization::ProtoCosIonization(XmlParameters* pParam, ublas::vector<double>* vAltGridKm):mpParams(pParam),mpAltGridKm(vAltGridKm)
{
}


ProtoCosIonization::~ProtoCosIonization()
{
}

void ProtoCosIonization::ComputeCosmoionization(std::deque<Specie*>& vpSp,
		     std::deque<Specie*>& rResult,
		     EFlux& rResultFlux)
{
	mpElecC=rResultFlux.mpElecCentEeV;
	mpElecB=rResultFlux.mpElecBotEeV;
	mpElecD=rResultFlux.mpElecDdengeV;
	mpGAngle=rResultFlux.ReturnAngle();
	mNbang=mpGAngle->mNbAngles;
	//mNbang2=mpGAngle->mNbAngles/2;
	mNben=mpElecC->size();
	mNbalt=mpAltGridKm->size();


	for(unsigned sp=0;sp<vpSp.size();++sp)
	{
		(vpSp)[sp]->InitCosmoImpact(mNbalt);
	}


	std::vector<TiXmlNode*> cosnodes= mpParams->GetNodes("/aero_main/proton_cosmic");

	std::vector<TiXmlNode*>::iterator it;
	for(it=cosnodes.begin();it!=cosnodes.end();++it)
	{ // We work with several nodes!

		CosmoIonize(*it,vpSp,rResultFlux,vpSp);
	}
	
	SpecieUtils::SpeciesToResu(vpSp,rResult);
	
}

void ProtoCosIonization::CosmoIonize(TiXmlNode* node,
		     std::deque<Specie*>& rResult,
		     EFlux& rResultFlux,
		     std::deque<Specie*>& vpSp
		     )
{
	string name=mpParams->Elem(node,"//name");
	//Log::SetPriority(Log::INFO);
	Log::mI<<" Proton-Cosmic ionisation process -> "<<name<<endl;
	double ionization_energy=0;
	mpParams->GetValue(node,"//ionization_energy",ionization_energy);

	string energydepofile=mpParams->GetSubFileName(node,"//energy_deposition_file");

	XmlParameters enedepo(energydepofile);

	ublas::matrix<double> endep;

	enedepo.Get2DArray("/analyse_xml/data",endep);

	ublas::vector<double> alt, ene;

	ublas::matrix_row< ublas::matrix<double> > acol(endep,0);
	alt=acol;

	ublas::matrix_row< ublas::matrix<double> > ecol(endep,1);
	ene=ecol;

	for(ublas::vector<double>::iterator et=ene.begin();et!=ene.end();++et)
	{
		if(*et<1E-42)
		{
			*et=1E-42;
		}
	}
	if(enedepo.Exists("/analyse_xml/smooth"))
	{
		unsigned window = 42;
		enedepo.GetValueOrDefault("/analyse_xml/smooth", window);
		ublas::vector<double> ene2;
		MathFunction::MoyGliss(ene, window, ene2);
		ene = ene2;
	}





	// We get the energy in function of altitude which is not taken
	// into account in the flux: if our grid is smaller than the
	// model grid, some energy is loss: we put it inside the
	// deposition energy
	ublas::vector<double> energy_deposition=MathFunction::IntLog(alt,ene,*mpAltGridKm);
	if(mpParams->Exists(node,"//electron_flux_file"))
	{
		string efluxfile=mpParams->GetSubFileName(node,"//electron_flux_file");
		ublas::vector<double> rest_energy=CosmoElecFlux(efluxfile,rResultFlux,vpSp);
		energy_deposition+=rest_energy;// The rest energy has also this 1E9 difference!
	}
	
	ublas::vector<double> ionizationcm_3s_1=energy_deposition*1E-9/ionization_energy;// 1E9 : GeV to eV
	//	Log::mL<<"ionization:"<<ionizationcm_3s_1<<endl<<"and"<<energy_deposition*1E-9<<endl;
	//	Log::mL<<"ionization:"<<"and ENE"<<ene*1E-9<<endl;
	ublas::vector<double> totdensity(mNbalt),inverse_totdensity(mNbalt);
	totdensity.clear();

	std::deque<Specie*>::iterator it;
	for(it=rResult.begin();it!=rResult.end();++it)
	{
		totdensity+=(*it)->mTotDensitycm_3;
	}
	for(unsigned i=0;i<mNbalt;++i)
	{
		inverse_totdensity[i]=1./totdensity[i];
	}


	for(it=rResult.begin();it!=rResult.end();++it)
	{ // We fill the different ions 
		// The branching ratio is computed in another place.
		//ublas::vector<double> production=ublas::inner_prod(((*it)->mTotDensitycm_3),ublas::inner_prod(ionizationcm_3s_1,inverse_totdensity));
		//	ublas::vector<double> production2=ublas::inner_prod(ionizationcm_3s_1,inverse_totdensity);
		//ublas::vector<double> production=ionizationcm_3s_1;
		ublas::vector<double> production(mNbalt); //=ublas::inner_prod(((*it)->mTotDensitycm_3),ublas::inner_prod(ionizationcm_3s_1,inverse_totdensity));

		assert((*it)->mTotDensitycm_3.size() == inverse_totdensity.size());
		for(unsigned i=0;i<mNbalt;++i)
		{
			production[i] = ionizationcm_3s_1[i] * inverse_totdensity[i] * (*it)->mTotDensitycm_3[i];
		}
		(*it)->mSpeciesProductioncm_3s_1.at(0)+=production;

		for(unsigned i=0;i<(*it)->mCosmoNumberOfElectrons.size();++i)
		{
			(*it)->mCosmoElecProductioncm_3s_1+=production*(*it)->mCosmoNumberOfElectrons.at(i);// The fraction is integrated in the number of electrons
		}
	}

}


ublas::vector<double> ProtoCosIonization::CosmoElecFlux(std::string vEFluxFile,EFlux& rResultFlux,std::deque<Specie*>& vpSp)
{
	ublas::vector<double> resuenergy(mNbalt);
	resuenergy.clear();
	//cout<<vEFluxFile<<" "<<rResultFlux.mpElecCentEeV->size()<<endl;

	XmlParameters elecfile(vEFluxFile);
	elecfile.ExistsOrDie("/electron_xml/nbalt","Your electron flux file "+vEFluxFile+" is not valid");
	elecfile.ExistsOrDie("/electron_xml/nben","Your electron flux file "+vEFluxFile+" is not valid");
	elecfile.ExistsOrDie("/electron_xml/nbang","Your electron flux file "+vEFluxFile+" is not valid");

	elecfile.ExistsOrDie("/electron_xml/altgrid","Your electron flux file "+vEFluxFile+" is not valid");
	
	elecfile.ExistsOrDie("/electron_xml/egrid","Your electron flux file "+vEFluxFile+" is not valid");
	elecfile.ExistsOrDie("/electron_xml/flux","Your electron flux file "+vEFluxFile+" is not valid");
	unsigned nbalt=0,nben=0,nbang=0;
	ublas::vector<double> altgrid,egrid,flux;

	elecfile.GetValue("/electron_xml/nbalt",nbalt);
	elecfile.GetValue("/electron_xml/nben",nben);
	elecfile.GetValue("/electron_xml/nbang",nbang);
	
	elecfile.Get1DArray("/electron_xml/altgrid",altgrid);
	elecfile.Get1DArray("/electron_xml/egrid",egrid);
	elecfile.Get1DArray("/electron_xml/flux",flux);

	std::reverse(egrid.begin(),egrid.end());
	assert(flux.size()==nbalt*nben*nbang);

	ublas::vector< ublas::matrix<double> > totflux(nben);
	totflux.resize(nben);

	MathFunction::GaussianAngle mesangles(nbang);


	for(unsigned e=0;e<nben;++e)
	{
		ublas::matrix<double> matene(nbalt,mNbang);
		ublas::matrix<double> altene(mNbalt,mNbang);
		for(unsigned i=0;i<nbalt;++i)
		{
			ublas::vector<double> ang(nbang);
			ublas::matrix_row< ublas::matrix<double> > ro(matene,i);
			for(unsigned k=0;k<nbang;++k)
			{
				ang(k)=flux(k+i*nbang*nben+e*nbang);
			}

			if(mNbang==nbang)
			{
				ro=ang;
			}else
			{
				ro=MathFunction::IntLin(mesangles.mXmu,ang,mpGAngle->mXmu);
			}
		}
		for(unsigned k=0;k<mNbang;++k)
		{
			ublas::matrix_column< ublas::matrix<double> > oldco(matene,k);
			ublas::matrix_column< ublas::matrix<double> > newco(altene,k);
			newco=MathFunction::IntLin(altgrid,oldco,*mpAltGridKm);
		}
		totflux(nben-1-e)=altene;// we want reverted energy
		assert(altene.size1()==mNbalt);
		assert(altene.size2()==mNbang);
	}

	// Sniff easy interpolation is finished
	//
	// Now les trucs bourrins!
	//
	// compute old grid bottom and top
	ublas::vector<double> obote(nben),otope(nben);
	assert(*egrid.begin()>*(egrid.end()-1));
	
	double diffo=(egrid(0)-egrid(1))/2.;
	obote(0)=egrid(0)-diffo;
	otope(0)=egrid(0)+diffo;

	for(unsigned e=1;e<nben-1;++e)
	{
		double diff=(egrid(e)-egrid(e+1))/2.;
		obote(e)=egrid(e)-diff;
		otope(e)=obote(e-1);
	}
	otope(nben-1)=obote(nben-2);
	obote(nben-1)=0;


	// Redistribute
	ublas::vector< ublas::matrix<double> > newflux;
	ublas::vector<double> ntop=*mpElecB+*mpElecD;
	resuenergy=MathGrid::RedistributeEAAFlux(mpGAngle,obote,otope,totflux,*mpElecB,ntop,newflux);
	for(unsigned e=0;e<mNben;++e)
	{
		for(unsigned i=0;i<mNbalt;++i)
		{
			for(unsigned k=0;k<mNbang;++k)
			{
				newflux(e)(i,k)*=1E-9/(2*PI*(*mpElecD)[e]);// Correction from planeto format to our one
			}
		}
	}

	
	
	rResultFlux.AddAnisotropicFlux(newflux,vpSp);
	
	return resuenergy;
}

/*
 * For the documentation
 *
 * Inside main:
 *  proton_cosmic  : the main markup, can be used several times
 *  proton_cosmic/name : the name of the specific process (user defined name,
 *                       typically cosmic ray ionization, MeV protons,
 *                       but can be as specialized as Alpha ionization
 *                       High Z ionization...)
 *  proton_cosmic/ionization_energy : mean energy of ionization, typically
 *  				     35 eV, see ....
 *  proton_cosmic/energy_deposition_file: file containing the energy deposed
 *  					it is an xml file containing
 *				<analyse_xml>
 *					<data>
 *					altitudes    energy_depositions	
 *					</data>
 *				<analyse_xml>
 *				Where altitudes and energy_depositions are
 *				two columns of an array!
 *				The ionizations are computed from these two values. Sometimes, the electrons grid does not fit the \Prog one, and the rest energy is added to the latest value when computing ionization.
 *  proton_cosmic/electron_flux_file: file containing the electron flux
 *  	<electron_xml>
 *  		<nbalt>nbalt</nbalt>
 *  		<nben>nben</nben>
 *  		<nbang>nbang</nbang>
 *  		<altgrid><altgrid
 *  		<egrid><egrid>
 *		<flux>
 *		
 *		ang-> alt -> ener
 *
 *
 *		</flux>
 *	</electron_xml>
 * 
 *
 * Specific data in species definition:
 * the branching ratio for the ionization by cosmic rays.
 * Ionization depends on the density of species, it is assumed that
 * the ionization is equiprobable. Put the product can be specified
 * by a branching ratio: defined in the species_file.
 *
 * The branching ratio is defined as
 * <cosmic> 
 *  	<Specie name="CH3+" state="X" fraction="0.4"  electrons="1"/>
 *  	<Specie name="CH4+" state="X" fraction="0.6" electrons="1"/>
 * </cosmic>
 *
 *
 */



