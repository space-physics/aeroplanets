/**
 * \file atmo/test.cpp
 * \brief Test for the atmosphere class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: test.cpp 916 2010-03-05 00:08:07Z gronoff $
 */

#ifndef HAVE_CONFIG_H
#include "config.h"
#endif

#include "atmo.hpp"
using namespace std;



#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>


using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


void TestLoad()
{
	cout<<"initialisation"<<endl;
	Atmo machin("./test/test.xml");

	cout<<"On va imprimer l'atmosphere interpolee"<<endl;
	machin.mpModelAtmosphere->PrintNeutralAtmo("testoutput/output_neutralatmo.dat");
	machin.mpModelAtmosphere->PrintNeutralColdens("testoutput/output_colden_neutralatmo.dat");

	machin.mpPlanet->ShowData();
	machin.PrintIonosphere("testoutput/output_ionosphere");

	machin.PrintSpecies("testoutput/output_species_crossec_");
//	machin.model_atmosphere->atmo_species[0]->PhotoCrs->PrintCrs("CO2_photo_crosssection.dat");
	/*
	for(unsigned i=0;i<((machin.mpModelAtmosphere)->mAtmoSpecies).size();++i)
	{
		cout<<"specie number"<<i<<endl;
		if(((machin.mpModelAtmosphere)->mAtmoSpecies[i])->CheckElecCrs())
		{
			cout<<"Cin size"<<(((machin.mpModelAtmosphere)->mAtmoSpecies[i])->mpElecCrs)->mCrsCincm2.size()<<endl;
			cout<<"Elastic size"<<(((machin.mpModelAtmosphere)->mAtmoSpecies[i])->mpElecCrs)->mElasticCrscm2.size()<<endl;
	//		(((machin.mpModelAtmosphere)->mAtmoSpecies[i])->mpElecCrs)->PrintCrs("test/elec"+ntostr(i)+".dat","test/redist"+ntostr(i)+".dat");

		}else
		{
			cout<<" undefined"<<endl;
		}
	}
	*/
	machin.Compute();

	machin.ProceedEmissions();
	cout<<"CHECK CHEMICAL REACTION LIST"<<endl;
	BOOST_CHECK(machin.CheckChemEmit());
	cout<<"output production"<<endl;
	SpecieUtils::PrintProduction(machin.mPhotoionizationResu,machin.mAltGridKm,"testoutput/output_productions.dat",true);
	cout<<"bye"<<endl;

	cout<<"We initialize the output species"<<endl;
	deque< SpecieId > mes_id;
	SpecieId CO2p("CO2+","X");
	mes_id.push_back(CO2p);
	SpecieId COp("CO+","X");
	mes_id.push_back(COp);
	SpecieId Op("O+","X");
	mes_id.push_back(Op);
	SpecieId Cp("C+","X");
	mes_id.push_back(Cp);
	SpecieId Cpp("C++","X");
	mes_id.push_back(Cpp);
	SpecieId CO2pp("CO2++","X");
	mes_id.push_back(CO2pp);
	SpecieUtils::SelectedPrintProduction(mes_id,machin.mPhotoionizationResu,machin.mAltGridKm,"testoutput/output_selected_photoions.dat");
	SpecieUtils::SelectedPrintProduction(mes_id,machin.mElectronImpactResu,machin.mAltGridKm,"testoutput/output_selected_elecions.dat");
	SpecieUtils::SelectedPrintProduction(mes_id,machin.mTotalResu,machin.mAltGridKm,"testoutput/output_selected_ions.dat");

	
	
	BOOST_CHECK(3 == 3);
	
	cout<<"Bonne nuit"<<endl;

}



boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	Log::Init("SuperLog.txt",Log::DEBUGG,Log::DEBUGG);
	test=BOOST_TEST_SUITE("ATMO TEST");
	test->add(BOOST_TEST_CASE(&TestLoad));
	return test;
}


