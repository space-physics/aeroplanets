/**
 * \file  elecimpact/test/test.cpp
 * \brief  test the electroionization
 *	  Copyright G Gronoff Sept 2009
 *	  Last Modification : $Id: test.cpp 745 2009-11-23 14:41:54Z gronoff $
 */

#ifndef HAVE_CONFIG_H
#include "config.h"
#endif

#include "electronionization.hpp"
using namespace std;
using namespace boost::assign;
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


void TestElectronImpact()
{

	/*
		vector<double> eCentEeV;//=MathGrid::GridExp(0.1,1600,30);
		vector<double> eDdengeV;//=MathGrid::WidthGrid(eCentEeV,1);
		double mSpfactor=0;
		if(!MathGrid::GridPolo(30,0.1,1600,eCentEeV,eDdengeV,mSpfactor))
		{
			BOOST_CHECK(false);
		}
		std::reverse(eCentEeV.begin(),eCentEeV.end());
		std::reverse(eDdengeV.begin(),eDdengeV.end());
		ElecCrossSection cross(eCentEeV,eDdengeV);
		cout<<"Print cente :"<<endl;
		MathString::print1d(eCentEeV);
		cout<<"print engdd"<<endl;
		MathString::print1d(eDdengeV);
		try
		{
			cout<<"===================="<<endl;
			cout<<"======Elec Crs Loading======="<<endl;
			cout<<"===================="<<endl;
		cross.LoadCrs("./eleconespectest.xml","CO2",false);
			cout<<"===================="<<endl;
			cout<<"===================="<<endl;
			cout<<"===================="<<endl;
		cross.TestPrint();
		}catch(Error& err)
		{
			err.Affiche();
			BOOST_CHECK(false);
		}

*/

	XmlParameters *pparam=new XmlParameters("test/test.xml");
	ublas::vector<double> telec,delec;
	ublas::vector<double> alt_grid;
	ElectronImpactIonization e_ioni(pparam,&delec,&telec,&alt_grid);
	BOOST_CHECK(true);
}



boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("ELECTRON IMPACT IONIZATION TEST");
	test->add(BOOST_TEST_CASE(&TestElectronImpact));
	return test;
}





