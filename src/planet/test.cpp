/**
 * \file planet/test.cpp
 * \brief Test the  planets
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: test.cpp 745 2009-11-23 14:41:54Z gronoff $
 */

#ifndef HAVE_CONFIG_H
#include "config.h"
#endif

#include "allplanets.hpp"
using namespace std;
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

void TestVenus()
{
	XmlParameters* p_the_param=new XmlParameters("test/test.xml");
	Venus venus(p_the_param);

	ublas::vector<double> mes_altitudes(4);
	mes_altitudes[0]=(100);
	mes_altitudes[1]=(200);
	mes_altitudes[2]=(300);
	mes_altitudes[3]=400;

	try
	{
		cout<<"Machin"<<endl;
		ublas::vector<double> mes_densites_elec=venus.TheisModelNecm_3(mes_altitudes,0.);
		cout<<"The electron density"<<endl;
		//MathString::print1d(mes_densites_elec);
		cout<<mes_densites_elec<<endl;
	}
	catch(Error &err)
	{
		err.Affiche();
		BOOST_CHECK(false);
	}
	catch(...)
	{
		cout<<"blablabla"<<endl;
	}

	BOOST_CHECK(true);
}


void TestMars()
{
	cout<<MarsData::MarinerConditions()<<endl;
}

boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("PLANET TEST");
	test->add(BOOST_TEST_CASE(&TestVenus));
	test->add(BOOST_TEST_CASE(&TestMars));
	return test;
}





