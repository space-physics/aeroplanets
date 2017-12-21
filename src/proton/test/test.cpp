/**
 * \file  proton/test/test.cpp
 * \brief  test the protoionization
 *	  Copyright G Gronoff Nov 2011
 *	  Last Modification : $Id$
 */

#ifndef HAVE_CONFIG_H
#include "config.h"
#endif

#include "proton.hpp"
using namespace std;
using namespace boost::assign;
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


void TestProtonImpact()
{
	BOOST_CHECK(true);
}




boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("Proton/Hydrogen IMPACT IONIZATION TEST");
	test->add(BOOST_TEST_CASE(&TestProtonImpact));
	return test;
}

