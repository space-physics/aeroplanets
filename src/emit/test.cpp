/** 
 * \file emit/test.cpp
 * \brief Test the emission class
 * Copyright G Gronoff March 2010
 * Last Nodification $Id: test.cpp 904 2010-03-02 02:01:51Z gronoff $
 */

#include "emission.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace std;

using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

void TestEmit()
{
	BOOST_CHECK(true);
}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("EMIT TEST");
	test->add(BOOST_TEST_CASE(&TestEmit));
	return test;
}

