/** 
 * \file chem/test.cpp
 * \ingroup Chem
 * \brief Test the chemistry class
 * Copyright G Gronoff Feb 2010
 * Last Nodification $Id: test.cpp 888 2010-02-22 23:21:04Z gronoff $
 */

#include "chem.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace std;

using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

void TestChem()
{
	BOOST_CHECK(true);
}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("CHEM TEST");
	test->add(BOOST_TEST_CASE(&TestChem));
	return test;
}







