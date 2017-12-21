#include "protocos.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace std;

using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


/// Dummy test for now
void TestCosmic()
{
	BOOST_CHECK(true);
}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("SPECIES TEST");
	test->add(BOOST_TEST_CASE(&TestCosmic));
	return test;
}





