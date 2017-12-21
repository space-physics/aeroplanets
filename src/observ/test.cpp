/** 
 * \file observ/test.cpp
 * \ingroup GeoObserv
 * \brief Test the  observation system
 * Copyright G Gronoff Feb 2010
 * Last Nodification $Id: test.cpp 1595 2012-11-03 01:30:26Z gronoff $
 */

#include "geoobs.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace std;

using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

void TestObserv()
{
	Log::Init("Log.txt");

	ublas::vector<double> altGridKm = MathGrid::GridExp(100, 800, 21);
	cout<<"==================================="<<endl;
	cout<<"45 degree"<<endl<<"-----------"<<endl;
	GeoPath gp(altGridKm, 45, 6300);
	cout<<"==================================="<<endl;
	cout<<"85 degree"<<endl<<"-----------"<<endl;
	GeoPath gp2(altGridKm, 85, 6300);
	cout<<"==================================="<<endl;
	cout<<"90 degree"<<endl<<"-----------"<<endl;
	GeoPath gp3(altGridKm, 90, 6300);
	cout<<"==================================="<<endl;
	cout<<"95 degree"<<endl<<"-----------"<<endl;
	GeoPath gp4(altGridKm, 95, 6300);
	cout<<"==================================="<<endl;
	cout<<"100 degree"<<endl<<"-----------"<<endl;
	GeoPath gp45(altGridKm, 100, 6300);
	cout<<"==================================="<<endl;
	cout<<"102 degree"<<endl<<"-----------"<<endl;
	GeoPath gp46(altGridKm, 102, 6300);


	cout<<"==================================="<<endl;
	cout<<"105 degree"<<endl<<"-----------"<<endl;
	GeoPath gp5(altGridKm, 105, 6300);
	cout<<"==================================="<<endl;
	cout<<"120 degree"<<endl<<"-----------"<<endl;
	GeoPath gp6(altGridKm, 120, 6300);
	cout<<"==================================="<<endl;
	cout<<"180 degree"<<endl<<"-----------"<<endl;
	GeoPath gp7(altGridKm, 180, 6300);
	/*
	// We create a GeoPoint for a "Earth"
	GeoPoint go(6300); // Earth: R = 6300km
	go.SetGeoPosition(100, 0, 0); // Alt = 100km, lat =0,long =0
	cout << go.GetDist() << endl;
	GeoPoint end = go.ReturnAllSortieAtmo(0, 2, 800);
	cout << end.GetDist() << endl;
	cout<< (end - go).GetDist() << endl;
	GeoPath gp(go, end, 420);

	ublas::vector<double> altGridKm = MathGrid::GridExp(100, 800, 210);
	Log::mL<<altGridKm<<endl;
	gp.TestNewPath(altGridKm);

	gp.ResetGrid(altGridKm);
	*/
	BOOST_CHECK(true);
}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("OBSERV TEST");
	test->add(BOOST_TEST_CASE(&TestObserv));
	return test;
}


