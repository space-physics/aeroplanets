/** 
 * \file eflux/test/test.cpp
 * \brief Test the EFlux class 
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: test.cpp 745 2009-11-23 14:41:54Z gronoff $
 *
 */



#include "eflux.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
using namespace std;

using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;



void TestFlux()
{
	MathFunction::GaussianAngle angle(8);
	cout<<"Premiere grille exponentielle"<<endl;
	ublas::vector<double> autre_grid_min=MathGrid::GridExp(10,1000,100);
	std::reverse(autre_grid_min.begin(),autre_grid_min.end());
	vector<double> autre_grid_max,autre_grid;
	autre_grid_max.push_back(2*autre_grid_min[0]-autre_grid_min[1]);
	for(unsigned i=0;i<autre_grid_min.size()-1;++i)
	{
		autre_grid_max.push_back(autre_grid_min[i]);
	}
	for(unsigned i=0;i<autre_grid_min.size();++i)
	{
		autre_grid.push_back(autre_grid_min[i]+(autre_grid_max[i]-autre_grid_min[i])*0.5);
	}

	XmlParameters* params=new XmlParameters("./test/test.xml");


//	ublas::vector<double> autre_grid_min2(autre_grid_min.size());
//	std::copy(autre_grid_min.begin(),autre_grid_min.end(),autre_grid_min2.begin());
	ublas::vector<double> autre_grid2(autre_grid.size());
	std::copy(autre_grid.begin(),autre_grid.end(),autre_grid2.begin());

	ublas::vector<double> autre_grid_max2(autre_grid_max.size());
	std::copy(autre_grid_max.begin(),autre_grid_max.end(),autre_grid_max2.begin());

	EFlux flux(params,&angle,&autre_grid_min,&autre_grid2,&autre_grid_max2);
	Log::mL<<"fin"<<endl;


}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("SPECIES TEST");
	test->add(BOOST_TEST_CASE(&TestFlux));
	return test;
}








