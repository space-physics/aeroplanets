/**
 * \file leastsq/test/test.cpp
 * \brief Test for the least square system
 * Copyright G Gronoff July 2010
 * Last Modification : $Id$
 */


#ifndef HAVE_CONFIG_H
#include "config.h"
#endif
using namespace std;
#include "../least.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
#include <deque>
using namespace boost::assign;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;




class Testfnc : public LMlsq
{
	public:
		Testfnc(std::deque< ublas::vector<double> > vX, std::deque< ublas::vector<double> > vY):LMlsq(vX,vY)
		{
		}
		
		int Exec(int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag);


};

int Testfnc::Exec(int vDiffsize, int vParamsize, const double* x, double* fvec, int vFlag)
{ // Test in the std minpack
	int i;
	double tmp1, tmp2, tmp3;
	assert(vDiffsize==15);
	assert(vParamsize==3);

	if (vFlag == 0)
	{
		cout<<"Zut Zut et re-Zut"<<endl;
		/*      insert print statements here when nprint is positive. */
		return 0;
	}
	for (i = 1; i <= 15; i++)
	{
		tmp1 = i;
		tmp2 = 16 - i;
		tmp3 = tmp1;
		if (i > 8) tmp3 = tmp2;
		fvec[i-1] = mYmeasu[0][i-1] - (x[1-1] + tmp1/(x[2-1]*tmp2 + x[3-1]*tmp3));
	}
	return 0;
}



void TestLsq()
{
	vector<double> measu;
	measu+= 1.4E-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39;
	ublas::vector<double> meas=StdToUblas(measu);
	std::deque< ublas::vector<double> > mm;
	mm.push_back(meas);
	Testfnc* leastobj=new Testfnc(mm, mm);

	ublas::matrix<double> comat;
	std::deque<double> params;
	//params+= 1.,1.,1.;
	
	params.push_back(1.);
	params.push_back(1.);
	params.push_back(1.);



	int ret=LstSq::LMleast(leastobj,-1,params,comat);


	cout<<" Return value : "<<ret<<endl<<" output params : "<<params[0]<<" "<<params[1]<<" "<<params[2]<<endl<<"Covariant matrix : "<<comat<<endl;


	delete leastobj;
	BOOST_CHECK( 1 == 1);
}

boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	Log::Init("SuperLog.txt",Log::DEBUGG,Log::DEBUGG);
	test=BOOST_TEST_SUITE("ATMO TEST");
	test->add(BOOST_TEST_CASE(&TestLsq));
	return test;
}


