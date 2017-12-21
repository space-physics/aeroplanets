/** 
 * \file species/test.cpp
 * \brief Test the specie class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: test.cpp 1299 2011-08-25 18:49:38Z gronoff $
 *
 */



#include "crs.hpp"
#include "species.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>


using namespace std;
using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;



void TestCrs()
{// We check the crs loading and writing



	ublas::vector<double> grid=MathGrid::GridExp(1.,1000.,500);
	CrossSection cross(grid);


	cout<<"Load CO2"<<endl;
	//cross.LoadCrs("./crs_test.xml","CO2");
	BOOST_CHECK(cross.LoadCrs("./crs_test.xml","CO2",false));
	cout<<"Fin load CO2"<<endl;
	cross.PrintCrs("testcrosssecforpython.dat");
	cout<<"Fin test crs"<<endl;



}

void TestElecCrs()
{
	ublas::vector<double> eCentEeV=MathGrid::GridExp(.1,1600,30);
	ublas::vector<double> eDdengeV=MathGrid::WidthGrid(eCentEeV,1);


	try
	{
		ElecCrossSection cross(eCentEeV,eDdengeV);
		cout<<"===================="<<endl;
		cout<<"======Elec Crs Loading======="<<endl;
		cout<<"===================="<<endl;
		cross.LoadCrs("./eleccrs_test.xml","CO2",true);
		cout<<"===================="<<endl;
		cout<<"===================="<<endl;
		cout<<"===================="<<endl;
	}catch(Error& err)
	{
		err.Affiche();
		BOOST_CHECK(false);
	}
	BOOST_CHECK(true);
	cout<<"Fin test elec crs"<<endl;



}

void TestSpecie()
{
	ublas::vector<double> ph_grid=MathGrid::GridExp(10.,700.,30);
	ublas::vector<double> el_grid;
	ublas::vector<double> pr_grid;

	Specie *mysp = new Specie("CO2",&ph_grid,&el_grid,&el_grid,&pr_grid,"specie_test.xml",false);
	BOOST_CHECK(mysp->mName=="CO2");


	Specie *mysp2 = new Specie("O+",&ph_grid,&el_grid,&el_grid,&pr_grid,"specie_test.xml",false);
	cout<<"le nom important"<<mysp2->mName<<endl;
	BOOST_CHECK(mysp2->mName=="O+");


}


void TestElecCrsCompar()
{
	ublas::vector<double> eCentEeV;//=MathGrid::GridExp(0.1,1600,30);
	ublas::vector<double> eDdengeV;//=MathGrid::WidthGrid(eCentEeV,1);

	double mSpfactor=0;

	if(!MathGrid::GridPolo(30,0.1,1600,eCentEeV,eDdengeV,mSpfactor))
	{
		BOOST_CHECK(false);
	}
	std::reverse(eCentEeV.begin(),eCentEeV.end());
	std::reverse(eDdengeV.begin(),eDdengeV.end());

	ElecCrossSection cross(eCentEeV,eDdengeV);
	cout<<"Print cente :"<<endl;
	//MathString::print1d(eCentEeV);
	cout<<eCentEeV<<endl;
	cout<<"print engdd"<<endl;
	//MathString::print1d(eDdengeV);
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
	XmlParameters compar("./comparison_transfortran.xml");
	cout<<"Comparison mCrsCincm2"<<endl;
	ublas::vector<double>  crscin;
	compar.Get1DArray("/resutrans/cin",crscin);
	std::reverse(crscin.begin(),crscin.end());

	cout<<crscin<<endl;
	cout<<cross.mCrsCincm2<<endl;
//	BOOST_CHECK(MathFunction::VectorCompare(cross.mCrsCincm2,crscin,30));
	cout<<"Comparison OmDegrad (whaaaa)"<<endl;
	//vector< vector<double> > degrad;
	ublas::matrix<double> degrad;
	compar.Get2DArray("/resutrans/omdeg",degrad);
//BOOST_CHECK(MathFunction::VectorCompare2D(cross.mOmDegradcm2,degrad,16));
//	BOOST_CHECK(MathFunction::VectorCompare2D(cross.mOmDegradcm2,degrad,20));
//	cout<<"Here, ~12 values have a 15\% deviation with the fortran, the other <1%. All the values are in a 16\% deviation from the fortran  "<<endl;
	//cout<<"A FEW WARNING IS NORMAL. ON MY COMPUTER, THERE ARE ~30 WARNINGS (compare with 900 values, with two different codes!)."<<endl;
	cout<<"Comparison cel "<<endl;
	ublas::vector<double> crscel;
	compar.Get1DArray("/resutrans/cel",crscel);
	std::reverse(crscel.begin(),crscel.end());
	BOOST_CHECK(MathFunction::VectorCompare(cross.mElasticCrscm2,crscel,20));



	/*
	   cout<<"Print cin "<<endl;
	   MathString::print1d(cross.mCrsCincm2);
	   cout<<"Print omdeg"<<endl;
	   MathString::print2d(cross.mOmDegradcm2);
	   BOOST_CHECK(true);
	   cout<<"Fin test elec crs"<<endl;

*/

}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("SPECIES TEST");
	test->add(BOOST_TEST_CASE(&TestCrs));
	test->add(BOOST_TEST_CASE(&TestSpecie));
	test->add(BOOST_TEST_CASE(&TestElecCrs));
	test->add(BOOST_TEST_CASE(&TestElecCrsCompar));
	return test;
}



