/**
 * \file photo/test.cpp
 * \brief The  test of the solar flux
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: test.cpp 745 2009-11-23 14:41:54Z gronoff $
 *
 */




#ifndef HAVE_CONFIG_H
#include "config.h"
#endif
#include "photoionization.hpp"
using namespace std;
using namespace boost::assign;


#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>


using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


void print_resu(ublas::vector<double> resu)
{
	cout<<"=============================="<<endl;
	cout<<"=============================="<<endl;
	cout<<"=============================="<<endl;
	cout<<endl<<endl<<endl;
	for(unsigned i=0;i<resu.size();++i)
	{
		cout<<resu[i]<<"\t";
	}
	cout<<endl<<endl<<endl;
	cout<<"=============================="<<endl;
	cout<<"=============================="<<endl;
	cout<<"=============================="<<endl;
}

void save_resu(string name,ublas::vector<double> energy, ublas::vector<double> resu)
{
	assert(energy.size()==resu.size());

	ofstream of(name.c_str());
	of.precision(9);
	of.setf(ios::scientific);
	of<<"# Flux grid "<<endl;
	of<<"# Energy in eV"<<endl;
	of<<"# Flux"<<endl;
	for(unsigned i=0;i<energy.size();++i)
	{
		of<<energy[i]<<"\t\t"<<resu[i]<<endl;
	}
	of.close();

}


void TestFlux()
{

	XmlParameters* param=new XmlParameters("./test/test.xml");
	Solar39Boxes solflux(param);
	ublas::vector<double> resu;
	vector<double> modelGrid,photonGrideVmax,photonGrideVmin;

	photonGrideVmax+=652.548421,
		413.280667,
		247.9684,
		123.9842,
		82.6561333,
		61.9921,
		49.59368,
		48.3846391,//
		48.3746391,
	       	43.6333627 ,
		41.3280667,
		40.8870565,//
		40.8770565,
		40.8238126,//
		40.8138126,
		35.4240571,
		33.6949512,//
		33.6849512,
		30.99605,
		27.5520444,
		26.6606599,//
		26.6506599,
		24.79684,
		22.5425818,
		22.3748827,//
		22.3648827,
		21.2281815,//
		21.2181815,
		20.6640333,
		20.3332787,
		19.6984697,//
		19.6884697,
		19.0744923,
		17.7120286,
		17.6386701,//
		17.6286702,
		16.5312267,
		16.2139077,//
		16.2039077,
		16.103275,//
		16.093275,
		15.7169271,//
		15.7069271,
		15.498025,
		14.5863765,
		13.7760222,
		13.0509684,
		12.7000371,//
		12.6900371,
		12.39842,
		12.0975288,//
		12.0875288,
		12.0250207,// 
		12.0150207;
	photonGrideVmin+= 413.280667,247.9684,123.9842,82.6561333,61.9921,49.59368,48.3846391,48.3746391,43.6333627,41.3280667,40.8870565,40.8770565,40.8238126,40.8138126,35.4240571,33.6949512,33.6849512,30.99605,27.5520444,26.6606599,26.6506599,24.79684,22.5425818,22.3748827,22.3648827,21.2281815,21.2181815,20.6640333,20.3332787,19.6984697,19.6884697,19.0744923,17.7120286,17.6386701,17.6286702,16.5312267,16.2139077,16.2039077,16.103275,16.093275,15.7169271,15.7069271,15.498025,14.5863765,13.7760222,13.0509684,12.7000371,12.6900371,12.39842,12.0975288,12.0875288,12.0250207,12.0150207,11.8080191;

	for(unsigned i=0;i<photonGrideVmin.size();++i)
	{
		modelGrid.push_back((photonGrideVmax[i]+photonGrideVmin[i])*0.5);
	}
/*
		21.2281815,48.3846391,40.8870565,40.8238126,33.6949512,26.6606599,22.3748827,19.6984697,17.6386701,16.2139077,16.103275,15.7169271,12.7000371,12.0975288,12.0250207, 
		48.3746391, 40.8770565, 40.8138126, 33.6849512, 26.6506599, 22.3648827, 21.2181815, 19.6884697, 17.6286702, 16.2039077, 16.093275, 15.7069271, 12.6900371, 12.0875288, 12.0150207;
*/

	ublas::vector<double> photonGrideVmin2(photonGrideVmin.size());
	std::copy(photonGrideVmin.begin(),photonGrideVmin.end(),photonGrideVmin2.begin());
	ublas::vector<double> photonGrideVmax2(photonGrideVmax.size());
	std::copy(photonGrideVmax.begin(),photonGrideVmax.end(),photonGrideVmax2.begin());

	ublas::vector<double> modelGrid2(modelGrid.size());
	std::copy(modelGrid.begin(),modelGrid.end(),modelGrid2.begin());
	
	solflux.RetrieveFlux(1.,modelGrid2,photonGrideVmin2,photonGrideVmax2,resu);
	print_resu(resu);
	save_resu("standard_grid.dat",modelGrid2,resu);
	// Premiere grille exponentielle
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
	//resu.erase(resu.begin(),resu.end());
	resu.clear();
	ublas::vector<double> autre_grid2(autre_grid.size());
	std::copy(autre_grid.begin(),autre_grid.end(),autre_grid.begin());
	ublas::vector<double> autre_grid_max2(autre_grid_max.size());
	std::copy(autre_grid_max.begin(),autre_grid_max.end(),autre_grid_max2.begin());
	solflux.RetrieveFlux(1.,autre_grid2,autre_grid_min,autre_grid_max2,resu);
	save_resu("exp_grid.dat",autre_grid2,resu);

	cout<<"seconde grille exponentielle"<<endl;
	ublas::vector<double> seconde_grid_min=MathGrid::GridExp(10,1000,100);
	std::reverse(seconde_grid_min.begin(),seconde_grid_min.end());
	vector<double> seconde_grid_max,seconde_grid;
	seconde_grid_max.push_back(2*seconde_grid_min[0]-seconde_grid_min[1]);

	for(unsigned i=0;i<seconde_grid_min.size()-1;++i)
	{
		seconde_grid_max.push_back(seconde_grid_min[i]);
	}

	vector<double> linemin,linemax;
	linemin+=
48.3846391,40.8870565,40.8238126,33.6949512,26.6606599,22.3748827,21.2281815,19.6984697,17.6386701,16.2139077,16.103275,15.7169271,12.7000371,12.0975288,12.0250207;

	linemax+=48.3746391, 40.8770565, 40.8138126, 33.6849512, 26.6506599, 22.3648827, 21.2181815, 19.6884697, 17.6286702, 16.2039077, 16.093275, 15.7069271, 12.6900371, 12.0875288, 12.0150207;

	/*
	cout<<"test avec des lignes plus ecartees :"<<endl;
	for(vector<double>:: iterator it=linemin.begin();it!=linemin.end();++it)
	{
		*it+=0.1;
	}
	for(vector<double>:: iterator it=linemax.begin();it!=linemax.end();++it)
	{
		*it-=0.1;
	}
	*/



	cout<<"Vecteurs des lignes effectues"<<endl;
	assert(linemin.size()==linemax.size());
	cout<<"We add the grid lines"<<endl;


	ublas::vector<double> second_grid_max2(seconde_grid_max.size());
	std::copy(seconde_grid_max.begin(),seconde_grid_max.end(),second_grid_max2.begin());

	ublas::vector<double> linemin2(linemin.size());
	std::copy(linemin.begin(),linemin.end(),linemin2.begin());
	ublas::vector<double> linemax2(linemax.size());
	std::copy(linemax.begin(),linemax.end(),linemax2.begin());

	MathGrid::AddGridLines(linemin2,linemax2,seconde_grid_min,second_grid_max2);
	cout<<"We compute the median line"<<endl;




	for(unsigned i=0;i<seconde_grid_min.size();++i)
	{
		seconde_grid.push_back(seconde_grid_min[i]+(seconde_grid_max[i]-seconde_grid_min[i])*0.5);
	}
	//resu.erase(resu.begin(),resu.end());

	resu.clear();

	ublas::vector<double> second_grid2(seconde_grid.size());
	std::copy(seconde_grid.begin(),seconde_grid.end(),second_grid2.begin());
	solflux.RetrieveFlux(1.,second_grid2,seconde_grid_min,second_grid_max2,resu);
	save_resu("seconde_exp_grid.dat",second_grid2,resu);




}


boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("PHOTOIONIZATION TEST");
	test->add(BOOST_TEST_CASE(&TestFlux));
	return test;
}




