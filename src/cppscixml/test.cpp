/** 
 * \file cppscixml/test.cpp
 * \brief Test the scixml class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: test.cpp 1111 2010-08-12 19:43:33Z gronoff $
 *
 */



#ifndef HAVE_CONFIG_H
#include "config.h"
#endif
using namespace std;
#include "scixml.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>


using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;


void TestLoad()
{
	//Log::SetPriority(Log::DEBUGG);
	XmlParameters montest("./test.xml");
	Log::mD<<"on va maintenant faire le test du nombre"<<endl;
	unsigned nombredea=montest.Numbers("/test/a");
	cout<<"fin de ce test"<<endl;
	BOOST_CHECK(nombredea == 3);
	Log::mD<<"fin de check"<<endl;
	Log::mD<<"Le texte dans test:"<<montest.Elem("/test")<<endl;;
	Log::mD<<"Le texte dans test:"<<montest.Process("/test/text()")<<endl;;
	
	bool machin=montest.Exists("/test/b");
	BOOST_CHECK(!machin);
	
	for(unsigned i=0;i<nombredea;i++)
	{
		Log::mD<<"a["<<i<<"] :"<<montest.Elem("/test/a",i)<<endl;
	}
	
	
	
	Log::mL<<"Bonne nuit"<<endl;

}

void TestKey()
{
	//Log::SetPriority(Log::DEBUGG);
	XmlParameters montest("./test.xml");
	unsigned nbtestkey=montest.NbParams("/testkey","fact");
	BOOST_CHECK(nbtestkey == 1 );
	BOOST_CHECK(montest.KeyExists("/testkey","fact"));
	BOOST_CHECK(!montest.KeyExists("/testkey","prout"));
	nbtestkey=montest.NbParams("/testkey","fact",2);
	BOOST_CHECK(nbtestkey == 0 );
	Log::mL<<"on ressort la clef"<<endl;
	
	Log::mL<<montest.Process("/testkey/attribute::fact")<<endl;
	string machin= montest.GetKey("/testkey","fact");
	BOOST_CHECK(machin == "3");
	Log::mL<<"On travaille avec une clef numerique"<<endl;
	int key;
	montest.GetNKey("/testkey","fact",key);
	BOOST_CHECK(key == 3);
	Log::mL<<"On travaille avec la valeur attendue"<<endl;
	montest.GetValue("/testkey",key);
	BOOST_CHECK(key == 126);

//	Log::SetPriority(Log::CONFIG);

	Log::mL<<"On choppe les testk"<<endl;
	vector<TiXmlNode*> mestk=montest.GetNodes("/testk");
	BOOST_CHECK(mestk.size()==6);
	Log::mL<<"C est bon, il y en a 6"<<endl;
	for(unsigned i=0;i<mestk.size();++i)
	{
		Log::mL<<"numero "<<i<<endl;
		Log::mL<<montest.GetKey(mestk[i],"/","a")<<endl;
		unsigned j;
		montest.GetNKey(mestk[i],"/","a",j);
		BOOST_CHECK(j==i+1);
	}
	BOOST_CHECK(true);



	Log::mL<<"fin test"<<endl;

}


void TestArray()
{

	XmlParameters montest("./test.xml");

	//vector< vector<int> > montbl;
	ublas::matrix<int> montbl;
	montest.Get2DArray("/testarr",montbl);
	//MathString::print2d(montbl);
	cout<<montbl<<endl;
	ublas::vector<int> montbl2;
	montest.Get1DArray("/testa",montbl2);
	//MathString::print1d(montbl2);
	cout<<montbl2<<endl;

	//Log::SetPriority(Log::DEBUGG,"TestArray()");
	Log::mL<<"test_path : "<<JoinPath("./test","machin.xml")<<endl;
	Log::mL<<"test_path_2 : "<<FilenameToPath("./test/machin.xml")<<endl;

	Log::mL<<"Check if the present file exists"<<endl;
	BOOST_CHECK(FileExists("./test.cpp"));
// Pas tres malin comme test...
	//	BOOST_CHECK(file_exists("/home/gronoff"));
	BOOST_CHECK(!FileExists("./testanepasmettre/test.cpp"));
	Log::mL<<"testrootpath : "<<FullPath("/home/gronoff/supertest/test.cpp")<<endl; 
	Log::mL<<"testrootpath : "<<FullPath("supertest/test.cpp")<<endl; 
	Log::mL<<"testrootpath : "<<FullPath("test.cpp")<<endl; 

	Log::mL<<"fin test"<<endl;
}


void TestNode()
{
	//Log::SetPriority(Log::DEBUGG);
	XmlParameters* montest;
	montest=new XmlParameters("./test2.xml");
	
	vector<TiXmlNode*> node_b=montest->GetNodes("/a/b");
	Log::mL<<"Nombre de nodes b : "<<node_b.size()<<endl;

	BOOST_CHECK(node_b.size()==3);

//	for(unsigned i=0;i<node_b.size();++i)
//	{
	vector<TiXmlNode*>::iterator nb;
	for(nb=node_b.begin();nb!=node_b.end();++nb)
	{
		//vector<TiXmlNode*> node_c=montest->get_nodes(node_b[i],"c");
		vector<TiXmlNode*> node_c=montest->GetNodes(*nb,"//bp/c");

		Log::mL<<"Nombre de nodes c : "<<node_c.size()<<endl;
		for(vector<TiXmlNode*>::iterator it=node_c.begin();it!=node_c.end();++it)
		{
			Log::mL<<montest->Elem(*it,"/")<<endl;
		}
	}



	delete montest;
}



/**
 *
 */
void TestPlus()
{// Test if an ionized specie can be retrieved without problems (the + at the end could have a meaning)
	XmlParameters* montest;
	montest=new XmlParameters("./test3.xml");
	
	BOOST_CHECK(montest->Exists("/specie/a"));
	BOOST_CHECK(montest->Exists("/specie/b"));
	BOOST_CHECK(!montest->Exists("/specie/c"));
/*
	BOOST_CHECK(montest->Exists("/specie/a+"));
	BOOST_CHECK(!montest->Exists("/specie/b+"));
	BOOST_CHECK(montest->Exists("/specie/c+"));
*/
	delete montest;
}



void TestRandom()
{
	XmlParameters *montest;
	montest=new XmlParameters("./test4.xml");
	montest->SetMonteCarloActive();
//	Log::SetPriority(Log::INFO);
	Log::mL<<"RANDOM ACTIVATED"<<endl;
	Log::mL<<"ONEvalue 100 30%"<<endl;
	for(unsigned i=0;i<10;++i)
	{
		double val;
		montest->GetValue("/onedimension",val);
	//	Log::SetPriority(Log::WARNING,"RESULTAT");
		Log::mL<<val<<endl;
	}
	Log::mL<<"ONEvalue2 10 +/-10 (1 sigma)"<<endl;
	for(unsigned i=0;i<10;++i)
	{
		double val;
	//	Log::SetPriority(Log::WARNING,"RESULTAT");
		montest->GetValue("/onedimension2",val);
		Log::mL<<val<<endl;
	}
	Log::mL<<"ARRAY"<<endl;
	ublas::vector<double> vals;
	montest->Get1DArray("/array",vals);
	cout<<"SIZE of vals : "<<vals.size()<<endl;
	//MathString::print1d(vals);
	cout<<vals<<endl;



	Log::mL<<"Test du facteur non defini"<<endl;
	ublas::vector<double> arr2;
	montest->Get1DArray("/arr2",arr2);
	cout<<endl<<arr2<<endl;
	BOOST_CHECK(nabs(arr2[0]-arr2[1])<1E-10);

	Log::mL<<"Test du facteur defini"<<endl;
	ublas::vector<double> arr3;
	montest->Get1DArray("/arr3",arr3);
	cout<<endl<<arr3<<endl;
	BOOST_CHECK(nabs(arr3[0]-arr3[1])<1E-10);


	// Test with the one dimension parameter
/*	ofstream of("./random30percent");
	for(unsigned i=0;i<10000;++i)
	{
		double val;
		montest->GetValue("/onedimension",val);
		of<<val<<endl;;
	}
	of.close();
*/
}

// Test the pointeur thing
void TestAppend()
{

	XmlParameters* file1=XmlParameters::AttachFileParameter("test4.xml");
	XmlParameters* file2=XmlParameters::AttachFileParameter("test4.xml");
	XmlParameters* file3=XmlParameters::AttachFileParameter("test4.xml");

	file3->DetachFileParameter();
	file1->DetachFileParameter();
	file2->DetachFileParameter();
	XmlParameters* file4=XmlParameters::AttachFileParameter("test4.xml");
	file4->DetachFileParameter();
	cout<<"Ouf"<<endl;
	BOOST_CHECK(!XmlParameters::DoINeedAGarbageCollector());


}


/// On initialise les tests boost:
boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{

	// We want less messages on the screen, but more information if we have to debug
	Log::Init("TestLog.txt",Log::CONFIG,Log::DEBUGG);
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("SCIXML TEST");
	test->add(BOOST_TEST_CASE(&TestLoad));
	test->add(BOOST_TEST_CASE(&TestKey));
	test->add(BOOST_TEST_CASE(&TestArray));
	test->add(BOOST_TEST_CASE(&TestNode));
	test->add(BOOST_TEST_CASE(&TestPlus));
	test->add(BOOST_TEST_CASE(&TestRandom));

	test->add(BOOST_TEST_CASE(&TestAppend));



//	test->add(BOOST_TEST_CASE(&MathString_test_lit_string));
//	test->add(BOOST_TEST_CASE(&MathString_test_lit_tout));
//	test->add(BOOST_TEST_CASE(&MathString_test_lit_tableau));
	return test;
}


