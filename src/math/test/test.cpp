/**
 * \file math/test/test.cpp
 * \brief test math functions
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: test.cpp 1305 2011-09-09 18:37:19Z gronoff $
 */




#ifdef HAVE_CONFIG_H
#include "config.h"

#endif


using namespace std;

#include "../mathstring.hpp"
#include "../mathfunction.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>

using namespace boost::assign;
using namespace MathString;


using namespace boost::unit_test;
using boost::unit_test_framework::test_suite;

void MathString_test_lit_string()
{
	string karde="42";
	string pitq = "3.14159";
	string testfou="0x10";


	int kard;
	double pit;
	int testf;
	BOOST_CHECK(LitString(karde,kard));
	BOOST_CHECK(LitString(pitq,pit));
	BOOST_CHECK(LitString(testfou,testf));

	BOOST_CHECK_EQUAL(kard,42);//(true);
	BOOST_CHECK(kard == 42);
	BOOST_CHECK(pit == 3.14159);
	
	cout<<"Nous avons nos 3 valeurs : "<<kard<<" et "<<pit<<" et "<<testf<<endl;


}


void MathString_test_lit_tout()
{
	string vals="32 11 32.44 33.12 \n 42.42";
	ublas::vector<int> intvals;
	ublas::vector<double> doublevals;
	BOOST_CHECK(LitToutString(vals,intvals));
	BOOST_CHECK(LitToutString(vals,doublevals));

	BOOST_CHECK_EQUAL(intvals.size(),3); // Attention, le second est aussi compte
	BOOST_CHECK_EQUAL(doublevals.size(),5);
	for(ublas::vector<int>::iterator it=intvals.begin();it!=intvals.end();it++)
	{
		cout<<*it<<endl;
	}

	cout<<"Et la suite : "<<endl;


	for(ublas::vector<double>::iterator it=doublevals.begin();it!=doublevals.end();it++)
	{
		cout<<*it<<endl;
	}



}


void MathString_test_lit_tableau()
{
	cout<<"Lit tableau"<<endl;
	string tableau1 = "\n0.1 0.2 0.3\n0.2 0.3 0.4\n	0.3 0.4 0.5\n";

	string tableau0 = "0.1 0.2 0.3\n0.2 0.3 0.4\n0.3 0.4 0.5\n";

	ublas::vector<double> vector0,vector1;
	BOOST_CHECK(LitTableau(tableau0,vector0,1));
	BOOST_CHECK(LitTableau(tableau1,vector1,1));
	ublas::vector<double> verif(3);
	verif(0)=0.1;
	verif(1)=0.2;
	verif(2)=0.3;
	BOOST_CHECK_EQUAL(vector0[0],0.2);
	BOOST_CHECK_EQUAL(vector1[0],0.2);

	for(ublas::vector<double>::iterator it=vector0.begin();it!=vector0.end();it++)
	{
		cout<<*it<<endl;
	}


	ublas::vector<double> vector2,vector3,vector4,vector5;

	BOOST_CHECK(LitTableauDouble(tableau0,vector2,1,vector3,2));
	BOOST_CHECK_EQUAL(LitTableauDouble(tableau0,vector4,1,vector5,5),false);

	cout<<"Nouvelle version : "<<endl;
	for(ublas::vector<double>::iterator it=vector2.begin();it!=vector2.end();it++)
	{
		cout<<*it<<endl;
	}
}


void MathString_test_trim()
{
	string machin="   salut  ";
	cout<<"la chaine ["<<machin<<"]"<<endl;
	string machin2=trim(machin);
	cout<<"la chaine ["<<machin2<<"]"<<endl;
	BOOST_CHECK(machin2 == "salut");

}



void MS_test_grid()
{

	cout<<"****************gridpow*****************"<<endl;
	ublas::vector<double> a=MathGrid::GridPow(12.,42.,10);
	ublas::vector<double> b=MathGrid::GridPow(42.,12.,10);
	cout<<a<<endl<<b<<endl;
	//MathString::print1d(b);
//	MathString::print1d(a);


	cout<<"****************gridcst*****************"<<endl;
	ublas::vector<double> c,d;
	c=MathGrid::GridCst(10.,42.,10);
	d=MathGrid::GridCst(42.,10.,10);
	cout<<c<<endl<<d<<endl;
	cout<<"****************gridpolo*****************"<<endl;

	//bool gridpolo(int ntab,double tmin,double tmax,std::vector<double>& tab,std::vector<double>&widthtab,double&spfac);

	ublas::vector<double> aa,ab;
	double asp;


	BOOST_CHECK(!MathGrid::GridPolo(10,1000,0.1,aa,ab,asp));
	if(MathGrid::GridPolo(10,1000,0.1,aa,ab,asp))
	{
		cout<<"Bizarre : "<<endl;
		cout<<aa<<endl;
	}
	
	if(MathGrid::GridPolo(100,5,1000,aa,ab,asp))
	{
		cout<<"Premiere grille <-"<<asp<<endl;
		//MathString::print1d(aa);
		cout<<aa<<endl;
	}
	ublas::vector<double> ba,bb;
	double bsp;
	if(MathGrid::GridPolo(10,10,1000,ba,bb,bsp))
	{
		cout<<"Seconde grille <-"<<bsp<<endl;
		cout<<ba<<endl;
		BOOST_CHECK(false);
	}


	cout<<"****************widthgrid*****************"<<endl;
	ublas::vector<double> gg;
	gg=MathGrid::WidthGrid(MathGrid::GridExp(10,42,10),1);
	//MathString::print1d(a);
	//MathString::print1d(gg);
	cout<<gg<<endl;
	BOOST_CHECK(true);



}


void Gaussangle()
{
	MathFunction::GaussianAngle test(8);
	cout<<"--------------------xmu-------------------"<<endl;
	cout<<test.mXmu<<endl;
	cout<<"------------------weight------------------"<<endl;
	cout<<test.mWeight<<endl;
	cout<<"------------------angzb-------------------"<<endl;
	cout<<test.mAngzb<<endl;
	cout<<"fin gaussianangle"<<endl;
}


void Interpol()
{
	vector<double> ox,oy,onewx;
	ox.push_back(1);
	ox.push_back(2);
	ox.push_back(3);
	ox.push_back(4);
	oy.push_back(1);
	oy.push_back(2);
	oy.push_back(3);
	oy.push_back(4);
	onewx.push_back(0.1);
	onewx.push_back(0.5);
	onewx.push_back(0.7);
	onewx.push_back(1);
	onewx.push_back(1.5);
	onewx.push_back(2.5);
	onewx.push_back(3.5);
	onewx.push_back(4.5);
	onewx.push_back(5.5);
	ublas::vector<double> x=StdToUblas(ox);
	ublas::vector<double> y=StdToUblas(oy);
	ublas::vector<double> newx=StdToUblas(onewx);

	ublas::vector<double> liny=MathFunction::IntLin(x,y,newx);
	cout<<"____________________intlin___________________"<<endl;
	cout<<liny<<endl;
	cout<<"____________________intlog___________________"<<endl;
	ublas::vector<double> logy=MathFunction::IntLog(x,y,newx);
	cout<<logy<<endl;

	cout<<"____________________intloglog___________________"<<endl;
	ublas::vector<double> logy2=MathFunction::IntLogLog(x,y,newx);
	cout<<logy2<<endl;


	cout<<"Fin intlinog"<<endl;
}




void AddVector()
{
	vector<double> supertest;
	supertest+= 0,1,2,3,4,5,6;
	vector<double> supertost;
	supertost+= 0,1,2,3,4,5,6;
	cout<<supertest.size()<<endl;
	BOOST_CHECK( supertest.size()==7);

	ublas::vector<double> fintest=StdToUblas(supertest)+StdToUblas(supertost);
	for(ublas::vector<double>::iterator it=fintest.begin();it!=fintest.end();++it)
	{
		cout<<*it<<endl;
	}
}



// To check if the two double vectors are equal at 1E-3
bool EqualVectors(const ublas::vector<double> & a, const ublas::vector<double> & b, double precision=0.01)
{
	if(a.size()!=b.size())
	{
		return false;
	}

	for(unsigned i=0;i<a.size();++i)
	{
		if(nabs(a[i]-b[i])>precision)
		{
			cout<<"Position "<<i<<" -> "<<a[i]<<" != "<<b[i]<<endl;
			cout<<nabs(a[i]-b[i])<<endl;
			return false;
		}
	}
	return true;
}


void TestRedistribute()
{
	vector<double> min1,max1,val;
	min1+=1.5,2.5,3.5,3.9,4.5,5.5,5.9;
	max1+=2.5,3.5,4.5,4.1,5.9,5.6,6.9;
	val+=10.0,10.,10.,10.,10.,10.,10.;
	vector<double> min2,max2;
	min2+=1,2,3,4,5,6;
	max2+=2,3,4,5,6,7;

	double reste=0;
	ublas::vector<double> resu=MathGrid::RedistributeChaoticFlux(StdToUblas(val),StdToUblas(min1),StdToUblas(max1),StdToUblas(min2),StdToUblas(max2),reste);

	cout<<"Resu"<<endl;
	cout<<resu<<endl;
	cout<<endl<<"Reste : "<<reste<<endl;

	vector<double> verification;
	verification+=5.,10.,15.,13.5714285714285716 ,17.42857142857142838,9;

	BOOST_CHECK( EqualVectors(resu,StdToUblas(verification)) );

}




void TestAddLines()
{
	vector<double> testmin,testmax;
	vector<double> linemin,linemax;
	testmin+=10,20,30,40,50;
	testmax+=20,30,40,50,60;
	linemin+=15,29,32,49,62;
	linemax+=16,31,33,51,63;

	ublas::vector<double> ulinemin=StdToUblas(linemin);
	ublas::vector<double> ulinemax=StdToUblas(linemax);
	ublas::vector<double> utestmin=StdToUblas(testmin);
	ublas::vector<double> utestmax=StdToUblas(testmax);
	MathGrid::AddGridLines(ulinemin,ulinemax,utestmin,utestmax);
	cout<<utestmin<<endl;


	cout<<"AUTRE TEST"<<endl;






	vector<double> photonGrideVmax,photonGrideVmin;

	// Modifications 
	// 48.3746   -> 48.3836
	// 40.877056 -> 40.887055
	// 40.8138 -> 40.8238
	// 33.684 -> 33.694
	// 26.6606 -> 26.6706
// 22.3648826596 22.3648826596
//21.218181507 21.218181507
//20.333278667 20.333278667
//19.6884696616 19.6884696616
//17.6286701455 17.6286701455
//16.2039077305 16.2039077305
//16.0932750094 16.0932750094
//15.7069271308 15.7069271308
//12.6900370514 12.6900370514
//12.0875287603 12.0875287603
//12.0150206898 12.0150206898

	photonGrideVmax+=652.54842105,413.28066667,247.9684,123.9842,82.65613333,61.9921,48.38463909,43.63336266,49.59368,40.88705648,40.82381263,41.32806667,33.69495123,35.42405714,30.99605,26.6606599,27.55204444,24.79684,22.37488266,21.22818151,22.54258182,20.33327867,19.69846966,20.66403333,19.07449231,17.63867015,17.71202857,16.21390773,16.10327501,15.71692713,16.53122667,15.498025,14.58637647,13.77602222,12.70003705,13.05096842,12.09752876,12.02502069,12.39842;
	photonGrideVmin+=413.28066667,247.9684,123.9842,82.65613333,61.9921,49.59368,48.37463909,43.63336266,41.32806667,40.87705648,40.81381263,35.42405714,33.68495123,30.99605,27.55204444,26.6506599,24.79684,22.54258182,22.36488266,21.21818151,20.66403333,20.33327867,19.68846966,19.07449231,17.71202857,17.62867015,16.53122667,16.20390773,16.09327501,15.70692713,15.498025,14.58637647,13.77602222,13.05096842,12.69003705,12.39842,12.08752876,12.01502069,11.80801905;


	vector<double> grid_min,grid_max;
	grid_min+=10;
	grid_max+=700;

	ublas::vector<double> nphotonGrideVmin= StdToUblas(photonGrideVmin);
	ublas::vector<double> nphotonGrideVmax=StdToUblas(photonGrideVmax);
	ublas::vector<double> ngrid_min=StdToUblas(grid_min);
	ublas::vector<double> ngrid_max=StdToUblas(grid_max);

	MathGrid::AddGridLines(nphotonGrideVmin,nphotonGrideVmax,ngrid_min,ngrid_max);
	cout<<"Lines resu"<<endl;
	cout<<ngrid_min<<endl;
	cout<<"Lines max resu"<<endl;
	cout<<ngrid_max<<endl;


}




void TestStrToN()
{
	string test="42";

	int supertest=0;
	strton(test,supertest);

	cout<<"Supertest : "<<supertest<<endl;
	BOOST_CHECK(supertest==42);


}



void TestComparison()
{
	cout<<"Test comparison"<<endl;
	vector<double> alpha,beta,gamma;
	double truc=1.;
	alpha.push_back(1.);
	beta.push_back(1.0001);
	gamma.push_back(1.);
	for(unsigned i=0;i<50;++i)
	{
		alpha.push_back(truc);
		beta.push_back(truc);
		gamma.push_back(truc*truc);
		truc+=1.13;
	}
	BOOST_CHECK(MathFunction::VectorCompare(StdToUblas(alpha),StdToUblas(beta)));
	BOOST_CHECK(!MathFunction::VectorCompare(StdToUblas(alpha),StdToUblas(gamma)));
}

void TestRandom()
{
	// On initialize le bazard
	MathRandom::Initialize();

	cout<<"TEST RANDOM"<<endl;

	for(unsigned i=0;i<5;++i)
		cout<<MathRandom::GetNormal()<<endl;
	for(unsigned i=0;i<5;++i)
		cout<<MathRandom::GetUniformReal()<<endl;
	for(unsigned i=0;i<5;++i)
		cout<<MathRandom::GetNormal()<<endl;

	cout<<"FIN TEST RANDOM"<<endl;

	BOOST_CHECK(true);

}

void TestLog()
{
	// We initialize the log function
	Log::Init("TestLog.txt");
	//Log::SetPriority(Log::DEBUGG);
	Log::mD<<"Hi";
	Log::mD<<endl;;
	Log::mD<<"Hi, this is a debug message"<<endl;
//	Log::SetPriority(Log::CONFIG,"TestLog()");
	Log::mL<<"Hi, this is config message"<<endl;
//	Log::SetPriority(Log::INFO,"TestLog()");
	Log::mI<<"Hi, this is an informative message, please read ..."<<endl;
	Log::AddVersionInfo("This message was wrote for the first test...");
	Log::AddVersionInfo(static_cast<string>("He"));
//	Log::SetPriority(Log::WARNING,"TestLog()");
	Log::mW<<"Hi, I am a warning"<<endl;
//	Log::SetPriority(Log::DEBUGG,"TestLog()");
	Log::mW<<"Hi, we can't grow in the warnings like that"<<endl;
	Log::mW<<"He, I did't tried to make a lot of lines..."<<endl;
	Log::mW<<"Yet."<<endl;
	Log::mW<<"What an idea, an what if I try to do math : 2+3 = "<<(2+3)<<endl;
	Log::mW<<"Do you have some beer ???"<<endl;
	Log::AddVersionInfo(static_cast<string>("H, after He...."));
//	Log::SetPriority(Log::ERROR,"TestLog()");
	Log::mE<<"OOOOOPS, I am an error.... NO  JUST KIDDING, it is for the TEST!!!"<<endl<<endl;;
	cout<<"The final version info: "<<Log::msMessageLog<<endl;
	
	
	if(cerr == cout)
	{
		cout<<"Whaaaa extraordinaire"<<endl;
	}else
	{
		cout<<" Heu, est-ce extra ici????"<<endl;
	}
	
	BOOST_CHECK(true);
}

void TestUblas()
{
	vector<double> a(20);

	for(unsigned i=0;i<20;++i)
		a[i]=static_cast<double>(i);

	ublas::vector<double> blasa=StdToUblas(a);
	cout<<"Ublas vector : "<<blasa<<endl;

	vector< vector<double> > machin(3,a);

	ublas::matrix<double> mat=StdToUblas(machin);
	cout<<"Ublas matrix : "<<mat<<endl;
	vector< vector<double> > bidule(0);
	cout<<StdToUblas(bidule)<<endl;

	ublas::vector<double> premier(10),second(10),troisieme(10),quatrieme(10);



	for(unsigned i=0;i<10;++i)
	{
		premier(i)=2;
		second(i)=3;
	}
	troisieme=ublas::element_prod(premier,second);
	cout<<"multiplication 1"<<troisieme<<endl;
	for(unsigned i=0;i<10;++i)
	{
		premier(i)=(double)i;
		second(i)=(double)i;
	}
	troisieme=ublas::element_prod(premier,second);
	cout<<"multiplication 2"<<troisieme<<endl;
	cout<<"fin"<<endl;
	quatrieme=premier*second;

	BOOST_CHECK(troisieme==quatrieme);
	BOOST_CHECK(not (quatrieme==premier));




	BOOST_CHECK(true);
	Log::mL<<"Fin test Ublas"<<endl;
}



void TestPathInterpolation()
{
	ublas::vector<double> ialt(5),values(5),outalt(5),outsza(5);
	std::deque<double> isza;
	std::deque< ublas::vector<double>* > ivalues;

	ialt[0]=0;
	ialt[1]=10;
	ialt[2]=20;
	ialt[3]=30;
	ialt[4]=40;
	values[0]=0;
	values[1]=10;
	values[2]=20;
	values[3]=30;
	values[4]=40;
	isza.push_back(0);
	isza.push_back(10);
	isza.push_back(20);
	ivalues.push_back(&values);
	ivalues.push_back(&values);
	ivalues.push_back(&values);

	outalt[0]=5;
	outalt[1]=15;
	outalt[2]=25;
	outalt[3]=35;
	outalt[4]=45;

	outsza[0]=0;
	outsza[1]=2;
	outsza[2]=4;
	outsza[3]=6;
	outsza[4]=8;

	ublas::vector<double> resu=MathFunction::IntLinPath(ialt,isza,ivalues,outalt,outsza);

	BOOST_CHECK(MathFunction::VectorCompare(resu,outalt));
	Log::mL<<"Double interpolation : "<<resu<<endl;
	resu[0]=-42;
	resu[1]=2;
	resu[2]=42;
	resu[3]=-2;
	resu[4]=4;
	NoNegative(resu);
	BOOST_CHECK(!(resu[0]<0));
	Log::mL<<"nonegative : "<<resu<<endl<<"Fin test interpolation"<<endl;

}


void TestSpline()
{
	ublas::vector<double> ialt(5),values(5),outalt(5);
	Log::mL<<"Initialization"<<endl;
	ialt[0]=0;
	ialt[1]=10;
	ialt[2]=20;
	ialt[3]=30;
	ialt[4]=40;
	values[0]=0;
	values[1]=10;
	values[2]=20;
	values[3]=30;
	values[4]=40;

	outalt[0]=5;
	outalt[1]=15;
	outalt[2]=25;
	outalt[3]=35;
	outalt[4]=45;
	Log::mL<<"splineInterp"<<endl;
	ublas::vector<double> resu=MathFunction::SplineInterp(ialt,values,outalt);
	Log::mL<<resu<<endl;
	BOOST_CHECK(MathFunction::VectorCompare(resu,outalt));

	values[0]=0;
	values[1]=100;
	values[2]=400;
	values[3]=900;
	values[4]=1600;
	resu=MathFunction::SplineInterp(ialt,values,outalt);
	Log::mL<<resu<<endl;
	Log::mL<<"End of spline test"<<endl;
}

void TestPhysTime()
{
	Log::mL<<"Computation of the julian time at the origin of time!"<<endl;
	double time = PhysTime::CalToJul(1984,1,1,17,15,1.42);
	Log::mL<<" Time: "<<time<<endl;
	BOOST_CHECK( nabs(time - 2445701.218766435 ) < 1E-10 );
}
/// On initialise les tests boost:
boost::unit_test_framework::test_suite* init_unit_test_suite(int,char*[])
{
	boost::unit_test_framework::test_suite* test;
	test=BOOST_TEST_SUITE("Mathstring Namespace test");
	test->add(BOOST_TEST_CASE(&TestLog));
	test->add(BOOST_TEST_CASE(&MathString_test_lit_string));
	test->add(BOOST_TEST_CASE(&MathString_test_lit_tout));
	test->add(BOOST_TEST_CASE(&MathString_test_lit_tableau));
	test->add(BOOST_TEST_CASE(&MathString_test_trim));
	test->add(BOOST_TEST_CASE(&MS_test_grid));
	test->add(BOOST_TEST_CASE(&Gaussangle));
	test->add(BOOST_TEST_CASE(&Interpol));
	test->add(BOOST_TEST_CASE(&AddVector));
	test->add(BOOST_TEST_CASE(&TestAddLines));
	test->add(BOOST_TEST_CASE(&TestStrToN));
	test->add(BOOST_TEST_CASE(&TestComparison));
	test->add(BOOST_TEST_CASE(&TestRedistribute));
	test->add(BOOST_TEST_CASE(&TestRandom)); 
	test->add(BOOST_TEST_CASE(&TestUblas));
	test->add(BOOST_TEST_CASE(&TestPathInterpolation));
	test->add(BOOST_TEST_CASE(&TestSpline));
	test->add(BOOST_TEST_CASE(&TestPhysTime));
	Log::mL<<"Fin tests math"<<endl;
	return test;
}

