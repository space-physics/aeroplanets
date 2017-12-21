/**
 * \file calcdens.cpp
 * \brief Implements the density computation functions
 * Copyrigth G Gronoff Feb 2010
 * Last Modification : $Id: calcdens.cpp 1282 2011-08-10 20:29:43Z gronoff $
 */
#include "calcdens.hpp"


# define Te mTempElecK
# define Ti mTempIonK
# define Tn mTempNeutreK
# define Ne (*mpElecDenscm_3)
# define REAC( i )  mChemList[i]->GetReactionRate(Te,Ti,Tn)
# define EFFIC( i , j , k )  mChemList[i]->GetEfficiency( j , k )
#define NAME rSubResultsNames
using namespace std;
typedef ublas::vector<double> Vec; 
typedef ublas::matrix<double> Mat;
typedef ublas::matrix_row< Mat > Row;
typedef ublas::matrix_column< Mat > Col;


/*
Vec operator*(Vec vA,Vec vB)
{
	assert(vA.size()==vB.size());
	Vec resu(vA.size());
	for(size_t i=0;i<vA.size();++i)
	{
		resu[i]=vA[i]*vB[i];
	}
	return resu;
}
*/


CalcDensN2A3S::CalcDensN2A3S():CalcDensPhEq(SpecieId("N2","A3S"))
{
}



ublas::vector<double> CalcDensN2A3S::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of N2(A3S)"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(7,mSize);
	NAME.resize(0);
	rWarnings="";

	Vec p1=(GetTotProd("N2","A3S")+GetTotProd("N2","B3P")+GetTotProd("N2","W3D")+GetTotProd("N2","C3P"));

	Vec l1=REAC(35) ; //radiative losses (Vegard Kaplan)
	Vec l2=GetDens("O","")*REAC(36);
	Vec l3=GetDens("CO","")*REAC(37);
	Vec l4=GetDens("CO2","")*REAC(38);
	Vec l5=GetDens("O2","")*REAC(39);

	Vec totloss=l1+l2+l3+l4+l5;


	resu=p1/totloss;

	std::deque< ublas::vector<double> > tr; // temporary results
	tr.push_back(l1*resu);
	NAME.push_back("VegardKaplan emission");

//	Log::mL<<"Put VK  in the matrix"<<tr.size()<<endl;
	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}

	return resu;
}


CalcDensO1S::CalcDensO1S():CalcDensPhEq(SpecieId("O","1S"))
{
}


ublas::vector<double> CalcDensO1S::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of O(1S)"<<endl;
	ublas::vector<double> resu(mSize);
//	rSubResults.resize(7,mSize);
//	rSubResults.clear();
	NAME.resize(0);
	rWarnings="";

	std::deque< ublas::vector<double> > tr; // temporary results

//	Log::mL<<"Losses"<<endl;
	Vec l1= REAC(13) + REAC(14); // radiative losses
	Vec l2=GetDens("O2","")*REAC(15);
	Vec l3=GetDens("O","")*REAC(16);
	Vec l4=GetDens("CO","")*REAC(17);
	Vec l5=GetDens("CO2","")*REAC(18);
	Vec l6=GetDens("N2","")*REAC(19);
	Vec l7=Ne*REAC(20);

//	Log::mL<<"Productions"<<endl;
	Vec p1=GetTotProd("O","1S");
//	Log::mL<<"Chem production"<<endl;
	Vec p2=Ne*GetDens("O2+","")*REAC(21);
	
//	Log::mL<<"Loss CO2"<<l5<<endl;
//	Log::mL<<"productions 1 :"<<p1<<" et p2 "<<p2<<endl;

	// Here, modification with the computation of N2(A3S)
/*	Vec q=2.8E-11*GetDens("O2","");
	for(Vec::iterator it=q.begin();it!=q.end();++it)
	{
		*it=0.36/(1+0.38/(*it));
	}
	Vec p3=(GetTotProd("N2","A3S")+GetTotProd("N2","B3P")+GetTotProd("N2","W3D")+GetTotProd("N2","C3P"));
	p3=p3*q;*/

	Vec p3=GetDens("N2","A3S")*GetDens("O","")*REAC(36)*EFFIC(36,"O","1S");
	Vec p4=GetDens("O2+","")*GetDens("N","")*REAC(22);

	Vec totprod=(p1+p2+p3+p4);
	Vec totloss=(l1+l2+l3+l4+l5+l6+l7);
	resu=totprod/totloss;


	tr.push_back( (p1/totloss));
	NAME.push_back("Physical productions");
	tr.push_back( (p2)/totloss);
	NAME.push_back("e + O2+ production");
	tr.push_back( (p3)/totloss);
	NAME.push_back("N2 deactivation production");
	tr.push_back( (p4)/totloss);
	NAME.push_back("O2+ + N (Frederick/Kopp)");
	tr.push_back(l1*resu);
	NAME.push_back("Total emission");
	tr.push_back(REAC(13)*resu);
	NAME.push_back("5577 A");
	tr.push_back(REAC(14)*resu);
	NAME.push_back("2972 A");

	// Losses
	tr.push_back(l2*resu); // O2 quenching
	NAME.push_back("O2 quenching");
	tr.push_back(l3*resu); // O quenching
	NAME.push_back("O quenching");
	tr.push_back(l4*resu); // CO quenching
	NAME.push_back("CO quenching");
	tr.push_back(l5*resu); // CO2 quenching
	NAME.push_back("CO2 quenching");
	tr.push_back(l6*resu); // N2 quenching
	NAME.push_back("Ne quenching");
	tr.push_back(l7*resu); // Ne quenching
	NAME.push_back("e- quenching");
	tr.push_back(totloss*resu); // total quenching
	NAME.push_back("Total quenching");

	rSubResults.resize(tr.size(),mSize);
	rSubResults.clear();
//	Log::mL<<"Put in the matrix"<<tr.size()<<endl;
	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}

//	Log::mL<<"End of the O1S density computation"<<endl;
	return resu;
}


CalcDensN2D::CalcDensN2D():CalcDensPhEq(SpecieId("N","2D"))
{
}


ublas::vector<double> CalcDensN2D::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of N(2D)"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(6,mSize);
	NAME.resize(0);
	rWarnings="";

	std::deque< ublas::vector<double> > tr; // temporary results

//	Log::mL<<"Losses"<<endl;
	Vec l1= REAC(31); // radiative losses
	Vec l2=GetDens("CO2","")*REAC(34);
	Vec l3=GetDens("O+","")*REAC(33);
	Vec l4=GetDens("O2+","")*REAC(32);
	Vec l5=Ne*REAC(30);
	Vec l6=GetDens("O","")*REAC(29);
	Vec l7=GetDens("O2","")*REAC(28);

	Vec p1=GetTotProd("N","2D");
	Vec p2=GetDens("N2+","")*GetDens("O","")*REAC(24);
	Vec p3=GetDens("N2+","")*Ne*REAC(25);
	Vec p4=GetDens("NO+","")*Ne*REAC(26);
	Vec p5=GetDens("N+","")*GetDens("O","")*REAC(27);

	Vec totprod=p1+p2+p3+p4+p5;
	Vec totloss=l1+l2+l3+l4+l5+l6+l7;
	resu=totprod/totloss;

	NAME.push_back("Total emission");
	tr.push_back(REAC(31)*resu);
	NAME.push_back("Physical productions");
	tr.push_back(p1);
	NAME.push_back("N2+ + O prod");
	tr.push_back(p2);
	NAME.push_back("N2+ + e prod");
	tr.push_back(p3);
	NAME.push_back("NO+ + e prod");
	tr.push_back(p4);
	NAME.push_back("N+ + O prod");
	tr.push_back(p5);


//	Log::mL<<"Put in the matrix"<<tr.size()<<endl;
	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}

//	Log::mL<<"End of the N2D density computation"<<endl;
	return resu;
}




CalcDensO1D::CalcDensO1D():CalcDensPhEq(SpecieId("O","1D"))
{
}


ublas::vector<double> CalcDensO1D::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of O(1D)"<<endl;
	ublas::vector<double> resu(mSize);
//	rSubResults.resize(11,mSize);
	NAME.resize(0);
	rWarnings="";

	std::deque< ublas::vector<double> > tr; // temporary results

//	Log::mL<<"Losses"<<endl;
	Vec l1= REAC(0) + REAC(1)+REAC(2); // radiative losses
	Vec l2=GetDens("O2","")*REAC(3);
	Vec l3=GetDens("O","")*REAC(4);
	Vec l4=GetDens("N2","")*REAC(5);
	Vec l5=Ne*REAC(6);
	Vec l6=GetDens("CO2","")*REAC(7);
	Vec l7=GetDens("CO","")*REAC(8);

	Vec p1=GetTotProd("O","1D");
	Vec p2=GetDens("O2+","")*Ne*REAC(9);
	Vec p3=Ne*GetDens("O","")*REAC(10);
	Vec p4=GetDens("O2","")*GetDens("N+","")*REAC(11);
	Vec p5=GetDens("CO+","")*Ne*REAC(12);
	Vec p6=GetDens("N","2D")*GetDens("O2","")*REAC(28)*EFFIC(28,"O","1D");

	Vec p7=(
			REAC(13)+
			REAC(15)*GetDens("O2","")*EFFIC(15,"O","1D")+
			REAC(18)*GetDens("CO2","")*EFFIC(18,"O","1D")+
			REAC(23)*Ne
	       );
	p7=GetDens("O","1S")*p7;
	Vec totprod=p1+p2+p3+p4+p5+p6+p7;
	Vec totloss=l1+l2+l3+l4+l5+l6+l7;
	resu=totprod/totloss;

	NAME.push_back("Total emission");
	tr.push_back(l1*resu);
	NAME.push_back("630");
	tr.push_back(REAC(0)*resu);
	NAME.push_back("636");
	tr.push_back(REAC(1)*resu);
	NAME.push_back("639");
	tr.push_back(REAC(2)*resu);

	NAME.push_back("Physical productions");
	tr.push_back(p1);
	NAME.push_back("O2+ + e");
	tr.push_back(p2);
	NAME.push_back("O + eth");
	tr.push_back(p3);
	NAME.push_back("O2 + N+");
	tr.push_back(p4);
	NAME.push_back("CO+ + e");
	tr.push_back(p5);
	NAME.push_back("O2 + N(2D) (?) ");
	tr.push_back(p6);
	NAME.push_back("O(1S) cascades");
	tr.push_back(p7);


	// Losses
	tr.push_back(l2*resu); // O2 quenching
	NAME.push_back("O2 quenching");
	tr.push_back(l3*resu); // O quenching
	NAME.push_back("O quenching");
	tr.push_back(l4*resu); // CO quenching
	NAME.push_back("N2 quenching");
	tr.push_back(l5*resu); // CO2 quenching
	NAME.push_back("e- quenching");
	tr.push_back(l6*resu); // N2 quenching
	NAME.push_back("CO2 quenching");
	tr.push_back(l7*resu); // Ne quenching
	NAME.push_back("CO quenching");
	tr.push_back(totloss*resu); // total quenching
	NAME.push_back("Total quenching");


	rSubResults.resize(tr.size(),mSize);
	rSubResults.clear();

//	Log::mL<<"Put in the matrix"<<tr.size()<<endl;
	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}

//	Log::mL<<"End of the O1D density computation"<<endl;
	return resu;




}




CalcDensCOa3Pi::CalcDensCOa3Pi():CalcDensPhEq(SpecieId("CO","a3Pi"))
{
}
ublas::vector<double> CalcDensCOa3Pi::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of CO(a3Pi)"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(4,mSize);
	NAME.resize(0);
	rWarnings="";

	std::deque< ublas::vector<double> > tr; // temporary results

//	Log::mL<<"Losses"<<endl;
	Vec l1= REAC(41); // radiative losses
	Vec l2=GetDens("O2","")*REAC(46);
	Vec l3=GetDens("N2","")*REAC(45);
	Vec l4=GetDens("CO2","")*REAC(42);
	Vec l5=GetDens("CO","")*REAC(43);
	Vec l6=GetDens("NO","")*REAC(44);

	Vec p1=GetTotProd("CO","a3Pi");
	Vec p2=GetDens("CO2+","")*Ne*REAC(40)*EFFIC(40,"CO","a3Pi");
	Vec totprod=p1+p2;
	Vec totloss=l1+l2+l3+l4+l5+l6;
	resu=totprod/totloss;


	NAME.push_back("Total emission");
	tr.push_back(l1*resu);
	NAME.push_back("Physical productions");
	tr.push_back(p1);
	NAME.push_back("CO2+ + e  productions");
	tr.push_back(p2);
	NAME.push_back("Total productions");
	tr.push_back(totprod);
//	Log::mL<<"Put in the matrix"<<tr.size()<<endl;

	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}
//	Log::mL<<"End of the COa3Pi density computation"<<endl;
	return resu;
}



CalcDensOplus2P::CalcDensOplus2P():CalcDensPhEq(SpecieId("O+","2P"))
{
}
ublas::vector<double> CalcDensOplus2P::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mL<<"We compute the density of O+(2P)"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(7,mSize);
	NAME.resize(0);
	rWarnings="";

	std::deque< ublas::vector<double> > tr; // temporary results
	Vec l1= REAC(50) + REAC(51) + REAC(52); // radiative losses
	Vec l2=GetDens("O","")*REAC(47);
	Vec l3=GetDens("N2","")*REAC(48);
	Vec l4=Ne*REAC(49);
	Vec l5=GetDens("CO2","")*REAC(53);
	Vec doublet = (REAC(50) + REAC(51));

	Vec p1=GetTotProd("O+","2P");
	Vec totprod=p1;
	Vec totloss=l1+l2+l3+l4+l5;
	resu=totprod/totloss;


	NAME.push_back("Total emission");
	tr.push_back(l1*resu);
	NAME.push_back("Doublet emission");
	tr.push_back(doublet  * resu);
	NAME.push_back("7320A");
	tr.push_back(REAC(50)*resu);
	NAME.push_back("7330A");
	tr.push_back(REAC(51)*resu);
	NAME.push_back("2470A");
	tr.push_back(REAC(52)*resu);

	NAME.push_back("Physical productions");
	tr.push_back(p1);
	NAME.push_back("Total productions");
	tr.push_back(totprod);
	for(unsigned i=0;i<tr.size();++i)
	{
		Row ro( rSubResults, i );
		ro=(tr[i]);
	}
	return resu;
}
//////////////////////////////////
CalcDensOpp::CalcDensOpp():CalcDensPhEq(SpecieId("O++","X"))
{
}
ublas::vector<double> CalcDensOpp::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of O++"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(7,mSize);
	NAME.resize(0);
	rWarnings="";

	Vec l1 = GetDens("N2","") * REAC(54);
	Vec l2 = GetDens("CO2","") * REAC(55);
	Vec l3 = GetDens("O","") * REAC(56);
	Vec l4 = GetDens("CO","") * REAC(57);
	Vec l5 = GetDens("He","") * REAC(58);
	Vec l6 = GetDens("O2","") * REAC(59);
	Vec l7 = Ne * REAC(60);
	Vec l8 = GetDens("H","") * REAC(61);


	Vec p1 = GetTotProd("O++","X");
	Vec totprod = p1;
	Vec totloss = l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8;
	resu = totprod / totloss;


	return resu;
}

//////////////////////////////////
CalcDensN2pp::CalcDensN2pp():CalcDensPhEq(SpecieId("N2++","X"))
{
}
ublas::vector<double> CalcDensN2pp::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of N2++"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(7,mSize);
	NAME.resize(0);
	rWarnings="";

	Vec l1 = GetDens("N2","") * REAC(62);
	Vec l2 = GetDens("CO2","") * REAC(64);
	Vec l3 = GetDens("O","") * REAC(63);
	Vec l4 = Ne * REAC(65);
	Vec l5 =  REAC(66);


	Vec p1 = GetTotProd("N2++","X");
	Vec totprod = p1;
	Vec totloss = l1 + l2 + l3 + l4 + l5;
	resu = totprod / totloss;


	return resu;
}




//////////////////////////////////
CalcDensCO2pp::CalcDensCO2pp():CalcDensPhEq(SpecieId("CO2++","X"))
{
}
ublas::vector<double> CalcDensCO2pp::GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)
{
	Log::mI<<"We compute the density of CO2++"<<endl;
	ublas::vector<double> resu(mSize);
	rSubResults.resize(7,mSize);
	NAME.resize(0);
	rWarnings="";
	Vec l1 = Ne * REAC(68);
	Vec l2 =  REAC(67);
	Vec l3 = GetDens("CO2","") * REAC(69);
	Vec l4 = GetDens("O","") * REAC(70);

	Vec p1 = GetTotProd("CO2++","X");
	Vec totprod = p1;
	Vec totloss = l1 + l2 + l3 + l4;
	resu = totprod / totloss;

	return resu;
}




