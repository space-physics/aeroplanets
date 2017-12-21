/**
 * \defgroup Math Math
 * This module aim is to do all the mathematical-related functions. As reading a file to get the values, or computing a grid.
 * \file mathstring.hpp
 * \brief defines mathematical functions, mainly associated with string retrieval.
 * E.G. if you  want to read an array from a file...
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathstring.hpp 1679 2013-02-15 23:31:18Z gronoff $
 *
 *
 */


#ifndef MATHSTRING_HPP
#define MATHSTRING_HPP


// Allows to define is debug or not
#include <config.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <boost/regex.hpp> 
#include <boost/assign/std/vector.hpp>
#include <cassert>
#include <fstream>
#include "mathrandom.hpp"



template< class T > ublas::vector<T> operator*(ublas::vector<T> vA,ublas::vector<T> vB)
{
	assert(vA.size()==vB.size());
/*	ublas::vector<T> resu(vA.size());
	for(size_t i=0;i<vA.size();++i)
	{
		resu[i]=vA[i]*vB[i];
	}
	return resu;*/
	return element_prod(vA,vB);
}

template< class T > ublas::vector<T> operator/(ublas::vector<T> vA,ublas::vector<T> vB)
{
	assert(vA.size()==vB.size());
/*	ublas::vector<T> resu(vA.size());
	for(size_t i=0;i<vA.size();++i)
	{
		resu[i]=vA[i]*vB[i];
	}
	return resu;*/
	return element_div(vA,vB);
}
template< class T > bool operator==(ublas::vector<T> vA,ublas::vector<T> vB)
{
	if(vA.size()!=vB.size())
		return false;
	bool resu=true;
	for(size_t i=0;i<vA.size();++i)
	{
		resu=resu&&(vA[i]==vB[i]);
	}
	return resu;
}



namespace mathessentials
{
	class MathError
	{
		private:
			std::string mMessage;
		public:
			MathError(std::string vMessage): mMessage(vMessage)
		{
		}
			void Affiche()
			{
				std::cout<<"================================"<<std::endl;
				std::cout<<"================================"<<std::endl;
				std::cout<<"Error In A Mathematical function"<<std::endl;
				std::cout<<"================================"<<std::endl;
				std::cout<<mMessage<<std::endl;
				std::cout<<"================================"<<std::endl;
				std::cout<<"================================"<<std::endl;
			}
	};




	/** Add two vector of same size 
	 * \ingroup Math
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 */
	template<class T> std::vector<T> AddVec(std::vector<T> vA,const std::vector<T>& vB)
	{
		unsigned size=vA.size();
		if(size!=vB.size())
		{
			MathError err("Error : vector size mismatch in AddVec");
			throw err;
		}
		for(unsigned i=0;i<size;++i)
		{
			vA[i]+=vB[i];
		}
		return vA;
	}

	/** Add two deque of same size 
	 * \ingroup Math
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 */
	template<class T> std::deque<T> AddDeQue(std::deque<T> vA,const std::deque<T>& vB)
	{
		unsigned size=vA.size();
		if(size!=vB.size())
		{
			MathError err("Error : vector size mismatch in AddVec");
			throw err;
		}
		for(unsigned i=0;i<size;++i)
		{
			vA[i]+=vB[i];
		}
		return vA;
	}

	/** Concatenate two deque 
	 * \ingroup Math
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return nothing: all is in the first vector
	 */
	template<class T> void ConcatenateDeQue(std::deque<T>& vA,const std::deque<T>& vB)
	{
		vA.insert(vA.end(),vB.begin(),vB.end());
	}


	/*
	 * Erase the whole content of a std::vector
	 * \param vVec : the vector
	 */
/*
	template<class T> void EraseVector(std::vector<T>& vVec)
	{
		vVec.erase(vVec.begin(),vVec.end());
	}
*/

	/** 
	 * allows to transform a std::vector<t> into a ublas::vector<t>
	 *
	 * \param va : the  vector
	 * \return the ublas::vector<t>
	 */
	template<class T> ublas::vector<T> StdToUblas(std::vector<T>& va)
	{
		ublas::vector<T> tmp(va.size());
		std::copy(va.begin(),va.end(),tmp.begin());
		return tmp;
	}

	/** 
	 * allows to transform a std::DeQue<t> into a ublas::vector<t>
	 *
	 * \param va : the  vector
	 * \return the ublas::vector<t>
	 */
	template<class T> ublas::vector<T> StdToUblas(std::deque<T>& va)
	{
		ublas::vector<T> tmp(va.size());
		std::copy(va.begin(),va.end(),tmp.begin());
		return tmp;
	}



	template<class T> ublas::matrix<T> PointerToMatrix(T vpP[],unsigned vNbLine,unsigned vNbCol)
	{
		ublas::matrix<T> resu(vNbLine,vNbCol);
		for(unsigned i=0;i<vNbLine;++i)
		{
			for(unsigned j=0;j<vNbCol;++j)
			{
				resu(i,j)=vpP[i*vNbCol+j];
			}
		}
		return resu;
	}

	/** 
	 * Allows to transform a std::vector< std::vector<T> > into a ublas::matrix<T>
	 *
	 * \param vA : the double vector
	 * \return the ublas::vector<T>
	 */
	template<class T> ublas::matrix<T> StdToUblas(std::vector< std::vector<T> >& vA)
	{
		unsigned lig=vA.size();
		if(lig==0)
		{
			return ublas::matrix<T>(0,0);
		}

		unsigned col=vA[0].size();

		ublas::matrix<T> tmp(lig,col);
		for(unsigned i=0;i<lig;++i)
		{
			ublas::matrix_row< ublas::matrix<T> > ro(tmp,i);
			std::copy(vA[i].begin(),vA[i].end(),ro.begin());
			//	ro=StdToUblas(vA[i]);
		}
		return tmp;



	}
	/*
	 * The inverse of the last function
	 * \param vA : the ublas matrix to pass to the function
	 * \return : the std vector 
	 *//*
	template<class T> std::vector< std::vector<T> > UblasToStd( ublas::matrix<T>& vA)
	{
		std::vector< std::vector<T> > resu;
		for(unsigned i=0;i<vA.size1();++i)
		{
			ublas::matrix_row< ublas::matrix<T> > ro(vA,i);
			std::vector<T> tmp(ro.size());
			for(unsigned j=0;j<ro.size();++j)
			{
				tmp[j]=ro(j);
			}
			resu.push_back(tmp);

			//resu.push_back(UblasToStd(ro));// UblasToStd(ro));
		}
		return resu;
	}*/	

	/** 
	 * allows to transform a ublas::vector<t> into a std::vector<t>
	 *
	 * \param va : the ublas::vector
	 * \return the std::vector<t>
	 */

	template<class T> std::vector<T> UblasToStd(const ublas::vector<T>& va)
	{
		std::vector<T> tmp(va.size());
		std::copy(va.begin(),va.end(),tmp.begin());
		//for(unsigned i=0;i<va.size();++i)
		//	tmp[i]=va(i);
		return tmp;
	}





	/**
	 * Takes two ublas::vectors, if one is full, and the other null, returns the full one
	 * if the two are filled, return the sum.
	 * Throw an error if they are of different non-null sizes.
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 *
	 */
	template<class T> ublas::vector<T>  AddOrFullVec(const ublas::vector<T>& vA,const ublas::vector<T>& vB)
	{
		if(0==vA.size())
		{
			return vB;
		}	
		if(0==vB.size())
		{
			return vA;
		}
		//return AddVec(vA,vB);
		return vA+vB;
	}


	/** Add two ublas::matrix of same size 
	 * \ingroup Math
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 *
	 *
	 */
	template<class T> ublas::matrix<T>  AddVec2(ublas::matrix<T> vA,const ublas::matrix<T> & vB)
	{
		unsigned size=vA.size1();
		unsigned size2=vA.size2();
		if(size!=vB.size1() || size2 != vB.size2() )
		{
			MathError err("Error : vector size mismatch in AddVec2");
			throw err;
		}
		return vA+vB;
	}


	/**
	 * Takes two vectors of vectors, if one is full, and the other null, returns the full one
	 * if the two are filled, return the sum.
	 * Throw an error if they are of different non-null sizes.
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 *
	 */
	template<class T>  ublas::matrix<T> AddOrFullVec2(const ublas::matrix<T> & vA,const ublas::matrix<T> & vB)
	{
		if(0==vA.size1())
		{
			return vB;
		}	
		if(0==vB.size1())
		{
			return vA;
		}
		return AddVec2(vA,vB);
	}

	/** Add two vector of matrix of same size 
	 * \ingroup Math
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 *
	 *
	 */

	template<class T> ublas::vector< ublas::matrix<T> > AddVec3(ublas::vector< ublas::matrix<T> > vA, const ublas::vector< ublas::matrix<T> >& vB)
	{
		unsigned size=vA.size();
		if(size!=vB.size())
		{
			MathError err("Error : vector size mismatch in AddVec3");
			throw err;
		}
		for(unsigned i=0;i<size;++i)
		{
			vA[i]=AddVec2(vA[i],vB[i]);
		}
		return vA;
	}


	/**
	 * Takes two vectors of matrix, if one is full, and the other null, returns the full one
	 * if the two are filled, return the sum.
	 * Throw an error if they are of different non-null sizes.
	 * \param vA : the first vector
	 * \param vB : the second vector
	 * \return the added vector
	 *
	 */
	template<class T> ublas::vector< ublas::matrix<T> > AddOrFullVec3(const ublas::vector< ublas::matrix<T> >& vA,const ublas::vector< ublas::matrix<T>  >& vB)
	{
		if(0==vA.size())
		{
			return vB;
		}	
		if(0==vB.size())
		{
			return vA;
		}
		return AddVec3(vA,vB);
	}




	/* Multiply a vector by a constant
	 * \ingroup Math
	 * \param vA : the vector
	 * \param vK : the constant
	 *
	 * \return the multiplied vector
	template<class T> std::vector<T> MultVec(std::vector<T> vA, T vK)
	{
		for(unsigned i=0;i<vA.size();++i)
		{
			vA[i]*=vK;
		}
		return vA;
	}
	 */



	/**
	 * Returns the position of the closest value in vBotteDeFoin
	 * \ingroup Math
	 * \param vAiguille : the constant to search
	 * \param vBotteDeFoin : the vector where the constant is hidden
	 * \return the position of the closest (always defined!)
	 */

	template<class T> unsigned CloseInVector(T vAiguille,const ublas::vector<T>& vBotteDeFoin)
	{
		if(vBotteDeFoin.size()<2)
			return 0;
		unsigned resu=0;
		double ecart=pow((vAiguille-vBotteDeFoin[0]),2);

		for(unsigned i=1;i<vBotteDeFoin.size();++i)
		{
			double necart=pow((vAiguille-vBotteDeFoin[i]),2);
			if(necart<ecart)
			{
				ecart=necart;
				resu=i;
			}
		}
		return resu;
	}


	/* Returns the numerical position in a vector, send -1 else
	 * \ingroup Math
	 * \param vAiguille : the constant to search
	 * \param vBotteDeFoin : the vector where the constant is hidden
	 * \return the position or -1
	 *

	template<class T> int PosInVector(T vAiguille,const std::vector<T>&  vBotteDeFoin)
	{
		for(unsigned i=0;i<vBotteDeFoin.size();++i)
		{
			if(vAiguille==vBotteDeFoin[i])
			{
				return (int)i;
			}
		}
		return -1;
	}

	 */



	/**
	 * Check if value in in the vector a. If so, return true, and position contains the value of the position in the vector else, return false, position=a.size()+1
	 * \ingroup Math
	 * \param vA : the vector
	 * \param vValue : the constant to search
	 * \param rPosition : the reference to the position
	 * \return true if it is found.
	 */
	template<class T> bool PosInVector(ublas::vector<T> vA, T vValue,unsigned& rPosition)
	{
		unsigned size=vA.size();
		rPosition = size+1;
		for(unsigned i=0;i<size;++i)
		{
			if(vA[i]==vValue)
			{
				rPosition=i;
				return true;
			}
		}
		return false;
	}



	/**
	 * Check if value in in the vector a. If so, return true, and position contains the value of the position in the vector else, return false, position=a.size()+1
	 * \ingroup Math
	 * \param vA : the vector
	 * \param vValue : the constant to search
	 * \param rPosition : the reference to the position
	 * \return true if it is found.
	 */
	template<class T> bool PosInVector(std::deque<T> vA, T vValue,unsigned& rPosition)
	{
		unsigned size=vA.size();
		rPosition = size+1;
		for(unsigned i=0;i<size;++i)
		{
			if(vA[i]==vValue)
			{
				rPosition=i;
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns the position of the maximum of a ublas::vector
	 * \ingroup Math
	 * \param vA: the vector
	 * \return size_t : the position of the maximum value
	 */
	template<class T> size_t PosOfMax(ublas::vector<T> vA)
	{
		size_t maxpos=0;
		size_t size=vA.size();
		if(size<=1)
		{
			return maxpos;
		}
		T maxval=vA[0];
		for(size_t i=1;i<size;++i)
		{
			if(vA[i]>maxval)
			{	
				maxval=vA[i];
				maxpos=i;
			}
		}
		return maxpos;
	}

	/** Fonction abs, return the absolute
	 * \ingroup Math
	 * \param a : the value we want the absolute
	 * \return the abs(a)
	 */
	template<class T> T nabs(T a)
	{// abs and max are the rare libraries where this form ( use of ?) is accepted!!!
		return ((a>0) ? (a): (-a));
		/*	if(a>0)
			return a;
			else
			return -a;
			*/
	}

	/** Fonction min
	 * \ingroup Math
	 * \param a : the first value
	 * \param b : the second value
	 * \return the min
	 *
	 */
	template<class T> T nmin(T a,T b)
	{
		return ( ( a > b ) ?  (b) : (a)  );
		/*
		   if(a>b)
		   return b;
		   else
		   return a;
		   */
	}

	/** Fonction max
	 * \ingroup Math
	 * \param a : the first value
	 * \param b : the second value
	 * \return the max
	 *
	 */ 
	template<class T> T nmax(T a,T b)
	{
		return ( ( a < b ) ?  (b) : (a)  );
		/*	if(a<b)
			return b;
			else
			return a;
			*/
	}



	/** Returns the modulo, always positive contrarly to fmod
	 * Here, the idea is to have a modulo that takes into account negative numbers!
	 * \ingroup math
	 * \param a : the number we want to put into the modulo part
	 * \param b : the modulo
	 * \return the modulo
	 *
	 */
	double modulo(double a, double b);

	/**
	 * Transforms a value into a string
	 * \ingroup Math
	 * \param vValue : the value to be transformed
	 * \return  : the string
	 *
	 *
	 */
	template<class T> std::string ntostr(T vValue)
	{
		std::ostringstream ostr;
		ostr<<vValue;
		return ostr.str();
	}

	/**
	 * Transforms a string into a value
	 * \ingroup Math
	 * \param value : the string
	 * \param machin : the value to return
	 */
	template<class T> void strton(std::string value,T& machin)
	{
		std::istringstream istr;
		istr.str(value);
		//	T machin;
		istr>>machin;
		//	return machin;
	}



	/** Allows to suppress the beginning and ending whitespaces
	 * \ingroup Math
	 * \param vText : the text to be transformed
	 * \param vChars : the leading and trailing characters to suppress,
	 * 	it is a regular expression. The predefined value suppress
	 * 	the spaces. \warning Do not modify it unless you known what you are doing
	 *
	 */

	std::string trim(std::string vText,std::string vChars="\\s*");


	/**
	 * Allows to replace a character, or a substring, into another substring in the whole string
	 * \param vStr : the string to modify (note this is not a reference....)
	 * \param vAig : the substring to replace
	 * \param vRep : the replacement substring
	 *
	 */
	std::string StrReplace(std::string vStr,std::string vAig,std::string vRep);

	/**
	 * Puts the negative values at zero in the matrix passed in argument
	 * \param vMat  : the matrix which needs to have values >=0
	 */
	void NoNegative(ublas::vector<double>& vMat);

	/**
	 * Puts the values below vMin at vMin
	 * \param vMat  : the matrix which needs to have values >=vMin
	 * \param vMin : the minimum value
	 */
	void MinValue(ublas::vector<double>& vMat,double vMin);
};
using namespace mathessentials;

/**
 * \ingroup Math
 * \~english Namespace to work with arrays of numerical values embedded in strings
 *
 * \~francais Namespace contenant les opérations nécessaires  pour travailler avec des tableaux de valeurs (doubles)   enfermées dans des strings
  */
namespace MathString
{




	/**
	 * \~english Reads a value inside a string.
	 * It is a template function  (you can read int, double...)
	 * therefore, the definition has to be in the header.
	 * \param vTruc : the string where the value is
	 * \param rValue : the searched value.
	 *
	 *
	 * \~francais  permet de lire une valeur dans la chaine de caractère.  
	 * Fonction template : le corps ne peut être mis dans le fichier cpp 
	 * \param vTruc : la chaine dans laquelle il y a le fichier 
	 * \param rValue : la valeur à trouver (paramètre) 
	 * \return bool : true si trouvé, false sinon
	 */
	template<class T> bool LitString(std::string vTruc,T &rValue)
	{
		std::istringstream istr;
		istr.str(vTruc);
		if(istr>>rValue)
		{
			return true;
		}
		return false;
	}

	/** 
	 * \~english Read a string to give a vector
	 * \param vTruc: the string
	 * \param rValue : the vector of values
	 * \return bool : true if found
	 *
	 *
	 * \~francais Permet de remplir un vecteur avec ce qu'il trouve dans une chaine de caractere
	 * \param  vTruc : la chaine de caractere
	 * \param  rValue : le vecteur des valeurs a trouver
	 * \return bool : true si trouve, false sinon
	 */
	template<class T> bool LitToutString(std::string vTruc,ublas::vector<T> & rValue)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::vector<T> tmp;
		T tmpvalue;
		while(istr>>tmpvalue)
		{
			tmp.push_back(tmpvalue);
		}
		rValue.resize(tmp.size());
		std::copy(tmp.begin(),tmp.end(),rValue.begin());
		if(rValue.size()>0)
		{
			return true;
		}
		return false;
	}

	/** 
	 * \~english Read a string to give a vector
	 * \param vTruc: the string
	 * \param rValue : the vector of values
	 * \return bool : true if found
	 *
	 *
	 * \~francais Permet de remplir un vecteur avec ce qu'il trouve dans une chaine de caractere
	 * \param  vTruc : la chaine de caractere
	 * \param  rValue : le vecteur des valeurs a trouver
	 * \return bool : true si trouve, false sinon
	 */
	template<class T> bool LitToutString(std::string vTruc,std::vector<T> & rValue)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::vector<T> tmp;
		T tmpvalue;
		while(istr>>tmpvalue)
		{
			tmp.push_back(tmpvalue);
		}
		rValue.resize(tmp.size());
		std::copy(tmp.begin(),tmp.end(),rValue.begin());
		if(rValue.size()>0)
		{
			return true;
		}
		return false;
	}

	/**
	 * \~english Fills a vector by the column of an array
	 * \param vTruc : the array (string)
	 * \param rValue : the vector to fill
	 * \param vPosition : which column we want to put in the vector (0=start)
	 *
	 *
	 * \~francais Permet de remplir un vecteur par une colonne de tableau
	 * \param  vTruc : le tableau sous forme de string
	 * \param  rValue : le vecteur à remplir
	 * \param  vPosition : la position de la colonne du tableau (0 = début)
	 */

	template<class T> bool LitTableau(std::string vTruc,ublas::vector<T> & rValue,unsigned int vPosition)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::string tmpstring;
		std::vector<T> resultat;
		while(getline(istr,tmpstring))
		{
			ublas::vector<T> tmpvec;
			if(LitToutString(tmpstring,tmpvec))
			{
				if(tmpvec.size()<(vPosition+1))
				{
					return false;
				}
				resultat.push_back(tmpvec[vPosition]);
			}
		}
		rValue.resize(resultat.size());
		std::copy(resultat.begin(),resultat.end(),rValue.begin());
	
		return true;
	}
	/**
	 * \~english Allows to fill a vector from an array (event in 2D)
	 * the array is read like
	 * 1 2 3
	 * 4 5 6
	 * 7 8 9
	 * \param vTruc : the array (string)
	 * \param rValue : the vector to fill
	 *
	 *
	 * \~francais Permet de remplir un vecteur  a partir d'un tableau.
	 * Celui-ci peut avoir plusieurs valeurs sur la même ligne
	 * il est alors lu de la suite:
	 * 1	2	3
	 * 4	5	6
	 * 7	8	9
	  \param  vTruc : le tableau sous forme de string
	  \param  rValue : le vecteur à remplir
	  */

	template<class T> bool LitTableau(std::string vTruc,ublas::vector<T> & rValue)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::string tmpstring;
		std::vector<T> resultat;
		while(getline(istr,tmpstring))
		{
			ublas::vector<T> tmpvec;
			if(LitToutString(tmpstring,tmpvec))
			{
				resultat.insert(resultat.end(),tmpvec.begin(),tmpvec.end());
			}
		}
		rValue.resize(resultat.size());
		std::copy(resultat.begin(),resultat.end(),rValue.begin());
		return true;
	}







	/**
	 * \~english fills a vector of vector by a 2D array (string)
	 * \param vTruc : the array (string)
	 * \param rValue : the vector of vector
	 *
	 *
	 * \~francais Permet de remplir vecteur de vecteur
	 * \param  vTruc : le tableau sous forme de string
	 * \param  rValue : le vecteur à remplir
	 */

	template<class T> bool Lit2DString(std::string vTruc,ublas::matrix<T> & rValue)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::string tmpstring;
		std::vector< ublas::vector<T> > resultat;
		size_t size1=0;
		size_t size2=0;

		unsigned lect=0;
		while(getline(istr,tmpstring))
		{
			ublas::vector<T> tmpvec;
			if(LitToutString(tmpstring,tmpvec))
			{
				size_t nsize=tmpvec.size();
				if(lect==0)
				{
					size2=nsize;
				}else
				{
					if(size2!=nsize)
					{
						MathError err("size mismatch in your 2d array");
						throw err;
					}
				}

				resultat.push_back(tmpvec);
			}
		}
		size1=resultat.size();
		rValue.resize(size1,size2);
		for(size_t i=0;i<size1;++i)
		{
			ublas::matrix_row< ublas::matrix<T> > ro(rValue,i);
			ro=resultat[i];

		}

		return true;

	}


	/**
	 * \~english Fills two vector by two columns of the array
	 * \param vTruc : the array (string)
	 * \param rValue1 : the first vector
	 * \param vPosition1 : its position in truc
	 * \param rValue2 : the second vector
	 * \param vPosition2 : its position in truc
	 *
	 * \~francais Permet de remplir deux vecteur par deux colonne de tableau
	 * \param vTruc : le tableau sous forme de string
	 * \param rValue1 : le vecteur à remplir
	 * \param vPosition1 : la position de la colonne du tableau (0 = début)
	 * \param rValue2 : le vecteur à remplir
	 * \param vPosition2 : la position de la colonne du tableau (0 = début)
	 */
	template<class T> bool LitTableauDouble(std::string vTruc, ublas::vector<T> & rValue1,unsigned int vPosition1,ublas::vector<T> & rValue2,unsigned int vPosition2)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::string tmpstring;
		std::vector<T> resultat1;
		std::vector<T> resultat2;
		while(getline(istr,tmpstring))
		{
			ublas::vector<T> tmpvec;
			if(LitToutString(tmpstring,tmpvec))
			{
				if(tmpvec.size()<(vPosition1+1))
				{
					return false;
				}
				resultat1.push_back(tmpvec[vPosition1]);

				if(tmpvec.size()<(vPosition2+1))
				{
					return false;
				}
				resultat2.push_back(tmpvec[vPosition2]);
			}
		}
		rValue1.resize(resultat1.size());
		rValue2.resize(resultat2.size());
		std::copy(resultat1.begin(),resultat1.end(),rValue1.begin());
		std::copy(resultat2.begin(),resultat2.end(),rValue2.begin());

		return true;

	}

	/**
	 * \~english Allows to fill 3 vectors by 3 columns of an array
	 *
	 * \param vTruc : the array (string)
	 * \param rValue1 : the first vector
	 * \param vPosition1 : its position in truc
	 * \param rValue2 : the second vector
	 * \param vPosition2 : its position in truc
	 * \param rValue3 : the first vector
	 * \param vPosition3 : its position in truc
	 *
	 *
	 * \~francais Permet de remplir trois vecteurs par trois colonnes de tableau
	 * \param vTruc : le tableau sous forme de string
	 * \param rValue1 : le vecteur à remplir
	 * \param vPosition1 : la position de la colonne du tableau (0 = début)
	 * \param rValue2 : le vecteur à remplir
	 * \param vPosition2 : la position de la colonne du tableau (0 = début)
	 * \param rValue3 : le vecteur à remplir
	 * \param vPosition3 : la position de la colonne du tableau (0 = début)
	 */
	template<class T> bool LitTableauTriple(std::string vTruc,ublas::vector<T> & rValue1,unsigned int vPosition1,ublas::vector<T> & rValue2,unsigned int vPosition2,ublas::vector<T> & rValue3,unsigned int vPosition3)
	{
		std::istringstream istr;
		istr.str(vTruc);
		std::string tmpstring;
		std::vector<T> resultat1;
		std::vector<T> resultat2;
		std::vector<T> resultat3;
		while(getline(istr,tmpstring))
		{
			std::vector<T> tmpvec;
			if(LitToutString(tmpstring,tmpvec))
			{
				if(tmpvec.size()<(vPosition1+1))
				{
					return false;
				}
				resultat1.push_back(tmpvec[vPosition1]);

				if(tmpvec.size()<(vPosition2+1))
				{
					return false;
				}
				resultat2.push_back(tmpvec[vPosition2]);

				if(tmpvec.size()<(vPosition3+1))
				{
					return false;
				}
				resultat3.push_back(tmpvec[vPosition3]);
			}
		}

		rValue1.resize(resultat1.size());
		rValue2.resize(resultat2.size());
		rValue3.resize(resultat3.size());
		std::copy(resultat1.begin(),resultat1.end(),rValue1.begin());
		std::copy(resultat2.begin(),resultat2.end(),rValue2.begin());
		std::copy(resultat3.begin(),resultat3.end(),rValue3.begin());
		return true;

	}








	/**
	 * Print an 1d vector
	 * \param  vec : the vector to print
	 */

	template<class T> void print1d(const std::vector<T>& vec)
	{
		std::cout<<"================================================="<<std::endl;
		std::cout.precision(9);
		//std::cout.setf(std::ios::scientific);
		for(unsigned i=0;i<vec.size();i++)
		{
			std::cout<<vec[i]<<std::endl;
		}
		std::cout<<"================================================="<<std::endl;
	}
	/**
	 * Print a 2d vector
	 * \param  vec : the vector to print
	 */

	template<class T> void print2d(const std::vector< std::vector<T> >& vec)
	{
		std::cout<<"================================================="<<std::endl;
		if(vec.size()==0)
			return;
		
		for(unsigned l=0;l<vec[0].size();++l)
		{
			for(unsigned c=0;c<vec.size();++c)
			{
				std::cout<<vec[c][l]<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<"================================================="<<std::endl;
	}






};


#endif
