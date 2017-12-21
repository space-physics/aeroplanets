/**
 * \file chapman.hpp
 * \brief Defines the chapman function and subfunction based on Huestis function http://www-mpl.sri.com/software/chapman/chapman.html
 * Copyright G Gronoff August 2010
 * Last modification : $Id$
 */

#ifndef CHAPMAN_HPP
#define CHAPMAN_HPP

#include "mathstring.hpp"
#include "mathflux.hpp"
#include "constphys.hpp"


/**
 * Class for the chapman function and subroutines, based on Huestis work
 * The class allow to pass data between subroutines through private members: it bypass the fortran common SH*.
 */
class ChapFunction
{
	private:
	/**
	 * Subroutine
	 * \param vX
	 * \param vChi
	 */
	double ChapDeq(double vX,double vChi);

	/**
	 * Subroutine
	 * \param vX
	 * \param vChi
	 */
	double ChapAsy(double vX,double vChi);


	/// Function of n
	double mFn[4];
	/// Functions outputs
	double mYI0,mYK0,mXI1,mXK1;

	double CxK1(double vX);
	double CxI1(double vX);
	double CyI0(double vX,double vChi);
	double CyK0(double vX,double vChi);

	void Cgd3(double vYo,double* vDn);
	
	double cXerfc(double vX);
	double Cgg06(double vX,double vX1,double* vGg);

	void Cgq85(double vX,double* vQbeta);


	public:
	/**
	 * Returns the chapman profile for the vXkm 'altitudes', for the vChi solar zenith angle
	 * \param vX: reduced altitude : (R+z)/H, R: planet radius, H scale height
	 * \param vChi: the solar zenith angle
	 */
	ublas::vector<double> Chapman(ublas::vector<double> vX, double vChi);


};




#endif

