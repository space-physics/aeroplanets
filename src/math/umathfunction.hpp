/**
 * \file umathfunction.hpp
 * \brief Defines basic functions, like integration...
 * contrarly to mathfunction, it works with ublas vectors
 * Copyright G Gronoff Nov 2009
 * Last modification : $Id: umathfunction.hpp 1440 2012-03-05 22:24:56Z gronoff $
 */

#ifndef UMATHFUNCTION_HPP
#define UMATHFUNCTION_HPP

#include "chapman.hpp"



/**
 *
 * Misc mathematical functions for computing mean, gaussian angle, interpolating...
 * \ingroup Math
 *
 */
namespace MathFunction
{

	/**
	 * \~english  Moving average
	 *  Allows to make a mean with n parameter to take into account
	 *  \param vData : the vector to mean
	 *  \param nMoy : the number of points to mean
	 *  \param rResult : the smoothed vector.
	 *
	 *
	 * \~francais Moyenne glissante: permet de faire une moyenne glissante avec n paramètres
	 * \param  vData : le vecteur à moyenner
	 * \param  nMoy : le nombre de points
	 * \param  rResult : le vecteur contenant les résultats
	 */

	void MoyGliss(const ublas::vector<double>& vData, int nMoy, ublas::vector<double>& rResult);



	/**
	 * \~english 
	 * Computes an exponential distribution of nb_pts points between min and max
	 * \param vMin : the min value
	 * \param vMax : the max value
	 * \param vNbPts : the number of points
	 * \return the distribution
	 *
	 * \warning Be carful, the resulting function is decreasing
	 *
	 * \~francais Permet de calculer une distribution exponentielle de nb_pts entre min et max
	 *	\param  vMin: la hauteur minimale
	 *	\param  vMax: la hauteur maximale
	 *	\param vNbPts : le nombre de points
	 *	\return vector<double> : la distribution
	 */
	ublas::vector<double> ExponentialDistrib(double vMin, double vMax,int vNbPts);



	/**
	 * Computes and store gaussian angles: allows to make efficient integration over a few number of angles.
	 */
	class GaussianAngle
	{
		private:
			/// The half number of angles: gaussian angle should be even
			int mNbmid;
			/// Init the different vectors
			void Gaussa();
			/// Temporary gaussian abscissa
			ublas::vector<double> mGmu;
			/// Temporary gaussian weight
			ublas::vector<double> mGwt;
		public:
			/// The number of angles
			int mNbAngles;

			/**
			 * Initializes the vectors, computes the angle abscissa, weight
			 * \param vNbAngle : number of angles
			 */
			GaussianAngle(int vNbAngle);

			/// Gaussian abscissa
			ublas::vector<double> mXmu;
			/// Gaussian weight
			ublas::vector<double> mWeight;
			/// Gaussian angle in degree
			ublas::vector<double> mAngzb;
			/// Difference in Gaussian angles: matrix, in radians; Corresponds to dtheta in the fortran proton code, but in the matrix format
			ublas::matrix<double> mDThetaR;
			/// For the tests with Gwt
			ublas::vector<double> returnGwt(){return mGwt;}
	};




	/*------------------------------------------------------
	 *
	 *
	 *  Cubic splines
	 *
	 *
	 * ---------------------------------------------------*/

	/**
	 * Derives the function X-Y, with the first and last derivative forced with D0 and Dn
	 * if D0 or Dn > 1E33, we take the natural way
	 * \param vX : the abscissa of the function
	 * \param vY : the value of the function
	 * \param vD0 : forcing of the derivative at the first point
	 * \param vDn : forcing of the derivative at the last point
	 */
	ublas::vector<double> SplineDeriv(ublas::vector<double> vX, ublas::vector<double> vY,double vD0, double vDn);
	/**
	 * Cubic spline interpolation at the point vPt
	 * \param vX : the abscissa of the function
	 * \param vY : the value of the function
	 * \param vDy : the derivative (typically computed with SplineDeriv)
	 * \param vPt : the point where the spline is computed
	 */
	double SplineInterpP(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vDy, double vPt);
	/**
	 * Cubic spline interpolation for the points of vNewx; with the derivative forced at the border, and computed with splinederiv.
	 * \param vX : the abscissa of the function
	 * \param vY : the value of the function
	 * \param vD0 : forcing of the derivative at the first point
	 * \param vDn : forcing of the derivative at the last point
	 * \param vNewx : the vector where the spline interpolation is computed
	 */
	ublas::vector<double> SplineInterpForce(ublas::vector<double> vX, ublas::vector<double> vY,double vD0, double vDn,ublas::vector<double> vNewx);
	/**
	 * Cubic spline interpolation for the points of vNewx; with the derivative computed at the border as the derivative for the extreme pts, and computed with splinederiv.
	 * \param vX : the abscissa of the function
	 * \param vY : the value of the function
	 * \param vNewx : the vector where the spline interpolation is computed
	 */
	ublas::vector<double> SplineInterp(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vNewx);

	/**
	 * Cubic spline interpolation for the points of vNewx; with the derivative computed at the border as the derivative for the extreme pts, and computed with splinederiv.
	 * The specificity of that function is that we take the exponent of the output (better suited for planetary atmospheres)
	 * \param vX : the abscissa of the function
	 * \param vY : the log of the value of the function
	 * \param vNewx : the vector where the spline interpolation is computed
	 */
	ublas::vector<double> SplineInterpExp(ublas::vector<double> vX, ublas::vector<double> vY,ublas::vector<double> vNewx);






	/*------------------------------------------------------
	 *
	 *
	 * Interpolations
	 *
	 *
	 * ---------------------------------------------------*/



	/** 
	 * Linear interpolation
	 * interpolate goldy on gnewx
	 * can be used for extrapolation
	 * see also intlog
	 * \param vOldx : x for the initial value
	 * \param vOldy : y for the initial value
	 * \param vNewx : x for the new value
	 * \return y for the new value
	 */
	ublas::vector<double> IntLin(ublas::vector<double>vOldx,ublas::vector<double>vOldy,ublas::vector<double>vNewx);
	
	/** 
	 * Logarithmic interpolation on log grid : the x grids are put at log too
	 * interpolate goldy on gnewx
	 * can be used for extrapolation
	 * see also intlin
	 *
	 * \warning null and negative values-> the interpolation
	 * must give 0 at the end. The problem is to optimize that function.
	 * (it is widely used).
	 * So, only an assert is done, the user of the function should
	 * test it before!
	 *
	 * \param vOldx : x for the initial value
	 * \param vOldy : y for the initial value
	 * \param vNewx : x for the new value
	 * \return y for the new value
	 */
	ublas::vector<double> IntLog(const ublas::vector<double>& vOldx,const ublas::vector<double>& vOldy,const ublas::vector<double>& vNewx);
	
	
	/** 
	 * Logarithmic interpolation
	 * interpolate goldy on gnewx
	 * can be used for extrapolation
	 * see also intlin
	 *
	 * \warning null and negative values-> the interpolation
	 * must give 0 at the end. The problem is to optimize that function.
	 * (it is widely used).
	 * So, only an assert is done, the user of the function should
	 * test it before!
	 *
	 * \param vOldx : x for the initial value
	 * \param vOldy : y for the initial value
	 * \param vNewx : x for the new value
	 * \return y for the new value
	 */
	ublas::vector<double> IntLogLog(const ublas::vector<double>& vOldx,const ublas::vector<double>& vOldy,const ublas::vector<double>& vNewx);




	/**
	 * Path interpolation
	 * Interpolate a series of values, with corresponding factor A and B (altitude and SZA) to a unique distribution of A and B (altitude and SZA).
	 * Typically, for a set of computed values (ex: densities) for a set of altitudes and SZA, computes these values on a path defined by an unique altitude and SZA for each point.
	 * \warning The output value can be positive or negative, if you want only positive value (eg densities), you can use the NoNegative function on the result.
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vInputVals : the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \param vNoNegative : true by default : all the negative values are set to 0.
	 * \return the values for each points
	 */
	ublas::vector<double> IntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector<double>* >& vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative=true);

	/**
	 * Path interpolation
	 * Same as IntLinPath, except that the interpolation is Log
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vInputVals : the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \return the values for each points
	 */
	ublas::vector<double> IntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector<double>* >& vInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA);

	/**
	 * Path interpolation for a set of input vals: we have a set of vInputVals (see IntLinPath). In this function, we can make more interpolations at the same time.
	 * Interpolate a series of values, with corresponding factor A and B (altitude and SZA) to a unique distribution of A and B (altitude and SZA).
	 * Typically, for a set of computed values (ex: densities) for a set of altitudes and SZA, computes these values on a path defined by an unique altitude and SZA for each point.
	 * \warning The output value can be positive or negative, if you want only positive value (eg densities), you can use the NoNegative function on the result.
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vSetInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \param vNoNegative : true by default : all the negative values are set to 0.
	 * \return the values for each points
	 */
	std::deque< ublas::vector<double> >  SetIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque<  std::deque< ublas::vector<double> > >& vSetInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative=true);

	/**
	 * Log implementation of SetIntLinPath
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vSetInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \return the values for each points
	 */
	std::deque< ublas::vector<double> >  SetIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque<  std::deque< ublas::vector<double> > >& vSetInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA);
	/**
	 * Path interpolation for a set of input matrices
	 * Interpolate a series of values, with corresponding factor A and B (altitude and SZA) to a unique distribution of A and B (altitude and SZA).
	 * Typically, for a set of computed values (ex: densities) for a set of altitudes and SZA, computes these values on a path defined by an unique altitude and SZA for each point.
	 * \warning The output value can be positive or negative, if you want only positive value (eg densities), you can use the NoNegative function on the result.
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vMatInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \param vNoNegative : true by default : all the negative values are set to 0.
	 * \return the values for each points
	 */
	ublas::matrix<double> MatIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::matrix<double>* > & vMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative=true);
		/**
	 * Path interpolation for a set of input matrices. Log implementation 
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vMatInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \return the values for each points
	 */
	ublas::matrix<double> MatIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::matrix<double>* > & vMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA);
	/**
	 * Path interpolation for a set of input  vector ofmatrices
	 * Interpolate a series of values, with corresponding factor A and B (altitude and SZA) to a unique distribution of A and B (altitude and SZA).
	 * Typically, for a set of computed values (ex: densities) for a set of altitudes and SZA, computes these values on a path defined by an unique altitude and SZA for each point.
	 * \warning The output value can be positive or negative, if you want only positive value (eg densities), you can use the NoNegative function on the result.
	 *
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vVecMatInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \param vNoNegative : true by default : all the negative values are set to 0.
	 * \return the values for each points
	 */
	ublas::vector< ublas::matrix<double> > VecMatIntLinPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector< ublas::matrix<double> >* > & vVecMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA,bool vNoNegative=true);
	/**
	 * Log implementation if VecMatIntLinPath
	 * \param vInputAlt : the input A factor (altitude). The A grid must be the same for every input values (each densities must be computed on the same altitude grid for each SZA)
	 * \param vInputSZA : the input B factor (list of SZA)
	 * \param vVecMatInputVals : deque containing  the values ( depending on B and A, vInputSZA[B][A])
	 * \param vOutputAlt : the output A (altitudes)
	 * \param vOutputSZA : the output B (for each A)
	 * \return the values for each points
	 */

	ublas::vector< ublas::matrix<double> > VecMatIntLogPath(const ublas::vector<double>& vInputAlt, const std::deque<double>& vInputSZA, const std::deque< ublas::vector< ublas::matrix<double> >* > & vVecMatInputVals,const ublas::vector<double>& vOutputAlt,const ublas::vector<double>& vOutputSZA);

	/**
	 *
	 * Integration
	 *        /
	 * return=| y(x) dx
	 *        /
	 *
	 * \param x the value dx
	 * \param y the value to integrate over x
	 * \param nb : -1 (default) if we integrate over the whole vectors (in that case, raise an error if sizes different). Else, integrate over the given number of parameters (raise an error...)
	 *
	 * old hint
	 */

	double TrapzInt(const ublas::vector<double>& x,const ublas::vector<double>& y,int nb=-1);


	/**
	 * Allows to compare two vectors of double : the values of the vectors can be very close, for example when comparing results with two
	 * different architecture, but not equal. This function is very important for unit testing.
	 * Returns false of it is out of the accepted error
	 *
	 * \param vVec1 the first vector to consider
	 * \param vVec2 the second vector to consider
	 * \param vPercent the percentage of error  2*abs(vVec[i]-vVec[2])/(abs(vVec1[i])+abs(vVec2[i]))*100 < vPercent by default at 1%
	 * \param vMin : minimum value for comparison. Below that value, errors above 100% are allowed if vPercent<100%! By default at 1E-20. Launch a warning if difference > 100%.
	 */

	bool VectorCompare(const ublas::vector<double>& vVec1,const ublas::vector<double> & vVec2, double vPercent=1,double vMin=1E-20);


	/**
	 * 2 dimensional vector comparison, the options are the same as VectorCompare
	 * Returns false of it is out of the accepted error
	 *
	 * \param vVec1 the first vector to consider
	 * \param vVec2 the second vector to consider
	 * \param vPercent the percentage of error  2*abs(vVec[i]-vVec[2])/(abs(vVec1[i])+abs(vVec2[i]))*100 < vPercent by default at 1%
	 * \param vMin : minimum value for comparison. Below that value, errors above 100% are allowed if vPercent<100%! By default at 1E-20. Launch a warning if difference > 100%.
	 */

	bool VectorCompare2D(const ublas::matrix<double> & vVec1,const ublas::matrix<double>  & vVec2, double vPercent=1,double vMin=1E-20);

        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo----------------------------------------------------+
         |                                                                   |
         |Function for the definition of new profiles (gaussian, chapman,...)|
         |                                                                   |
         +------------------------------------------------------------------*/

	/**
	 * Returns a profile following the Bates-Walker algorithm
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude 0, necessary for initialization
	 * \param vN0cm_3 : the density at altitude 0
	 * \param vT0K : the temperature at altitude 0
	 * \param vTexoK : the exospheric temperature
	 * \param vS : the shape of the function
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particles in AMU
	 */
	ublas::vector<double> BatesWalkerProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vT0K, double vTexoK, double vS, double vRkm, double vGoms_2, double vMassamu);



	/**
	 * Returns a profile following  an exponential law (Scale height)
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude 0, necessary for initialization
	 * \param vN0cm_3 : the density at altitude 0
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particles in AMU
	 */
	ublas::vector<double> ExpProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vRkm, double vGoms_2, double vMassamu);
	
	/**
	 * Returns a profile following  an exponential law (Scale height)
	 * for the species, with the addition of another exponential law
	 * for the mixed atmosphere
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude 0, necessary for initialization (mesospheric influence should be negligible)
	 * \param vAlt1Km : the altitude 1, mesosphere, necessary for initialization (thermospheric influence should be negligible)
	 * \param vN0cm_3 : the density at altitude 0 (thermosphere)
	 * \param vN1cm_3 : the density at altitude 1 (mesosphere)
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vTmesoK : the mesospheric temperature (for the scale height)
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particle in AMU
	 * \param vMixMassamu : the mean particle mass in mesosphere in AMU
	 */
	ublas::vector<double> DoubleExpProfile(ublas::vector<double> vAltKm, double vAlt0Km, double vAlt1Km, double vN0cm_3, double vN1cm_3, double vTexoK, double vTmesoK, double vRkm, double vGoms_2, double vMassamu, double vMixMassamu);
	/**
	 * Returns a profile following  an exponential of an hyperbolic law
	 * the asymptotes of the hyperbolic law are the barometric profiles
	 * for  the specie (with Texo) at high altitude, and the mixed mesospheric
	 * at the mesosphere. The vC parameter is characteristic of the transition
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAltmesoKm : the altitude 0, necessary for initialization (mesospheric influence should be negligible)
	 * \param vAltexoKm : the altitude 1, mesosphere, necessary for initialization (thermospheric influence should be negligible)
	 * \param vNmesocm_3 : the density at altitude 1 (mesosphere)
	 * \param vNexocm_3 : the density at altitude 0 (thermosphere)
	 * \param vTmesoK : the mesospheric temperature (for the scale height)
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vC : the parameter for the hyperbola, characteristic 
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particle in AMU
	 * \param vMixMassamu : the mean particle mass in mesosphere in AMU
	 */
	ublas::vector<double> ExpHyperbolaProfile(ublas::vector<double> vAltKm, double vAltmesoKm, double vAltexoKm, double vNmesocm_3, double vNexocm_3, double vTmesoK, double vTexoK, double vC, double vRkm, double vGoms_2, double vMassamu, double vMixMassamu);



	/**
	 * Makes a gaussian profile
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAltPeakKm : the altitude of the peak
	 * \param vN0cm_3 : the density at the peak
	 * \param vDevKm : the deviation in Km
	 */

	ublas::vector<double> GaussianProfile(ublas::vector<double> vAltKm, double vAltPeakKm, double vN0cm_3, double vDevKm);


	/**
	 * Makes a Chapman profile
	 *
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAltPeakKm : the altitude of the peak
	 * \param vN0cm_3 : the density at the peak
	 * \param vSHkm : the SH, shape+scale height in Km
	 */
	ublas::vector<double> ChapmanProfile(ublas::vector<double> vAltKm, double vAltPeakKm, double vN0cm_3, double vSHkm);

	/**
	 * Makes a Chapman profile, cos corrected
	 *
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude of the peak
	 * \param vN0cm_3 : the density at the peak
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vSZAdeg : the SZA in degree
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particle in AMU
	 */
	ublas::vector<double> ChapmanCos(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vSZAdeg, double vRkm, double vGoms_2, double vMassamu);

	/**
	 * Makes a Chapman profile, with a multiplicative factor in the exp
	 *
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude of the peak
	 * \param vN0cm_3 : the density at the peak
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vC : half of the multiplicative factor (must be 1 for beta chapman, 0.5 for alpha chapman)
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particle in AMU
	 */
	ublas::vector<double> ChapmanVar(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vC, double vRkm, double vGoms_2, double vMassamu);

	/**
	 * Makes an Epstein profile
	 *
	 * \param vAltKm : the altitudes where we compute the profile
	 * \param vAlt0Km : the altitude of the peak
	 * \param vN0cm_3 : the density at the peak
	 * \param vTexoK : the exospheric temperature (for the scale height)
	 * \param vRkm : the radius of the planet in Km
	 * \param vGoms_2 : the acceleration of gravity at the surface in  m/s2
	 * \param vMassamu : the mass of the particle in AMU
	 */
	ublas::vector<double> Epstein(ublas::vector<double> vAltKm, double vAlt0Km, double vN0cm_3, double vTexoK, double vRkm, double vGoms_2, double vMassamu);


};




#endif
