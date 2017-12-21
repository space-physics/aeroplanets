/**
 * \file path.hpp
 * \brief Defines the path class, used to have path in the atmosphere 
 * (for the magnetic field lines)
 * Copyright G Gronoff January 2010
 * Last Modification : $Id: path.hpp 1336 2011-11-08 22:11:28Z gronoff $
 */

#ifndef PATH_HPP
#define PATH_HPP


#include "eleccrs.hpp"
#include "protoncrs.hpp"


/**
 * Class to create a path in the atmosphere of a planet.
 * Allows to make transformations between the different geometries
 */
class Path
{
	public:
		/// The vector of altitudes: altitudes of each point of the path in km
		ublas::vector<double> mAltitudeKm;
		/// The vector of length: distance to the bottom point (remember that the paradigm is to have decreasing arrays)
		ublas::vector<double> mLengthKm;
		/// The vector of SZA
		ublas::vector<double> mSZADegree;

		/// Simple constructor
		Path()
		{}
		/**
		 * Constructor initialing all the arrays
		 * \param vAltKm : the altitude vector
		 * \param vLenKm : the length vector
		 * \param vSZADeg : the SZA vector
		 */
		Path(ublas::vector<double> vAltKm,ublas::vector<double> vLenKm,ublas::vector<double> vSZADeg);

		/**
		 * Init the path with the Altitude and length arrays
		 * If the two have not the same size, raise an error.
		 * Else put the SZA vector at the default value.
		 * \param vAltitudeKm : the altitude vector
		 * \param vPathKm : the path length vector
		 * \param vSZADegree : the default value SZA 
		 */
		void InitAltLen(ublas::vector<double> vAltitudeKm, ublas::vector<double> vPathKm,double vSZADegree);

		/**
		 * Init the path with the Altitude, length and SZA arrays
		 * If the tree have not the same size, raise an error.
		 * \param vAltitudeKm : the altitude vector
		 * \param vPathKm : the path length vector
		 * \param vSZADegree : the  SZA vector
		 */
		void InitAltLenSZA(ublas::vector<double> vAltitudeKm, ublas::vector<double> vPathKm,ublas::vector<double> vSZADegree);



		/**
		 * Creates the  vectors by taking into account 3 cartesian coordinates (vXRp toward sun). Each coordinate is given in unit of planet radius, the planet radius is given in parameter)
		 * \param vXRp: Coordinates of the points with respect to X, sunward direction. In unit of planet radius.
		 * \param vYRp : Coordinates with respect to Y. In unit of planet radius.
		 * \param vZRp : Coordinates with respect to Z. In unit of planet radius.
		 * \param vRpKm : Radius of the planet in Km
		 */
		void CartesianToPath(ublas::vector<double> vXRp,ublas::vector<double> vYRp,ublas::vector<double> vZRp,double vRpKm);

		/// get the altitude vector
		ublas::vector<double> GetAltKm()
		{
			return mAltitudeKm;
		}
		/// get the length vector
		ublas::vector<double> GetLenKm()
		{
			return mLengthKm;
		}
		/// get the SZA vector
		ublas::vector<double> GetSZADeg()
		{
			return mSZADegree;
		}


};



#endif
