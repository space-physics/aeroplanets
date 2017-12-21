/**
 * \file geopath.hpp
 * \brief Defines the path of a line of sight
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geopath.hpp 1595 2012-11-03 01:30:26Z gronoff $
 */

#ifndef GEO_PATH_HPP
#define GEO_PATH_HPP


#include "geopoint.hpp"


/**
 * \ingroup GeoObserv
 * Allows to define a path in the atmosphere, create the points on this path, and integrate some values on this path
 */
class GeoPath
{
	protected: 
		/// The vector of altitudes: altitudes of each point of the path in km
		ublas::vector<double> mAltitudeKm;
		/// The vector of length: distance to the bottom point (remember that the paradigm is to have decreasing arrays)
		ublas::vector<double> mLengthKm;
		/// The vector of SZA
		ublas::vector<double> mSZADegree;


		/// The path, in terms of GeoPoint
		std::deque<GeoPoint> mPointsPath;

		/// To transform the GeoPoint path into the arrays, it can work with non equidistants points in the path
		void PointsPathToArrays();

		/// The start point for the computation
		GeoPoint mStartPoint;
		/// The end point for the computation
		GeoPoint mEndPoint;
		/// Number of points in the Path
		unsigned mNbPoints;
		/** 
		 * Create the geo points path between mStartPoint and mEndPoint
		 * given the number of points you whish between the two.
		 * Then fill the arrays (without calling pointspath...)
		 * \param vNbPoints : the number of points in the path
		 */
		void StraightExtremas(unsigned vNbPoints);

		/// Check if the path is correctly defined
		bool mbIsPathDefined;

		
	public:
		/**
		 * Initialize the path with two points, and ask for making a straight line between them with vNbPoints
		 * \param vStartPoint : the observer point
		 * \param vEndPoint : the out of the atmosphere point
		 * \param vNbPoints : the number of points in the Path
		 */
		GeoPath(GeoPoint vStartPoint, GeoPoint vEndPoint,unsigned vNbPoints);
		
		/** 
		 * Initialize the path for ionization computation
		 * \param vAltGridKm: the altitude grid in km in decreasing order
		 * \param vSZAdeg: the solar zenith angle in degree
		 * \param vRplanetKm: the radius of the planet in Km
		 */
		GeoPath(ublas::vector<double> vAltGridKm, double vSZAdeg, double vRplanetKm);

		/**
		 * Initialize the path with the observer point
		 * In that case, the path functions should be used
		 * (it is impossible in the constructor to make a difference
		 * between the declination and the altitude)
		 * \param vStartPoint : the position of the observer.
		 */
		GeoPath(GeoPoint vStartPoint);

		/** 
		 * Creates a GeoPath as wanted
		 * uses 
		
GeoPoint GeoPoint::ReturnAllSortieAtmo(double vAzimuthDegree, double vDecDegree, double vAltOutKm)
		 */

		/**
		 * Creates the Path with the azimuth of the path, the altitude of the tangent point, and the out altitude of the atmosphere
		 * \param vAzimuthDegree : the azimuth
		 * \param vTanAltKm : the altitude of the tangent point
		 * \param vOutAltKm : the maximum altitude of the atmosphere
		 * \param vNbPoints : the number of points in the path
		 */
		void InitAzTanPath(double vAzimuthDegree,double vTanAltKm,double vOutAltKm,unsigned vNbPoints);

		/**
		 * Creates the Path with the azimuth of the path,  the declination of the path (no restrictions positive or negative), and the out altitude of the atmosphere.
		 * \param vAzimuthDegree : the azimuth
		 * \param vDecDegree : the declination
		 * \param vOutAltKm : the maximum altitude of the atmosphere
		 * \param vNbPoints : the number of points in the path
		 */
		void InitAzDecPath(double vAzimuthDegree,double vDecDegree,double vOutAltKm,unsigned vNbPoints);


		/**
		 * Integrates a value defined by its altitude and 'intensity', on the Path. An interpolation-extrapolation is made between this alt-value  on the path altitudes, and then integrated
		 * \param vAltKm : the values altitudes
		 * \param vVals : the values to integrate on the Path
		 * \return the value of the integration, in value_units*KM. \warning If you want to work in m or cm, you have to make an unit change!
		 */

		double PathIntegrate(ublas::vector<double> vAltKm, ublas::vector<double> vVals);
		/**
		 * Integrates a value defined by its altitude and 'intensity', on the Path. An interpolation-extrapolation is made between this alt-value  on the path altitudes, and then integrated. In that case, an absorption exists, in terms of exp(- int vAbsorption_Km dkm).
		 * \param vAltKm : the values altitudes
		 * \param vVals : the values to integrate on the Path
		 * \param vAbsorption_Km : the absorption in unit of X*km-1 (the integration is directly performed!) X should be coherent with value_units to compute the absorption.
		 * \return the value of the integration, in value_units*KM. \warning If you want to work in m or cm, you have to make an unit change!
		 */

		double PathIntegrateAbsorbed(ublas::vector<double> vAltKm, ublas::vector<double> vVals, ublas::vector<double> vAbsorption_Km);



		/* Test the technique to interpolate the point from a list of altitudes */
		void TestNewPath(ublas::vector<double> vAltListKm);

		/**
		 * Changes the separation of the points in a path 
		 * to follow the given altitude grid
		 * \warning When the path goes below the lowest altitude point in the grid
		 * there is some problem in the interpolation, creating a bad value.
		 * \todo put the value for these point at a high default!
		 * \param vAltListKm : the altitude grid in Km
		 */
		void ResetGrid(ublas::vector<double> vAltListKm);


		/** Returns the altitude grid
		 * \return altitude grid in km for the path 
		 */
		ublas::vector<double> ReturnAltGrid()
		{
			return mAltitudeKm;
		}


		/** Returns the length grid
		 * \return length grid in km for the path 
		 */
		ublas::vector<double> ReturnLenGrid()
		{
			return mLengthKm;
		}
};





#endif
