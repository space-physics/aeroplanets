/**
 * \defgroup GeoObserv Observation point
 * \file geopoint.hpp
 * \brief Defines the geographical point 
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geopoint.hpp 892 2010-02-24 19:40:56Z gronoff $
 */

#ifndef GEO_POINT_HPP
#define GEO_POINT_HPP



#include <planet/allplanets.hpp>




/**
 * \ingroup GeoObserv
 * Defines a point on a sphere (planet are considered as spheres).
 * Allows to work with lat-lon coordinates on that sphere (subsolar coordinates!)
 */

class GeoPoint
{
	protected:
		/// The latitude in radians
		double mLatrad;
		/// The longitude in radians
		double mLonrad;
		/// The altitude in Km
		double mAltKm;
		/// The radius of the planet in Km
		double mRpKm;

		/// The distance of the point to the center of the planet (sum of Alt and Rp)
		double mRKm;
		/// The SZA in radians
		double mSZArad;

		/// Cartesian coordinates: The direction toward sun, X, in km
		double mXKm;
		/// Cartesian coordinates: Y in km
		double mYKm;
		/// Cartesian coordinates:  the direction toward the 'pole' in km (subsolar coordinates, it can actually be really different of the planet pole)
		double mZKm;


		/// Converts from geographic coordinates into cartesian coordinates
		void GeoToCar();
		/// Converts from cartesian coordinates into geographic (subsolar) coordinates
		void CarToGeo();
		/// Converts the coordinates (geo actually) into SZA
		void PositionToSZA();
		/// Returns a deque with the cartesian coordinates, reduced on the unit sphere
		std::deque<double> GetCartesianReduced()
		{
			std::deque<double> resu;
			resu.push_back(mXKm/mRpKm);
			resu.push_back(mYKm/mRpKm);
			resu.push_back(mZKm/mRpKm);
			return resu;
		}
		

	public:

		/// Returns the Latitude in degree
		double GetLatDegree()
		{
			return mLatrad*180./PI;
		}
		/// Returns the Longitude in degree
		double GetLonDegree()
		{
			return mLonrad*180./PI;
		}

		/// Returns the SZA in degree
		double GetSZADegree()
		{
			return mSZArad*180./PI;
		}
		/// Returns the altitude
		double GetAlt()
		{
			return mAltKm;
		}
		/// Returns the distance to the center
		double GetDist()
		{
	//		return mRpKm+mAltKm;
			return mRKm;
		}
		/// Returns the planet radius
		double GetRp()
		{
			return mRpKm;
		}
		/// Returns a deque with the cartesian coordinates
		std::deque<double> GetCartesian()
		{
			std::deque<double> resu;
			resu.push_back(mXKm);
			resu.push_back(mYKm);
			resu.push_back(mZKm);
			return resu;
		}


		/** The constructor, initialize the  radius of the planet
		 * \param vRpKm the radius of the planet in Km
		 */
		GeoPoint(double vRpKm);


		/**
		 * Returns a copy of the present point
		 */
		GeoPoint Copy();

		/**
		 * Set the position of the point, in subsolar coordinates
		 * \param vAltKm : the altitude of the point
		 * \param vLatDegree : the latitude in degree
		 * \param vLonDegree : the longitude in degree
		 */
		void SetGeoPosition(double vAltKm,double vLatDegree,double vLonDegree);

		/**
		 * Set the position of the point, in cartesian coordinates ( X towards the sun)
		 * \param vXKm
		 * \param vYKm
		 * \param vZKm
		 */
		void SetCartPosition(double vXKm, double vYKm, double vZKm);

		
		/**
		 * Gives the angle between the present point, and the point in parameter, with respect to the planet center, in Degree
		 * \param vC : the second point to compute the angle
		 * \return : the angle between the present point, the center, and vC
		 */
		double EcartAngleDegree(GeoPoint vC);


		/**
		 * Returns the altitude of the tangent point at the given altitude and longitude, if we suppose a line of sight starting from the present point
		 * \param vLatDegree : the latitude of the tangent point
		 * \param vLongDegree : the longitude of the tangent point
		 * \return the altitude of the tangent point (if < 0; it means it is inside the planet. This point is always defined since mRKm=0 is always a solution).
		 */
		double ReturnTangentPointAltitude(double vLatDegree, double vLongDegree);


		/**
		 * Given an Azimuth (0 towards north), returns a geopoint on the great circle, defined by the present point and the azimuth direction, at 90 Degree from the present point.
		 * \param vAzimuthDegree: the azimuth of the great circle
		 * \return The GeoPoint
		 */
		GeoPoint ReturnTrajPoint(double vAzimuthDegree);


		/**
		 * Returns the GeoPoint in the azimuth direction, at altitude vHoKm, so that the line of sight from the present point to this point is tangent at that point.
		 * \param vAzimuthDegree : The azimuth
		 * \param vHoKm The altitude of the tangent point
		 * \return The GeoPoint
		 */
		GeoPoint ReturnSubPoint(double vAzimuthDegree,double vHoKm);


		/**
		 * Returns the GeoPoint in the azimuth direction, at Declination direction (down!) so that the line of sigth is tangent at that point
		 * \param vAzimuthDegree : the direction azimuth
		 * \param vDecDegree : the direction declination
		 * \return The tangent point
		 */

		GeoPoint ReturnSubViz(double vAzimuthDegree, double vDecDegree);

		
		/**
		 * Returns the point where the los goes out of the atmosphere, by using the azimuth and declination of the los, and the limit altitude of the atmosphere
		 * \param vAzimuthDegree the LOS azimuth
		 * \param vDecDegree the LOS declination, NEGATIVE
		 * \param vAltOutKm : the limit altitude of the atmosphere
		 * \return the out point
		 */
		GeoPoint ReturnSortieAtmo(double vAzimuthDegree, double vDecDegree, double vAltOutKm);

		
		/**
		 * Returns the point where the los goes out of the atmosphere, by using the azimuth and declination of the los, and the limit altitude of the atmosphere. This one has no restrictions on the declination, wich can be positive.
		 * \param vAzimuthDegree the LOS azimuth
		 * \param vDecDegree the LOS declination, without restrictions
		 * \param vAltOutKm : the limit altitude of the atmosphere
		 * \return the out point
		 */
		GeoPoint ReturnAllSortieAtmo(double vAzimuthDegree, double vDecDegree, double vAltOutKm);




		/**
		 * Returns the point where the los goes out of the atmosphere, by using the azimuth the los, the altitude of the tangent point and the limit altitude of the atmosphere
		 * \param vAzimuthDegree the LOS azimuth
		 * \param vHoKm : the tangent point altitude
		 * \param vAltOutKm : the limit altitude of the atmosphere
		 * \return the out point
		 */
		GeoPoint ReturnSortieAtmoHo(double vAzimuthDegree, double vHoKm, double vAltOutKm);
		/**
		 * Add two geopoint
		 * \param vA : the first point
		 * \param vB : the second point
		 */
		friend GeoPoint operator+(const GeoPoint& vA,const GeoPoint& vB)
		{

			assert(vA.mRpKm==vB.mRpKm);
			GeoPoint resu(vA.mRpKm);
			resu.SetCartPosition(vA.mXKm+vB.mXKm,vA.mYKm+vB.mYKm,vA.mZKm+vB.mZKm);
			return resu;
		}
		/**
		 * Substract two geopoint
		 * \param vA : the first point
		 * \param vB : the second point
		 */
		friend GeoPoint operator-(const GeoPoint& vA,const GeoPoint& vB)
		{
			assert(vA.mRpKm==vB.mRpKm);
			GeoPoint resu(vA.mRpKm);
			resu.SetCartPosition(vA.mXKm-vB.mXKm,vA.mYKm-vB.mYKm,vA.mZKm-vB.mZKm);
			return resu;
		}	
		/**
		 * Multiply the geopoint by a constant
		 * \param vA : the first point
		 * \param vMult : the constant
		 */
		friend GeoPoint operator*(const GeoPoint& vA,double vMult)
		{

			GeoPoint resu(vA.mRpKm);
			resu.SetCartPosition(vA.mXKm*vMult,vA.mYKm*vMult,vA.mZKm*vMult);
			return resu;
		}
		/**
		 * Multiply the geopoint by a constant
		 * \param vA : the first point
		 * \param vMult : the constant
		 */
		friend GeoPoint operator*(double vMult, const GeoPoint& vA)
		{
			return vA*vMult;
		}
};





#endif

