/**
 * \file geoobs.hpp
 * \brief  Define the experiment : the set of observations. Can be seen as a collection of line of sight
 * Copyright G Gronoff 2010
 * Last Modification : $Id: geoobs.hpp 1594 2012-10-30 08:20:29Z gronoff $
 */


#ifndef GEO_OBS_HPP
#define GEO_OBS_HPP
#include "geopath.hpp"



/**
 * \ingroup GeoObserv
 * Reads the definition of the GeoPaths in the XML file, create and store the
 * list of theses Paths.
 * This object allows to make a comparison between the simulation and the experiment by making a realistic comparison.
 *
 *
 */
class GeoObs
{
	private:
		/// Main parameter file
		XmlParameters* mpParam;

		/// Reads the parameter file
		void ReadParam();

		/// If the satellite is fixed
		bool mbIsSatFixed;
		/// The altitude of the satellite
		std::deque<double> mSatAltitudes;
		/// The latitude of the satellite
		std::deque<double> mSatLatDegree;
		/// The longitude of the satellite
		std::deque<double> mSatLonDegree;
		/// The azimuth of the path
		std::deque<double> mAzimut;
		/// The tangent altitude of the path (either mDec or mTanAltKm defined, not the two)
		std::deque<double> mTanAltKm;
		/// The declination of the path (either mDec or mTanAltKm defined, not the two)
		std::deque<double> mDec;

		/// The list of Paths which wants in our experiment
		std::deque<GeoPath> mPaths;
		
		/// True is the tangent altitude is defined (therefore declination is undefined)
		bool mbIsTanDefined;
		/// True if the declination is defined (therefore tangent altitude in undefined)
		bool mbIsDecDefined;

		/**
		 * Reads a node (satellite altitude, latitude, longitude, azimuth, tangent altitude or declination), to find if a constant, a list or a range is defined. Returns the defined value (1 value inside the deque if it is constant, the corresponding else). 
		 * \param vNode : the node where to search the value
		 */
		std::deque<double> ReadNode(TiXmlNode* vNode);
		/**
		 * Reads a range node.
		 * This function is simple now because it is only possible to have a constant range of data. But could be extended in the future
		 * \param vNode : the node of the range
		 * \return the list of values, defined by the node.
		 *
		 */
		std::deque<double> ReadRange(TiXmlNode* vNode);

		/**
		 * Creates the paths by using the data read inside the xml file
		 */
		 void CreatePath();

		 /**
		  * The radius of the planet in km
		  */
		 double mRpKm;

		 /**
		  * The out point of the planet in Km
		  */
		 double mOutAtmoKm;

		 /**
		  * The number of points in the paths
		  */
		 unsigned mNbPoints;

		 /**
		  * True if we use the geometry
		  */
		 bool mbIsGeoUsed;


	public:
		/**
		 * The constructor. Reads the parameters and creates the path
		 * \param vpParams : the pointer to the XmlParameters file
		 * \param vRpKm : the radius of the planet in Km
		 * \param vOutAtmoKm : the output of the planet (can be overloaded by the xml
		 */
		GeoObs(XmlParameters* vpParams,double vRpKm,double vOutAtmoKm);

		/**
		 * To check if the GeoObs class should be used
		 *
		 * \return true if the path is used
		 */
		bool IsPath()
		{
			return mbIsGeoUsed;
		}

		/**
		 * Allows to integrate a value on the list of paths
		 * If the Satellite altitude is above the maximum altitude of the atmosphere, the data are extrapolated.
		 * \param vAltGridKm : the altitude grid of the values to integrates
		 * \param vValue : the values corresponding to the altitude
		 * \return the integrated list.
		 */
		ublas::vector<double> PathIntegrate(ublas::vector<double> vAltGridKm, ublas::vector<double> vValue);

		/**
		 * Allows to integrate a value on the list of paths, with an absorption
		 * If the Satellite altitude is above the maximum altitude of the atmosphere, the data are extrapolated.
		 * \param vAltGridKm : the altitude grid of the values to integrates
		 * \param vValue : the values corresponding to the altitude
		 * \param vAbsorb_Km : the value absorption coefficient units of value/Km
		 * \return the integrated list.
		 */
		ublas::vector<double> PathIntegrateAbsorbed(ublas::vector<double> vAltGridKm,ublas::vector<double> vValue,ublas::vector<double> vAbsorb_Km);


		/// Return true if the tangent altitude are defined
		bool IsTanDefined()
		{
			return mbIsTanDefined;
		}

		/// Return true if the observation declination are defined
		bool IsDecDefined()
		{
			return mbIsDecDefined;
		}
		/// Return the tangent altitude list in Km
		std::deque<double> GetTanList()
		{
			return mTanAltKm;
		}

		/// Return the declination list in degree
		std::deque<double> GetDecList()
		{
			return mDec;
		}

};






#endif
