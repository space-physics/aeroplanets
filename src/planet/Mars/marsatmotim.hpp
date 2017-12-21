/**
 * \file marsatmotim.hpp
 * \brief Defines the martian atmosphere from MarTim lecture
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: marsatmotim.hpp 831 2010-01-08 17:58:49Z gronoff $
 */

#ifndef MARS_ATMO_TIM_HPP
#define MARS_ATMO_TIM_HPP
#include "../planet.hpp"


const unsigned MARTIM_NEUTRAL_SPECIES = 7;
const unsigned MARTIM_ION_SPECIES = 6;
const unsigned MARTIM_POSITION_EDENS = 5;
typedef boost::multi_array<double, 3> bcube;
typedef bcube::index bcindex;
//typedef bcube::extent_gen bcextents;

typedef boost::multi_array<double, 4> bfour;
typedef bfour::index bqindex;
//typedef bfour::extent_gen bfextents;

class MarsAtmoTim
{
	private:

		/// The composition filename
		std::string mCompoFilename;

		/// The temperature filename
		std::string mTempFilename;

		/// Read the composition and fills the vector
		void ReadCompo();

		/// Read the temperature and fills the vector
		void ReadTemp();

		/// number of lat 
		unsigned mLatDim;
		/// number of long 
		unsigned mLonDim;
		/// number of altitudes 
		unsigned mHDim;

		/** the neutral atmo composition [specie][lat][long][alt]
		 *
		 * N2 O CO2 Ar CO O2 NO
		 *
		 */
		bfour mAtmosphere;

		/** the ionospheric composition [specie][lat][long][alt]
		 *
		 * CO2+ N2+ O+ O2+ CO+ e
		 *
		 */
		bfour mIonosphere;

		/// the altitude cube : dependance of the alt because computation on isobares [lat][long][alt]
		bcube mAltitude;

		/// the neutral temperature [lat][long][alt]
		bcube mTempNeutre;

		/// the mass of the atmo [lat][long][alt]
		bcube mMassAtmo;




		/// The solar longitude of the sun
		double mLsDeg;
		/// The longitude in respect to planet coordinates
		double mLocalLongDeg;
		/// The latitude in respect to planet coordinates
		double mLocalLatDeg;
		
		/// The longitude in respect to the subsolar point
		double mLongDeg;
		/// The latitude in respect to the subsolar point
		double mLatDeg;

		/// The solar declination angle
		double mSdaDeg;

		/// The solar zenith angle
		double mSZADeg;


		/**
		 * Find the solar declination angle in respect to the Ls
		 * \param vLsDeg : the solar longitude
		 */
		double LsToSda(double vLsDeg);

		/**
		 * Modification of the coordinate system. The passage from the two systems is given by the  vPsiDeg parameter
		 * if vPsiDeg= mSdaDeg this code changes the subsolar coordinate to local coordinate 
		 * if vPsiDeg= -mSdaDeg this code changes the local coordinates to subsolar coordinate
		 * if vPsiDeg = 90-latitude this codes changes the Hour coordinates (proportional to alpha-deg, but with a time dependance) into azimuth-dec
		 * if vPsiDeg = -(90-latitude) this codes changes the azimuth-deg into Hour coordinates 
		 * \param vLatDeg : the latitude in degree in the first frame
		 * \param vLonDeg : the longitude in degree in the first frame
		 * \param vPsiDeg : the Difference in degree between the two poles
		 * \param rLatDeg : the latitude in degree in the final frame -> result, reference
		 * \param rLongDeg : the longitude in degree in the final frame -> result, reference
		 */
		void ChgtCoord(double vLatDeg,double vLonDeg, double vPsiDeg, double& rLatDeg,double& rLongDeg);


		/**
		 * Use mLocalLatDeg and mLocalLongDeg to find the closest position in the martim array
		 * it fills the mMartimLat and mMartimLong values
		 *
		 */
		void InitRelPosition();

		/// The corresponding latitude in the martim matrix. computed when initlocalcoords or initsubsolarcoords
		unsigned mMartimLat;
		/// The corresponding longitude in the martim matrix. computed when initlocalcoords or initsubsolarcoords
		unsigned mMartimLong;

	public:

		/**
		 * Read the composition and the temperature filenames
		 * \param vCompoFilename : the composition filename
		 * \param vTempFilename : the temperature filename
		 */
		MarsAtmoTim(std::string vCompoFilename, std::string vTempFilename);

		void InitLocalCoords(double vLsDeg,double vLocalLatDeg,double vLocalLongDeg);

		void InitSubsolarCoords(double vLsDeg, double vLatDeg, double vLongDeg);


		/// gives the local latitude in degree
		double ReturnLocalLatDeg()
		{
			return mLocalLatDeg;
		}

		/// gives the  local longitude in degree
		double ReturnLocalLongDeg()
		{
			return mLocalLongDeg;
		}
		/// gives the subsolar latitude in degree
		double ReturnLatDeg()
		{
			return mLatDeg;
		}

		/// gives the subsolar  longitude in degree
		double ReturnLongDeg()
		{
			return mLongDeg;
		}
		/// gives the ls in degree
		double ReturnLsDeg()
		{
			return mLsDeg;
		}
		/// gives the sza in degree

		double ReturnSZADeg()
		{
			return mSZADeg;
		}



		/**
		 * Returns the Martim neutral atmosphere
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return vector of vector : the neutral atmosphere
		 */
		std::deque< ublas::vector<double> > ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);
		/**
		 * Returns the Martim  ionospheric composition and densities
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * The last value is the electron density!
		 * \return vector of vector : the ionosphere
		 */
		std::deque< ublas::vector<double> > ReturnIono(const ublas::vector<double>& vAltitudeGridKm);
		
		
		/**
		 * Returns the Martim  neutral temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm);
		
};







#endif
