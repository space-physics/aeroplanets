/**
 * \defgroup Planetes Planetes
 * \file planet.hpp
 * to define the planet, and the position on the planet if we do a vertical computation.
 * It also defines the atmosphere model.
 *
 * \brief Defines the planet abstract class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: planet.hpp 1356 2011-11-16 21:20:46Z gronoff $
 *
 */


#ifndef GEOPHYS_HPP
#define GEOPHYS_HPP

#include <species/species.hpp>
#include <math/mathfunction.hpp>





/**
 * \ingroup Planetes
 *
 * Class to compute planet specific parameters:
 * Distance to the sun
 *
 *
 *
 *
 */
class Planete
{
	protected:
		/// The parameters object
		XmlParameters* mpParameter;

		/** To store a possible error in loadcoords : the error can be
		 * neglected when it is overloaded by something else!
		 */
		bool mIsLoadCoordError;

	public:


		/** Void initialization: 
		 * the UA need to be initialized, as for the other parameters
		 * \param pParam : the xmlparameters object -> Planet can read these
		 */
		Planete(XmlParameters* pParam);

		/**
		 * The distance to the sun is initialized
		 *
		 * \param pParam : the xmlparameters object -> Planet can read these
		 * \param vDistSun : the distance to the sun
		 */
		Planete(XmlParameters* pParam,double vDistSun);

		/**
		 * Destruction of the object
		 */
		virtual ~Planete();
		
		/**
		 * Display the distance of the planet to the sun
		 */
		void ShowDistance();

		/**
		 * Display the name of the planet
		 */ 
		void ShowName();


		/**
		 * Display the important facts about the planet
		 *
		 */
		void ShowData();

		/**
		 * Set the latitude/longitude. On this point, the 0:0 refers to the subsolar point.
		 * \param vLatitudeDegree: the latitude
		 * \param vLongitudeDegree : the longitude
		 */
		void SetCoord(double vLatitudeDegree,double vLongitudeDegree);


		// "Typical usage of the function. For very important, and probably redefined - improved later"
		/**
		 * Returns the atmospheric model for the considered species and the the considered altitude grid.
		 *
 		* \xrefitem usage "usage"  "Usage" Is called from the neutral atmosphere class when the atmosphere is not an interpolation of the parameters: typically, it loads an atmosphere model, and use the parameters stored in XmlParameters to read the other parameters.
		 * \param vAltitudeGridKm : the altitude grid where the densities are computed
		 * \param vSpNames : vector containing the (string) name of the species searched. 
		 * \param vType :  the type of the model, starting at 1 (because the type parameter of the atmosphere part considers 0 as the interpolation of the data.)
		 * \return : map  containing the name of the specie against its density array (cm-3)
		 */
		virtual std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string >  vSpNames,int vType)=0;



		/**
		 * Defines the electron density with respect to the mode
		 * \param vAltGridKm : the altitude grid used
		 * \param vType : the type of the model
		 * \return a model of electron density cm-3
		 */
		virtual ublas::vector<double> ElectronDensity(const ublas::vector<double> & vAltGridKm,const int & vType)=0;


		/**
		 * Defines the electron density with respect to the mode
		 * \param vAltGridKm : the altitude grid used
		 * \param vType : the type of the model
		 * \return the electron temperature grid in K
		 */
		virtual ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)=0;

		/**
		 * Defines the electron density with respect to the mode
		 * \param vAltGridKm : the altitude grid used
		 * \param vType : the type of the model
		 * \return the ion temperature grid in K
		 */
		virtual ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType)=0;



		/// If the atmosphere defined the temperature
		bool mIsTemperatureDefined;
		/// The temperature grid if defined
		ublas::vector<double> mTemperatureModelGridK;



	protected:
		/**
		 * return true if SZA, hrloc, lat were successfully computed
		 * false else: does not kill the app because inherited 
		 * functions could use other parameters to compute this.
		 */
		bool LoadCoords();




	public:
                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo--------------------------------------+
                 |The values typically used when calculations are made:|
                 |        g value                                      |
                 |        UA distance                                  |
                 |        ...                                          |
                 +----------------------------------------------------*/


		/// The name of the planet
		std::string mName;

		/** The distance of the planet to the sun in the sun model
 		* \xrefitem usage "usage"  "Usage"  important for the sun model
		 */
		double mUA;

		/// The gravity at the surface of the planet (point considered). In m/s2
		double mGms_2;

		/// The radius of the planet (point considered). In km
		double mRKm;

		/// The latitude consiedered: \warning it is relative to the subsolar point. For planet with important inclination, the computation must be made inside the specialized part. (we need an inclination, and a reference date for a given equinoxe or solstice. When the solar longitude is known (ref: the eq) the computation, in python is done by return asin(sin(inclination*pi/180)*sin(ls*pi/180.)) )
		double mLatDegree;

		/// The longitude consiedered: \warning it is relative to the subsolar point
		double mLoDegree;

		/// The solar zenith angle.
		double mSZADegree;

		/// The local hour
		double mHrLoc;


		/*
		 *
		 * For the magnetic field
		 *
		 *
		 */

		/** Returns the vector for dB/B : the divergence of the magnetic field used for the mirror effect
		 */
		virtual ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm) = 0;

	protected:
		/**
		 * Reads the parameters in search for dB_B, or the type
		 * \param vAltGridKm : the altitude grid
		 * \param vrdB_B : the result if readed
		 * \return The dB_B model type; 0 if it is what was just read, if nothing is read, just give 0;
		 */
		int ReaddB_B(const ublas::vector<double>& vAltGridKm, ublas::vector<double>& vrdB_B);



};



#endif


