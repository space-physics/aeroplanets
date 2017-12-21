/**
 * \file venus.hpp
 * \brief Defines the Venus planet
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: venus.hpp 1356 2011-11-16 21:20:46Z gronoff $
 */


#ifndef VENUS_PLANET_HPP
#define VENUS_PLANET_HPP
#include "planet.hpp"
#include <cstring>
extern "C"
{
/** 
 * \ingroup Planetes
 * Loader for the Fortran Venus VTS3 model Hedin83
 * The venus model needs the vts3 model
 * \param altitude : the altitude considered
 * \param yrd : not used
 * \param xlat : the latitude considered
 * \param xloc : the local hour considered
 * \param f107av : the average over the past 3 monthes of the f107 parameter
 * \param f107 : the f107 parameter for the considered day
 * \param ma : mass. 48 in fact, not to modify unless you're an expert of this code
 * \param densite : array of output densities. Size is 7.  0-> total mass density (g/cm3); 1-> CO2 density (cm-3); 2->O; 3->CO; 4->He ; 5-> N ; 6 ->N2
 * \param temp : the output neutral temperature. Size is 2. 0-> exospheric temperature, 1-> local temperature
 */ 

	extern void vts3_(float* altitude,float*yrd,float*xlat,float*xloc,float*f107av,float*f107,int*ma,float*densite,float*temp);
/** 
 * \ingroup Planetes
 * Loader for the Fortran Venus VenGram model 
 * \param i : integer, should be 0
 * \param chgt : the altitude for the computation
 * \param clat : the latitude for the computation
 * \param clon : the longitude for the computation
 */ 

	extern void datastep_v05_(int* i, // Position: used if there is a change. 0 here
			double* chgt, // altitude
			double* clat, // latitude
			double* clon, // longitude
			double* csec, // variation in time: put 0
			double* day0, // the julian day - use the function
			double* RH0d, // For the perturbation of the density : put 0
			double* RHOu, // For the perturbation of the density : put 0
			double* RHOv, // For the perturbation of the density : put 0
			int* eof, // Put 0, is 1 if a problems appears with the files
			double* DELHGT, // Modification in altitude : put 0
			double* DELLAT,// Modification in latitude : put 0
			double* DELLON,// Modification in longitude : put 0
			double* DELTIME,// Modification in time : put 0
			double* TEMP,  // The temperature OUTPUT
			double* PRES, // The pressure OUTPUT
			double* DENSLO, // OUTPUT
			double* DENS,// The density OUTPUT
			double* DENSHI,// OUTPUT
			double* DENSP, // OUTPUT
			double* EWWIND, // OUTPUT
			double* EWpert, // OUTPUT
			double* NSWIND, //OUTPUT
			double* NSpert, // OUTPUT
			double* Hrho, // OUTPUT
			double* HSCALE, // OUTPUT
			double* dsunlat, // Forces the latitude of the Sun if dradau>0
			double* dsunlon, // Forces the longitude of the Sun if dradau>0
			double* dsunLs, // Forces the Ls of the Sun if dradau>0
			double* dradau, // Forces the distance to the Sun, and the 3 previous parameters
			double* dowlt, // Venus-earth one way light time (minutes) if dradau>0
			int* LonEast, // 1 for east longitudes positives
			double* corlim, //OUTPUT
			double* DENSTOT, //OUTPUT
			int* IERT, // 1 for time input as Earth-receive time, or 0 (recommended) Venus-event time
			int* IUTC, // 1 for time input as UTC, or 0 for Terrestrial dynamical time
			double* fmolCO2, // Molar fraction CO2 OUTPUT
			double* fmolN2,// Molar fraction N2 OUTPUT
			double* fmolO,// Molar fraction O OUTPUT
			double* fmolCO,// Molar fraction CO OUTPUT
			double* fmolHe,// Molar fraction He OUTPUT
			double* fmolN,// Molar fraction N OUTPUT
			double* fmolH,// Molar fraction H OUTPUT
			double* AMz, // OUTPUT
			double* pertstep, // perturbation step (0)
			double* corlmin, // minimum relative step size for perturbation (0)
			int* iupdate, // Storage for warnings
			double* ALS, // OUTPUT
			double* SZA, // OUTPUT
			double* owlt, //OUTPUT 
			double* sunlat, // OUTPUT
			double* sunlon, // OUTPUT
			double* VenusAU, // OUTPUT
			double* TLOCAL, // OUTPUT
			double* profnear, // Lat-lon radius within which weight for auxiliary profile is 1.0
			double* proffar, // Lat-lon radius beyond which weight for auxiliary profile is 0.0
			int* nprof, // Number of profiles (1)
			char* DATADIR); // Array of 300 char for the directory where the data are!

			


}

/**
 * \ingroup Planetes
 * To work with the planet Venus
 *
 */
class Venus : public Planete
{
	protected:
	
		/** Launches the VenGram model 
		 * \param vAltGridKm : the altitude grid
		 * \param rAtmocm_3 : reference to the atmosphere (output), species with density in cm_3
		 * \param rTempK : reference to the temperature of the atmosphere (output), in K 
		 * */
		void LaunchVenGram(ublas::vector<double> vAltGridKm, ublas::matrix<double>& rAtmocm_3 , ublas::vector<double>& rTempK);



		/**
		 *
		 * Scanned curves for the dayside Ne density
		 * (Fox & Sung 2001 -> used to compute the density at low altitude)
		 * \param vAltGridKm : the altitude grid
		 * \param vF107 : the F107 parameter.
		 * \return The ne densites in cm-3
		 *
		 */
		ublas::vector<double> NeDensityScannedModel(ublas::vector<double> vAltGridKm, double vF107);

		/**
		 *
		 * Scanned curves for the Ionic Temperature (dayside)
		 * (Fox & Sung 2001 -> used to compute the temperature at low altitude)
		 * \param vAltGridKm : the altitude grid
		 * \param vF107 : the F107 parameter.
		 * \return The ion temperature in K
		 *
		 */

		ublas::vector<double> IonKTempScannedModel(ublas::vector<double> vAltGridKm,double vF107);


		/**
		 *
		 * Scanned curves for the Electron Temperature (dayside)
		 * (Fox & Sung 2001 -> used to compute the density at low altitude)
		 * \param vAltGridKm : the altitude grid
		 * \return The electron temperature in K
		 *
		 */


		ublas::vector<double> ElectronKTempScannedModel(ublas::vector<double> vAltGridKm);
		/**
		 * Modele Theis JGR 1980 1984
		 * \param vAltKm the altitude in km
		 * \param vSZADegree the SZA in degree
		 * \return the electron density (cm-3) at the altitude for the SZA
		 * this function is valid above 150 km.
		 * Therefore, you need another model below.
		 */ 
		double FsModNecm_3(double vAltKm,double vSZADegree);
	
		/**
		 * Modele Theis JGR 1980 1984
		 * \param vAltKm the altitude in km
		 * \param vSZADegree the SZA in degree
		 * \return the electron density (cm-3) at the altitude for the SZA
		 * this function is valid above 150 km.
		 * Therefore, you need another model below.
		 */ 
		double FsModTeK(double vAltKm,double vSZADegree);
		/**
		 * Option 1 for the electron temperature model : theis and al 80 84
		 * with Fox 2001 (corrected) for the lower altitudes.
		 * \param vAltGridKm : the altitude grid
		 * \param vSZADegree : the solar zenith angle
		 * \return the electron temperature in K
		 */
		ublas::vector<double> TheisModelTeK(ublas::vector<double> vAltGridKm,double vSZADegree);



	public:

		/**
		 * Option 1 for the electron density model : theis and al 80 84
		 * with Fox 2001 (corrected) for the lower altitudes.
		 * \param vAltGridKm : the altitude grid
		 * \param vSZADegree : the solar zenith angle
		 * \return the electron density in cm-3
		 * 
		 * This function is public only to unit testing.
		 * I know this is not the best, but I had to do these testings.
		 *
		 */
		ublas::vector<double> TheisModelNecm_3(ublas::vector<double> vAltGridKm,double vSZADegree);


		Venus(XmlParameters* pParam);
		Venus(XmlParameters* pParam,double vUA);
		~Venus();
		std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType);
		ublas::vector<double> ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm);
};










#endif
