/**
 * \file mars.hpp
 * \brief Defines Mars
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: mars.hpp 1356 2011-11-16 21:20:46Z gronoff $
 */
#ifndef MARS_PLANET_HPP
#define  MARS_PLANET_HPP
#include "planet.hpp"
#include "Mars/marsbinary.hpp"
/**
 * \ingroup Planetes
 */
class Mars : public Planete
{
	protected:

		/**
		 * Retrieve the distance thanks to the 
		 * solar longitude parameter
		 * \param vLs the solar longitude
		 *
		 */
		void LsInit(double vLs);
		/// The solar longitude
		double mLsDegree;


		/// Check if the mars tim atmosphere is loaded
		bool mbIsMartimLoaded;
		/// Check if the binary file is loaded
		bool mbIsBinaryLoaded;

		/// The Martim object
		MarsAtmoTim* mMarTim;

		/// The Binary object
		DumpData* mBin;

		/// To load the Martim object
		void LoadMarTim();

		/// To load the dumpdata object
		void LoadDump();





	public:
		Mars(XmlParameters* pParam);
		Mars(XmlParameters* pParam,double vUA);
		~Mars();
		std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType);

		ublas::vector<double> ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int & vType);
		ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm);

};




#endif
