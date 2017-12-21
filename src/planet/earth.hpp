/**
 * \file earth.hpp
 * \brief Defines the earth planet
 * Copyright G Gronoff Nov 2009
 * Last Modifcation $Id: earth.hpp 1356 2011-11-16 21:20:46Z gronoff $
 */

#ifndef EARTH_PLANET_HPP
#define EARTH_PLANET_HPP
#include "Earth/irimsis.hpp"
/**
 * \ingroup Planetes
 * To work with the earth...
 * This is not an original thing...
 */
class Earth : public Planete
{
	protected:
		bool mbIsModelLoaded;
		EarthIriMsis* mEarth;
		void LoadEarth();

	public:
		Earth(XmlParameters* pParam);
		Earth(XmlParameters* pParam,double vUA);
		~Earth();
		std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType);
		ublas::vector<double> ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm);
		
};




#endif
