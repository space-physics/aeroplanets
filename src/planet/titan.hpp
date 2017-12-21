/**
 * \file titan.hpp
 * \brief Defines the Titan planet
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: titan.hpp 1356 2011-11-16 21:20:46Z gronoff $
 */
#ifndef TITAN_PLANET_HPP
#define TITAN_PLANET_HPP
#include "planet.hpp"
#include "Titan/titandata.hpp"

/**
 * \ingroup Planetes
 * To work with the satellite of Saturn Titan
 */
class Titan : public Planete
{
	protected:
	public:
		Titan(XmlParameters* pParam);
		Titan(XmlParameters* pParam,double vUA);
		~Titan();
		std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType);
		ublas::vector<double> ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm);
		
};




#endif
