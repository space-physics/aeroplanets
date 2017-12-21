/**
 * \file anonym.hpp
 * \brief Defines the anonym planet : the user must define all the parameter of the planet
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: anonym.hpp 1356 2011-11-16 21:20:46Z gronoff $
 */

#ifndef ANONYM_PLANET_HPP
#define ANONYM_PLANET_HPP

#include "planet.hpp"


/**
 * \ingroup Planetes
 * To work with an anonymous planet: we must define all the parameters inside the xml parameter file
 */

class AnonymousPlanet : public Planete
{
	protected:

	public:
		AnonymousPlanet(XmlParameters* pParam);
		AnonymousPlanet(XmlParameters* pParam,double vUA);
		~AnonymousPlanet();

		std::map< std::string, ublas::vector<double> > AtmoModel(const ublas::vector<double>& vAltitudeGridKm,std::deque< std::string > vSpNames,int vType);
		ublas::vector<double> ElectronDensity(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ElectronTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> IonTemperature(const ublas::vector<double>& vAltGridKm,const int& vType);
		ublas::vector<double> ReturndB_B(const ublas::vector<double>& vAltGridKm);
};




#endif
