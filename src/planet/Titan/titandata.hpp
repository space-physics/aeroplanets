/**
 * \file titandata.hpp
 * \brief Defines the standards atmospheres for Titan (See Gronoff 2009)
 * Copyright G Gronoff Nov 2009
 * Last modification $Id: titandata.hpp 773 2009-12-04 14:43:05Z gronoff $
 * \author G Gronoff
 *
 */

#ifndef TITAN_DATA_HPP
#define TITAN_DATA_HPP
#include "../planet.hpp"

namespace TitanData
{

	/**
	 * Return the Wodarg 2000 electron density
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	ublas::vector<double> WodargDensE(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the Wodarg 2000 electron temperature
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	ublas::vector<double> WodargTempE(const ublas::vector<double>& vAltitudeGridKm);
	

	/**
	 * Return the Wodarg 2000 ion temperature
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	ublas::vector<double> WodargTempI(const ublas::vector<double>& vAltitudeGridKm);
	

	/**
	 * Return the neutral atmosphere (N2,CH4) of wodarg 2000
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	std::deque< ublas::vector<double> > WodargNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);


	/**
	 * Return the neutral temperature for the Wodarg 2000 (ok, it is mainly 175K
	 * \param vAltitudeGridKm : the altitude grid
	 */
	ublas::vector<double>  WodargTn(const ublas::vector<double>& vAltitudeGridKm);


	/**
	 * Return the neutral atmosphere (N2,CH4) of waite 2005
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	std::deque< ublas::vector<double> > WaiteNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the neutral temperature for the Wodarg 2000 (ok, it is mainly 175K
	 * \param vAltitudeGridKm : the altitude grid
	 */
	ublas::vector<double>  WaiteTn(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the neutral atmosphere (N2,CH4) of Cui 2009
	 * \param vAltitudeGridKm : the altitude Grid
	 */
	std::deque< ublas::vector<double> > CuiNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);


	/**
	 * Return the neutral temperature for the Wodarg 2000 (ok, it is mainly 175K
	 * \param vAltitudeGridKm : the altitude grid
	 */
	ublas::vector<double>  CuiTn(const ublas::vector<double>& vAltitudeGridKm);


};




#endif
