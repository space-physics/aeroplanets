/**
 *
 * \file datamarinerviking.hpp
 * \brief Defines the martian atmosphere during mariner and viking conditions
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: datamarinerviking.hpp 831 2010-01-08 17:58:49Z gronoff $
 */

#ifndef DATAMARINERVIKING_HPP

#define DATAMARINERVIKING_HPP
#include "marsatmotim.hpp"

namespace MarsData
{
	/*             CONDITIONS MARINER  iflag=2            
		       Alt    |    n_CO2     |      n_N2     |      n_CO     |      n_O      |     n_O2      |      n_H      |      Tn      |
		       (km)    |   (cm-3)     |     (cm-3)    |     (cm-3)    |     (cm-3)    |    (cm-3)     |     (cm-3)     |      (K)    |#endif
		       ---------------------------------------------------------------------------------------------------------------------------
		       */
	/**
	
	  Returns a matrix giving the Viking1 conditions (altitudes, CO2, N2,CO, O, O2,  H, Tn
Viking 1 : LS=100
July 20 08:51 UT. 1976 
f 107Mars= 28
(f 107Earth= 72), the Sun–Mars distance d is 1.61 AU, and the
solar zenith angle  reaches 44
	*/
	ublas::matrix<double> MarinerConditions();
	/**
	
	  Returns a matrix giving the Viking conditions (altitudes, CO2, N2,CO, O, O2,  H, Tn

	  Mariner 6 :  LS=200  M
	TGCM (Bougher et al., 1988, 1990, 1999, 2000, 2006) calculated
	in the conditions of Mariner 6 (high solar activity f107mars=98   F107 Earth=200, d = 1.45 AU, solar zenith angle 57 )
	*/
	ublas::matrix<double> VikingConditions();



	/**
	 * Return the Mariner neutral atmospherer grid in function of the altitude given in parameter.
	 * Performs a log interpolation
	 * \param vAltitudeGridKm : the final altitude grid.
	 *
	 */
	std::deque< ublas::vector<double> > MarinerNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the Mariner temperature grid in function of the altitude given in parameter.
	 * Performs a log interpolation
	 * \param vAltitudeGridKm : the final altitude grid.
	 *
	 */
	ublas::vector<double> MarinerNTemp(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the Viking neutral atmospherer grid in function of the altitude given in parameter.
	 * Performs a log interpolation
	 * \param vAltitudeGridKm : the final altitude grid.
	 *
	 */
	std::deque< ublas::vector<double> > VikingNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);

	/**
	 * Return the Mariner temperature grid in function of the altitude given in parameter.
	 * Performs a log interpolation
	 * \param vAltitudeGridKm : the final altitude grid.
	 *
	 */
	ublas::vector<double>  VikingNTemp(const ublas::vector<double>& vAltitudeGridKm);











}


#endif
