/**
 * \file mathflux.hpp
 * \brief Defines the mathematical function to initialize fluxes
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathflux.hpp 1679 2013-02-15 23:31:18Z gronoff $
 *
 */


#ifndef MATHFLUX_HPP
#define MATHFLUX_HPP
#include "mathstring.hpp"


/**
 * \ingroup Math
 * Namespace to work with the fluxes of particle.
 * It allows to put the fluxes inside energy grids
 * as the one computed with the MathGrid functions
 */
namespace MathFlux
{
	/**
	 * Computes a dirac energy distribution for the vCentEeV grid, with vNang angles. 
	 * The results are stored in fxdown and fxup parameters. fxup[energy][angle]
	 * \param vEnTotErg : total energy of the distribution
	 * \param vEzeroeV : energy  where the dirac is found
	 * \param vPowLaw : if true  add a power law centered at 10, with a total energy of 6.25E11*vEnTotErg
	 * \param vNang : the number of angles
	 * \param vCentEeV : the energy grid
	 * \param vIsotro :  2 -> forward 0 -> non vIsotrope 1 -> vIsotrope. Be careful : the flux is by cm3 by sr by s, that explain why there is no dependance on angle on the isotropic case
	 * \param vXmu : gaussian abscissa of the angles (computed with the MathFun::GaussianAngle)
	 * \param rFxDown : OUTPUT the downward energy distribution 
	 * \param rFxUp : OUTPUT the upward energy distribution 
	 */
	void InDirac(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,const ublas::vector<double>& vXmu,   ublas::matrix<double>  & rFxDown, ublas::matrix<double> &rFxUp);
	
	
	/**
	 * Computes a gaussian energy distribution for the vCentEeV grid, with vNang angles. 
	 * The results are stored in fxdown and fxup parameters. fxup[energy][angle]
	 * \param vEnTotErg : total energy of the distribution
	 * \param vEzeroeV : energy  where the dirac is found
	 * \param vPowLaw : if true  add a power law centered at 10, with a total energy of 6.25E11*vEnTotErg
	 * \param vNang : the number of angles
	 * \param vCentEeV : the energy grid
	 * \param vIsotro :  2 -> forward 0 -> non vIsotrope 1 -> vIsotrope . Be careful : the flux is by cm3 by sr by s, that explain why there is no dependance on angle on the isotropic case
	 * \param vXmu : gaussian abscissa of the angles (computed with the MathFun::GaussianAngle)
	 * \param rFxDown : OUTPUT the downward energy distribution 
	 * \param rFxUp : OUTPUT the upward energy distribution 
	 */
	
	
	void InGauss(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,const ublas::vector<double>& vXmu,            ublas::matrix<double> & rFxDown, ublas::matrix<double> & rFxUp);



	/**
	 * Computes a Maxwellian energy distribution for the vCentEeV grid, with vNang angles. 
	 * The results are stored in fxdown and fxup parameters. fxup[energy][angle]
	 * \param vEnTotErg : total energy of the distribution
	 * \param vEzeroeV : energy  where the dirac is found
	 * \param vPowLaw : if true  add a power law centered at 10, with a total energy of 6.25E11*vEnTotErg
	 * \param vNang : the number of angles
	 * \param vCentEeV : the energy grid
	 * \param vIsotro :  2 -> forward 0 -> non vIsotrope 1 -> vIsotrope . Be careful : the flux is by cm3 by sr by s, that explain why there is no dependance on angle on the isotropic case
	 * \param vXmu : gaussian abscissa of the angles (computed with the MathFun::GaussianAngle)
	 * \param rFxDown : OUTPUT the downward energy distribution 
	 * \param rFxUp : OUTPUT the upward energy distribution 
	 */
	
	void InMaxw(double vEnTotErg,double vEzeroeV,bool vPowLaw,int vNang,const ublas::vector<double>& vCentEeV,unsigned vIsotro,ublas::vector<double> vXmu,           ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp);



	/**
	 * Normalize the energy in precipitation grid
	 * \param vEnTotErg : the total energy
	 * \param rFxDown : the downward energy flux
	 * \param rFxUp : the upward energy flux
	 * \param vWeight : the gaussian weight
	 * \param vXmu : the gaussian abscissa
	 * \param vEeV : the energy grid
	 * \param vDdengeV : the energy grid width
	 */
	void NormFlux(double vEnTotErg, ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp,ublas::vector<double> vWeight,ublas::vector<double> vXmu,ublas::vector<double> vEeV,ublas::vector<double> vDdengeV);



	/**
	 * Fills the matrix (ene, ang) up and down with a flux per steradian (isotropic)
	 * \param rFxDown : the downward energy flux
	 * \param rFxUp : the upward energy flux
	 * \param vFx : the energy flux to fill with
	 */

	void FillFlux(ublas::matrix<double>& rFxDown, ublas::matrix<double>& rFxUp, ublas::vector<double> vFx);


};





#endif
