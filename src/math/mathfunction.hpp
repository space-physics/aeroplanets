/**
 * \file mathfunction.hpp
 * \brief Defines basic functions, like integration...
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathfunction.hpp 1309 2011-09-14 16:29:09Z gronoff $
 */



#ifndef MATHFUNCTION_HPP
#define MATHFUNCTION_HPP
#include "mathstring.hpp"
#include "mathflux.hpp"
#include "constphys.hpp"
#include "umathfunction.hpp"



/**
 * To compute physical/time properties
 * \ingroup Math
 */

namespace PhysTime
{
	/**
	 * Computes the julian day (from subroutine in VenGram) by method of Meeus, Astronomical
C       Algorithms, 2nd Edition, 1998, page 61
	 *
	 * \param vYear : the year
	 * \param vMonth : the month
	 * \param vDay : the day
	 * \param vHour : the hour
	 * \param vMin : the minute
	 * \param vSec : the seconds
	 *
	 */
	double CalToJul(int vYear, int vMonth, int vDay, int vHour, int vMin, double vSec);
	/**
	 * Computes the julian day, whith UT 
	 * \param vYear : the year
	 * \param vMonth : the month
	 * \param vDay : the day
	 * \param vHour : the hour
	 *
	 */
	double CalToJul(int vYear, int vMonth, int vDay, double vHour);


};


/**
 * To compute grids
 * \ingroup Math
 *
 */
namespace MathGrid
{
	/**
		Permet de calculer une distribution exponentielle de nb_pts entre min et max, alias de exponential_distrib.
		Correspond a une grille initiale.


		\param  vMin: la hauteur minimale
		\param  vMax: la hauteur maximale
		\param vNbPts : le nombre de points
		\return vector<double> : la distribution

	*/
	ublas::vector<double> GridExp(double vMin, double vMax,int vNbPts);


	/** gridpow : defines a power law decrease grid (dirk89,lilen...)
	  it can work with either begin>end or end>begin
	  \param begin : the start value of the grid
	  \param end : the end value of the grid
	  \param nb_pts : the size of the grid
	  \return the grid
	  */
	ublas::vector<double> GridPow(double begin,double end,int nb_pts);
	
	

	/** 
	 * \~english gridcst : defines a constant width grid (dirk89,lilen...)
	 * it can work with either begin>end or end>begin
	 * \warning, in this implementation, begin is at the start!!! 
	 * It is inverted compared to gridcst fortran. If you invert
	 * begin and end, the grids are inverted!
	 * \param begin : the start value of the grid
	 * \param end : the end value of the grid
	 * \param nb_pts : the size of the grid
	 * \return the grid
	 *
	 *
	 * \~francais \warning: dans cette implementation, begin est le debut du tableau
	 * de sortie. La grille est calculée comme dans le prog fortran, mais 
	 * elle est ensuite inversee. 
	 * Les grilles crées en inversant begin et end sont mirroir.
	 */
	ublas::vector<double> GridCst(double begin,double end,int nb_pts);
	
	
	//gridpolo

	/**
	entab : subroutine for gridpolo to avoid overflow...
		\param ntab : the size of the table
		\param spfac : the growth
		\param tmin : ?
	  */
	double EnTab(int ntab,double spfac,double tmin);
	/**
		gridpolo:power law decreasing grid.
		\warning: the min value SHOULD be inferior.
		the values should not be negative
		the tmax/tmin should be superior to 20!
		\param ntab : the size of the grid
		\param tmin : the min value
		\param tmax : the max value
		\param tab : the output grid, given as a reference.
		\param widthtab : the output widthgrid, given as a reference.
		\param spfac: the outpout growth parameter, given as a reference
		\return false if an error occured in the computation.
	*/
	bool GridPolo(int ntab,double tmin,double tmax,ublas::vector<double>& tab,ublas::vector<double>&widthtab,double&spfac);

	// widthexpgrid
	//widthnonexpgrid
	// widthgrid

	/**
	 * Computes the width of a grid.
	 * If the grid has an exp type, the type should be 1
	 * else, it is computed as a power decrease.
	 * \warning the size of the grid should be superior to 3
	 * \param tab : the grid
	 * \param itype : the type of the grid (1 exp)
	 * \return the widthgrid
	 */
	ublas::vector<double> WidthGrid(ublas::vector<double> tab,int itype);
	
	/**
	 * Computes the width of an exponential grid
	 * \param ttab : the grid
	 * \return the widthgrid
	 */
	ublas::vector<double> WidthExpGrid(ublas::vector<double> ttab);
	/**
	 * Computes the width of a power law decrease grid
	 * \param ttab : the grid
	 * \return the widthgrid
	 */
	ublas::vector<double> WidthNonExpGrid(ublas::vector<double> ttab);

	/**
	 * Interpolate a grid on another one. Both grids have data per box (hence a simple interpolation is not possible)
	 * Redistibute the flux (vFlux1) given for a grid, defined by its upper and lower values (vGridMax1, vGridMin1)
	 * on the second grid (vFlux2,vGridMax2,vGridMin2).
	 * Hypothesis : The grid are increasing or decreasing perfectly (the interpolation is possible on the average value for the grid).
	 *
	 * \warning This function is correct if we consider that it is a number of things inside the boundaries that we want to redistribute:
	 * e.g. 10 photons/cm3/s inside the 10-15eV box.
	 * 
	 * \warning A null width boxe is a problem.
	 *
	 * \param vFlux1 : the flux to redistribute
	 * \param vGridMin1 : the vector of the minimum values of the flux1 grid
	 * \param vGridMax1 : the vector of the maximum values of the flux1 grid
	 * \param vGridMin2 : the vector of the minimum values of the output grid
	 * \param vGridMax2 : the vector of the maximum values of the output grid
	 * \param rSupplement : reference to a value of the flux not taken into account.
	 * \return The redistributed grid.
	 *
	 *
	 */
	ublas::vector<double> InterpolateNonChaoticFlux(const ublas::vector<double>& vFlux1,const ublas::vector<double>& vGridMin1,const ublas::vector<double>& vGridMax1,
			const ublas::vector<double>& vGridMin2,const ublas::vector<double>& vGridMax2,double &rSupplement);



	/**
	 * Redistribution of a grid on another one.
	 * Redistibute the flux (vFlux1) given for a grid, defined by its upper and lower values (vGridMax1, vGridMin1)
	 * on the second grid (vFlux2,vGridMax2,vGridMin2).
	 * Hypothesis : the second grid must be increasing perfectly (ie can be given by one vector, ie vGridMax(i)==vGridMin(1+i), or decreasing perfectly.
	 *              but the first one can be a s... one, with lines added... 
	 *              The hypothesis on the first one is that gridmin[i] < gridmax [i] 
	 * \warning This is not an efficient function, because, consisidering the possibility of a first grid chaotic, I a to
	 * compute the things a lot
	 *
	 * \warning This function is correct if we consider that it is a number of things inside the boundaries that we want to redistribute:
	 * e.g. 10 photons/cm3/s inside the 10-15eV box.
	 * 
	 * \warning A null width boxe is a problem.
	 *
	 * \param vFlux1 : the flux to redistribute
	 * \param vGridMin1 : the vector of the minimum values of the flux1 grid
	 * \param vGridMax1 : the vector of the maximum values of the flux1 grid
	 * \param vGridMin2 : the vector of the minimum values of the output grid
	 * \param vGridMax2 : the vector of the maximum values of the output grid
	 * \param rSupplement : reference to a value of the flux not taken into account.
	 * \return The redistributed grid.
	 *
	 *
	 */
	ublas::vector<double> RedistributeChaoticFlux(const ublas::vector<double>& vFlux1,const ublas::vector<double>& vGridMin1,const ublas::vector<double>& vGridMax1,
			const ublas::vector<double>& vGridMin2,const ublas::vector<double>& vGridMax2,double &rSupplement);

	/**
	 * Add lines in the grid : the lines are defined as small portions which, in the case of a flux, can be very strong: it is important for the 
	 * accuracy of photoionization models that lines are taken into account as separate points.
	 *
	 * There is no hypothesis on the width of the line for this function (therefore, we can add something larger than a line). More precisely
	 * if the line is on a gap, it will be splited eg : A | B  => A' | L1 | L2 | B'
	 * \param vNewGridMin : the lower value of the lines grid
	 * \param vNewGridMax : the upper value of the lines grid
	 * \param rGridMin : the lower values of the old grid. This parameter is UPDATED (reference)
	 * \param rGridMax : the upper values of the old grid. This parameter is UPDATED (reference)
	 *
	 *
	 * The created grids are increasing or decreasing perfectly if it is the case for rGridMin, rGridMax, else, it is chaotic.
	 *
	 */

	void AddGridLines(const ublas::vector<double>& vNewGridMin, const ublas::vector<double>& vNewGridMax, ublas::vector<double>& rGridMin, ublas::vector<double>& rGridMax);

/**
 * Function needed by AddGridLines : insert cut, with the value vCut inside the grids
 * \param vCut : the cut
 *
 *  \param rGridMin : the lower values of the old grid. This parameter is UPDATED (reference)
 * \param rGridMax : the upper values of the old grid. This parameter is UPDATED (reference)
 */
	void InsertGridCut(const double & vCut, ublas::vector<double>& rGridMin,ublas::vector<double>&rGridMax);

	/**
	 * implemets the equivalent of the std::vector  insert function -> very slow!!!
	 *
	 */
	template<class T> ublas::vector<T> UblasInsert(const ublas::vector<T>& vVec,size_t position,T vVal)
	{
		ublas::vector<T> resu(vVec.size()+1);
		if(position!=0)
			std::copy(vVec.begin(),vVec.begin()+position,resu.begin());
		if(position!=vVec.size())
			std::copy(vVec.begin()+position,vVec.end(),resu.begin()+position+1);
		resu[position]=vVal;
		return resu;
	}

	/**
	 * Allows to redistribute an electron distro depending on altitude and angle
	 * (Energy, Altitude, Angle) on a new energy grid.
	 * This function does not divide the distribution by the grid width, so it cannot
	 * use a linear interpolation.
	 * It returns the loss of energy in the redistribution, depending on the altitude
	 * \param vpGangle : the Gaussian angle corresponding to the angle distribution in the flux
	 * \param vBotE1 : the bottom energy grid for the old layer
	 * \param vTopE1 : the top energy grid for the old layer
	 * \param vOldElecDistro : the old energy distribution (energy)(altitude,angle)
	 * \param vBotE2  : the bottom energy grid for the old layer
	 * \param vTopE2 : the top energy grid for the old layer
	 * \param vNewElecDistro : the new energy distribution (energy)(altitude,angle)
		
	 *
	 */
	ublas::vector<double> RedistributeEAAFlux(MathFunction::GaussianAngle* vpGangle,
			ublas::vector<double> vBotE1,
			ublas::vector<double> vTopE1,
			ublas::vector< ublas::matrix<double> > vOldElecDistro,
			ublas::vector<double> vBotE2,
			ublas::vector<double> vTopE2,
			ublas::vector< ublas::matrix<double> >& vNewElecDistro
			);


};




#endif
