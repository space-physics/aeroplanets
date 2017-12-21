/**
 * \ingroup Chem
 * \file calcdens.hpp
 * \brief Defines the implementation of the different density computation functions. Inheriting from CalcDensPhEq.
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: calcdens.hpp 1170 2010-10-27 21:53:28Z gronoff $
 */

#ifndef CALC_DENS_HPP
#define CALC_DENS_HPP
#include "reaclist.hpp"




/**
 * Computes the N2(A3S) density (responsible for the VK emissions
 */

class CalcDensN2A3S : public CalcDensPhEq
{

	public:
		/// The constructor
		CalcDensN2A3S();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};



/**
 * Computes the O1S density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 *
 */

class CalcDensO1S : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensO1S();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};



/**
 * Computes the N(2D) density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 * It is important for the computation of O(1D) for the earth.
 * \todo find C
 */

class CalcDensN2D : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensN2D();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};


/**
 * Computes the O1D density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 *
 */

class CalcDensO1D : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensO1D();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};



/**
 * Computes the COa3Pi density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 * The density was not taken into account before, so we test direct emissions and the influence of the quenching on the process...
 *
 */

class CalcDensCOa3Pi : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensCOa3Pi();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};




/**
 * Computes the O+(2P) density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 * The density was not taken into account before, so we test direct emissions and the influence of the quenching on the process...
 *
 */

class CalcDensOplus2P : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensOplus2P();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};



/**
 * Computes the O++ density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 *
 */
class CalcDensOpp : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensOpp();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};

/**
 * Computes the N2++ density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 *
 */
class CalcDensN2pp : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensN2pp();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};



/**
 * Computes the CO2++ density. 
 * For that species, the photochemical equilibrium is a very good approximation.
 *
 */
class CalcDensCO2pp : public CalcDensPhEq
{
	public:
		/// The constructor
		CalcDensCO2pp();
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */	
		 ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings);
};








#endif
