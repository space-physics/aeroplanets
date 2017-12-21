/**
 * \file neutralatmo.hpp
 * \brief Defines the neutral atmo classes
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: neutralatmo.hpp 1109 2010-08-06 00:04:38Z gronoff $
 */

#ifndef NEUTRALATMO_HPP
#define NEUTRALATMO_HPP

#include <species/species.hpp>
#include <math/mathfunction.hpp>
#include <planet/allplanets.hpp> // allows to use the planet parameters for the atmosphere.
#include <chem/chem.hpp>


/**
 * Reads the neutral atmosphere data and initializes species 
 * and their cross sections
 * \ingroup atmospheres_models
 */
class NeutralAtmo
{
	protected:
		/// The parameters object
		XmlParameters* mpParameter;
		/// The altitude grid
		ublas::vector<double> mAltGridKm;

		/// The considered planet
		Planete* mpPlanet;

		/** Initializes the species vector, read the file
		 * Read the atmosphere
		 *       
		 */
		void ReadParameters();

		/// The proton grid
		ublas::vector<double>* mpPhotonGrideV;
		/// The electron grid
		ublas::vector<double>* mpElectronGrideV;
		/// The electron grid width
		ublas::vector<double>* mpElectronGridWidtheV;
		/// The photon grid
		ublas::vector<double>* mpProtonGrideV;


	public:

		/**
		 * Initializes the neutral atmo class
		 * \param pParam	 : the xmlparameters, to read the atmospheric definition
		 * \param vAltGridKm     : the altitude grid
		 * \param pPlanet	 : the planete object. Allows to load planet specific functions as temperature, atmosphere models, physical properties...
		 * \param pPhGrideV	 : photon grid. These grid are considered for the interpolation of the cross sections.
		 * \param pElGrideV	 : electron grid
		 * \param pElGridEngddeV : the electron grid width
		 * \param pPrGrideV	 : proton grid
		 *
		 */
		NeutralAtmo(XmlParameters* pParam,ublas::vector<double> vAltGridKm,Planete* pPlanet,ublas::vector<double> *pPhGrideV,ublas::vector<double> *pElGrideV, ublas::vector<double>* pElGridEngddeV,ublas::vector<double> *pPrGrideV);
		/**
		 * The destructor
		 */
		~NeutralAtmo();
		/// Vector containing the atmospheric species used for computation.
		std::deque<Specie*> mAtmoSpecies;

		/// Returns  a vector containing the neutral temperature model
		ublas::vector<double> Temperature();


		/// Prints the neutral atmosphere in a file \param vFilename the file to write
		void PrintNeutralAtmo(std::string vFilename);

		/// Prints the neutral atmosphere column densities in a file \param vFilename the file to write
		void PrintNeutralColdens(std::string vFilename);


		/*
		 * Print the bent neutral atmosphere (in parameter)
		 * \param vAtmo : the neutral atmosphere to print
		 * \param vKm : the length in Km
		 * \param vFilename : the filename
		 */
		void PrintBentNeutralAtmo(std::deque<Specie*> vAtmo,ublas::vector<double> vKm,std::string vFilename);

		/*
		 * Print the bent neutral column density (in parameter)
		 * \param vAtmo : the neutral atmosphere to print
		 * \param vKm : the length in Km
		 * \param vFilename : the filename
		 */
		void PrintBentNeutralColdens(std::deque<Specie*> vAtmo,ublas::vector<double> vKm,std::string vFilename);

		/**
		 * Computes the column density for the different species.
		 * This function suppose a vertical atmosphere
		 * that was the case for the old Trans*, but the objective
		 * of this new version is also to have much more complex
		 * geometries. One of the point is that we can approximate
		 * a very complex magnetic field geometry for electron
		 * by creating a vertical atmosphere which is the atmosphere 
		 * along the line. It is not compatible with the photoionization
		 * but it works (Gronoff et al 2009 titan ionization II)
		 *
		 * In that case, this function needs to be modified.
		 *
		 * \warning this function needs to be called after the neutral temperature
		 * initialization -> it depends on it
		 * \param vNeutralT -> the neutral temperature vector
		 */

		void ComputeVerticalColumnDensity(const ublas::vector<double>& vNeutralT);


		/**
		 * Clear the production of the different species
		 * important task to perform between two 
		 * ionization-excitation processes!
		 */

		void ClearProductions();



		/**
		 * Modifies the density of a species, and also the column density (only for non-bent).
		 * All of this based on a model of density variation
		 * \param vName : the name of the species to modify (not excited state yet)
		 * \param vModel : the model of density (gaussian, chapman,...)
		 * \param vParams : the necessary parameters for the model
	 	* \param vAddParams : additionary parameters (non variable for adjustement)
		 * \param vNeutralT -> the neutral temperature vector
		 */
		void ResetSpecieDens(std::string vName,int vModel, std::deque<double> vParams, std::deque<double> vAddParams, const ublas::vector<double>& vNeutralT );

		/**
		 * Modifies the density of a species, based on data
		 * \param vName : the name of the species to modify (not excited state yet)
		 * \param vDenscm_3 : the values of the density at the altitudes of the grid
		 * \param vNeutralT -> the neutral temperature vector
		 */
		void ResetSpecieDensInterp(std::string vName, ublas::vector<double> vDenscm_3,const ublas::vector<double>& vNeutralT);




};





#endif
