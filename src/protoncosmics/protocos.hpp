/**
 * \file protocos.hpp
 * \brief Definition for the proton-cosmics impact.
 * The effects of protons and cosmics rays is computed here with
 * planetocosmics. We use the output of the code averaged here.
 * Copyright G Gronoff Dec 2009
 * Last Modification : $Id: protocos.hpp 1342 2011-11-09 22:39:36Z gronoff $
 */

#ifndef PROTO_COS_HPP
#define PROTO_COS_HPP
#include <species/species.hpp>
#include <eflux/eflux.hpp>


/**
 *
 * \ingroup ionization_processes
 * Allows to read the output of planetocosmics (modification G Gronoff
 * for previous work on Titan (Gronoff 2009). To be included in 
 * the ionization processes.
 * 
 */
class ProtoCosIonization
{
	private:
		/// The Xmlparameter object
		XmlParameters* mpParams;
		/// The altitude grid
		ublas::vector<double>* mpAltGridKm;

		/// The center of the elec grid
		ublas::vector<double>* mpElecC;
		/// The bottom of the elec grid
		ublas::vector<double>* mpElecB;
		/// The width of the elec grid
		ublas::vector<double>* mpElecD;

		/// The number of energies
		unsigned mNben;
		/// The number of angles
		unsigned mNbang;
		/// The number of altitudes
		unsigned mNbalt;
		/// Direct link to the Gaussian angle. It is not defined in that function
		MathFunction::GaussianAngle* mpGAngle; 


		/**
		 * Compute the ionization and the resulting flux (does not overload the result flux but add : the computation can be added several times
		 */
		void CosmoIonize(TiXmlNode* node,std::deque<Specie*>& rResult,EFlux& rResultFlux,std::deque<Specie*>& vpSp);

		/**
		 * Reads the electron flux files, fills the eflux file (by performing interpolations!
		 * and returns the rest of the energy at each altitudes
		 */
		ublas::vector<double> CosmoElecFlux(std::string vEFluxFile,EFlux& rResultFlux,std::deque<Specie*>& vpSp);

	public:
		/** The constructor: inits the xmlparameters object
		 * \param pParam : the object parameter
		 * \param vAltGridKm : the altitude grid (pointer)
		 */
		ProtoCosIonization(XmlParameters* pParam, ublas::vector<double>* vAltGridKm);
		/// The destructor
		virtual ~ProtoCosIonization();
		/**
		 * Computes the cosmic rays impact  with the densities - column densities
		 * embedded in the species files 
		 *
		 * \param vpSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section
		 * \param rResult : the result vector of species -> were the computations are written
		 * \param rResultFlux : flux of secondary electron, needed to compute the secondary production.
		 *
		 *
		 */
		void ComputeCosmoionization(std::deque<Specie*>& vpSp,
					     std::deque<Specie*>& rResult,
					     EFlux& rResultFlux);


};



#endif
