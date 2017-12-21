/**
 * \file protoncrs.hpp
 * \brief Defines the proton cross section class
 * Copyright G Gronoff November 2011
 * Last Modification : $Id$
 */


/**
 * As for the class ElecCrossSection, we inherit from Cross Section, 
 * and we add more data to account for the physics:
 * - the presence of proton and hydrogen cross sections
 * - the elastic, inelastic, ionizations
 * - the charge exchange
 */

#ifndef PROTON_CRS_HPP
#define PROTON_CRS_HPP

#include "crs.hpp"
#include "documentation.hpp"


class ProtCrossSection:public CrossSection
{
	private:
		/// Computes the elastic cross section for Hydrogen, following Kozelov and Ivanov 92. As explained in the Fortran code
		void HydrogenElasticFunction();
		/// Computes the elastic cross section for Protons, following Kozelov and Ivanov 92. As explained in the Fortran code
		void ProtonElasticFunction();
	protected:
		/// The maximum defined energy in the elastic cross section (computed by ExtractCrs or ExtractCrsLog)
		double mMaxDefinedEnergyElasticCrseV;
	public:
		/// The total ionization cross section
		ublas::vector<double> mTotIonizationCrscm2;
		/// The total excitation cross section
		ublas::vector<double> mTotExcitationCrscm2;
		/// The total exchange cross section
		ublas::vector<double> mTotExchangeCrscm2;

		/// The name of the impactor (proton or hydrogen typically)
		std::string mParticleName;


		/// The  Ionization energy loss: vioni  -- Wio
		ublas::vector<double> mElossIoni;

		/// The minimum ionization energy
		double mVion;


		/// The  electron pulling energy loss: varra -- W01
		ublas::vector<double> mElossHioni;
		/// The elastic energy loss: velas_(species) ----  Welas 
		ublas::vector< ublas::matrix<double> > mElossElas;

		/// True if mElossElas is computed
		bool mbIsElossElasComputed;


		/// The excitation energy -- Wex
		double mWexcitation;
		/// The capture energy (only for protons) -- W10
		double mW10;

		

	public:
		/**
		 * Initialize the Grid and the cross section
		 * \param vGreV: the proton energy grid in eV
		 */
		ProtCrossSection(ublas::vector<double> vGreV, std::string vParticleName);
		/**
		 * Load the Cross section.
		 *
		 * All of that from the given object parameter.
		 * \param pParams : object to read the cross section file
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : true if MonteCarlo active
		 * \param bLogInterp : true if you want log interpolation, the default HERE
		 * \return true if loaded.
		 *
		 * When this function has been called, all the cross section are loaded.
		 * 
		 */
		bool LoadCrs(XmlParameters *pParams,std::string vName,bool bIsMonteCarlo,bool bLogInterp=true);

		/**
		 * Load the Cross section.
		 * All of that from the given object file -> calls LoadCrs(pParams afterwards)
		 *
		 * \param vFilename:  the filename for the crs
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : true if MonteCarlo active
		 * \param bLogInterp : true if you want log interpolation, the default HERE
		 * \return true if loaded.
		 */
		bool LoadCrs(std::string vFilename,std::string vName,bool bIsMonteCarlo,bool bLogInterp=true);

		/**
		 * Print the cross section with the new energy grid
		 * into a file. This file can be read with a python-pylab load function
		 * \param vFilename : name of the file. It is created.
		 * \param vRedisFilename : name of the file for the redistribution function (optional)
		 */

		void PrintCrs(std::string vFilename,std::string vRedisFilename="");

		/**
		 * Compute the elastic energy loss, Welas in the fortran code
		 * \param vSpeciesMass : the neutral species mass
		 * \param vGa : the Gaussian Angle object for the proton transport (computed while initializing the proton transport system)
		 * \param vbIsRedis : True if the redistribution in angle system is used (should be the default).
		 * \param vbForce : false by default. To prevent the costly computation of the elastic loss at each call, the result is protected by the variable mbIsElossElasComputed. Anyway, if it is necessary to force the computation, the vbForce can be set to true.
		 */
		void ComputeElasticLoss(const double& vSpeciesMass,const MathFunction::GaussianAngle& vGa, bool vbIsRedis, bool vbForce = false);

};




#endif 

