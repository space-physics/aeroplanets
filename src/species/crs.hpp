/**
 * \defgroup specie_parameters Species
 * 
 * \file crs.hpp
 * \brief Defines the cross section class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: crs.hpp 1340 2011-11-09 21:26:42Z gronoff $
 *
 */


#ifndef CRS_HPP
#define CRS_HPP
#include <cppscixml/scixml.hpp>
#include "shirai.hpp"



/**
 * \ingroup specie_parameters
 * Class to define a species plus its energy state
 * It allows to verify that two species are identical, but with 
 * different energy states...
 *
 */
class SpecieId
{
	public:
		/// The name of the species
		std::string mName;
		/// The name of the species' state
		std::string mState;

		/// Constructor that does nothing
		SpecieId()
		{
		}
	
		/** 
		 * Returns the name of the species, with its excitation
		 */
		std::string StandardName()
		{
			return mName+"("+mState+")";
		}

	
		/**
		 * Constructor that initializes name and state
		 * \param vNa : the name of the specie
		 * \param vSt : the state of the specie
		 *
		 */
		SpecieId(std::string vNa,std::string vSt)
		{
			// nb : we trim to be sure
			mName=trim(vNa);
			mState=trim(vSt);
		}


		/**
		 * Allows to compare two species.
		 * \param A :the first species
		 * \param B :the second species
		 * \return  true if the names and states are equal
		 */
		friend bool operator == (const SpecieId&A,const SpecieId&B)
		{
			return (A.mName==B.mName)&&(A.mState==B.mState);
		}

		/**
		 * Allows to fill the binary tree
		 * \param A :the first species
		 * \param B :the second species
		 * \return  true if the names and states of A < B
		 */
		friend bool operator < (const SpecieId&A,const SpecieId&B)
		{
			return (A.mName+A.mState< B.mName+B.mState);
		}

		/**
		  Check if it is the same specie (same name, the state can be different)
		 * \param A :the first species
		 * \param B :the second species
		 * \return  true if the names are equal
		  */

		bool Same(const SpecieId&A,const SpecieId&B)
		{
			return (A.mName==B.mName);
		}


};



/**
 * \ingroup specie_parameters
 * Class to load and interpolate cross sections against an energy grid.
 * Normalises the search for excited species.
 * To compute fluxes, the TotalCrs cross section should be used
 * while to compute electron, ion, and excited species production
 * Crs should be used.
 *
 */
class CrossSection
{
	protected:
		/// The xml parameter!!!
		XmlParameters* mpParam;
		/// True if the param is loaded here
		bool mParamLoaded;
		/// The cross section energy grid
		ublas::vector<double> mGrideV;
		/// The name of the specie
		std::string mSpecies;


		/// The exprapolation is in loglog if this parameter is true; false by default, but put at true when LoadCrs bLogInterp=True
		bool mExtrapolationLogLog;



		/// To check if created_species and created_species_multiplicator have the same dimensions
		bool CheckCreatedSpecies();

		/**
		 *
		 * Get the process or total crs node
		 * and returns the interpolated cross section.
		 * It verfies if the cross section is not negative after the linear interpolation.
		 * \param pParams : the xml parameters object where the node is extracted
		 * \param pNode : the node
		 * \param  vThreshold : the threshold for the process: all values with energies below will be put to 0
		 * \param rMaxEnergy : reference parameter. It allows to store the maximum energy defined. See also mMaxDefinedEnergyTotalCrs and mMaxDefinedEnergyProcess.
		 * \return : the cross section
		 *
		 * \warning It performs a linear interpolation, because the log interpolation is not valid yet for the  values equal to 0.
		 */
		ublas::vector<double> ExtractCrs(XmlParameters* pParams,TiXmlNode* pNode,double vThreshold,double& rMaxEnergy);



		/**
		 *
		 * Get the process or total crs node
		 * and returns the interpolated cross section.
		 * It verfies if the cross section is not negative after the linear interpolation.
		 * \param pParams : the xml parameters object where the node is extracted
		 * \param pNode : the node
		 * \param  vThreshold : the threshold for the process: all values with energies below will be put to 0
		 * \param rMaxEnergy : reference parameter. It allows to store the maximum energy defined. See also mMaxDefinedEnergyTotalCrs and mMaxDefinedEnergyProcess.
		 * \return : the cross section
		 *
		 * \warning It performs a log interpolation, all the values eqal to 0 are put to 1E-42.
		 */
		ublas::vector<double> ExtractCrsLog(XmlParameters* pParams,TiXmlNode* pNode,double vThreshold,double& rMaxEnergy);

		/**
		 * True if total crs is defined by a cross section
		 * This point is very useful for  the electron part
		 * it allows to refine the interpolation at high 
		 * energies.
		 *
		 *
		 */
		bool mDefinedTotalCrsDirectly;

		/**
		 * If the total cross section is defined direcly
		 * (ie mDefinedTotalCrsDirectly is true)
		 * so we need the maximum energy at which this 
		 * cross section is defined.
		 * This parameter store it.
		 */
		double mMaxDefinedEnergyTotalCrs;

		/**
		 * We store the maximum defined energy for each processes
		 * (mProcessNames)
		 */
		std::deque<double> mMaxDefinedEnergyProcess;


		/// If the total cross section is defined for an electron impact process: add the inelastic (ionization+excitation) and the elastic cross sections!
		bool mbTotalSumIneel;



	public:

		// Public : because needed to compute the
		// energy deposition

		/** To store the threshold for ionization processes
		 * useful for ElecCrs
		 */
		std::deque<double> mIonizationThresholdeV;

		/** To store the threshold for excitation processes
		 * useful for ElecCrs
		 */
		std::deque<double> mExcitationThresholdeV;

		/** To store the threshold for exchange processes
		 * useful for ProtCrs (Hydrogen and Protons)
		 */
		std::deque<double> mExchangeThresholdeV;




		/** To store the ionization processes cross section
		 * useful for ElecCrs
		 */
	//	std::vector< std::vector<double>* > mIonizationCrscm2;
		std::deque< ublas::vector<double>* > mIonizationCrscm2;	
		/// Allows to fill mIonizationCrscm2: we can only add the address when the mCrscm_2 vector is full (is it still necessary with deque???)
		std::deque<unsigned> mIonizationCrsPosition;
		
		/** To store the excitation processes cross section
		 * useful for ElecCrs
		 */
		std::deque< ublas::vector<double>* > mExcitationCrscm2;
//		ublas::matrix<double> mExcitationCrscm2;
		/// Allows to fill mExcitationCrscm2: we can only add the address when the mCrscm_2 vector is full
		std::deque<unsigned> mExcitationCrsPosition;

		/** To store the exchange processes cross section
		 * useful for ProtCrs
		 */
		std::deque< ublas::vector<double>* > mExchangeCrscm2;
		/// Allows to fill mExchangeCrscm2: we can only add the address when the mCrscm_2 vector is full
		std::deque<unsigned> mExchangeCrsPosition;


		/// Deque for defining if there is an Auger process. If so, use the other Deque for the Auger electron energy
		std::deque<bool> mIsAuger;

		/// Deque for defining the Auger electron energy. It is a deque of deque in case of double K-Shell ionization (very high energies!)
		std::deque< std::deque<double> > mAugerEnergy;
		
		/// Deque for defining the Auger electron efficiency (theoretically 1, but can be different if double K-Shell)
		std::deque< std::deque<double> > mAugerEfficiency;



	public:

		/**
		 * Initialize the grid
		 * \param  vGreV : the default energy grid
		 */
		CrossSection(ublas::vector<double> vGreV);
		virtual ~CrossSection();

		/// The total cross section in cm2: used for computing the degradation -> flux
		ublas::vector<double> mTotalCrscm2;

		/// True if the elastic cross section is defined
		bool mbIsDefinedElastic;
		/// Cross section for the elastic processes
		ublas::vector<double> mElasticCrscm2;

		/**
		 *
		 * True is the cross section is defined, false else
		 *
		 */
		bool mIsDefinedCrs;

		/** True if the Total cross section is defined.
		 * This function allows to raise a warning if not defined
		 * and to allow the code to compute that total crs
		 * by adding the different processes cross sections.
		 * 
		 * \attention This flag is added by getting in mind the fact that the 
		 * cross section for electron impact is unclear in 
		 * the fortran version of Trans*.
		 * The main problem with that total cross section for the electrons
		 * is that it is parametrised at high energy!
		 */
		bool mIsDefinedTotalCrs;

		/*=====================================================
		 * Here, the detail of the different cross sections
		 * in order to compute the different products
		 *=====================================================*/


		/// The name of the different processes to ionize or excite. Vector : different processes
		std::deque<std::string> mProcessNames;

		/// The name of the corresponding created species. Vector : different processes. Vector vector : many species can be created
		std::deque< std::deque< SpecieId > > mCreatedSpecies;
		
		
		/// If the species is created in double, or if a branching ratio is used, the multiplicator is needed
		std::deque< std::deque<double> > mCreatedSpeciesMultiplicator;

		/// The number of electron created through each processes. It is possible to put 0 even if an ion is created, to account for a sub-process. (In that case, do not add two time the ion....). This number of electrons is used to compute the energy per electron produced
		std::deque<double> mNumberOfElectrons;

		/// The number of ions created through each processes. It is possible to put 0 even if an electron is created, to account for a sub-process. (In that case, do not add two time the electron produced....). If the key ions is not set, the number of ions created through that process is automatically set to 1. You need to overload it. This number of ions is used to compute the energy per ion produced.
		std::deque<double> mNumberOfIons;
		/// The threshold for the process
		std::deque<double> mThresholdseV;

		/// The corresponding cross sections in cm2
		std::deque< ublas::vector<double> > mCrscm2;


	
		/*=====================================================
		 * Here the functions to init and print cross sections
		 *=====================================================*/



		/**
		 * Load the cross sections from the file given
		 * \param pParams : object to read the cross section file
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : if monte carlo activated
		 * \param bLogInterp : true if you want log interpolation, false else.
		 * \return true if loaded.
		 * load the sub processes
		 * performs the interpolation on the grid
		 * be careful with negative interpolation. 
		 * \warning Make log interpol???
		 */
		bool LoadCrs(XmlParameters *pParams,std::string vName,bool bIsMonteCarlo, bool bLogInterp=false);

		/**
		 * Load the cross sections from the file given
		 * \param vFilename : name of the file where crs can be loaded
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : if monte carlo is active
		 * \param bLogInterp : true if you want log interpolation, false else.
		 * \return true if loaded.
		 */
		bool LoadCrs(std::string vFilename,std::string vName, bool bIsMonteCarlo,bool bLogInterp=false);

		/**
		 * Print the cross section with the new energy grid
		 * into a file. This file can be read with a python-pylab load function
		 * \param vFilename : name of the file. It is created.
		 */

		void PrintCrs(std::string vFilename);



		/// Set up the monte carlo
		void SetMonteCarloActive();


};






#endif
