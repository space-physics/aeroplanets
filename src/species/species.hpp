 /** 
 * \file species.hpp
 * \brief Defines the specie class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: species.hpp 1922 2014-02-26 22:38:52Z gronoff $
 *
 */


#ifndef SPECIES_HPP
#define SPECIES_HPP
#include "path.hpp"

/// Forward declaration in order to compile
class Specie;
/**
 * \ingroup specie_parameters
 * Namespace for useful tools to use with specie
 *
 *
 */
namespace SpecieUtils
{
	/**
	 * Print the productions computed inside a species vector.
	 * It prints the TotProduction, and if the option is checked, also the state production
	 * \param vSpecies : the result species vector
	 * \param vAltGridKm : the altitude grid
	 * \param vFilename : the file were to write
	 * \param vPrintStateProduction : if true, print also the state production, initialized at false
	 * \param vIsLength : true if we work with the length in a bent frame
	 */
	void PrintProduction(std::deque<Specie*> vSpecies,ublas::vector<double> vAltGridKm,std::string vFilename,bool vPrintStateProduction=false,bool vIsLength=false);


	/**
	 * Print the electron production.
	 * The production comes from the species vector of electron production. Therefore, it is not possible to use the result species to compute the total electron production. You should use the neutral atmosphere that has been ionized
	 *
	 * \param vSpecies : the result species vector
	 * \param vAltGridKm : the altitude grid
	 * \param vFilename : the file were to write
	 * \param vIsLength : true if we work with the length in a bent frame
	 *
	 */
	void PrintElectronProduction(std::deque<Specie*> vSpecies,ublas::vector<double> vAltGridKm,std::string vFilename,bool vIsLength=false);

	/**
	 * Print the production of the chosen species
	 * More precisely, we give a list of species-state (specieID)
	 * and the code search into species if this state exists
	 * if the state in specieID is void, we plot the total production
	 * of the species (sum of states).
	 * \param rStates : the states to print
	 * \param rSpecies : the vector of the Specie
	 * \param vAltKmGrid : the grid
	 * \param rFilename : the name of the file to print
	 * \param vIsLength : true if we work with the length in a bent frame
	 *
	 */
	void SelectedPrintProduction(std::deque<SpecieId>& rStates, std::deque<Specie*> rSpecies, ublas::vector<double>  vAltKmGrid, std::string rFilename,bool vIsLength=false);

	/**
	 * Takes a species vector and create the result species vector
	 * This function merges the different productions, for example O+ can come from photoionization (or electron impact...) of CO O O2 CO2...
	 * Therefore, the O+ production is stored in CO O O2 CO2 Species object
	 * Thanks to that function, a unique O+ Species is created, and all the different sources are added.
	 * \param rInitspec : the initial species vector
	 * \param rResuspec : the resulting species vector
	 */

	void SpeciesToResu(std::deque<Specie*>& rInitspec,std::deque<Specie*>& rResuspec);

	/**
	 * Merge two result species vector.
	 * If a vector cames from photoionization, and the other from electron impact ionization, it can be very interesting to add them to compute 
	 * total ionization.
	 * \param vS1 :the first vector to merge
	 * \param vS2 : the second vector to merge
	 * \return the merged vector
	 */
	std::deque<Specie*> MergeResult(std::deque<Specie*> vS1,std::deque<Specie*> vS2);


	/**
	 * Find the position of the specie named by its string in the vector of species pointer.
	 * \xrefitem nb "Nb" "Nb" I think this function  can be improved
	 * \param name : the name of the searched species. Be careful, this is just the name, not the excitation value
	 * \param sp : the vector of species where to find it. There must be one occurence of the specie in that vector (else, it find the first)
	 * \return the position
	 */

	int PosOfSp(std::string name,const std::deque<Specie*>& sp);

	/**
	 * Get the production of one specie in the list of the species
	 * \param vName : the name of the species
	 * \param vState : the state of the species
	 * \param vSp : the list of species where to search for
	 * \param rDensity : the output production
	 * \return bool : true if the production was found
	 */
	bool GetSpecieProduction(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double> & rDensity);
	/**
	 * Get the density of one specie in the list of the species
	 * \param vName : the name of the species
	 * \param vState : the state of the species
	 * \param vSp : the list of species where to search for
	 * \param rDensity : the output density
	 * \return bool : true if the density was found
	 */
	bool GetSpecieDensity(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double> & rDensity);
	
	/**
	 * Put the density of one specie in the list of the species
	 * \param vName : the name of the species
	 * \param vState : the state of the species
	 * \param vSp : the list of species where to search for
	 * \param vDensity : the intput density
	 * \return bool : true if the species was found and the density addsd
	 */
	bool PutSpecieDensity(std::string vName,std::string vState,const std::deque<Specie*>& vSp,ublas::vector<double>  vDensity);




	/**
	 * Bent the atmosphere, to have a set of species along a magnetic field line, but with the electron and species production considered :D
	 * \param vAtmoSpecies : the set of set of atmosphere Species. set of (set of atmospheric species usually considered as the atmosphere)
	 * \param vAltitudesKm : the altitude corresponding to the vSpecies grids
	 * \param vSZADegree : SZA for each vSpecies
	 * \param vPath : The path for interpolating the new species
	 *
	 * \return OutSpecies : the set of species considered as the output atmosphere (data along the field line, and especially the column density).
	 */
	std::deque<Specie*> BendAtmosphere(std::deque< std::deque<Specie*> > vAtmoSpecies,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath);


	/**
	 * Copy an atmosphere:
	 * allows to modify the specificities of each atosphere without computing again all the cross sections
	 * \param vAtmo: the deque of the inital species (the atmo)
	 */
	std::deque<Specie*>  CopyAtmo(std::deque<Specie*> vAtmo);
};






/**
 *  The Specie: base of the system
 *  It contains all the parameters: density...
 *  And store the output (production)...
 *
 *  It is also used when computing the density of ions......
 * \ingroup specie_parameters
 */
class Specie
{
	private:
		/**
		 * True if the PhotCrs is loaded. It allows to delete it when the program stops
		 */
		mutable bool mIsPhotoCrsLoaded;
		/// True if PhotoCrs copied
		mutable bool mIsPhotoCrsCopied;
		/// True if ElecCrs loaded
		mutable bool mIsElecCrsLoaded;
		/// True if ElecCrs copied
		mutable bool mIsElecCrsCopied;
		/// True if ProtCrs loaded
		mutable bool mIsProtCrsLoaded;
		/// True if ProtCrs copied
		mutable bool mIsProtCrsCopied;


		/// True if HCrs loaded
		mutable bool mIsHCrsLoaded;
		/// True if HCrs copied
		mutable bool mIsHCrsCopied;



		/// True if spparams is read only in that object: in that case, it is destroyed at the end
		bool mIsmpParamsHere;

		/** Merges the production of different species
		  increase the values if exists
		  */
		void MergingProd();


		/**
		 * Load the different parameters for the specie
		 * by reading the corresponding xml file
		 * \param bIsMonteCarlo : true if monte carlo is active.
		 */
		void LoadSpecie(bool bIsMonteCarlo);

		/**
		 * Load the different parameters, but not
		 * the cross sections: this is not useful for result species
		 *
		 */
		void LoadResuSpecie();

		// The xmlparameters file to read the initialized...
		XmlParameters* mpParams;


		/// The photon energy grid : needed for the interpolation of crs
		ublas::vector<double>* mpPhotonEgrideV;

		/// The electron energy grid : needed for the interpolation of crs
		ublas::vector<double>* mpElectronEgrideV;

		// The electron energy grid width : needed for the redistribution of the electron cross section
		ublas::vector<double>* mpElectronEngddeV;

		/// The proton energy grid : needed for the interpolation of crs; the proton grid is also used for the Hydrogen 
		ublas::vector<double>* mpProtonEgrideV;

		/**
		 * Put a state production. Put the production in the total production of the species (except if "-NOTOT" is added to the name of the species")
		 * \param vName : the name of the state
		 * \param vProductioncm_3s_1 : the vector production in cm-3s-1
		 *
		 */

		void PutStateProduction(std::string vName,ublas::vector<double> vProductioncm_3s_1);

		/**
		 * Fills the TotProduction vector thanks to the different state production
		 * \deprecated Now the production is filled in the putstateproduction
		 * 
		 */
		void StatesToTotProd();

		/// Check if porter parameters are loaded
		bool mbPorterLoaded;

		/// The porter energy
		ublas::vector<double>* mpEPortereV;
		/// The porter Gamma
		ublas::vector<double>* mpGPorter;
		/// The porter Beta
		ublas::vector<double>* mpBPorter;
		/// The porter Alpha
		ublas::vector<double>* mpAPorter;

	public:
		/// For the proton redistribution: the sigma parameter
		ublas::vector<double> mProtonSigma;
		/// For the proton redistribution: the redistribution option
		ublas::vector<int> mProtonRedistributionFunction;

		/// The phase function redistribution matrix for the protons
		ublas::vector< ublas::matrix<double> > mProtonPhase;


	public:
		/// The number of electron in each ionization for the cosmic rays
		std::deque<double> mCosmoNumberOfElectrons;


		/**
		 * Stores the position of the pic of the production computed by
		 * the ComputePicProd() function.
		 * This parameter can be user to compute the altitude of the production peak.
		 */
		unsigned mProdPicPosition;

		/**
		 * Computes the position of the maximum of the mTotProductioncm_3s_1 and stores it in mProdPicPosition.
		 */
		void ComputePicProd();


	public:




		/**
		 * Definition of the specie:
		 * \param vName : the name of the specie (to read its properties 
		 * 				in the file)
		 * \param pPhgrideV : photon energy grid
		 * \param pElgrideV : electron energy grid
		 * \param pElEngddeV : electron energy grid width
		 * \param pPrgrideV : proton energy grid
		 * \param spParams : the xml parameters object. Link to the
		 * 		      species properties (see example)
		 * \param bIsMonteCarlo : true if the MC is activated
		 */
		Specie(std::string vName,
				ublas::vector<double>* pPhgrideV,
				ublas::vector<double>* pElgrideV,
				ublas::vector<double>* pElEngddeV,
				ublas::vector<double>* pPrgrideV,
				XmlParameters* spParams,
				bool bIsMonteCarlo);

		/**
		 * Definition of the specie, but with a file link instead of 
		 * an xmlparameter object
		 *
		 * \param vName : the name of the specie (to read its properties 
		 * 				in the file)
		 * \param pPhgrideV : photon energy grid
		 * \param pElgrideV : electron energy grid
		 * \param pElEngddeV : electron energy grid width
		 * \param pPrgrideV : proton energy grid
		 * \param vParamfilename : name of the specie property file
		 * \param bIsMonteCarlo : true if the MC is activated
		 */
		Specie(std::string vName,
				ublas::vector<double>* pPhgrideV,
				ublas::vector<double>* pElgrideV,
				ublas::vector<double>* pElEngddeV,
				ublas::vector<double>* pPrgrideV,
				std::string vParamfilename,
				bool bIsMonteCarlo);

		/** 
		 * Definition of a specie, with a limited initalization
		 * stands for the result species
		 * \param vName : the name of the species
		 * \param spParams : the specie property parameter
		 */
		Specie(std::string vName,XmlParameters* spParams);

		/** 
		 * Definition of a specie, with a limited initalization
		 * stands for the result species
		 * \param vName : the name of the species
		 * \param vParamfilename : name of the specie property file
		 */
		Specie(std::string vName,std::string vParamfilename);
		/**
		 *
		 * Definition of the transformation of a deque of species, with a path
		 * to a species
		 * \param vSpecies : the list of species, with inside values corresponding to the different parameters (nb, all the vectors inside will be transformed, therefore, a solution of an ionization can be transformed thanks to that function)
		 * \param vAltitudesKm : the altitude corresponding to the vSpecies grids
		 * \param vSZADegree : SZA for each vSpecies
		 * \param vPath : The path for interpolating the new species
		 */
		Specie(std::deque< Specie* > vSpecies,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath);


		/// Destructeur
		~Specie();

		/// To clear the different productions when computing a new one
		void ClearProd();

		/**
		 * Puts the density of a state (add the state in the system if it does not exists; resize the density list if necessary)
		 * \param vName the name of the state we are filling
		 * \param vDensitycm_3 the density to put
		 */
		void PutStateDensity(std::string vName,ublas::vector<double> vDensitycm_3);

		/**
		 * Returns the density of a state, returns false if the density does not exists
		 * \param vName the name of the state we are search for
		 * \param rDensitycm_3 the density to save
		 * \return true if the density was found
		 */
		bool GetStateDensity(std::string vName,ublas::vector<double>& rDensitycm_3);

		/**
		 * Returns the production of a state, returns false if the production does not exists
		 * \param vName the name of the state we are search for
		 * \param rProductioncm_3s_1 the production to save
		 * \return true if the production was found
		 */
		bool GetStateProduction(std::string vName,ublas::vector<double>& rProductioncm_3s_1);
		/*               _\|/_
				 (o o)
				 +----oOO-{_}-OOo----+
				 |Physical properties|
				 +------------------*/




		/// The name of the species
		std::string mName;
		/// The mass of the species
		double mMass;

		/// States of the species, beggining with the ground state!
		std::deque< std::string > mStates;
		/*************************/
		/*                       */
		/* The porter parameters */
		/*                       */
		/*************************/

		/// Loads the porter parameters
		void LoadPorter();
		/// Returns the porter energy parameter
		ublas::vector<double>* ReturnEPorter();   
		/// Returns the porter gamm parameter
		ublas::vector<double>* ReturnGPorter();   
		/// Returns the porter beta parameter
		ublas::vector<double>* ReturnBPorter();   
		/// Returns the porter alpha parameter
		ublas::vector<double>* ReturnAPorter();   
		/*               _\|/_
				 (o o)
				 +----oOO-{_}-OOo----------+
				 |Total production and loss|
				 +------------------------*/

		/// The total density of the species (including excited states) (function of altitude)
		ublas::vector<double> mTotDensitycm_3;

		/// The total column density cm-1
		ublas::vector<double> mColDenscm_2;

		/// The scaleheight of the specie in function of the altitude \warning used to compute the photoionization.  \xrefitem usage "usage" "Usage" Used by the photoionization in vertical case. To work along magnetic field line, you mustnot use this one. It is in cm
		ublas::vector<double> mScaleHcm;


		/// The total production of the species (including excited states) (function of altitude)
		ublas::vector<double> mTotProductioncm_3s_1;
		/// The total loss of the species (including excited states) (function of altitude)
		ublas::vector<double> mTotLosscm_3s_1;

		/// For working with the escape
		double mTotEscapecm_2s_1;


		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo-------------+
				 |Density, production and loss|
				 |depending on the state      |
				 +---------------------------*/

		/// The  density of the species (including excited states) (function of state and altitude)
		std::deque< ublas::vector<double> > mStateDensitycm_3;
		/// The  production of the species (including excited states) (function of state and  altitude)
		std::deque< ublas::vector<double> > mStateProductioncm_3s_1;
		/// The  loss of the species (including excited states) (function of state and  altitude)
		std::deque< ublas::vector<double> > mStateLosscm_3s_1;
		/// The escape of each state
		std::deque<double> mStateEscapecm_2s_1;



		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo-----------------------------+
				 |The Cross sections, depending on the process|
				 +-------------------------------------------*/

		/// Photoionization cross section. Pointer, to load it correctly after
		CrossSection* mpPhotoCrs;

		/// Electron impact cross section. Pointer, to load it correctly after
		ElecCrossSection* mpElecCrs;

		///  Proton impact cross section. Pointer, to load it correctly after
		ProtCrossSection* mpProtCrs;

		///  Hydrogen impact cross section; necessary for the coupled transport H/p. This is still a ProtCrossSection object (there are not enough differences between protons and hydrogen physical properties for the cross sections to need another class). Pointer, to load it correctly after
		ProtCrossSection* mpHCrs;

		/*    		 _\|/_
				 (o o)
				 +----oOO-{_}-OOo-------------------------------------+
				 |Result of the computing                             |
				 |The electron production, depending on each processes|
				 +---------------------------------------------------*/

		/// Electron production due to photoionization of the species
		ublas::vector<double> mPhotoElecProductioncm_3s_1;
		/// Electron production due to electroionization of the species
		ublas::vector<double> mElecElecProductioncm_3s_1;
		/// Electron production due to protoionization of the species
		ublas::vector<double> mProtElecProductioncm_3s_1;
		
		/// Electron production due to Hionization of the species
		ublas::vector<double> mHElecProductioncm_3s_1;

		/// Electron production due to cosmic ray ionization
		ublas::vector<double> mCosmoElecProductioncm_3s_1;


		/// Total electron production
		ublas::vector<double> mElecProductioncm_3s_1;

		/*    		 _\|/_
				 (o o)
				 +----oOO-{_}-OOo------------------------------------+
				 |Production of Species, depending upon the processus|
				 +--------------------------------------------------*/

		/**
		 * compute the production of the subspecies
		 * calls merging_prod. This function is called by 
		 * species_to_resu.  That is why it is public
		 *
		 */
		void FinishProduction()
		{
			MergingProd();
		}


		/**
		 *
		 * Init the 3 following vectors for the photoionization processes
		 * \param vNbAlt : the number of altitudes considered.
		 * can be different from the size of the ColDens if we work on a curved atmosphere.
		 */
		void InitPhotoIonization(unsigned vNbAlt);

		/**
		 *
		 * Init the 3 following vectors for the electron impact ionization processes
		 * \param vNbAlt : the number of altitudes considered.
		 * can be different from the size of the ColDens if we work on a curved atmosphere.
		 */
		void InitElectronImpact(unsigned vNbAlt);

		/**
		 * Init the vectors for the proton/hydrogen ionization.
		 * Reads the ponderation array for that
		 * \param vNbAlt : the number of altitudes in the grid
		 * \param vGa : smart pointer to the gaussian angle object for the proton transport
		 * \param vIsRedis : boolean to check if the redistribution must take place
		 * \param vBeamSpreading : correction factor for the beam Spreading (see proton.cpp)
		 */
		void InitProtonImpact(unsigned vNbAlt, boost::shared_ptr<MathFunction::GaussianAngle> vGa, bool vIsRedis, double vBeamSpreading);

	protected:
		/**
		 * Init the phase function for the proton/hydrogen ionization.
		 * \param vGa : smart pointer to the gaussian angle object for the proton transport
		 * \param vIsRedis : boolean to check if the redistribution must take place
		 * \param vBeamSpreading : correction factor for the beam Spreading (see proton.cpp)
		 */
		void InitProtonPhase(boost::shared_ptr<MathFunction::GaussianAngle> vGa, bool vIsRedis, double vBeamSpreading);
	public:

		/**
		 * Init the vectors for the cosmic ray ionization.
		 * Reads the ponderation array for that
		 * \param vNbAlt : the number of altitudes considered.
		 */
		void InitCosmoImpact(unsigned vNbAlt);

		/** 
		 * Allows to define at 0 the vector containing the state production, density...
		 * Necessary for merging neutral atmosphere with productions!
		 */
		void InitDummyStates(unsigned vNbAlt);

		/// The name of the different processes to ionize or excite. Vector : different processes
		std::deque<std::string> mProcessNames;

		/// The name of the corresponding created species. Vector : different processes. Vector vector : many species can be created
		std::deque< std::deque< SpecieId > > mCreatedSpecies;

		/// The multiplicative factors
		std::deque< std::deque<double> > mMultiplicativeFactor;

		/// The corresponding production
		std::deque< ublas::vector<double> > mSpeciesProductioncm_3s_1;

		// Probability of a production
		std::deque<double> mSpeciesProductionProbabilitys_1;
		double mPhotoElecProductionProbabilitys_1;
		double mPhotoDissociationProbabilitys_1;


		std::string ReturnPhotoProbabilityInfo();

		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo------+
				 |Production of Species|
				 +--------------------*/


		/// The different species, comes from a merging of created_species. The merging is a private function
		std::deque<SpecieId> mTotSpecies;
		/// The different species production
		std::deque< ublas::vector<double> >  mTotSpeciesProductioncm_3s_1;




		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo----+
				 |                   |
				 |Print useful things|
				 |                   |
				 +------------------*/

		/** \todo function to print the densities in a file
		 * \param vFilename : the name of the file
		 */
		void PrintDensities(std::string vFilename);

		/** 
		 * \param vFilename : the name of the file
		 * \param vAltGridKm : the altitude grid
		 */
		void PrintProductions(std::string vFilename,const ublas::vector<double>& vAltGridKm);

		/**
		 * Print the cross section of the species thanks to the prefix.
		 * The name of the species is added, and the specific cross section (electron, photon...) too!
		 * \param vFileprefix : the prefix of the filename
		 */
		void PrintCrossSections(std::string vFileprefix);

		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo----------------------------+
				 |                                           |
				 |                                           |
				 |To make the merging easy : friend operators|
				 |                                           |
				 |                                           |
				 +------------------------------------------*/



		/**
		 * Allows to add a species to the present species.
		 * \param s2 : the second specie
		 */
		Specie& operator+=(const Specie& s2);

		/**
		 * Copy operator
		 * \param s1 : the specie to be copied
		 */
		// Specie operator=(const Specie& s1);
		Specie(const Specie&s1);


		/// We need to declare that function friend to access the mainParameter. Used by the new species
		friend void SpecieUtils::SpeciesToResu(std::deque<Specie*>& rInitspec,std::deque<Specie*>& rResuspec);





		/*       	 _\|/_
				 (o o)
				 +----oOO-{_}-OOo--------------------------------------------------+
				 |                                                                 |
				 |Accessors : allows to check if we have the  cross sections loaded|
				 |                                                                 |
				 +----------------------------------------------------------------*/

		/// True if the electron cross section is available
		bool CheckElecCrs() const
		{
			return mIsElecCrsLoaded||mIsElecCrsCopied;
		}
		/// True if the photon cross section is available
		bool CheckPhotCrs() const 
		{
			return mIsPhotoCrsLoaded||mIsPhotoCrsCopied;
		}
		/// True if the proton/H cross section is available
		bool CheckProtCrs() const
		{
			return (mIsProtCrsLoaded||mIsProtCrsCopied) and (mIsHCrsLoaded||mIsHCrsCopied);
		}



};






#endif
