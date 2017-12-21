/**
 * \defgroup atmospheres_models Atmospheres
 * \file atmo.hpp
 * \brief Defines the Atmo class : allows to define the grids and start the different computation for ionization and excitation
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: atmo.hpp 1921 2014-02-26 05:37:14Z gronoff $
 *
 */

#ifndef ATMO_HPP
#define ATMO_HPP
#include <species/species.hpp>
#include <photo/photoionization.hpp>
#include <math/mathfunction.hpp>
#include "neutralatmo.hpp"
#include <planet/allplanets.hpp>
#include <elecimpact/electronionization.hpp>
#include <proton/proton.hpp>
#include <emit/emission.hpp>
/**
 * \ingroup atmospheres_models
 * Class to define the atmosphere and the energy grids.
 * Stores the different species, including their densities.
 *
 * It is one of the MAIN classes
 */

class  Atmo
{
	private:
		/// The name of the parameter file
		std::string mParameterFile;
		/// Main parameter file
		XmlParameters* mpParameter;


		/// The chemistry object
		boost::shared_ptr<Chem> mChem;

		/// The emission object
		boost::shared_ptr<Emission> mEmit;
		/// To check if the parameter is loaded
		bool mIsMainParameter;
		
		/// To check if the model atmosphere is loaded
		bool mIsModelAtmosphere;

		/// To check if the photoionization model is up
		bool mIsPhotoionizationModel;

		/// To check if the electron impact model is up
		bool mIsElectroionizationModel;
		
		/// To check if the planet is loaded
		bool mIsPlanet;

		/// Init the chemical list \TODO mettre les protons dans cette partie a la FIN
		void InitChem();


		/// Init the altitude, energy grids
		void InitGrids();
		/** Init the altitude grid
		  This code reads the xml file /atmosphere/alt_grid/use_model
		  if the type is 0: it reads the /atmosphere/alt_grid/st_grid parameters
		  if the type is 1: it directly reads the altitude grid
		  */
		void InitAltGrid();

		/** Init the electrons grid
		   Reads the xml file /electron/grid/st_grid
		   if type 0: exp decrease
		   if type 1 : power law decrease

		  */
		void InitElecGrids();

		/** Init the electrons grid
		   Reads the xml file /proton/grid/st_grid
		   if type 0: exp decrease
		   if type 1 : power law decrease

		  */
		void InitProtGrid();


		/* Returns the standard grid
		  \param number : size of the grid
		  \param type : type of the grid 0 -> exp decrease 1 ->  power law decrease  2 -> constant width
		  \param min : min value
		  \param max : max value
		  \param cent : center value of the  grid
		  \param ddeng : width of the grid
		  \param spfac : increasing factor, only for type 1 grid
		  */
	//	void StandardGrid(int number,int type,double min,double max,std::vector<double>& cent,std::vector<double> ddeng,double& spfac);


/*
 *
 * BUG CORRECTED 14 SEPT 2009 -> new grids
 * 		       The Torr and Torr standard 39 box grid is a total sh...
		       in fact, if you use it to interpole the cross sections
		       you have a problem because the mean of this grid is not
		       increasing or decreasing!
		       The smartest thing to do now, with the ease of C++
		       is to cut the continuum boxes to add the lines.
		       The redistribution of the flux function will do the most
		       difficult job : redistribute the values of the solar
		       flux model inside these new boxes.
		       Another point will be to automatically add the line-size
		       boxes inside the photon Grids.
 

 */




		/**
		  Inits the photon grid.
		  Reads /sun/grid/use_model
		  
		  if 0 -> Torr and Torr grid 
		  if 1 -> standard grid 
		  bug CORRECTED line-size boxes not added. Please use the standard Torr and Torr grid
		  	0-> exp decrease
			1-> pow law decrease
			2-> cst width
		   for the power law grid, add the lines. It already works in the tests, and it puts the error in energy redistribution a 1% instead of 5%! -> done 27 sept 09
		  */
		void InitPhotGrid();


		/**
		 *
		 * The gaussian angle class, initialized for the number
		 * of angle found in the xmlparameter file.
		 *
		 */
		MathFunction::GaussianAngle* mpGAngle;

		/**
		 * true if the gaussian angle is initialized
		 */
		bool mIsGaussianAngle;

		/// True if the chemistry is initialized
		bool mbIsChem;



		boost::shared_ptr<MathFunction::GaussianAngle> mpProtonGaussianAngle;

		/// Number of angle for the Proton transport...
		int mNbProtonAngle;



	private:
                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo----------------------------------+
                 |                                                 |
                 |Multiplicators, mainly for Monte-Carlo simulation|
                 |which can also be used inside the chemistry      |
                 |and thus are put inside the model!               |
                 |                                                 |
                 +------------------------------------------------*/


		/// The electron density multiplicator
		double mElecDensMult;
		/// The electron temperature multiplicator
		double mElecTempMult;
		/// The ion temperature multiplicator
		double mIonTempMult;
		/// The neutral temperature multiplicator
		double mTempMult;



	public:
		/**
		 * The constructor, need parameter
		 * \warning  the vertical column density is computed here
		 * BUG CORRECTED jan 2010 if you want to use a modified atmosphere for bent lines -> the modification is made inside another part
		 * this atmosphere is bugged because the vertical column density is computed
		 * \param vParameterFile the file containing all the informations for the run
		 */
		Atmo(std::string vParameterFile);
		/// The destructor
		~Atmo();


		/// Check the chemistry and emissions objects
		bool CheckChemEmit()
		{
			bool resu=true;
			if(mbIsChem)
			{
				resu=resu&&mChem->CheckReacId();
				resu=resu&& mEmit->CheckEmitList();
			}
			return resu;
		}



		/// Class planete: to know the physical parameters and define the atmosphere model
		Planete* mpPlanet;

/*       _\|/_
         (o o)
 +----oOO-{_}-OOo------------+
 |Initialization of the grids|
 |hv grid                    |
 |electron grid              |
 |altitude grid              |
 +--------------------------*/

		/// Number of angle for the electron transport...
		int mNbAngle;


		/// Altitude grid in km
		ublas::vector<double> mAltGridKm;
		

		/// The dB_B parameter
		ublas::vector<double> mdB_B;

		// In a near future
		// we will need the length grid
		// decoupled from the alt grid
		// length grid in km
		//vector<double> length;


		/// Photon grid mean in eV. This should be used for the interpolations. This grid is decreasing
		ublas::vector<double> mPhotonGrideV;
		/// Photon grid in eV : min 
		ublas::vector<double> mPhotonGrideVmin;
		/// Photon grid eV : max
		ublas::vector<double> mPhotonGrideVmax;

		/// Electron bottom grid, the grid is decreasing
		ublas::vector<double> mElecBotEeV;
		/// Electron center grid, the grid is decreasing
		ublas::vector<double> mElecCentEeV;
		/// Electron width grid, the grid is decreasing
		ublas::vector<double> mElecDdengeV;

		/// Proton center grid
		ublas::vector<double> mProtonGrideV;

		/// Proton Width grid
		ublas::vector<double> mProtonWidthGrideV;

		/// Sp factor: not used. But still here to load the grids!
		double mSpfactor;




/*       _\|/_
         (o o)
 +----oOO-{_}-OOo--------------+
 |Grid initialization functions|
 +----------------------------*/




/*       _\|/_
         (o o)
 +----oOO-{_}-OOo---------------------------------------------------------+
 |Model species : the different species for ionization, electron impact...|
 |Result species : the different species coming from the ionization ...   |
 +-----------------------------------------------------------------------*/
// The model species are initialized in the class NeutralAtmo
// The result species are stored in the vector ResultSpecies
//

	/// The object containing the model species
	NeutralAtmo* mpModelAtmosphere;
	/// The vector containing the result species
	std::deque<Specie*> mResultSpecies;
	/// The vector containing the electron density in cm-3
	ublas::vector<double> mElectronDensity;
// Also : temperatures

	/// The vector containing the electron temperature in K
	ublas::vector<double> mElectronTemperature;
	/// The vector containing the neutral temperature in K
	ublas::vector<double> mNeutralTemperature;
	/** The vector containing the ion temperature
	 * All the ions are considered having the same temperature, in  K
	 */
	ublas::vector<double> mIonTemperature;


/*       _\|/_
         (o o)
 +----oOO-{_}-OOo--------------------+
 |Atmosphere initialization functions|
 +----------------------------------*/

	/**
	 * Initializes the planet
	 * \xrefitem usage "usage"  "Usage"  To add a new planet, you have to make an heritage of your planet object. Then, you have to add it in this function
	 */
	void InitPlanet();

	/**
	 * Initializes the neutral atmosphere
	 */
	void InitNeutralAtmo();

	/**
	 * Initializes the neutral temperature
	 * called by InitNeutralAtmo.
	 * \xrefitem nb "Nb" "Nb" now, it is a very simple function: it calls the model_atmosphere
	 * temperature() function, but it can evolve e.g. by calling an empiric function
	 */
	void InitNeutralTemp();

	/**
	 * Initializes the ionosphere:
	 * electron density and temperature
	 * ion temperature
	 */
	void InitIonosphere();



        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo----------------------------------------------------------------------+
         |                                                                                     |
         |Atmosphere modification function (for fitting data, see multiple computation section)|
         |                                                                                     |
         +------------------------------------------------------------------------------------*/

	/**
	 * Modify the electron precipitation spectra
	 * \param vModel : the model of the precipitation (gaussian, Dirac,...)
	 * \param vParams : the necessary parameters for the model
	 * \param vAddParams : additionary parameters (non variable for adjustement)
	 */
	void ResetElecPrecip(int vModel, std::deque<double> vParams, std::deque<double> vAddParams);
	
	/**
	 * Modifies the density of a species, and also the column density (only for non-bent).
	 * All of this based on a model of density variation
	 * \param vName : the name of the species to modify (not excited state yet)
	 * \param vModel : the model of density (gaussian, chapman,...)
	 * \param vParams : the necessary parameters for the model
	 * \param vAddParams : additionary parameters (non variable for adjustement)
	 */
	void ResetSpecieDens(std::string vName,int vModel, std::deque<double> vParams, std::deque<double> vAddParams);
	
	/**
	 * Modifies the density of a species, based on data
	 * \param vName : the name of the species to modify (not excited state yet)
	 * \param vAltKm : the altitudes for the interpolation
	 * \param vDenscm_3 : the values of the density at these altitudes
	 */
	void ResetSpecieDensInterp(std::string vName, ublas::vector<double> vAltKm, ublas::vector<double> vDenscm_3);

	/**
	 * Modifies the density of the electron (themal)
	 * \param vModel : the model of density (gaussian, chapman...)
				0 : Exponential (Bates-Walker profile)
				1 : Chapman function
				2 : Gaussian function
	 * \param vParams : the necessary parameters for the model
	 * \param vAddParams : additionary parameters (non variable for adjustement)
	 */
	void ResetElectronDens(int vModel, std::deque<double> vParams, std::deque<double> vAddParams);


	/**
	 * Modifies the electron density, based on data
	 * \param vAltKm : the altitudes for the interpolation
	 * \param vDenscm_3 : the values of the density at these altitudes
	 */
	void ResetElectronDensInterp(ublas::vector<double> vAltKm, ublas::vector<double> vDenscm_3);



	/**
	 * For retrieving the main computed parameters: VER, production, density,...
	 * This function is necessary for the fitting capabilities
	 *
	 * \param vId: the id of the species to retrieve
	 * \param vParams: Identifiant of the data to retrieve: the first parameter gives the general tyoe:
	 * 		- 0 : production (/volume)
	 * 		- 1 : volume emission rate (needs a second parameter for the wavelength)
	 * 		- 2 : limb flux  (needs a second parameter for the wavelength)
	 * 		- 3 : density
	 * \param vComparisonScale : the scale where to interpolate the results. If length = 0; the results are not modified
	 * \return the vector with the selected data
	 */
	ublas::vector<double> RetrieveObservable(SpecieId vId, std::deque<int> vParams,ublas::vector<double> vComparisonScale );

 	/*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-------+
         |                      |
         |The bending parameters|
         |                      |
         +---------------------*/

	// True if the bending is activated!
	bool mIsBendingActivated;

	/**
	 *  The atmosphere Path: the geometry of the electron
	 *  precipitation if it is not vertical.
	 *
	 */
	Path mAtmoPath;


	/// The SZA where the parameters must be computed
	std::deque<double> mComputeSZA;
	//ublas::vector<double> mComputeSZA;


	/// The bent electron temperature
	ublas::vector<double> mBentElectronTemperature;
	/// The bent electron density
	ublas::vector<double> mBentElectronDensity;
	/// The bent ion temperature
	ublas::vector<double> mBentIonTemperature;
	/// The bent neutral temperature
	ublas::vector<double> mBentNeutralTemperature;

	/// The species neutral atmosphere in the bent frame
	std::deque<Specie*> mAtmoBentSpecies;


	/// List the neutral atmospheres in the SZA where it is computed
	std::deque< std::deque<Specie*> > mAtmoSpeciesSZAs;
	/// List the neutral temperatures...
	std::deque< ublas::vector<double>* >  mNeutralTempSZAs;
	/// List the electron temperatures...
	std::deque< ublas::vector<double>* >  mElecTempSZAs;
	/// List the ion temperatures...
	std::deque< ublas::vector<double>* >  mIonTempSZAs;
	/// List the electron densities...
	std::deque< ublas::vector<double>* >  mElectronDensitySZAs;
	

	private:
	/// Bend the neutral atmosphere and the ionosphere... \todo modify the atmosphere according to the SZA
	void BendMyAtmosphere();
	
	public:
	/**
	 * Compute the photoionization in the case
	 * of a bent atmosphere
	 */
	void PhotoIonizeBent();
	/// The list of the photoelectron fluxes, related to mComputeSZA
	std::deque< EFlux* > mPhotoFluxSZAs;
	/// The list of the species ionization for the different SZAs
	std::deque< std::deque< Specie* > > mResuPhotoSZAs; 
	/// If the photon flux is set for the sza
	bool mbIsPhotoFluxSZAs;

	// The photoionization electrons fluxes will need a deque

	/**
	 * Compute the cosmic ray ionization in the case of a bent atmosphere
	 */
	void CosmoIonizeBent();

	/// The list of the cosmoelectron fluxes, related to mComputeSZA
	std::deque< EFlux* > mCosmoFluxSZAs;
	/// The list of the species ionization for the different SZAs
	std::deque< std::deque< Specie* > > mResuCosmoSZAs; 
	/// If the cosmo flux is set for the sza
	bool mbIsCosmoFluxSZAs;


	/**
	 * Compute the proton/hydrogen ionization in the case of a bent atmosphere
	 */
	void ProtonIonizeBent();

	/// The list of the protonelectron fluxes, related to mComputeSZA
	std::deque< EFlux* > mProtonEFluxSZAs;
	/// The list of the species ionization for the different SZAs
	std::deque< std::deque< Specie* > > mResuProtonSZAs; 
	/// If the proton flux is set for the sza
	bool mbIsProtonFluxSZAs;





	/**
	 * And the actual electron impact ionization computation.
	 * Not so far from the initial one...
	 */
	void ElectroIonizeBent();



// Model to Result


        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-----------+
         |                          |
         |TO DO THE MAIN COMPUTATION|
         |                          |
         +-------------------------*/

	/// The main computation function, if a bent atmosphere is found, the ComputeBent function is called.
	void Compute();


	/// If we do not want to compute the production but to only read it in a file, this is the function
	void ReadProductions();

	/** The computation in the case of a bent atmosphere 
	*/
	void ComputeBent();

	/// The main species result
	std::deque<Specie*> mTotalResu;


	/// The computation of the photodissociation branching ratio
	void ComputePhotodissociation(std::string suffix="");

        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo------------+
         |                           |
         |Do the photoionization work|
         |                           |
         +--------------------------*/


	/// The photoionization model. only active for the PhotoIonize function. But it could be interesting later...
	Photoionization* mpPhotomodel;


	/**
	 * Performs the photoionization computation
	 * in the vertical atmosphere
	 */
	void PhotoIonize();
	/**
	 * Computes the absorption by the different species:
	 * along the solar path if nothing is defined, or along a specified 
	 * orientation
	 */
	void PhotoAbsorb();

	/**
	 * Stores the photoionization  resulting species
	 *
	 */
	std::deque<Specie*> mPhotoionizationResu;

	/**
	 * Stores the photoelectron flux
	 */
	EFlux* mpPhotoFlux;

	/// True if the mpPhotoFlux is defined
	bool mbIsPhotoFluxDefined;
        
	
	/*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-------+
         |                      |
         |                      |
         |Do the cosmic ray work|
         |                      |
         +---------------------*/

	/// The cosmic ionization function
	
	void CosmoIonize();
	/**
	 * Stores the cosmoelectron flux
	 */
	EFlux* mpCosmoFlux;

	/// True if the mpCosmoFlux is defined
	bool mbIsCosmoFluxDefined;

	/// True if the CosmoModel is defined
	bool mbIsCosmoModel;
	/// The model for ionization by cosmic rays
	ProtoCosIonization* mpCosmoModel;

	/**
	 * Stores the cosmoionization  resulting species
	 *
	 */
	std::deque<Specie*> mCosmoionizationResu;
	/*       _\|/_
                 (o o)
         +----oOO-{_}-OOo-------+
         |                      |
         |                      |
         | Do the Proton/H work |
         |                      |
         +---------------------*/
	/// The proton ionization function
	void ProtonIonize();
	/**
	 * Stores the protolectron flux
	 */
	EFlux* mpProtonEFlux;

	/// True if the mpCosmoEFlux is defined
	bool mbIsProtonEFluxDefined;

	/// True if the ProtonModel is defined
	bool mbIsProtonModel;
	/// The model for ionization by proton/hydrogen
	ProtonHydrogenTransport* mpProtonModel;




	/// Stores the proton ionization
	std::deque<Specie*> mProtoionizationResu;



        /*       _\|/_
                 (o o)
         +----oOO-{_}-OOo------------+
         |                           |
         |Do the electron import work|
         |                           |
         +--------------------------*/
	/**
	 * Performs the electron impact computation
	 * in the vertical atmosphere.
	 * Load the electron precipitation flux model
	 * add it with the different electron fluxes
	 */
	void ElectroIonize();

	/// As for the photoionization model, only active in ElectroIonize()...
	ElectronImpactIonization* mpElecModel;

	/// True if the elec model defined
	bool mbIsElecModelDefined;
	/**
	 * Stores the electron precipitation flux
	 */
	EFlux* mpElecPrecipFlux;
	
	
	/** 
	 * Stores the total precipitation flux
	 */
	EFlux* mpTotalFlux;
	
	/// True if the mpElecPrecipFlux is defined
	bool mbIsElecPrecipFluxDefined;
	/**
	 * Store the electron impact resulting species
	 *
	 */
	std::deque<Specie*> mElectronImpactResu;

/*

	We initialise the objects so that several computations can be made with the same ones. 
	This is necessary for working with least square methods

*/

	/** 
	 * Initialize the several computation objects for multiple calculations
	 * \todo put that compatible with bent atmosphere
	 */
	void InitMultipleComp();

	/**
	 * Start one computation compatible with multiple computations
	 */
	void MCompute();

	/**
	 * Delete the objects initialized for multiple calculations
	 */
	void FinishMultipleComp();


	private:
		/// True if the photoionization is launched in MCompute
		bool mMultiplePhoto;
		/// True if the cosmoionization is launched in MCompute
		bool mMultipleCosmo;
		/// True if the protonionization is launched in MCompute
		bool mMultipleProton;
		/// True if the electroionization  is launched in MCompute
		bool mMultipleElec;
		/// True if the Chemistry  is launched in MCompute
		bool mMultipleChem;
		/// True if the emissions are launched in MCompute
		bool mMultipleEmit;

	public:




/*       _\|/_
         (o o)
 +----oOO-{_}-OOo-------------------+
 |Print the different inputs/outputs|
 +---------------------------------*/

	/**
	 * Creates two files, one for the ionosphere density model
	 * ie the electron density model
	 * and one for the temperature model.
	 * The two files have the same prefix, given in parameter
	 * and they give _density.dat and _temp.dat as suffixes
	 * \xrefitem nb "Nb" "Nb" now, it only print the electron
	 * density model in _density. In the future, it will be possible
	 * to improve the function by adding the results
	 *
	 * \param vFilenamePrefix : the prefix for the two files created.
	 *
	 *
	 */
	void  PrintIonosphere(std::string vFilenamePrefix);
	
	/**
	 * Calls the NeutralAtmo::PrintNeutralAtmo function. Used to check LM convergence
	 * \param vFile : the name of the output file
	 */
	void PrintNeutralAtmo(std::string vFile);
	/* Same functionm except it works only with the bent atmosphere
	 *
	 * \param vFilenamePrefix : the prefix for the two files created.
	 */
	void  PrintBentIonosphere(std::string vFilenamePrefix);
	
	
	
	/**
	 * Creates the cross section files for every 
	 * neutral atmosphere species
	 *
	 */
	void  PrintSpecies(std::string vFilenamePrefix);

	/**
	 * Compute and writes the chemistry and emissions processes
	 * It should be called before ProceedOutputsm so that the outputs can write suppl parameters
	 * \param suffix : allows to add a suffix at the end of the filenames
	 *  The suffix is needed when you want to do monte-carlo simulations.
	 */
	void ProceedEmissions(std::string suffix="");


	/**
	 * Read the parameters file for the output neededs.
	 * \param suffix : allows to add a suffix at the end of the filenames
	 *  The suffix is needed when you want to do monte-carlo simulations.
	 *  \todo make option to print the SZA results when a bent atmosphere is selected
	 *  \todo print electron production when bent atmosphere
	 */
	void ProceedOutputs(std::string suffix="");
};

	



#endif
