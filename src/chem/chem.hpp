/**
 * \file chem.hpp
 * \brief Allows to set up the chemistry, initialize the different reactions and density computation, and allows to retrieve the additional atmospheric species necessary for the chemistry.
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: chem.hpp 1099 2010-08-03 22:06:10Z gronoff $
 *
 */


#ifndef CHEM_HPP
#define CHEM_HPP

#include "calcdens.hpp"



/**
 * The chemistry class:
 * - is used to initialize all chemical reactions
 * - is used to initialize the photochemical equilibrium density computation
 * - is used to read the atmosphere for the chemistry computation 
 * (and store them)
 *
 * It is important to note that the neutral densities for the chemistry
 * can be totally different than the densities for the other computation
 * it can also have a minimum altitude smaller than the std one.
 */


class Chem
{
	protected:
		/// Pointer to the planet: to compute more densities for the chemistry
		Planete* mpPlanet;

		/// The Xml parameter, needed for ReadParameters()
		XmlParameters* mpParam;
		/// The altitude grid in Km
		ublas::vector<double> mAltGridKm;


		/// List of the species defined by Chem
		std::deque<Specie*> mDefinedHere;


		/// List of pointers to the chemical densities species, used for chemistry (some are from the vNeutralDenscm_3, some are read from the xml file, some are added from the planet section.
		std::deque<Specie*> mChemDensComputcm_3;
		/// List of pointers to the total production species
		std::deque<Specie*> mTotProdcm_3s_1;
		/// List of pointers to the photo production species
		std::deque<Specie*> mPhotProdcm_3s_1;
		/// List of pointers to the elec production species
		std::deque<Specie*> mElecProdcm_3s_1;
		/// List of pointers to the proton production species
		std::deque<Specie*> mProtProdcm_3s_1;
		/// List of pointers to the cosmic production species
		std::deque<Specie*> mCosmoProdcm_3s_1;
		/// Neutral temperatures
		ublas::vector<double> mTempNeutreK;
		/// Electron temperatures
		ublas::vector<double> mTempElecK;
		/// Ion temperatures
		ublas::vector<double> mTempIonK;
		/// Electron density
		ublas::vector<double> mElecDenscm_3;
		/// The list of chemical reactions
		std::deque< boost::shared_ptr<ChemReact> >  mChemList;

		/// The list of photoequilibrium density classes
		std::deque< boost::shared_ptr<CalcDensPhEq> >  mCalcPhEqList;
		/**
		 * To set up the different densities used
		 * \param vNeutralDenscm_3 the densities in the ionization computation
		 */
		void InitDensities(std::deque<Specie*> vNeutralDenscm_3);
		/**
		 * To set up the different chemical reactions
		 */
		void InitChemReact();
		/**
		 * To set up the different photochemical density calculations
		 */
		void InitCalcDensPhEq();

		/**
		 * To init the densities  (inside xml; old and new mixing)
		 */
		void InitNewDensities();


		/**
		 * To init the densities declared inside the xml file
		 */
		void ReadDensities();

		std::map<std::string,std::string> mDensityWarnings;



		/// The electron density multiplicator: can be used to modify a main ion density in relationship with the modification of the electron density
		double mElecDensMult;

		/**
		 * Get the multiplicative factor of the species defined in the chemistry
		 * it can take into account the existence of the electron density multiplicator to make a modification proportional to this one
		 * (mainly used when we have a species which is proportional to the e dens, eg O+ or O2+ in terristrial planets
		 * \param vName : the name of the species which could be multiplied
		 */
		double GetMultiplicator(std::string vName);



	public:

		/**
		 * The constructor
		 * \param vpParam : pointeur to the xmlparameters object for reading the reaction rate modification of the error...
		 * \param vpPlanet : pointer to the planet (used for adding more data about the ion densities...)
		 * \param vAltGridKm : the altitude grid in km
		 * \param vNeutralDenscm_3  List of pointers to the neutral densities species
		 * \param vTotProdcm_3s_1  List of pointers to the total production species
		 * \param vPhotProdcm_3s_1  List of pointers to the photo production species
		 * \param vElecProdcm_3s_1  List of pointers to the elec production species
		 * \param vProtProdcm_3s_1  List of pointers to the proton production species
		 * \param vCosmoProdcm_3s_1  List of pointers to the cosmic production species
		 * \param vTempNeutreK  Neutral temperatures
		 * \param vTempElecK  Electron temperatures
		 * \param vTempIonK Ion temperatures
		 * \param vElecDenscm_3 Electron density
		 * \param vElecDensMult : the multiplicator of the electron density, used to correct possible additional parameters here

*/
		Chem(XmlParameters* vpParam,
				Planete* vpPlanet,
				ublas::vector<double> vAltGridKm,
				std::deque<Specie*> vNeutralDenscm_3,
				std::deque<Specie*> vTotProdcm_3s_1,
				std::deque<Specie*> vPhotProdcm_3s_1,
				std::deque<Specie*> vElecProdcm_3s_1,
				std::deque<Specie*> vProtProdcm_3s_1,
				std::deque<Specie*> vCosmoProdcm_3s_1,
				ublas::vector<double> vTempNeutreK,
				ublas::vector<double> vTempElecK,
				ublas::vector<double> vTempIonK,
				ublas::vector<double> vElecDenscm_3,
				double vElecDensMult);

		/// Destroy the new species
		~Chem();

		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * Some extra storages have been prepared for the user, like rSubResults where the VER can be put...
		 *
		 * \param vId : the SpecieId of the density we want
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */
		ublas::vector<double> GetDensity(SpecieId vId,ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string& rWarnings);



		/**
		 * Prints the selected densities
		 */
		void SelectedPrintDensities(std::deque<SpecieId>& rStates, ublas::vector<double>  vAltKmGrid, std::string rFilename);

		/**
		 * Prints the list of chemical reactions into a file
		 * \param vFilename : the file where we write the chemical reaction
		 *
		 */
		void PrintChemList(std::string vFilename);


		/**
		 * Check if the chemical reactions have the good Id
		 * \return false if there was an error
		 */
		bool CheckReacId();

		/**
		 * Get the reaction rate
		 * \param vId : the id of the position of the reaction
		 * \return the reaction rate
		 */
		ublas::vector<double> GetRate(unsigned vId)
		{
			return mChemList[vId]->GetReactionRate(mTempElecK,mTempIonK,mTempNeutreK);
		}
		/**
		 * Get the reaction efficiency
		 * \param vId : the id of the position of the reaction
		 * \param vName : the name of the product
		 * \param vState : the state of the product
		 * \return the efficiency rate
		 */
		double GetEfficiency(unsigned vId,std::string vName,std::string vState)
		{
			return mChemList[vId]->GetEfficiency(vName,vState);
		}

		/**
		 * Prints the list of density computation warnings into a file
		 * \param vFilename : the name of the file to write
		 */

		void PrintDensWarnings(std::string vFilename);

		/**
		 * Prints the atmosphere used in the chemistry model
		 * \param vFilename : the name of the file to write
		 */

		void PrintChemDensModels(std::string vFilename);




		/**
		 * Returns the density of your species.
		 * It also put the density of your species in the chemistry if the main species exists.
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the density
		 */
		ublas::vector<double> GetDens(std::string vSpecieName,std::string vSpecieState);


		/**
		 * Put the density of your species-state if the species exists in the chemistry list.
		 * The "" state should not be used
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \param vDensitycm_3 : the density to add
		 */
		void PutDens(std::string vSpecieName,std::string vSpecieState,ublas::vector<double> vDensitycm_3);
		/**
		 * Returns the total production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production
		 */
		ublas::vector<double> GetTotProd(std::string vSpecieName,std::string vSpecieState);
		/**
		 * Returns the Photon production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production
		 */
		ublas::vector<double> GetPhotProd(std::string vSpecieName,std::string vSpecieState);
		/**
		 * Returns the Electron production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production
		 */
		ublas::vector<double> GetElecProd(std::string vSpecieName,std::string vSpecieState);
		/**
		 * Returns the Proton production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production
		 */
		ublas::vector<double> GetProtProd(std::string vSpecieName,std::string vSpecieState);
		/**
		 * Returns the cosmic production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production
		 */
		ublas::vector<double> GetCosmoProd(std::string vSpecieName,std::string vSpecieState);

		/**
		 * Reset the electron density
		 * \param vElecDenscm_3: the electron density in cm_3
		 */
		void ResetEdens(ublas::vector<double> vElecDenscm_3)
		{
			mElecDenscm_3=vElecDenscm_3;
		}
/*
		void ResetDensity(std::deque<Specie*> vNeutralDenscm_3)
		{
			InitDensities(vNeutralDenscm_3);
		}
		void ResetProduction(std::deque<Specie*> vTotProdcm_3s_1,
				std::deque<Specie*> vPhotProdcm_3s_1,
				std::deque<Specie*> vElecProdcm_3s_1,
				std::deque<Specie*> vProtProdcm_3s_1,
				std::deque<Specie*> vCosmoProdcm_3s_1
				)
		{
			mTotProdcm_3s_1=(vTotProdcm_3s_1);
			mPhotProdcm_3s_1=(vPhotProdcm_3s_1);
			mElecProdcm_3s_1=(vElecProdcm_3s_1);
			mProtProdcm_3s_1=(vProtProdcm_3s_1);
			mCosmoProdcm_3s_1=(vCosmoProdcm_3s_1);
		}
		*/
};



#endif

