/**
 * \defgroup Chem Chemistry
 * \file reaction.hpp
 * \brief Defines the chemical reaction class. A list of this class is used in the chemistry to compute the densities...
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: reaction.hpp 1106 2010-08-04 18:16:56Z gronoff $
 */









#ifndef REACTION_CHEM_HPP
#define REACTION_CHEM_HPP

#include <planet/allplanets.hpp>

/**
 * \ingroup Chem
 * Defines the chemical reactions. This class, meant to be derived, defines the basis of the chemical reactions.
 * It allows to modify the main value of the chemical reaction (for the derived class!), and its uncertainty, through the configuration file.
 * One of the advantage of this definition is the fact that if the chemical reaction is used several times, the study of the uncertainty  will be paussible because in the different sub-utilisationm, the same constant will be applied.
 * In this class, we also store the Einstein coefficients for the emission.
 */
class ChemReact
{
	protected:


		/** Reads the parameter in the xml file, to check if we must use 
		 * the modified cross section, or if we must use the monte carlo system.
		 */
		void ReadParameters();

		/// The Xml parameter, needed for ReadParameters()
		XmlParameters* mpParam;

		/// The main reaction rate. This value can be totally modified at the beginning, by overload, or by MC modification for sensitivity studies.
		double mMainValue;

		/// The main reaction rate, unmodified by the montecarlo
		double mOrigMainValue;

		/// The uncertainty, in the for of string: allows to put a % or not
		std::string mUncertainty;

		/// Complementary information about the reaction
		std::string mSupplInfo;
		/// Abstract function to set the id of the reaction computed by this (inherited) class
	public:


		/// The Id of the reaction (must match the position in the Chem vector)
		unsigned mId;

		/// Returns the id \return id of the chemical reaction
		unsigned GetId()
		{
			return mId;
		}

		double ReturnMainValue()
		{
			return mMainValue;
		}

		/// True if the monte carlo system is activated
		bool mbIsMC;

		/// If the system is used for a einstein coefficient value
		bool mbIsEmit;
		/// The corresponding frequencies in that case
		double mEmitFreqnm;

		/// List of the reacting species
		std::deque< SpecieId > mReactant;
		/// List of the catalysers
		std::deque< SpecieId > mCatalys;
		/// List of the products
		std::deque< SpecieId > mProducts;
		/// Efficiency of the products production
		std::map< SpecieId, double > mEfficiency;

		/**
		 * Returns the efficiency of the reaction for the yield of the 
		 * species-state.
		 * \param vName : the name of the species
		 * \param vState : the state of the species
		 * \return the efficiency
		 *
		 */
		double GetEfficiency(std::string vName,std::string vState);



		/// Return the main informations of the reaction.
		std::string GetInfo();

		/**
		 * The constructor
		 * \param vpParam : pointeur to the xmlparameters object for reading the reaction rate modification of the error...
		 * \param vId : the identifiant of the reaction
		 */
		ChemReact(XmlParameters* vpParam,unsigned vId): mpParam(vpParam),mId(vId)
		{
			mbIsMC=false;
			mbIsEmit=false;
	
		}


		/// The destructor
		virtual ~ChemReact();

		/** Virtual function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 * Nb: the implementation of this GetValue in the inherited class must start by something like assert(vTe.size()==vTi.size); assert(vTe.size()==vTn.size()); so that there is no error  of unused variables when the reaction rate is constant
		 */
		virtual ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn)=0;



};




/**
 * \ingroup Chem
 * Defines the  class for the computation of the photochemical equilibrium densities
 * It is an abstract class : it must be inherited to give the density.
 * Some useful function to extend it is to use the 
 * SpecieUtils::GetSpecieDensity and SpecieUtils::GetSpecieProduction functions.
 */

class CalcDensPhEq
{
	protected:
		/// List of pointers to the neutral densities species
		std::deque<Specie*> mNeutralDenscm_3;
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
		ublas::vector<double>* mpElecDenscm_3;

		/// Abstract function to set the id of the density computed by this (inherited) class
		/// The id of the density computed
		SpecieId mId;
		/// The list of chemical reactions
		std::deque< boost::shared_ptr<ChemReact> >  mChemList;

		/// The size of the arrays.  Necessary to get the productions and densities when they are not defined (returns a vector of mSize, filled with 0)
		unsigned mSize;
		/**
		 * Returns the density of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the density (with coherent size)
		 */
		ublas::vector<double> GetDens(std::string vSpecieName,std::string vSpecieState);

		/**
		 * Returns the total production of your species
		 * \param vSpecieName : the name of the species you are searching for
		 * \param vSpecieState : the state of the species you are searching for
		 * \return the production (with coherent size)
		 */
		ublas::vector<double> GetTotProd(std::string vSpecieName,std::string vSpecieState);


	public:
		/// Constructor: it calls the Id creator (abstract)
		CalcDensPhEq(SpecieId vId):mId(vId)
	{
		mSize=0;
		// The set Id will be called for every inherited functions
	}
		/// The destructor
		virtual ~CalcDensPhEq(){}
		/// Returns the Id
		SpecieId GetId()
		{
			return mId;
		}
		/**
		 * Initialize the object: it is not inside the constructor because 
		 * there are a lot of parameters. Moreover, as these reactions
		 * will be created inside the Chem object, it will be able to 
		 * call it one by one.
		 *
		 * \param vNeutralDenscm_3  List of pointers to the neutral densities species
		 * \param vTotProdcm_3s_1  List of pointers to the total production species
		 * \param vPhotProdcm_3s_1  List of pointers to the photo production species
		 * \param vElecProdcm_3s_1  List of pointers to the elec production species
		 * \param vProtProdcm_3s_1  List of pointers to the proton production species
		 * \param vCosmoProdcm_3s_1  List of pointers to the cosmic production species
		 * \param vTempNeutreK  Neutral temperatures
		 * \param vTempElecK  Electron temperatures
		 * \param vTempIonK Ion temperatures
		 * \param vpElecDenscm_3 Electron density
		 * \param vChemList The list of chemical reactions
		 */
		void Init( std::deque<Specie*> vNeutralDenscm_3,
				std::deque<Specie*> vTotProdcm_3s_1,
				std::deque<Specie*> vPhotProdcm_3s_1,
				std::deque<Specie*> vElecProdcm_3s_1,
				std::deque<Specie*> vProtProdcm_3s_1,
				std::deque<Specie*> vCosmoProdcm_3s_1,
				ublas::vector<double> vTempNeutreK,
				ublas::vector<double> vTempElecK,
				ublas::vector<double> vTempIonK,
				ublas::vector<double>* vpElecDenscm_3,
				std::deque< boost::shared_ptr<ChemReact> >  vChemList
			 );
		/**
		 * The main function : compute the density of the element thanks to the different parameters (productions, densities, chemical reactions...).
		 * Some extra storages have been prepared for the user:
		 * \param rSubResults: matrix containing some extra information (like the influence of a specific process on the total density). The data in that matrix depends on what the creator of the inherited class wants
		 * \param rSubResultsNames : names for the matrix.
		 * \param rWarnings : warnings (as non existing neutral species...) 
		 * \return the density as a ublas::vector
		 */
		virtual ublas::vector<double> GetDensity(ublas::matrix<double>& rSubResults, std::deque<std::string>& rSubResultsNames,std::string rWarnings)=0;

};




#endif
