/**
 * \file emission.hpp
 * \brief Defines the emission object : it stores the different emits, calls the geometry system, creates the synthetic spectrum 
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emission.hpp 1191 2010-11-23 21:22:19Z gronoff $
 * \ingroup Emit
 */

#ifndef EMISSION_HPP
#define EMISSION_HPP
#include "emitlist.hpp"

/**
 * Store the emissions, and read the configuration file
 * \ingroup Emit
 *
 */
class Emission
{
	private:

		/// The Xml parameter
		XmlParameters* mpParam;
		/// The geometry object,  needed for limb observations
		boost::shared_ptr<GeoObs> mPath;

		/// The chemistry
		boost::shared_ptr<Chem> mChem;
		
		/// Pointer to the planet: to compute more densities for the chemistry
		Planete* mpPlanet;
		/// The altitude grid in Km
		ublas::vector<double> mAltGridKm;
		
		/// The Emit list
		std::deque< boost::shared_ptr<Emit> > mEmitList;


		/// List of filenames for the column emission
		std::deque< std::string >  mColFileName;

		/// List of filenames for the profile emission
		std::deque< std::string >  mProfFileName;

		/// List of filenames for the extra informations
		std::deque< std::string >  mExtraFileName;

		/// List of bool for using absorption
		std::deque< bool >  mUseAbsorptionList;

		/**
		 * Deque containing the absorbing species, the frequency, and the absorption value for the specific frequency
		 */
		std::deque< std::map< std::string, std::map<double,double> > > mAbsorptionValscm_2;


		/// List of id or the profile emission
		std::deque< bool >  mbCompute;
		
		/// Initialization of the Path
		void InitPath();
		/// Initialization of the emission list
		void InitEmission();

		bool mbComputeAll;///< if everything is computed in multiplecompute

		/// The integrated against the wavelength, in R
		std::map<double,double> mEmitSpectrumIntegratedR;
	public:


		/**
		 * The constructor
		 * \param vpParam pointer to the xml file
		 * \param vChem pointer to the chemical reactions
		 * \param vpPlanet pointer to the planet
		 * \param vAltGridKm the altitude grid
		 *
		 */
		Emission(XmlParameters* vpParam,boost::shared_ptr<Chem> vChem,
				Planete* vpPlanet,
				ublas::vector<double> vAltGridKm
				);

		/**
		 * Computes the emissions, and write the corresponding files
		 * \param vSuffix the suffix for the output files.
		 */
		void ComputeEmissions(std::string vSuffix="");


		/**
		 * Init for multiple call of emissions : basically init the system for computeemissions, but 
		 * without print
		 */
		void InitMultipleEmit();
		/**
		 *
		 * Compute one call for multiple call compatible emissions
		 */
		void CompMultipleEmit();



		/// Check if the emission list is integrated correctly
		bool CheckEmitList();

		/**
		 * To print the entire spectrum computed in a file
		 * \param vFilename : the file where to write the spectrum
		 */
		void PrintEntireSpectrum(std::string vFilename);

		/**
		 * To print the informations on the different emissions
		 * \param vFilename : the file where to write 
		 */
		void PrintInfo(std::string vFilename);


		/**
		 * Returns the VER of the emitting specie vSpec, at the position vPos
		 * \param vSpec : id of the emitting species
		 * \param vPos : position of the emission in the emission list of the species. If <0, it is the total emission
		 */
		ublas::vector<double> GetVER(SpecieId vSpec,int vPos);

		/**
		 * Returns the limb emission of the emitting specie vSpec, at the position vPos
		 * \param vSpec : id of the emitting species
		 * \param vPos : position of the emission in the emission list of the species. If <0, it is the total emission
		 */
		ublas::vector<double> GetLimb(SpecieId vSpec,int vPos);

		/**
		 * Returns the variable for the limb observations
		 */
		ublas::vector<double> GetLimbParameter()
		{
			std::deque<double> parameter;
			if(mPath->IsTanDefined())
				parameter=mPath->GetTanList();
			if(mPath->IsDecDefined())
				parameter=mPath->GetDecList();
			ublas::vector<double> resu(parameter.size());
			for(size_t i=0;i<parameter.size();++i)
				resu[i]=parameter.at(i);
			return resu;
		}
};

#endif
