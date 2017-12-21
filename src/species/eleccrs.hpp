 /** 
 * \file eleccrs.hpp
 * \brief Defines the electron cross section class
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: eleccrs.hpp 1336 2011-11-08 22:11:28Z gronoff $
 *
 */



/**
 * We inherit from Cross section, to make the redistribution!
 * In fact, the electron cross section is much more problematic
 * than the photon one:
 * - the electron cross section is extrapolated to quite high energies
 *   this can be dangerous when we are outside the region where 
 *   this physic is valid
 * - there is a separation between elastic and inelastic cross sections
 * - the original Trans* code used to add the inelastic cross sections
 *   in our code, we will create the possibility to put directly
 *   a total inelastic cross section for the electrons.
 *   The problem with that approach, is that we lose some
 *   points concerning the threshold energy. 
 *   I (GG 19 sept 09) think that it can be solved thank to
 *   a mean threshold energy, which is close to the main ionization energy.
 *
 *
 */

#ifndef ELEC_CRS_HPP
#define ELEC_CRS_HPP
#include "crs.hpp"
#include "documentation.hpp"


/**
 * \ingroup specie_parameters
 * Class to extend the CrossSection class to the electrons
 *
 * Thanks to that class, we will be able to compute
 * the redistribution factors, needed to compute 
 * the multi-stream electron impact 
 *
 */
class ElecCrossSection:public CrossSection
{
	private:
		/// If the ionization are defined with respect to the position: true (false if there is one ionization cross section defined by  /IonizationCrs)
		bool mIonizationPosition;

		/// the center grid for searching the position of the energy
		ublas::vector<double> mCenterGrideV;

		/// the number of energies
		unsigned mNben;
		/**
		 * Returns the position of vEnereV in the energy grid
		 * with a minimum position (maximum energy) at vPosMin
		 * \param vEnereV : the energy considered
		 * \param vPosMin : the minimum position (max energy)
		 * \return the position of vEnereV  in the grid
		 *
		 *
		 */
		unsigned PosEner(double vEnereV,unsigned vPosMin);



		/**
		 * Average for the energy of the primary
		 * \warning if you understand something, please mail me guillaume.gronoff@free.fr
		 * \param vN the pointer to the target energy
		 * \param vNp the pointer to the primary energy
		 * \param vW the Threshold corrected
		 * \param vEmaxeV the maximum energy
		 */
		double AverageEnergySecondary(unsigned vN, unsigned vNp, double vW, double vEmaxeV);
		/**
		 * Average for the energy of the secondary
		 * \warning if you understand something, please mail me guillaume.gronoff@free.fr
		 * This computes the total cross section correction for 
		 * the electron mean energy loss...
		 * \param vN the pointer to the target energy
		 * \param vNp the pointer to the primary energy
		 * \param vW the Threshold corrected
		 * \param vEmaxeV the maximum energy
		 */

		double AverageEnergyPrimary(unsigned vN, unsigned vNp, double vW, double vEmaxeV);
		/**
		 * Normalized energy for the secondary electron
		 * Computes the differential oscillator strenght (Rees 69)
		 * normalized.
		 * It is integrated over the energies vEmineV and vEmaxeV
		 *
		 * \warning if you better understand this thing, please mail me guillaume.gronoff@free.fr
		 *
		 * \param  vEmineV : 
		 * \param  vEmaxeV : 
		 * \param  vYy : Energy loss of the primary
		 * \param  vZ : Energy of the secondary electron (Es)
		 */

		double NormalizedSecondaryElectron(double vEmineV,double vEmaxeV,double vYy,double vZ);
		/**
		 * Computes the secondary production energies with the rees method.
		 * In fact this subfunction computes
		 * f'/A*ln(...) (with Rees 69 notation)
		 * \param  vX : The threshold in eV
		 * \param  vY : Energy of the incident  (primary...) electron
		 * \param  vZ : energy of the secondary electron in eV
		 *
		 * \bug Possible error there : e_limit for computation
		 * of secondary at 5600 eV!!!!
		 *
		 */	
		double SecondaryElectronFunction(double vX,double vY,double vZ);




	protected:


		// Maximum energy defined for elastic process
		double mMaxDefinedEnergyElasticCrseV;

		/**
		 * Extrapolates the cross section (rCrscm2) for the high energies
		 * when the maximum cross section is defined for vEnerMaxeV.
		 * This one is specific for the elastic cross section
		 * Nb: rCrscm2 is a reference
		 * \param vEnerMaxeV : the max energy where the cross section
		 * 			is defined in the data
		 * \param rCrscm2 : the cross section vector, which needs extrpolation
		 * 			it is a reference
		 */


		void ExtrapolateElastic(double vEnerMaxeV,ublas::vector<double>& rCrscm2);


		/**
		 * Do the redistribution: computes the redistribution factors for the ionization process(es).
		 * \param vIonProcessCrscm2 : cross section of the process(es) that ionize. They can be one (for example: total ionization), or many.
		 * \param vIonProcessThresholdeV : ionization threshold
		 *
		 * \warning : the main problem when we work with doubly charged ions, it the fact that 2 electrons are not created here.
		 * It is not really a problem now.
		 *
		 */
		void RedistributeIonization(std::deque< ublas::vector<double>* >  vIonProcessCrscm2,std::deque<double> vIonProcessThresholdeV);


		/**
		 * Do the redistribution: computes the redistribution factors for the excitation process(es).
		 * \param vExcProcessCrscm2 : cross section of the process(es) that excite.
		 * \param vExcProcessThresholdeV : cross sections for the thresholds.
		 * \warning doesnt performs a modified adding in cin as the old version. It is in fact very unlikely that the energy width is larger than the threshold for these low energies! That must be the reason why the .lsec. is not added for that part. So the cin is only the sum, without modification of the different processes. For the C++ code, this is the default (only way 20 sept 09).
		 *
		 */
		void RedistributeExcitation(std::deque< ublas::vector<double>* >  vExcProcessCrscm2,std::deque<double> vExcProcessThresholdeV);


		/// The width grid
		ublas::vector<double> mGrEngddeV;
		/// When the excitation crs is specifically defined
		ublas::vector<double> mExcitationUniqueCrscm2;
		/// When the ionization crs is specifically defined
		ublas::vector<double> mIonizationUniqueCrscm2;



		/// If the Opal double differential cross section system is used. False by default
		bool mbIsOpal;
		
		/// The opal parameter Ebar
		double mEbareV;
		
		/// The Rees alpha parameter (see Gronoff et al. 2011) (orig = 1/31.5)
		double mReesAlpha;
		/// The Rees Beta parameter (see Gronoff et al. 2011) (orig = 339.)
		double mReesBeta;
		/// The Rees Gamma parameter (see Gronoff et al. 2011) (orig = 1/2.49)
		double mReesGamma;
		

	public:

		/** The redistribution matrix,
		* created by LoadCrs
		* and filled by RedistIonization
		* and RedistExcitation
		* Please note that mOmDegradcm2(n1,n2)
		* is the cross section for inelastic scattering
		* from n1 to n2
		* and mOmDegradcm2[n2][n1] is the cross section
		* for the production of a secondary with 
		* energy n2.
		* Therefore, this is important to consider
		* that ionization cross section is different
		* than excitation cross section.
		*/

		ublas::matrix<double>  mOmDegradcm2;
		
		/**
		 * Matrix to store the creation fluxes
		 * it is not used for the computation of the degradation in energy
		 * but is important to check for each species in comparison with the 
		 * double differential cross sections in the litterature
		 */
		ublas::matrix<double>  mOutSecondarycm2;


		/**
		 * The corrected cross section for ionization and
		 * excitation:
		 * this cross sections allows better computation
		 * of the secondary and degradation redistribution
		 * as explained in Swartz 85
		 */
		ublas::vector<double> mCrsCincm2;



         // Put inside the crs : for printing purposes
	//	// Cross section for the elastic processes
	//	ublas::vector<double> mElasticCrscm2;


		/**
		 * Initialize the Grid and the CrossSection class
		 * \param vGreV : the default energy grid
		 * \param vGrEngddeV : the energy grid width
		 */
		ElecCrossSection(ublas::vector<double> vGreV,ublas::vector<double> vGrEngddeV);

		/**
		 * Load the Cross section.
		 * Performs the electron interpolation
		 * do the redistribution
		 *
		 * All of that from the given object parameter.
		 *
		 * \param pParams : object to read the cross section file
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : true if MonteCarlo active
		 * \param bLogInterp : true if you want log interpolation, the default HERE
		 * \param vReesAlpha : value of the Alpha parameter for the Rees evaluation of the secondary electron cross section (doubly differential cross section). See G. Gronoff 2011/2 Uncertainty of ionization-airglow models. 
		 * \param vReesBeta : value of the Beta parameter for the Rees evaluation.
		 * \param vReesGamma : value of the Beta parameter for the Rees evaluation.
		 * \param bIsOpal : if the Opal parameterization is used instead of the Rees one (false by default)
		 * \param vOpalEbareV : the Opal parameter for the species in case of Opal.
		 * \return true if loaded.
		 * load the sub processes
		 * performs the interpolation on the grid
		 * be careful with negative interpolation. 
		 */
		bool LoadCrs(XmlParameters *pParams,std::string vName,bool bIsMonteCarlo,bool bLogInterp=true,double vReesAlpha = 1/31.5, double vReesBeta = 339., double vReesGamma =1/2.49, bool bIsOpal = false, double vOpalEbareV=0. );

		/**
		 * Load the Cross section.
		 * Performs the electron interpolation
		 * do the redistribution
		 *
		 * All of that from the given object file.
		 *
		 * \param vFilename:  the filename for the crs
		 * \param  vName : name of the species
		 * \param bIsMonteCarlo : true if MonteCarlo active
		 * \param bLogInterp : true if you want log interpolation, the default HERE
		 * \param vReesAlpha : value of the Alpha parameter for the Rees evaluation of the secondary electron cross section (doubly differential cross section). See G. Gronoff 2011/2 Uncertainty of ionization-airglow models. 
		 * \param vReesBeta : value of the Beta parameter for the Rees evaluation.
		 * \param vReesGamma : value of the Beta parameter for the Rees evaluation.
		 * \param bIsOpal : if the Opal parameterization is used instead of the Rees one (false by default)
		 * \param vOpalEbareV : the Opal parameter for the species in case of Opal.
		 * \return true if loaded.
		 * load the sub processes
		 * performs the interpolation on the grid
		 * be careful with negative interpolation. 
		 */
		bool LoadCrs(std::string vFilename,std::string vName,bool bIsMonteCarlo,bool bLogInterp=true,double vReesAlpha = 1/31.5, double vReesBeta = 339., double vReesGamma =1/2.49, bool bIsOpal = false, double vOpalEbareV=0. );

		/**
		 * Print the cross section with the new energy grid
		 * into a file. This file can be read with a python-pylab load function
		 * \param vFilename : name of the file. It is created.
		 * \param vRedisFilename : name of the file for the redistribution function
		 */

		void PrintCrs(std::string vFilename,std::string vRedisFilename);


		/**
		 * Allows to print necessary outputs for the tests
		 *
		 *
		 */
		void TestPrint();







};






#endif
