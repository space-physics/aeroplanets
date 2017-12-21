/**
 * \file electronionization.hpp
 * \brief Definition for the electron impact ionization and excitation,
 *	  and electron transport classe.
 *	  Depends on some fortran functions...
 *	  Copyright G Gronoff Sept 2009
 *	  Last Modification : $Id: electronionization.hpp 1899 2014-01-16 20:34:04Z gronoff $
 */



#ifndef ELECTRON_IONIZATION_HPP
#define ELECTRON_IONIZATION_HPP
#include <species/species.hpp>
#include <eflux/eflux.hpp>
#include <protoncosmics/protocos.hpp>
#include <proton/proton.hpp>



// The fortran source

extern "C"
{

	/**
	 * The disort function : computes the  electron transfert
	 * \param vPosEner : position in the energy grid     <--
	 * \param vNbLayer : number of layers in the problem <--
	 * \param vDiffTau : diff in tau                     <--
	 * \param vSSALB : ???                               <--------
	 * \param vPmom : the phase moments -> matrix        <--------
	 * \param vUsrTau : logical, by default = 1
	 * \param vNbTauLayers : the number of layers in tau (2*nalt-1) <-
	 * \param vUTau : the tau diff                       <--------
	 * \param vNbStr : the number of angles              <--
	 * \param vUsrAng : logical, by default 0
	 * \param vUmu : integer, by default 0
	 * \param vNumu : double, array of size (vNbStr) not used
	 * \param vNphi : integer, 0
	 * \param vPhi : double*, size 1, 0
	 * \param vFbeam : double 0
	 * \param vUmu0 : double 1
	 * \param vSrc : double* the photoelectron source     <------
	 * \param vSrcu : void, same size as vSrc
	 * \param vPhi0 : double, 1
	 * \param vPrecipFlux : the precipitating flux        <------
	 * \param vLamber : default 1
	 * \param vAlbedo : the albedo at the bottom          <------
	 * \param vHl : array double, 2*nbang+1        
	 * \param vBsrc : 0. double
	 * \param vTsrc : 0. double
	 * \param vLdeltam : integer, default
	 * \param  vExorc : integer, default 
	 * \param vOnlyfl : integer, default
	 * \param vAccur : double default 0.
	 * \param vPrtn : integer*7, default...
	 * \param vNaltm1 : integer nalt-1
	 * \param v2Naltm1 : integer 2*nalt-1
	 * \param vNAng1 : integer nbang
	 * \param vNAng2 : integer nbang
	 * \param vMaxphi : int 1
	 * \param vRfldir : double* 2*nalt-1
	 * \param rFlDown : the resulting down flux          -------->
	 * \param rFlUp : the resulting up flux              -------->
	 * \param vDfdt : double 2*nalt-1 0.
	 * \param rUavg : double 2*nalt-1                    -------->
	 * \param vUu   : double 0(2*nbrango2,2*nbralt-1,maxphi)
	 * \param rUou : (2*nbrango2,2*nbralt-1)	     -------->
	 * \param vNbAng_2 : nbrango2
	 * \param vLinear : the linearity
	 * \param vChecking : true (1) if disort is checking the entries.
	 */

	extern void disort_(int* vPosEner, int* vNbLayer,double* vDiffTau, double* vSSALB,double* vPmom,int* vUsrTau, int* vNbTauLayers,double* vUTau,int* vNbStr, int* vUsrAng, int* vUmu,double* vNumu,int* vNphi,double* vPhi,double* vFbeam,double* vUmu0,double* vSrc,double *vSrcu,double* vPhi0,double* vPrecipFlux,int* vLamber,double* vAlbedo,double* vHl,double* vBsrc,double* vTsrc, int* vLdeltam, int*vExorc,int* vOnlyfl,double* vAccur,int* vPrtn,int* vNaltm1,int*v2Naltm1,int*vNAng1,int*vNAng2,int* vMaxphi,double* vRfldir, double * rFlDown, double * rFlUp,double* vDfdt,double*rUavg,double* vUu, double* rUou,int* vNbAng_2,int* vLinear,int* vChecking);



};




/**
 * \ingroup electron_impact_ionization_process
 * \ingroup ionization_processes
 * Allows to compute the electro-impact ionization.
 */

/// Maximum energy for the porter approximation of the secondary electron phase angle
const double RUTHERFORD_PORTER_THRESHOLD=447.;

class ElectronImpactIonization
{
	protected:
		/// The parameters object
		XmlParameters* mpParameter;

		/// The flux object
		EFlux* mpElecFlux;

		/// the center grid for searching the position of the energy
		ublas::vector<double> mCenterGrideV;

		/// The center of the elec grid
		ublas::vector<double>* mpElecC;
		/// The bottom of the elec grid
		ublas::vector<double>* mpElecB;
		/// The width of the elec grid
		ublas::vector<double>* mpElecD;

		/// Direct link to the Gaussian angle. It is not defined in that function
		MathFunction::GaussianAngle* mpGAngle; 
		
		/** Returns the position of the energy in the grid
		 * This is exactly the same function as
		 * ElecCrossSection::PosEner() with a pre-selected posmin
		 * \param vEnergyeV : the energy
		 * \param vPosmin : the minimal position
		 * \return the position in the energy grid
		 */
		unsigned PosEner(double vEnergyeV,unsigned vPosmin=0);

		/**
		 * Compute the integration of the flux over its energy
		 * \param vEne1eV : the min energy for the integration
		 * \param vEne2eV : the max energy for the integration
		 * \param vFluxeV : the matrix of the flux vFluxeV(energy,angle)
		 * \param vPosAngle : against which angle we compute the flux
		 * \param rSum : the first result flux \f$ \sum vFlux(E) dE \f$
		 * \param rSumE : the second result flux \f$ \sum vFlux(E) E dE \f$
		 *
		 */
		void IntFlux(double vEne1eV,double vEne2eV,ublas::matrix<double> vFluxeV,unsigned vPosAngle,double& rSum, double& rSumE);

		/**
		 * Computes the integration of the flux over its energy. The parameter here is the position
		 * \param vPos1 : position of the min energy
		 * \param vPos2 : position of the max energy
		 * \param vFlux : the flux
		 * \param vPos3 : the position in the flux (angle or altitude)
		 * \param rSum : the first result flux \f$ \sum vFlux(E) dE \f$
		 * \param rSumE : the second result flux \f$ \sum vFlux(E) E dE \f$
		 *
		 */

		void IntFluxP(unsigned vPos1,unsigned vPos2, ublas::matrix<double> vFlux,unsigned vPos3,double& rSum, double& rSumE);	
//define ANALYSE_ELEC_FLUX	
#ifdef ANALYSE_ELEC_FLUX		
		/**
		 * Computes the integration of the flux over its energy. 
		 * But with carefully checking that every step is positive
		 * The parameter here is the position
		 * \param vPos1 : position of the min energy
		 * \param vPos2 : position of the max energy
		 * \param vFlux : the flux
		 * \param vPos3 : the position in the flux (angle or altitude)
		 * \param rSum : the first result flux \f$ \sum vFlux(E) dE \f$
		 * \param rSumE : the second result flux \f$ \sum vFlux(E) E dE \f$
		 *
		 */

		void IntFluxPa(unsigned vPos1,unsigned vPos2, ublas::matrix<double> vFlux,unsigned vPos3,double& rSum, double& rSumE);
		/// Test du flux en electrons
		double mTestFlux;
#endif
		/*
		 * Computes the reflected thermal energy
		 * \param vNtherm : vector containing the position of the crossover energy grid -thermal electrons distribution. 
		 * \return the reflected energy value
		 */
		double ReflectedEnergy(ublas::vector<unsigned> vNtherm);
		/*
		 * Computes the transmitted thermal energy
		 * \param vNtherm : vector containing the position of the crossover energy grid -thermal electrons distribution. 
		 * \return the transmitted energy value
		 */
		double TransmittedEnergy(ublas::vector<unsigned> vNtherm);


		/**
		 * Computes the energy absorbed by heating (coulombian interactions)
		 * \param vNtherm : vector containing the position of the crossover energy grid -thermal electrons distribution. 
		 * \param vQIntensity : the flux of electron
		 * \return the absorbed energy value
		 */

		double Heating(ublas::vector<unsigned> vNtherm,const ublas::vector< ublas::matrix<double> >& vQIntensity );


		/// The heating by electron coulombian collisions
		ublas::vector<double> mQe;
		/**
		 * Calls the disort function for the energy given by posener
		 * \param vPosEner : the position in the energy vector
		 * \param rQInt : the energy of the upper level
		 * \param rQIntensity : reference to the intensity : main thing computed by mstream
		 * not be used as a const...
		 */

		void Mstream(unsigned vPosEner,const ublas::matrix<double>& rQInt,ublas::vector< ublas::matrix<double> >& rQIntensity );

		/**
		 * Computes the rQint : energy for the next loop
		 * \param vPosEner : the position in the energy vector
		 * \param rQInt : reference to the energy of the upper level
		 * \param rQIntensity : reference to the intensity : main thing computed by mstream,used to compute qint
		 */
		void QMstream(unsigned vPosEner,ublas::matrix<double>& rQInt
				,ublas::vector< ublas::matrix<double> >& rQIntensity );



		/**
		 * Computes the ionization and excitation of the species. Returns the computation of the deposited elastic and inelastic energy
		 * \param vQIntensity : intensity of the electron flux (computed by mstream)
		 * \param vPositionThermic : position of the crossover energy : the separation between thermic electrons and suprathermic energy
		 * \param rInelasticDepeV  : reference to the third result : the deposition of energy into inelastic chocs.
		 * \param rEnrate : the inelastic deposition for the species at each altitudes
		 */
		void IonizeSpecies(ublas::vector< ublas::matrix<double> > vQIntensity,ublas::vector<unsigned> vPositionThermic,  ublas::vector<double>& rInelasticDepeV,ublas::matrix<double>& rEnrate);


		/** The albedo for the electrons at the bottom of the line
		 * Should not be necessary in theory
		 * in practice, for vertical precipitation 1.
		 *
		 * This parameter MUST BE 0 when you work with BENT lines 
		 * which GOES OUTSIDE the atmosphere (kind of U turn :-))
		 *
		 */
		double mAlbedo;




		// Electron density
		ublas::vector<double>* mpElectronDensitycm_3;
		/// Electron temperature pointer
		ublas::vector<double>* mpElectronTemperatureK;

		/**
		 * Compute the electron loss function, called by Init.
		 *(continuous slowing down approximation) due to 
		 * electron-electron and Coulomb interaction, (See Stamnes,
		 * Rees: Heating of thermal ionospheric electrons...., JRL,
		 * 10, 309-312,1983) and the backscatter ratio for e-e interaction                                                
		 * The variable mEnergyLosseV contains the loss function a:
		 *  ([e]=e-density)  ELOSSE = [e] * L(E)                   
		 *
		 * \warning This function returns an energy loss of 0 if E>150eV
		 * moreover, if we are at Thermal energy -> 0! T=k*E
		 */
		void ComputeEnergyLoss();

		/**
		 *
		 * Compute the phase functions for all altitudes 
		 * and energies.
		 * The outer loop is in energy, because when we do not call
		 * porter, the function is the same for all energies!
		 * 
		 * \warning The output is in [energy][angle][layer]
		 * (where nblayer=nbalt-1 the altitude at the top is 
		 * not taken into account).
		 * This particularity is to optimize the transformation
		 * into pmp(angle,layer), which is used by disort!
		 *
		 */
		void ComputePhaseFunction();


		/**
		 * The Phase function in [energy][angle][layer]
		 * Modification 9 oct 09 : energy layer angle
		 * for matrix usage and fortran call (inversion row column
		 * for the matrices!!!)
		 * for optimization of Disort!
		 *
		 */
		//	std::vector< std::vector< std::vector<double> > > mPhaseFunction;
		ublas::vector< ublas::matrix<double> > mPhaseFunction;

		/** The energy loss computation result
		 * mEnergyLosseV[alt][energy]
		 */
		ublas::matrix<double> mEnergyLosseV;


		/// To give the specie vector to other functions.
		std::deque<Specie*>* mpSp;


		/**
		 * We initialize the Rutherford function.
		 *
		 */

		void InitRuther();

		// The Rutherford matrix in terms of energy,angle
		//std::vector< std::vector<double> > mRutherford;
		ublas::matrix<double> mRutherford;

		/**
		 * We initialize the Porter Matrix.
		 * To be used in real life, this matrix has to be pondered by the species densities!
		 */
		void InitPorter();
		/// The porter vector, in term of specie, energy, angle
		ublas::vector< ublas::matrix<double> > mPorter;

		/**
		 * Returns the Porter value for energy and altitude
		 * \param vAlt : the altitude number
		 * \param vEner : the energy number
		 */

		ublas::vector<double> Porter(unsigned vAlt,unsigned vEner);

		/**
		 * Subroutine for Porter and Rutherford
		 * \param vNbAngle : the number of angles + 1 (0 is used to check the isotropy...)
		 * \param vAlpha : the parameter of Porter et al
		 * \param vBeta: the parameter of Porter et al
		 * \param vGamma: the parameter of Porter et al
		 * \param vBG :  false if we call from Rutherford
		 */
		ublas::vector<double> RecUp(unsigned vNbAngle,double vAlpha,double vBeta, double vGamma,bool vBG=true);

		/**
		 * Subroutine for Porter and Rutherford
		 * \param vNbAngle : the number of angles + 1 (0 is used to check the isotropy...)
		 * \param vAlpha : the parameter of Porter et al
		 * \param vBeta: the parameter of Porter et al
		 * \param vGamma: the parameter of Porter et al
		 * \param vBG :  false if we call from Rutherford
		 */
		ublas::vector<double> RecDown(unsigned vNbAngle,double vAlpha,double vBeta, double vGamma,bool vBG=true);

		/*       _\|/_
			 (o o)
			 +----oOO-{_}-OOo---------------------+
			 |                                    |
			 |                                    |
			 |Necessary arrays for the energy loop|
			 |                                    |
			 |                                    |
			 +-----------------------------------*/


		/// Loss by cm! mCTot(Energy,altitude)
		ublas::matrix<double> mCTot;
		//	std::vector<double> mCTot;

		/**
		 * Fills the mCTot vector
		 */
		void ComputeCTot();


		/** The vectors for working with PWOM 
		 */
		ublas::vector<double> mNumberDensity1, mNumberDensity2, mNumberFlux1, mNumberFlux2;

		/**
		 * Returns the inverse of the sinus of the magnetic dip angle (vertical = 1)
		 * Not bent safe if a magnetic dip angle is defined, but why defining it when bent???
		 */
		ublas::vector<double> ReturnInverseSinMagneticDipAngle();

		/// The standard altitude grid
		ublas::vector<double>* mpAltGridKm;

		/// The vector of inverse of sinus magnetic dip angle
		ublas::vector<double> mIsmgdpa;

		/// number of altitudes 
		unsigned mNbalt;
		/// number of energies
		unsigned mNben;
		/// number of angles
		unsigned mNbang;
		/// number of angles divided by 2
		unsigned mNbang2;

		/// number of species
		unsigned mNbsp;


		/// True if disort is checking
		int mIsDisortChecking;


		/// Factor to modify the L(E) main parameter. It is mainly used to compute uncertainties on that parameter.
		double mFactUncertLE;

		// HERE, we define 3 parameters to check in depth the energy conservation
		// through the redistribution. We will see the total influence of the 
		// secondary electrons, the primary degraded electrons, and the
		// electrons wich undergo a coulombian interation
		// Please have in mind that a degraded electron car undergo a coulombian
		//  interaction and then ionize. The idea here is to have an idea if
		//  there is a problem of quantity when the energy conservation explodes

		/// The total energy in the secondary
		double mSecondarySumeV;
		/// The total energy in the degraded
		double mDegradedSumeV;
		/// The total energy in the coulombian degraded
		double mCoulombianSumeV;



	public:
		/**
		 * Initializes the class
		 * \param pParam : the object parameter
		 * \param vpElectronDensitycm_3 : the electron density in cm_3
		 * \param vpElectronTemperatureK : the electron temperature in K
		 * \param vpAltGridKm : the altitude grid in km
		 */
		ElectronImpactIonization(XmlParameters* pParam,ublas::vector<double>* vpElectronDensitycm_3,ublas::vector<double>* vpElectronTemperatureK,ublas::vector<double>* vpAltGridKm);



		/**
		 * Init the different parameter in the vertical approximation
		 */

		void Init();
		/**
		 * Compute the electron impact:
		 * performs the transport of electrons
		 * and computes the resulting species.
		 * \param vSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section \warning HERE it is a vector of species : too difficult, an not effective to make a pointer of pointer when there is so few things in the vector (the number of species is a small thing! I think this statement can be valid even if hundreds of species).
		 * \param vFlux: The input flux parameter.
		 * \param rResult : the result vector of species -> were the computations are written
		 * Note: to avoid a lot of difficulties, the rResult is also a reference instead of a pointer
		 */

		void ComputeElectronImpact(std::deque<Specie*> vSp,
				const EFlux& vFlux,
				std::deque<Specie*>& rResult);




		/*       _\|/_
			 (o o)
			 +----oOO-{_}-OOo---------------------------------+
			 |                                                |
			 |These internal vector and matrices are in public|
			 |        because they can be interesting         |
			 |                                                |
			 +-----------------------------------------------*/

		/// Layer diff of tau (collision depth) nalt
		ublas::matrix<double> mDiffTau;
		/// Differential in tau (collision depth) 2*nalt-1
		ublas::matrix<double> mUTau;
		/// Phase function (no moments) : the absorption
		ublas::matrix<double> mOmegaFunction;
		/// The SSALB parameter it corresponds to the albedo of single scattering/ elastic scattering for electron impact.
		ublas::matrix<double> mSsalb;


		/*       _\|/_
			 (o o)
			 +----oOO-{_}-OOo-----+
			 |                    |
			 |The different fluxes|
			 |                    |
			 +-------------------*/


		/// The downward hemispheric flux
		ublas::matrix<double> mFluxHemisphericDown;
		/// The upward hemispheric flux
		ublas::matrix<double> mFluxHemisphericUp;
		/// The total hemispheric flux
		ublas::matrix<double> mFluxHemisphericTot;
		/// The net hemispheric flux
		ublas::matrix<double> mFluxHemisphericNet;

		/// The total electron energy interaction (flux input)
		double mTotalEnergy;
		/// The reflected electron energy
		double mReflectedEnergy;
		/// The transmitted electron energy
		double mTransmittedEnergy;

		/// The absorbed energy by heating
		double mAbsorbedEnergy;
		/// The absorbed energy by inelastic deposition
		double mInelasticEnergy;
		/// The resullt total energy (without reflected)
		double mResultTotalEnergy;

		/// Energy per elec production : relative to computed energy
		double mElos1;
		/// Energy per elec production : relative to total (input) energy - reflected
		double mElos2;
		/// Energy per elec production : relative to total (input) energy
		double mElos3;

		/// Energy per ion production : relative to computed energy
		double mIlos1;
		/// Energy per ion production : relative to total (input) energy - reflected
		double mIlos2;
		/// Energy per ion production : relative to total (input) energy
		double mIlos3;

		/// The relative error between computed energy and input energy
		double mError;
		/// Total electron production
		double mProdTot;
		/// Total ion production (can be different than electron production due to double ionisation processes)
		double mProdIonTot;

		/// Compute the energy total deposition, the relative error!
		void ComputeDeposition();



		/**
		 * To print the qintensity flux
		 * \param vFilename : the filename
		 * \param qintensity : the qintensity matrix
		 */
		void PrintQintensity(std::string vFilename,ublas::vector< ublas::matrix<double> > qintensity);


		/**
		 * To print the flux at the selected altitudes
		 * \param vFilename : the filename
		 * \param vAlts : vector of altitudes where to compute the flux
		 * \param option : option for the flux :
		 *         			- 0 : flux down
		 *         			- 1 : flux up
		 *         			- 2 : flux total
		 *         			- 3 : flux net
		 */
		void PrintFluxes(std::string vFilename,ublas::vector<double> vAlts,unsigned option);

		/**
		 * To print the flux at the selected altitudes
		 * \param vFilename : the filename
		 * \param vPath : the path, to get the correspondance between altitudes and lengths
		 * \param vAkm : vector of altitudes where to compute the flux
		 * \param vLkm : vector of length where to compute the flux
		 * \param option : option for the flux :
		 *         			- 0 : flux down
		 *         			- 1 : flux up
		 *         			- 2 : flux total
		 *         			- 3 : flux net
		 */
		void PrintBentFluxes(std::string vFilename,const Path & vPath,ublas::vector<double> vAkm, ublas::vector<double> vLkm,unsigned option);
		       


		/**
		 * To print the heating by the electrons
		 *  \param vFilename : the filename
		 */
		void PrintElectronHeating(std::string vFilename);

		/**
		 * To print the PWOM outputs
		 *  \param vFilename : the filename
		 */
		void PrintPWOM(std::string vFilename);



		/**
		 * To print the mean pair production energy
		 *  \param vFilename : the filename
		 */
		void PrintEnergyPairs(std::string vFilename);

		/**
		 * To print the mean pair production energy for every sub species produced
		 * \param vSpecies : the produced species.
		 *  \param vFilename : the filename
		 */
		void PrintEnergyPerProduction(std::deque<Specie*> vSpecies, std::string vFilename);


		/**
		 * To print the conservation of the energy, in an xml format
		 *  \param vFilename : the filename
		 */
		void PrintEnergyConservation(std::string vFilename);
};







#endif
