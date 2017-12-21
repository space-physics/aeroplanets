/**
 * \file photoionization.hpp
 * \brief The  definition of the photoionization class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: photoionization.hpp 1922 2014-02-26 22:38:52Z gronoff $
 */
 

/**
 * \defgroup photoionization_process Photoionization
 * \defgroup ionization_processes Particle impact processes
 *
 */

#ifndef PHOTOIONIZATION_HPP
#define PHOTOIONIZATION_HPP
#include "flux_model.hpp"
#include "lines.hpp"
/**
 * \ingroup photoionization_process
 * \ingroup ionization_processes 
 * Class Photoionization: do all the photoionization task:
 *  - set up the solar flux model (and computes the solar flux)
 *  - compute photoionization and photoexcitation
 *  - compute the photoelectron flux 
 *
 *
 */

class Photoionization
{
	protected:
		/// The parameters object
		XmlParameters* mpParameter;

		/// The photonGrid
		ublas::vector<double>* mpPhotonGreV;
		/// The photonGrid min
		ublas::vector<double>* mpPhotonGrMineV;
		/// The photonGrid max
		ublas::vector<double>* mpPhotonGrMaxeV;

		/// The electron bottom grid
		ublas::vector<double>* mpEBotEeV;
		/// The electron center grid
		ublas::vector<double>* mpECentEeV;
		/// The electron width grid
		ublas::vector<double>* mpEDdengeV;

		/// The altitude grid 
		ublas::vector<double>* mpAltGridKm;
		/// The solar flux model
		SolarFlux* mpSolflux;

		/// If the solflux is defined
		bool mbIsSolfluxDefined;


		/// The solar flux
		ublas::vector<double> mFluxPhcm_2s_1;



		/// The distance to the sun
		double mUA;
		/// The solar zenith angle
		double mSZADegree;

		/// The size of the planet in km
		double mRKm;


		/** To init the chapman function
		 *
		 * \param vpSp : the species vector to init atomic mass
		 *
		 */
		void ChapInit(std::deque<Specie*>& vpSp);

		/** The chapman function, used to init xchap
		 *
		 * \param vChiDeg : the chi angle in degree
		 * \param vAltKm : the vector giving the altitudes to consider
		 * \param vScaleHcm : the vector giving the scaleH to consider
		 *
		 *
		 */
		ublas::vector<double> Chap(double vChiDeg,const ublas::vector<double>& vAltKm,const ublas::vector<double>& vScaleHcm);

		/**
		 * Sperfc needed to compute the chapman function
		 * \param d: the parameter needed.
		 */
		double Sperfc(double d);


		/**
		 * Used to store the chapman value computed through chapinit
		 * this map is read by xchap(ialt,specie*)
		 *
		 */
		std::map< Specie*, ublas::vector<double> > mSpecieToChap;



		/**
		 * Used to store the density outside the specie
		 * is very important when the atmosphere is curved.
		 *
		 */

		std::map< Specie*,ublas::vector<double>* > mSpecieToDensity;


		/**
		 * Allows to retrieve the chapman function computed
		 * \param iAlt : the reference to the altitude
		 * \param vpSp : the reference to the species
		 *
		 *
		 */
		double Xchap(int iAlt,Specie* vpSp);



		/**
		 * Fills the photoionization - photoexcitation
		 * parameters in the species.
		 * \param vpSp  : the pointer to the species vector
		 * \param vEnergy : the energy of the flux
		 * \param iEv : the position on the energy loop
		 * \param iAlt : the position of the altitude in the loop
		 * \param vFl : the flux at this energy
		 * \param vExa : the absorption parameter
		 * \param rProelec : reference to the output : electron vs altitude
		 * \param rProelecE : reference to the output : electrons vs altitude vs energy
		 */
		void IonizeSpecies(std::deque<Specie*>& vpSp,double vEnergy, unsigned iEv,unsigned iAlt,double vFl,double vExa,ublas::vector<double>&rProelec,ublas::matrix<double> &rProelecE);


		/**
		 * Fills the photodissociation probability
		 * parameters in the species.
		 * \param vpSp  : the pointer to the species vector
		 * \param iEv : the position on the energy loop
		 * \param vEnergy : the energy of the flux
		 * \param vFl : the flux at this energy
		 */
		void PhotodissocieSpecies(std::deque<Specie*>& vSp,double vEnergy, unsigned iEv, double vFl);
		/**
		 * Search at which position put the electron of energy given
		 * in parameter in the electron grid.
		 * \param vEnergy : the energy of the electron
		 * \return the position
		 *
		 */
		unsigned Search(double vEnergy);



	public:

		/**
		 * Initializes the class 
		 * \param pParam : the object parameter
		 */
		Photoionization(XmlParameters* pParam);

		/** 
		 * The destructor
		 */
		virtual ~Photoionization();
		/**
		 * Initializes the main parameters
		 */
		void Init(ublas::vector<double>* pPhotonGrideV,ublas::vector<double>* pPhotonGridmineV,ublas::vector<double>* pPhotonGridmaxeV,ublas::vector<double>*pElecBotEeV,ublas::vector<double>* pElecCentEeV,ublas::vector<double>* pElecDdengeV,ublas::vector<double>* vAltGridKm, double vUA,double vRKm,double vSZADegree);

		/**
		 * Computes the solar flux impact  with the densities - column densities
		 * embedded in the species files 
		 *
		 * \param vpSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section
		 * \param rResult : the result vector of species -> were the computations are written
		 * \param rResultFlux : flux of secondary electron, needed to compute the secondary production.
		 *
		 *
		 */
		void ComputePhotoionization(std::deque<Specie*>& vpSp,
					     std::deque<Specie*>& rResult,
					     EFlux& rResultFlux);
		/**
		 * Computes the solar flux photodissociation probability of a particle, for the species in the vpSP
		 * writes a file
		 * \param vpSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section
		 * \param vSuffix : suffix for the written file
		 *
		 *
		 */
		void ComputePhotodissociationBranching(std::deque<Specie*>& vpSp, std::string vSuffix);


		/**
		 * Computes the transmission function in function of the wavelength box.
		 * \param vpSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section
		 * \param vOutFile : the output file, containing the boxes and the corresponding transmissions
		 */
		void ComputeTransmission(std::deque<Specie*>& vpSp,std::string vOutFile);


		/* =====================================
		 * No more useful with the new technique 
		 * =====================================
		 * Computes the solar flux impact with the given densities - column densities
		 *
		 * \param vpSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section
		 * \param rResult : the result vector of species -> were the computations are written
		 * \param rResultFlux : flux of secondary electron, needed to compute the secondary production.
		 * \param vAltGridKm : the altitude grid for the computation (vertical)
		 * \param vSpDensitycm_3 : the density for the species map specie->density grid
		 *
		 * \todo Implement the bent case!!!

		void ComputePhotoionization(std::deque<Specie*>& vSp,
					     ublas::vector<double> vAltGridKm,
					     std::map< std::string, ublas::vector<double> > vSpDensitycm_3,
					     std::deque<Specie*>& pResult,
					     EFlux& rResultFlux
					     );
		 */



};




#endif
