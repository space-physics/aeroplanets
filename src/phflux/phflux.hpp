/** 
 * \file phflux.hpp
 * \brief Defines the PHFlux class 
 * Copyright G Gronoff Nov 2011
 * Last Modification $Id$
 */
#ifndef PHFLUX_HPP
#define PHFLUX_HPP
#include <eflux/eflux.hpp>
/**
 * \ingroup flux_parameters
 *
 * class to compute the proton and hydrogen flux parameters
 * - it can compute the precipitating proton and hydrogen flux
 *   and can be used as an input for the proton hydrogen transport code
 *
 */


class PHFlux
{

	protected:



		/// Direct link to the Gaussian angle. It is not defined in that function
		boost::shared_ptr<MathFunction::GaussianAngle> mpGAngle; 

		/// The parameter, to read the precipitation flux if needed
		XmlParameters* mpParameter;

		/// The name of the precipitating particle to retrieve (proton, hydrogen...)
		std::string mParticleName;

		/// The name of the precipitating particle to retrieve (proton, hydrogen...)
		std::string mParticlePosition;

		/// Read the measured precipitation 
		void ReadMeasuredPrecipitation();
		
		/**
		 *
		 * Redistributes the measured and read particle precipitation
		 * on the standard energy grid.
		 * The angle grid is from vertical precipitation (toward planet) to horizontal (nang/2) to vertical toward space).
		 *
		 * \param vGrideV : the energy grid of the read precip
		 * \param vAnglesDegree : the read angles in degree
		 * \param vInputFluxcm_2sr_1s_1 : the read flux
		 * \param rFluxcm_2sr_1s_1 : the result flux
		 * \param vbDown : true if the flux is downward (default) false if upwards
		 * 
		 */

		void RedistributeFlux(ublas::vector<double> vGrideV,
				std::deque<double> vAnglesDegree,
				const std::deque< ublas::vector<double> > & vInputFluxcm_2sr_1s_1,
				ublas::matrix<double>& rFluxcm_2sr_1s_1,
				bool vbDown=true
				);


		/**
		 *
		 * Warning for non existing cross section
		 *
		 */
		bool mWarningMissingCrs;
		/**
		 * Messages with respect to the warning
		 */
		std::deque<std::string> mMissingCrsMsg;

	public:
                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo-------------------+
                 |                                  |
                 |Parameters passed as pointer      |
                 |(main parameter for the atmo class|
                 | in fact)                         |
                 |                                  |
                 +---------------------------------*/
		
		
		/// Particle center grid
		ublas::vector<double>* mpGrideV;
		/// Particle width grid
		ublas::vector<double>* mpWidthGrideV;



		/**
		 *
		 * Initialization
		 * \param pParam : the standard xmlparameter
		 * \param pAngle : the gaussian angle class
		 * \param vpGrideV : the particle energy grid
		 * \param vpWidthGrideV : the particle energy width grid
		 * \param vParticleName : the name of the precipitating particle (now, proton and hydrogen
		 */
		PHFlux(XmlParameters* pParam, boost::shared_ptr<MathFunction::GaussianAngle> pAngle,ublas::vector<double>* vpGrideV, ublas::vector<double>* vpWidthGrideV, std::string vParticleName);



		/**
		 * Initialization for a bend line, by using the photofluxes.
		 * \param vPHFlux : list of the different fluxes (computed by the transport function)
		 * \param vAltitudesKm : the altitude corresponding to the vSpecies grids
		 * \param vSZADegree : SZA for each vSpecies
		 * \param vPath : The path for interpolating the new species
		 */
		PHFlux(std::deque< PHFlux* > vPHFlux,const ublas::vector<double>& vAltitudesKm, const std::deque<double> & vSZADegree,const Path & vPath);

		/**
		 * Initialization of the array at 0
		 * \param vNbAlt : the number of altitudes
		 */
		void InitVoid(unsigned vNbAlt);

		/**
		 * Computes the flux when considered the productions by energy-altitude, and an isotropic distribution of the primary .srticle
		 * \param vProdcm_3s_1eV_1 stands for the input particle production by altitude and energy (cm-3.s-1.eV-1)
		 * \param vProdcm_3s_1 stands for input particle production by altitude (cm-3.s-1)
		 * \param vSp : vector of species
		 * \warning: the densities are obtained through the species it is a problem for non-vertical studies!
		 */

		void InitIsotropic(ublas::matrix<double> vProdcm_3s_1eV_1, ublas::vector<double> vProdcm_3s_1,std::deque<Specie*> vSp);


		/**
		 * Allows to initialize the flux with the anisotropic production (for example from cosmic ray input)
		 * be careful: no interpolation is made on energy angle.
		 * \param vFluxcm_3s_1eV_1sr_1 : the number of particle per phase space volume (energy, volume, angle)
		 * \param vSp : deque of the species of your atmosphere
		 */
		void AddAnisotropicFlux(ublas::vector< ublas::matrix<double> > vFluxcm_3s_1eV_1sr_1,std::deque<Specie*> vSp);


                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo------------------------------+
                 |                                             |
                 |                                             |
                 |Useful functions : total energy, print fluxes|
                 |                                             |
                 |                                             |
                 +--------------------------------------------*/


		/**
		 * Return the total energy in the precipitating flux 
		 * and in the altitude dependant flux.
		 *
		 * Called for each type of precipitation
		 * this function allows to compute the
		 * energy dependant on the type.
		 * When we call it on the sum of these 
		 * precipitations, it is the total energy
		 * =>  easier than fortran for the programmer!
		 *     Hum, lets say,once again!
		 *
		 * Needs mProdcm_3s_1eV_1 and/or PrecipFluxDowncm_2s_1eV_1
		 * to be defined
		 *
		 *
		 * \param vAltGridKm -> altitude grid for the integration
		 * (it can be a length grid. In that case, we integrate
		 * against that one). In km
		 * \return the energy, in eV, in the column of 1 cm2
		 */
		double FluxEnergy(const ublas::vector<double> & vAltGridKm );


		/**
		 * Computes the energy in the precipitating flux (downward)
		 */
		double FluxEnergyDown();
		/**
		 * Computes the energy in the precipitating flux (downward)
		 * but computes it from the first layer
		 */
		double FluxEnergyDown2();

		/**
		 * Computes the energy in the escaping flux (upward)
		 */
		double FluxEnergyUp();
		
		/**
		 * Computes the energy in the precipitating flux (upward)
		 * but computes it from the first layer
		 */
		double FluxEnergyUp2();



		/**
		 * Print the precipitating energy mPrecipFluxDowncm_2s_1eV_1
		 * \param vFilename : the name of the file to write
		 */
		void PrintPrecip(std::string vFilename);

		/**
		 * Print  mFluxcm_2s_1eV_1
		 * \param vFilename : the name of the file to write
		 */
		void PrintFluxInt(std::string vFilename, const ublas::vector<double>& vAltGridKm);


		/**
		 * Print the total energy flux (energy/cm2/s/km)
		 * \param vFilename : the name of the file to write
		 * \param vAltGridKm : the altitude grid
		 * \param vIsLength : true if we work in term of length instead of altitude
		 */
		void PrintEnergyFlux(std::string vFilename,const ublas::vector<double>& vAltGridKm,bool vIsLength=false);
		/**
		 * Print the particle profile at each altitude, in function of energy
		 */
		void PrintProfile(std::string vFilename,const ublas::vector<double>& vAltGridKm);


                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo----------------+
                 |                               |
                 |Parameters embedded when needed|
                 |(mainly for solar ionization)  |
                 |                               |
                 +------------------------------*/


		/**
		 * Particle production by energy
		 * cm-3.s-1.eV
		 * (be careful with eV: it is not a distribution in boxes!)
		 * (it is divided by the size of the energy box in photoionization)
		 * Prod[z][E]
		 *
		 */
		//std::vector< std::vector<double> > mProdcm_3s_1eV_1;
		ublas::matrix<double> mProdcm_3s_1eV_1;

		/**
		 * Particle production by altitude
		 * Prod[z] cm-3
		 */
		ublas::vector<double> mProdcm_3s_1;


                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo---------------------+
                 |                                    |
                 |Parameters computed thanks to PHFlux|
                 |                                    |
                 +-----------------------------------*/


		/**
		 * The particle flux, in function of altitude and energy
		 * flux[z][E]
		 * cm-2.s-1 
		 */
	//	std::vector< std::vector<double> > mFluxcm_2s_1eV_1;
		ublas::matrix<double> mFluxcm_2s_1eV_1;
		/**
		 * The particle flux, in function of altitude, energy, and angle
		 * q[z][E][A]
		 * \warning Unlike in Trans*, I have normalized the fluxes function
		 * therefore, it is ALTITUDE, ENERGY, ANGLE -not eza-
		 * cm-2.s-1.sr-1
		 */
		//std::vector< std::vector< std::vector<double> > > mFluxAcm_2s_1eV_1sr_1;
		ublas::vector< ublas::matrix<double> > mFluxAcm_2s_1eV_1sr_1;


		/**
		 * Allows to compute mFluxcm_2s_1eV_1 from  mFluxAcm_2s_1eV_1sr_1
		 */
		void AnisotropicFluxToAverage();



                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo------------------------+
                 |                                       |
                 |Computation of precipitation parameters|
                 |                                       |
                 +--------------------------------------*/


		/**
		 * Reads the mainParameter to find if particle
		 * precipitation is set. If so, the objects 
		 * read these precipitation and stores it.
		 * Because it is a Flux function, it can
		 * be directly put in degrad.
		 * \param vAltitudes: if we have an altitude for the precipitation, we read it in that vector
		 * \param vSp : vector of species
		 */

		void ReadPrecipitation(const ublas::vector<double> & vAltitudes,std::deque<Specie*> vSp);



		/// True if we have a precipitation flux at the top of  the atmosphere
		bool mPrecipitationDefined;

		/// The downward flux at the top of the atmosphere
//		std::vector< std::vector<double> > mPrecipFluxDowncm_2s_1sr_1;
		ublas::matrix<double> mPrecipFluxDowncm_2s_1sr_1;
		/// The upward flux, for comparison at the top of the atmosphere
//		std::vector< std::vector<double> > mPrecipFluxUpcm_2s_1sr_1;
		ublas::matrix<double>  mPrecipFluxUpcm_2s_1sr_1;

                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo-------------------------------------+
                 |                                                    |
                 |To allow a better communication with other functions|
                 |                                                    |
                 +---------------------------------------------------*/

		/// Return the number of angles
		unsigned NbAngles()
		{
			return mpGAngle->mNbAngles;
		}

		/**
		 * The addition operator:
		 * thanks to this one, you can add the different particle
		 * fluxes classes, and the sum can be passed in parameter
		 * inside the degrad function.
		 */
		friend  PHFlux  operator+(const PHFlux& vA, const PHFlux& vB);



		PHFlux& operator+=(const PHFlux& vB);

//		void printelecadress()
//		{
//			std::cout<<"electron adress : "<<mpElecCentEeV<<std::endl;;
//		}
		/**
		 * Return the angle weight
		 * \param position : the position of the angle in the std vector
		 * \return the angle weight
		 *
		 */
		double ReturnAngleWeight(unsigned position)
		{
			return mpGAngle->mWeight[position];

		}
		
		/**
		 * Returns the flux integrated for the altitude vAlt, the energy vEne
		 * and the option vOption
		 * \param vAlt: the considered altitude position
		 * \param vEne: the considered energy position
		 * \param vOption: the option from 0-down, 1 - up, 2 total, 3 net
		 */
		double FluxInt(unsigned vAlt, unsigned vEne, unsigned vOption);
		/**
		 * Return a reference to the angle object
		 * \return the reference to the angle object
		 */
		boost::shared_ptr<MathFunction::GaussianAngle> ReturnAngle()
		{
			return mpGAngle; 
		}
};




#endif
