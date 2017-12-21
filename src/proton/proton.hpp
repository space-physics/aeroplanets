/**
 * \file proton.hpp
 * \brief Definition for the proton/hydrogen impact ionization and excitation,
 *	  and transport
 *	  Copyright G Gronoff Nov 2011
 *	  Last Modification : $Id$
 */



#ifndef PROTON_IONIZATION_HPP
#define PROTON_IONIZATION_HPP
#include <species/species.hpp>
#include <phflux/phflux.hpp>


// The fortran solver for matrix exponential
// EXPOKIT function modified

extern "C"
{
	/**
	 * DGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
	 * of the nonhomegenous linear ODE problem w' = Aw + u, w(0) = v.
	 * phi(z) = (exp(z)-1)/z and A is a General matrix.
	 * See the fortran source for more information, or expokit
	 *
	 * \param vOrder: the order of the matrices
	 * \param vKrylovSize: the maximum size for the Krylov basis (10 by default in the code)
	 * \param vTime : time at wich  the solution is needed (fractional step max) (1. in the code)
	 * \param vU : input line matrix of size vOrder. Operand vector with respect to the phi function (forcing term of the ODE problem)
	 * \param vV : input line matrix of size vOrder. Operand vector with respect to the exp function (initial condition of the ODE problem)
	 * \param vOutput : output line matrix of size vOrder. Computed approximation of exp(t*A)v + t*phi(tA)u 
	 * \param vTol : requested accuracy tolerance on vOutput.
	 * \param vAnorm : input approximation of some norm of A
	 * \param vWsp : workspace 
	 * \param vWsplength : the length of the first workspace
	 * \param vIwsp : workspace
	 * \param vIwsplength : the length of the first workspace
	 * \param vItrace : running mode (0: silent) 
	 * \param vIflag :  Exit flag (0 : no problem!)
	 * \param vValmatrix : the input matrix values
	  \param vIindices : the pointer to the first value of vValmatrix rows
	 * \param vJindices : the pointer to the column of vValmatrix
	 */
	extern void dgphiv_(int* vOrder, int* vKrylovSize, double* vTime, double* vU, double* vV, double* vOutput, \
	double* vTol, double* vAnorm, double* vWsp, int* vWsplength, int* vIwsp,  int* vIwsplength, int* vItrace, \
	int* vIflag, double* vValmatrix, unsigned* vIindices, unsigned* vJindices);


	extern void coocsr_(int* nrow, int* notnull, double* mat, int* line, int*col, double* outmat, int* outcol, int* outl);
};
/*
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
*/



/**
 *
 * \ingroup ionization_processes
 *  Computes the transport of Hydrogen and Proton
 *  Computes the ionization, excitations... by H,p impact
 *
 */


class ProtonHydrogenTransport
{
	private:
		/// The switch to know if the redistribution is active
		bool mbIsRedis;
		/// The number of angles
		int mNbAngle;

		/// The magnetic mirror matrix
		ublas::matrix<double> mMirror;

		/// The magnetic mirror type
		int mMirrorType;

		/// The beam spreading parameter
		double mBeamSpreading;

		/** Initialization of the Mirror
		 *
		 * \param vpdB_B : the dB_B vector
		 * \param type : the mirror type (0 : nothing, -1: symmetric cf fortran proton code)
		 *
		 */
		void InitMirror(ublas::vector<double>* vpdB_B, int type);


		/**
		 * Initialization of the matrices mLoss, mBdiag and mBinf
		 * \param vSp : the species vector
		 */
		void InitMatrices(std::deque<Specie*> vSp);



		/**
		 * Returns the phase function for the species vK, for the angle vMu and the energy vE
		 * \param vSp the vector of species
		 * \param vK : the position of the species in the vSp
		 * \param vJ : the final angle for the phase
		 * \param vJJ : the initial angle for the phase
		 * \param vE : the energy
		 * \param vSt : the proton phase type to be selected (0: excitation-ionization; 1: elastic; 2: exchange)
		 */
		double PhaseR(std::deque<Specie*> vSp, unsigned vK, unsigned vJ, unsigned vJJ, unsigned vE, unsigned vSt);
		/**
		 * Multidimensional matrices (see protmat.f)
		 * [species][energy][FinalParticle, InitialParticle][FinalAngle, InitialAngle]
		 * (Particle: Proton or Hydrogen)
		 */
		ublas::vector< ublas::vector< ublas::matrix< ublas::matrix<double> > > > mLoss, mBdiag, mBinf;


		/**
		 * The matrix to be solved
		 * [altitude][energy][finalparticle,initialparticle][finalangle, initialangle]
		 */
		ublas::vector< ublas::vector< ublas::matrix< ublas::matrix<double> > > > mMatInf, mMatDiag;

		/** Transforms the multidim matrices into sparce matrices for the computation
		 * \param vAltPos : the position in altitude
		 * \param vItype : The bloc type, in term of proton/hydogen going (down) leaving another proton/hydrogen going (down) for 1. I.E. 1 down down, 2 up up, 3 down up, 4 up down 
		 *
		 * nb: in the fortran, this is a mix between coocsr and maxrow
		 *
		 * The submatrices down down, and up up, are the homogeneous matrices.
		 * The submatrices size are (nb-energy * nb-angle/2)**2
		 *
		 * \param vrVal : vector the values of the sparse matrix
		 * \param vrRowP : vector[row] = the position of the first value for this row in the sparce matrix
		 * \param vrCol : the column for each values
		 * the indices for RowP and Col are in fortran style [ 1 -> size ], and the rowp vector is taller: the number of value + 1 is the last element
		 *
		 *
		 * The matrix to be used in the fortran is in the form particle x energy x angle / 2 (for up or down). (This is squared)
		 * Therefore, the index for a specific column is particletype x nbE x nbA / 2 + energy x nbA / 2 + angle ( + 1 for Fortran indexing)
		 * For the inferior matrix, which defines the degradation in energy (redistribution), the energy is put towards the inferior energy => e - 1.
		 *
		 *                           P                                                          H
		 *  --------------------------------------------------- ---------------------------------------------------------------
		 *
		 *  E1              E2       ...                    En   E1             E2            ...                         En
		 *  ============   ===========                          ============    ============                              =====
		 *  Mu1 Mu2  Mui   Mu1 Mu2  Mui                         Mu1 Mu2  Mui    Mu1 Mu2  Mui                    
		 *
		 * (The same structure applies for both row and columns)
		 * (The column indexing corresponds to the initial particles, energy, and angle; the row indexing to the final state(particle, energy, and angle)).
		 *
		 *
		 * \return the maximum value of the matrix
		 */
		double BlocMatrix(unsigned vAltPos,unsigned vItype, ublas::vector<double>& vrVal, ublas::vector<unsigned>& vrRowP, ublas::vector<unsigned>& vrCol);


		/**
		 * multiply vX by the sparce matrix vVal, vRowP, vCol.
		 * \param vX : the multiplying vector
		 * \param vVal : vector the values of the sparse matrix
		 * \param vRowP : vector[row] = the position of the first value for this row in the sparce matrix
		 * \param vCol : the column for each values
		 * \return the multiplication vector
		 */
		ublas::vector<double> VectMultSparceMat(const ublas::vector<double>& vX, const ublas::vector<double> &vVal, const ublas::vector<unsigned> & vrRowP, const ublas::vector<unsigned> & vrCol);


		/**
		 * FluxToVec : reads the fluxes mpProtFlux and mpHFlux and puts them in vectorial form to be used in the DGPHIV function (and VectMultSparceMat).
		 * \param vAlt : the altitude at which the flux is read. If 0, it reads the precipitated flux (forced by the user).
		 * \param vbIsUp : true if the needed flux is upwards.
		 */
		ublas::vector<double> FluxToVec(unsigned vAlt, bool vbIsUp);

		/**
		 * Inits the Proton and H downward fluxes at altitude 0
		 */
		void FluxToVecAlt0();
		/**
		 * VecToFlux updates the fluxes in the objects mpProtFlux and mpHFlux using the data in vInputFlux. It is the reverse of FluxToVec (anyway, it does not update the user defined precipitated fluxes).
		 *
		 * \param vInputFlux : the flux to be put in the object (the structure is the same as the one described in BlocMatrix, but in vector form)
		 * \param vAlt : the altitude at which the flux is evaluated. 
		 * \param vbIsUp : true if the needed flux is upwards.
		 */
		void VecToFlux(const ublas::vector<double>& vInputFlux, unsigned vAlt, bool vbIsUp);


		/// The initial tolerance for the loops
		double mInitialTolerance;
		/// The coeff for the modification of the  tolerance 
		double mCoeffTol;
		/// The maximum loop number
		unsigned mLoopMax;
		/// Variable to force the reflection of the protons at the bottom TODO: put an albedo
		bool mbIsReflected; // True (default) if the bottom is reflected


		// Width grid for the protons, based on the difference between the energies (closer to the fortran results)
		ublas::vector<double> mProtonWidthGrideV;

	protected:
		/// The parameters object
		XmlParameters* mpParameter;

		/// The altitude grid
		ublas::vector<double>* mpAltGridKm;

		/// The dB?B parameter
		ublas::vector<double>* mpdB_B;

		/// The inverse sinus of the dip angle 
		ublas::vector<double> mInvSin;

		/// Retrieve the inverse sinus
		void RetrieveInverseSinus();

		/// The proton flux object
		boost::shared_ptr<PHFlux> mpProtFlux;

		/// The hydrogen flux object
		boost::shared_ptr<PHFlux> mpHFlux;

		/// The gaussian angle object for the proton transport
		boost::shared_ptr< MathFunction::GaussianAngle > mpGa;


		/// The input energy, including the up and down fluxes
		double mEinput;

		/*-------------
		 * Initialization and vectors
		 * for the phase functions
		 *
		 * -----------*/

		/**
		 * Initialization of the phase function
		 * Reads the parameters -- called by the constructor (?)
		 */
		void InitPhaseFunction();


		/**
		 * Performs the transport, and computes the hydrogen and proton fluxes
		 * using the continuous slowing down (CSD) approximation.
		 */
		void Transport();

		/**
		 * Computes the energy in the fluxes, and the energy conservation
		 * \param vSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section \warning HERE it is a vector of species : too difficult, an not effective to make a pointer of pointer when there is so few things in the vector (the number of species is a small thing! I think this statement can be valid even if hundreds of species).
		 */

		void EnergyConservation(std::deque<Specie*> vSp);

		/**
		 * Computes the ionization and the secondary electron flux
		 * \param vSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section \warning HERE it is a vector of species : too difficult, an not effective to make a pointer of pointer when there is so few things in the vector (the number of species is a small thing! I think this statement can be valid even if hundreds of species).
		 * \param rResultFlux: The resulting electron flux.
		 */

		void IonizeSpecies(std::deque<Specie*> vSp,EFlux& rResultFlux);

	public:

		/**
		 * Initializes the class
		 * \param pParam : the object parameter
		 * \param vAltGridKm : the altitude grid (pointer)
		 * \param vpGa : shared pointer to the Proton Gaussian angle grid
		 * \param vpdB_B : the dB/B magnetic field parameter
		 */
		ProtonHydrogenTransport(XmlParameters* pParam, ublas::vector<double>* vAltGridKm, boost::shared_ptr< MathFunction::GaussianAngle > vpGa, ublas::vector<double>* vdB_B);
		/**
		 * Compute the transprot of protons and hydrogen
		 * and computes the resulting species.
		 * \param vSp : the vector of species, which embeed the density and column density. And which point to the photoionization cross section \warning HERE it is a vector of species : too difficult, an not effective to make a pointer of pointer when there is so few things in the vector (the number of species is a small thing! I think this statement can be valid even if hundreds of species).
		 * \param vPFlux: The input proton flux parameter.
		 * \param vHFlux: The input hydrogen flux parameter.
		 * \param rResult : the result vector of species -> were the computations are written
		 * \param rResultFlux: The resulting electron flux.
		 * Note: to avoid a lot of difficulties, the rResult is also a reference instead of a pointer
		 */

		void ComputeProtonImpact(std::deque<Specie*> vSp,
				boost::shared_ptr<PHFlux> vPFlux, boost::shared_ptr<PHFlux> vHFlux,
				std::deque<Specie*>& rResult, EFlux& rResultFlux);


		// The number of electrons and ions produced at each altitudes
		ublas::vector<double> mElecProduction, mIonProduction;

		/**
		 * Print the fluxes
		 * \param vFilename : the file where to write
		 * \param vAlts : the positions of the altitudes in the altitude grid
		 * \param vOption: option for the kind of flux to print: 0 - down, 1 - up, 2 total, 3 net
		 * \param vIsProton : true if we print for the protons, false for the H
		 */
		void PrintFluxes(std::string vFilename,ublas::vector<double> vAlts,unsigned vOption, bool vIsProton);
		
		void PrintEnergyFluxes(std::string vFilename,ublas::vector<double> vAlts, bool vIsProton, bool vIsLength=false);


		void PrintEnergyConservation(std::string vFilename);
	private:
		//
		void PrintFluxTot(std::string vFilename);
		void ForceFluxTot();
		double mError1, mError2, mReflectedEnergy, mAbsorbedEnergy, mInelasticEnergy;

};


#endif

