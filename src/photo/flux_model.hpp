/**
 * \file flux_model.hpp
 * \brief The  definition of the solar flux
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: flux_model.hpp 1550 2012-08-22 19:29:22Z gronoff $
 *
 */




#ifndef FLUX_MODEL_HPP
#define FLUX_MODEL_HPP

#include <species/species.hpp>
#include <eflux/eflux.hpp>
#include <phflux/phflux.hpp>



/**
 * \ingroup photoionization_process
 * \ingroup flux_parameters
 * Main solar flux model class.
 * is intended to be derived.
 *
 */

class SolarFlux
{
	protected:
		/// The parameters object
		XmlParameters* mpParameter;
	public:
		/**
		 * 
		 * Initializes the solar flux model
		 * \param pParam : the parameter to read more inputs
		 */
		SolarFlux(XmlParameters* pParam);


		/**
		 * The destructor
		 *
		 */
		virtual ~SolarFlux();

		/**
		 * Retrieve the flux in the resuGrid (reference) parameter 
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid in eV, decreasing
		 * \param vMinGrideV : the min energy gridin eV, decreasing
		 * \param vMaxGrideV : the max energy grid in eV, decreasing
		 * \param rResuFluxPhcm_2s_1 : the result grid in eV cm-2 s-1, decreasing
		 *
		 */
		virtual void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1)=0;

		/**
		 * Retrieve the different multiplicators for the flux boxes, and applies the multiplication to 
		 * the boxes.
		 * The system reads the /aero_main/sun/model/BoxFactor markup, which is organised as an array of position - factor; starting at 0. Raise an error if we are out of the box.
		 * \param vFluxBoxescm_2s_1 : the flux boxes to be modified.
		 */
		void BoxMultiplication(ublas::vector<double>&  vFluxBoxescm_2s_1);

		/**
		 * Retrieve the different multiplicators for the flux in function of the range (nm or eV)
		 * The system reads the /aero_main/sun/model/RangeFactor markup
		 * \param vMainGrideV : the main energy grid in eV, decreasing
		 * \param vFluxBoxescm_2s_1 : the flux boxes to be modified.
		 */
		void RangeMultiplication(const ublas::vector<double>& vMainGrideV, ublas::vector<double>&  vFluxBoxescm_2s_1);


};




/**
 * The torr and torr fluxe
 * \ingroup flux_parameters
 * \ingroup photoionization_process
 */
class Solar39Boxes : public SolarFlux
{
	protected:
		/**
		 * The function to transform the f107 parameter, with the distance to the sun, in a standard flux
		 *
		 * \warning if you want to extend the grid, please extends the parts where we put rResuPhcm_2s_1 at 0!!!
		 *
		 * \param vF107 : the f107 flux
		 * \param vUA : the distance to the sun
		 * \param vFudge : the fudge factor when <25.7nm ie energie > 49 eV
		 * \param vMaxGrideV : the max of the grid on what we want the flux
		 * \param vMinGrideV : the min of the  grid on what we want the flux
		 * \param rResuPhcm_2s_1 : the reference to the result
		 */
		void F107ToFlux(double vF107, double vUA,double vFudge,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1);
	public:
		
		/**
		 * The constructor
		 * \param pParam : the needed XmlParameter
		 */
		Solar39Boxes(XmlParameters* pParam);

		/**
		 * The implementation of the function
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid
		 * \param vMinGrideV : the min energy grid
		 * \param vMaxGrideV : the max energy grid
		 * \param rResuFluxPhcm_2s_1 : the result grid 
		 *
		 */

		void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1);

};


/**
 * The User defined flux fluxe
 * \ingroup flux_parameters
 * \ingroup photoionization_process
 */

class SolarUserDefined : public SolarFlux
{
	private:
		/**
		 * Allows to transform the measurements into a flux for the considered planet (hence vUA)
		 * \param vUA : the distance of the planet to the sun
		 * \param vMaxGrideV : the max of the grid on what we want the flux
		 * \param vMinGrideV : the min of the  grid on what we want the flux
		 * \param rResuPhcm_2s_1 : the reference to the result
		 */
		void MeasurementToFlux(double vUA,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1);

	public:
		SolarUserDefined(XmlParameters* pParam);
/**
		 * The implementation of the function
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid
		 * \param vMinGrideV : the min energy grid
		 * \param vMaxGrideV : the max energy grid
		 * \param rResuFluxPhcm_2s_1 : the result grid 
		 *
		 */

		void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1);
};

/**
 * The torr and torr fluxes, with the OW additional flux to 175nm
 * This flux is basically what we have inside Trans*: it does not takes into account the O I lines. (Very bad inside the 120-140nm box; but takes into account Ly alpha). The main problem of this part is that Ly alpha needs a radiative transfer.
 * \ingroup flux_parameters
 * \ingroup photoionization_process
 */
class SolarTransBoxes : public SolarFlux
{
	protected:
		/**
		 * The function to transform the f107 parameter, with the distance to the sun, in a standard flux
		 *
		 * \warning if you want to extend the grid, please extends the parts where we put rResuPhcm_2s_1 at 0!!!
		 *
		 * \param vF107 : the f107 flux
		 * \param vUA : the distance to the sun
		 * \param vFudge : the fudge factor when <25.7nm ie energie > 49 eV
		 * \param vMaxGrideV : the max of the grid on what we want the flux
		 * \param vMinGrideV : the min of the  grid on what we want the flux
		 * \param rResuPhcm_2s_1 : the reference to the result
		 */
		void F107ToFlux(double vF107, double vUA,double vFudge,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1);
	public:
		
		/**
		 * The constructor
		 * \param pParam : the needed XmlParameter
		 */
		SolarTransBoxes(XmlParameters* pParam);

		/**
		 * The implementation of the function
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid
		 * \param vMinGrideV : the min energy grid
		 * \param vMaxGrideV : the max energy grid
		 * \param rResuFluxPhcm_2s_1 : the result grid 
		 *
		 */

		void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1);

};




/**
 * \ingroup photoionization_process
 * \ingroup flux_parameters
 * The HEUVAC model, to be used with the SolarFlux classes
 */

namespace HEUVAC
{

	/**
	 * Returns the wavelength and the flux in parameters (corrected with ModF74113)
	 * \param vWavelengthnm the wavelength bin in nm
	 * \param vFluxphcm_1s_1_bin : the flux in photon cm-1 s-1 per bin
	 */
	void GetF74113(ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin);
	/**
	 * Correction of the flux used by GetG74113 
	 * \param vWavelengthnm the wavelength bin in nm
	 * \param vFluxphcm_1s_1_bin : the flux in photon cm-1 s-1 per bin
	 */
	void ModF74113(const ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin);
	
	/**
	 * Find a value in an array with increasing values
	 * \param vArr : the array
	 * \param vAl : the value
	 */
	unsigned BiSplt(const std::vector<double>& vArr, double vAl);
	
	/**
	 * Returns the wavelength and the flux of HEUVAC, for the given f107 and f107 averaged
	 * \param vF107 the f 10.7 flux
	 * \param vF107a the f 10.7 flux averaged over 3 month
	 * \param vWavelengthnm the wavelength bin in nm
	 * \param vFluxphcm_1s_1_bin : the flux in photon cm-1 s-1 per bin
	 */
	void HeuvacFlux(double vF107, double vF107a, ublas::vector<double>& vWavelengthnm, ublas::vector<double>& vFluxphcm_1s_1_bin);

}



/**
 * The HEUVAC model
 * Problem: it does not go to the ly alpha
 * \ingroup flux_parameters
 * \ingroup photoionization_process
 */
class HeuvacFlux : public SolarFlux
{
	protected:
		/**
		 * The function to transform the f107 parameter, with the distance to the sun, in a standard flux
		 *
		 * \warning if you want to extend the grid, please extends the parts where we put rResuPhcm_2s_1 at 0!!!
		 *
		 * \param vF107 : the f107 flux
		 * \param vF107av : the f107 averaged over 3 month (see HEUVAC manual)
		 * \param vUA : the distance to the sun
		 * \param vMaxGrideV : the max of the grid on what we want the flux
		 * \param vMinGrideV : the min of the  grid on what we want the flux
		 * \param rResuPhcm_2s_1 : the reference to the result
		 */
		void HeuvacF107ToFlux(double vF107, double vF107av, double vUA, const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1);
	public:
		
		/**
		 * The constructor
		 * \param pParam : the needed XmlParameter
		 */
		HeuvacFlux(XmlParameters* pParam);

		/**
		 * The implementation of the function
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid
		 * \param vMinGrideV : the min energy grid
		 * \param vMaxGrideV : the max energy grid
		 * \param rResuFluxPhcm_2s_1 : the result grid 
		 *
		 */

		void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1);

};

/**
 * The HEUVAC model extended to low energies
 * Problem: it does not go to the ly alpha
 * \ingroup flux_parameters
 * \ingroup photoionization_process
 */
class HeuvacFluxLow : public SolarFlux
{
	protected:
		/**
		 * The function to transform the f107 parameter, with the distance to the sun, in a standard flux
		 *
		 * \warning if you want to extend the grid, please extends the parts where we put rResuPhcm_2s_1 at 0!!!
		 *
		 * \param vF107 : the f107 flux
		 * \param vF107av : the f107 averaged over 3 month (see HEUVAC manual)
		 * \param vUA : the distance to the sun
		 * \param vMaxGrideV : the max of the grid on what we want the flux
		 * \param vMinGrideV : the min of the  grid on what we want the flux
		 * \param rResuPhcm_2s_1 : the reference to the result
		 */
		void HeuvacLowF107ToFlux(double vF107, double vF107av, double vUA, const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuPhcm_2s_1);
	public:
		
		/**
		 * The constructor
		 * \param pParam : the needed XmlParameter
		 */
		HeuvacFluxLow(XmlParameters* pParam);

		/**
		 * The implementation of the function
		 * \param vUA : the distance of the planet to the sun
		 * \param vMainGrideV : the main energy grid
		 * \param vMinGrideV : the min energy grid
		 * \param vMaxGrideV : the max energy grid
		 * \param rResuFluxPhcm_2s_1 : the result grid 
		 *
		 */

		void RetrieveFlux(double vUA,const ublas::vector<double>& vMainGrideV,const ublas::vector<double>& vMinGrideV,const ublas::vector<double>& vMaxGrideV,ublas::vector<double>& rResuFluxPhcm_2s_1);

};



#endif
