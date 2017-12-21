/**
 * \file photo/documentation.hpp
 * \brief The documentation for the photo part
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: documentation.hpp 643 2009-10-01 12:56:29Z gronoff $
 *
 */





/**
 * \page photoionizat_model Solar fluxes and photoionization model
 * 
 *  In this library, we set up a solar flux, and we lauch the photoionization model.
 * 
 * The Photoionization model, inside the init function, reads the /sun/model/flux_model:type to know which type of
 * flux model to use. There is a switch between  these different fluxes. 
 *
 * \section photoionizat_model_exist The different flux models
 *
 * The solar fluxes inherit from the SolarFlux class. 
 * To get the flux, in ph/cm2/s (so there is a given number of photons inside
 * each boxes, that means that we have to use the RedistributeChaoticFlux to modify our grid).
 *
 * - Solar39Boxes is the first EUV-XUV model created for Trans++
 *		it is based on the Torr and Torr fluxes and Boxes
 *		unfortunately, it is a problem to redistribute the 
 *		fluxes into another Grid (problem that instaured the
 *		need for the RedistributeChaoticFlux function).
 *		Unfortunately, the original model had not line width
 *		which creates a bug in the RedistributeChaoticFlux.
 *		So a width of 0.01 nm was put for the lines
 *		(no effect on the results of course).
 *
 *		To have a better accuracy when we use that model, it
 *		is better to add lines in your standard photon grid.
 *		Or to use the modified Torr and Torr grid.		
 *
 *
 *
 *
 * 
 * \section photoinization_flux_model_create  How to create a flux model
 *
 *
 * To create a flux model, you have to inherit the SolarFlux
 * \code
 * class Solar39Boxes : public SolarFlux
 * \endcode
 * You can then add your private functions.
 *
 * The most important is the constructor:
 *
 * \code
 * Solar39Boxes(XmlParameters* pParam);
 * \endcode
 * which needs the pParam parameter to read things in the configuration.
 *
 * The second point is the RetrieveFlux function:
 *
 * \code
 * void RetrieveFlux(double vUA,const std::vector<double>& vMainGrideV,const std::vector<double>& vMinGrideV,const std::vector<double>& vMaxGrideV,std::vector<double>& rResuFluxPhcm_2s_1);
 * \endcode
 *
 * It is declared virtual in the SolarFlux class.
 * This function has to read the necessary inputs for the flux model (ex f107 for 
 * the Solar39Boxes).
 * And it has to use the given grids to interpolate (or Redistribute) on the
 * rResuFluxPhcm_2s_1 parameter.
 * Be careful with the units, it is Ph.cm-2.s-1 -> there is a count of the number
 * of photons inside the box. So, this is not *really* a flux. That's why the
 * Solar39Boxes uses the RedistributeChaoticFlux function.
 *
 *
 *
 *
 *
 *
 *
 * \section photoionization_model The photoionization model
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */





