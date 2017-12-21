/**
 * \file main/timeline.hpp
 * \brief The documentation of the timeline : allows to understand the different processes occuring 
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: timeline.hpp 643 2009-10-01 12:56:29Z gronoff $
 *
 */




#ifndef TIMELINE_HPP
#define TIMELINE

/**
 * 
 * \page architecture_timeline Timeline
 *
 * These timelines allows you to understand the fundamental processes wich happen during the run of Trans++.
 * It is also a good source of information to know where to start browsing the datas.
 * So, it is very important to hack Trans++!
 *
 *
 * \section photoionization PhotoIonization Timeline, with vertical atmosphere
 *
 * This version is valid for the 0.10 version of Trans++.
 * Future version will also work with the electron precipitation, and will compute electron impact ionization
 * wich is not the case here.
 *
 *
 * \subsection photoionization_mainline Main execution line
 *
 * 		- The main program runs, it uses the file given in parameter to
 * 		  initalize the atmosphere Atmo object.
 * 		- The atmosphere creates an XmlParameter object to read this file.
 * 		- The atmosphere :
 * 				* Inits the different Grids
 * 				* Inits the planet object
 * 				* Inits the neutral atmosphere
 * 				* Inits the neutral temperature
 * 				* Inits the Ionosphere model
 * 				* And launches the computation of the 
 * 				  vertical column density.
 * 		  These actions are explained below.
 * 		- In main, we ask atmosphere to Compute:
 * 				* Compute check if the photoionization is required
 * 				* If so, it launches the PhotoIonize function
 * 				* This function initializes a Photoionization object
 * 				* This function initalizes a EFlux object
 * 				* This function launches the ComputePhotoionization member function. This function gives mPhotoionizationResu, a vector of Specie* objects where the main productions are computed. This object can be added with another one thanks to the SpecieUtils::MergeResu function. It gives the mpPhotoFlux, which is the flux object for the photoelectrons. This object is needed to the electron impact ionization computations.
 * 		- In main, we ask atmosphere to ProceedOutputs: it read the parameters, and create the required output files.
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




#endif
