/**
 * \file lines.hpp
 * \brief  Definition of the class to add lines to solar flux grids
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: lines.hpp 928 2010-03-15 23:37:58Z gronoff $
 *
 */


#ifndef LINES_HPP
#define LINES_HPP
#include "flux_model.hpp"

/**
 * \ingroup photoionization_process
 *
 * Class SolarGridLines : add the lines vectors.
 * It allows to improve the quality of the solar
 * flux grids (the torr and torr grid includes
 * lines by default, but this grid is restrained
 * to EUV-XUV analysis, and we want to extend this!) 
 *
 */

class SolarGridLines
{
	protected:
		/// To read the parameter for the lines
	       	XmlParameters* mpParameters;
		/// The vector for the minimum lines
		ublas::vector<double> mMinLines;
		/// The vector for the maximum lines
		ublas::vector<double> mMaxLines;
		/// Initialization of the standard grid lines (for the std model, Torr et Torr)
		void InitStandard();

		/// Initialization of the standard extended grid lines (for the trans* model, Torr et Torr+ extention to the SR continuum)
		void InitStandardExtended();
		/// Initialization for the user defined lines
		void InitUserDefined();


	public:
		SolarGridLines(XmlParameters* vpParam);

		/**
		 * Add the solar lines to the extremum lines vector
		 * \param rPhGrMin : reference to the min vector (it has more energy...)
		 * \param rPhGrMax : reference to the min vector (it has less energy...)
		 */
		void AddSolarLines(ublas::vector<double>& rPhGrMin,ublas::vector<double>& rPhGrMax);
		




};


#endif
