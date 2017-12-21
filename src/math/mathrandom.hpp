/**
 * \ingroup Math 
 * \file mathrandom.hpp
 * \brief defines  random parameters -> allows to have a truly random system
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: mathrandom.hpp 1141 2010-09-08 22:06:18Z gronoff $
 *
 */
#ifndef MATH_RANDOM_HPP
#define MATH_RANDOM_HPP

#include <iostream>
#include <fstream>
#include <boost/limits.hpp>
#include <boost/random.hpp>
#include <ctime>
#include "logging.hpp"

/**
 * Class to perform random operations
 *
 */
class MathRandom
{
	private:
		/// The static random generator
		static boost::mt19937 msRng;
		/// The seed, to have a real generator
		static unsigned long msSeed;

		static bool mIsInitialized;
	

	public:
		/// Initialize the class : the mrRng and msSeed
		static void Initialize();

		/**
		 * Get a random number which follows
		 * the normal law, with sigma and mean
		 * \param vMean : the mean of the law
		 * \param vSigma : the sigma of the law (standard when error are given)
		 */
		static double GetNormal(double vMean=0.,double vSigma=1.);

		/**
		 * Get a random number which follows an uniform law
		 * \param vMin : the minimum possible
		 * \param vMax : the maximum possible
		 */
		static double GetUniformReal(double vMin=0.,double vMax=1.);





};




#endif
