/**
 * \file least.hpp
 * \brief Definition for the least square system
 * Copyright G Gronoff July 2010
 * Last Modification : $Id$
 */

#ifndef LEAST_SQUARE_HPP
#define LEAST_SQUARE_HPP

#include <math/mathfunction.hpp>



/**
 * The Levenberg Marquardt least square prototype
 * for an easier use of LM, just derive from that class, and everything
 * will be ok.
 *
 */

class LMlsq
{
	protected:
		/// The position of the point for the measure: will be used as parameter for the comparison with model
		std::deque< ublas::vector<double> > mXinit;
		/// The measurements at the points defined in mXinit: will be used as comparison with the output of the model
		std::deque< ublas::vector<double> > mYmeasu;
		/// The total size of the vectors
		unsigned mSize;
	public:
		/**
		 * Initialisation of the class: fill mXinit and mYmeasu
		 * \param vXinit : the position of the measurement points
		 * \param vYmeasu : the value of the measure
		 *
		 */
		LMlsq(std::deque< ublas::vector<double> > vXinit, std::deque< ublas::vector<double> > vYmeasu):mXinit(vXinit),mYmeasu(vYmeasu)
		{
			mSize = 0;
			assert(mXinit.size()==mYmeasu.size());
			for(size_t i = 0; i < mXinit.size(); ++i)
			{
				assert(mXinit[i].size()==mYmeasu[i].size());
				mSize += mXinit[i].size();
			}
		}

		/**
		 * Virtual function, called by LM. That function uses the parameters in 'pParam' to compute its output. The difference between the computed values and the measurement should be stored in 'pDiff'.
		 * \param vDiffsize : size of the pDiff array, it is equal to mYmeasu.size()
		 * \param vParamsize : size of the parameter array, ie the number of unknowns
		 * \param pParam : the parameters, to compute the model output
		 * \param pDiff : output array of double, will contain the differences between the model and the measure: used in the least sqr function
		 * \param vFlag: flag from the least sqr function.
		 * \return 0 if everything is ok
		 */
		virtual int Exec(int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag)=0;

		/**
		 * Returns the size of measured parameters
		 */
		unsigned Psize()
		{
			return mSize;
		}
		virtual ~LMlsq(){}
};


/**
 * Namespace to store the leastsquare functions
 */
namespace LstSq
{
	/**
	 * Performs a least square interpolation using the Levenberg Marquardt method, while computing the Jacobian.
	 * \param pLm : pointer to an object of type LMlsq, which contains the function and the data to be fitted
	 * \param vTol : the tolerance, defined by the experimental data
	 * \param rParams: deque containing the first guess of the parameters, when the computation is finished, it contains the result of the LM.
	 * \param rCovar : matrix containing the covariant at the end of the computation. This covariant matrix is used to compute the uncertainty of the inversion.
	 * \return the information from the LM algorithm
	 */
	int LMleast(LMlsq* pLm, double vTol, std::deque<double>& rParams, ublas::matrix<double> & rCovar);


	/**
	 * Subroutine for the LMleast function. Should not be used elsewhere.
	 * It uses the parameters in LMleast to call the object (inside void*).
	 *
	 * \param pP : pointer containing the pointer to the object
	 * \param vDiffsize : size of the pDiff array, it is equal to mYmeasu.size()
	 * \param vParamsize : size of the parameter array, ie the number of unknowns
	 * \param pParam : the parameters, to compute the model output
	 * \param pDiff : output array of double, will contain the differences between the model and the measure: used in the least sqr function
	 * \param vFlag: flag from the least sqr function.
	 * \return 0 if everything is ok
	 */
	int Least(void* pP, int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag);
};



#endif
