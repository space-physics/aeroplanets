#ifndef RUNLSQ_HPP
#define RUNLSQ_HPP

#include <leastsq/least.hpp>
#include <atmo/atmo.hpp>
/**
 * The least square function to call
 *
 */
class RunAtmoLsq: public LMlsq
{
	private:
		/// The error of the measured
		std::deque< ublas::vector<double> > mErrors;
		/// The name of the species measured (for adjustements)
		std::deque< std::string > mAdjustSpecie;
		/// The state of the species measured (for adjustements)
		std::deque< std::string > mAdjustState;
		/// The Id of the species measured to adjust
		std::deque< SpecieId > mAdjustSpecieId;
		/// The parameters of the species to adjust (type: production, VER, limb, density... + ...)
		std::deque< std::deque<int> > mAdjustParams;
		/// The atmosphere object to run
		Atmo* mpAtmo;
		/// The species name to adjust
		std::deque< std::string > mSpecieName;
		/// The species model to adjust
		std::deque<int> mSpecieModel;
		/// The additional; non-adjusted; parameters
		std::deque< std::deque<double> > mAdditionalParameters;
		/// If true, the calibration of the data can be adjusted
		bool mbUseCalib;

		/////// To perform the fitting of the electron precipitation

		/// To check if we use the electron precipitation fitting
		bool mbIsElecPrecipAdjusted;
		/// The electron precipitation model used 
		unsigned mElecPrecipModel;
		/// The additional; non-adjusted; parameters for the electron precipitation model
		std::deque<double> mElecPrecipAdditionalParameters;




		/// Number of calls of exec
		unsigned mNbCalls;
		/// The tolerance to errors
		std::deque<double> mTol;

		/// The Chi square reduced (divided by nb of data point - nb of degree of freedom)
		double mChiSquarev;

		/// True if we print the ionosphere at each steps
		bool mbPrintIonosphere;
		/// The Printionosphere fileprefix
		std::string mPrintIonoFile;
		/// True if we print the thermosphere at each steps
		bool mbPrintThermosphere;
		/// The PrintThermosphere fileprefix
		std::string mPrintThermoFile;
		/// True if we print the list of all the parameters at the end
		//	bool mbPrintParamsSteps;
		/// The filename for printing the parameters
		//	std::string mParamsStepsFile;

		std::deque< std::deque< std::deque<double> > > mStepsPars;


		/// True if we print the comparison between data and model at each steps
		bool mbPrintDiffplot;
		/// The PrintDiff fileprefix
		std::string mPrintDiffFile;

		/** The option for the pDiff calculation:
		 * pDiff is the parameter showing the progression of the fit
		 * different options have to be available
		 */
		unsigned mDiffOption;
	public:
		/**
		 * Initialize the  class
		 * \param vX : the position of the measurement points
		 * \param vY : the value of the measure
		 * \param vZ : the errors of the measure (ponderation of the tolerance, 1 if nothing)
		 * \param vAdjustSpecie : the specie measured
		 * \param vAdjustState : the state of the specie measured
		 * \param vAdjustParams : some more informations about the state of the measured species
		 * \param vpAtmo : the aeroplanets model
		 * \param vSpecieName : name of the species that are ajusted
		 * \param vSpecieModel : type of the model for the ajusted species
		 * \param vAdditionalParameters : the non adjusted parameters
		 * \param vbIsElecPrecipAdjusted : if the electron precipitation is adjusted
		 * \param vElecPrecipModel : the electron precipitation model
		 * \param vElecPrecipAdditionalParameters : the additional parameters for the electron precipitation model
		 * \param vTol : the tolerance to errors
		 */
		RunAtmoLsq(std::deque< ublas::vector<double> > vX,
				std::deque< ublas::vector<double> > vY,
				std::deque< ublas::vector<double> > vZ,
				std::deque< std::string > vAdjustSpecie,
				std::deque< std::string> vAdjustState,
				std::deque< std::deque<int> > vAdjustParams,
				Atmo* vpAtmo,
				std::deque< std::string > vSpecieName,
				std::deque<int> vSpecieModel,
				std::deque< std::deque<double> > vAdditionalParameters,
				bool vbIsElecPrecipAdjusted,
				unsigned vElecPrecipModel,
				std::deque<double> vElecPrecipAdditionalParameters,
				std::deque<double> vTol
				):LMlsq(vX,vY),mErrors(vZ),mAdjustSpecie(vAdjustSpecie),mAdjustState(vAdjustState),mAdjustParams(vAdjustParams),mpAtmo(vpAtmo),mSpecieName(vSpecieName),mSpecieModel(vSpecieModel),mAdditionalParameters(vAdditionalParameters), mbIsElecPrecipAdjusted(vbIsElecPrecipAdjusted), mElecPrecipModel(vElecPrecipModel), mElecPrecipAdditionalParameters(vElecPrecipAdditionalParameters), mTol(vTol)
	{
		mbUseCalib=false;
		mNbCalls=0;
	
		for(size_t j=0; j < vX.size(); ++j)
		{
			if(mErrors[j].size()==0)
			{
				mErrors[j].resize(vY[j].size());
				for(size_t i=0; i<mErrors[j].size();++i)
					mErrors[j][i]=1.;
			}
		}
		assert(mErrors.size()==vY.size());
		mbPrintIonosphere=false;
		mbPrintThermosphere=false;
		//mbPrintParamsSteps=false;
		mbPrintDiffplot=false;
		mDiffOption = 0;
		for(size_t i=0; i < vAdjustSpecie.size(); ++i)
		{
			SpecieId tmp(vAdjustSpecie[i],vAdjustState[i]);
			mAdjustSpecieId.push_back(tmp); 
		}
	}
		/**
		 * Sets the diff option
		 * \param vDiffOption: the parameter changing the convergence criteria
		 *
		 * 0 and default: (log(measure) - log(data)) / error
		 * 1 : CHi2 minimizing (bad)
		 * 2: percentage (measure - data) / (measure + data) / error
		 */
		void SetDiffOption(unsigned vDiffOption)
		{
			mDiffOption = vDiffOption;
		}

		/**
		 * Set the print of each step plot active
		 * \param vFilen : the prefix of the output file
		 */
		void SetDiffPlot(std::string vFilen)
		{
			mbPrintDiffplot=true;
			mPrintDiffFile=vFilen;
		}

		/**
		 * Set the print of each step thermosphere active
		 * \param vFilen : the prefix of the output file
		 */
		void SetThermo(std::string vFilen)
		{
			mbPrintThermosphere=true;
			mPrintThermoFile=vFilen;
		}
		/**
		 * Set the print of each step ionosphere active
		 * \param vFilen : the prefix of the output file
		 */
		void SetIono(std::string vFilen)
		{
			mbPrintIonosphere=true;
			mPrintIonoFile=vFilen;
		}

		/**
		 * Prints the convergence list
		 * \param vFilename : name of the file where to print the list
		 */
		void PrintConvergence(std::string vFilename);

		/**
		 * Allows the use of the modification of the calibration
		 */
		void SetUseCalib()
		{
			mbUseCalib=true;
		}
		/**
		 * Returns the number of calls to exec
		 */
		unsigned GetNbCalls()
		{
			return mNbCalls;
		}

		/**
		 * Returns the ChiSquarev value
		 */
		double GetChiSquarev()
		{
			return mChiSquarev;
		}

		/**
		 * Returns the observable value, computed with the model, corrected with the multiplicator
		 * \param vTmpMult : the last value of the parameters, it can be the multiplicator
		 * \return  the data
		 */
		std::deque< ublas::vector<double> > RetrieveData(double vTmpMult);


		/**
		 * Function, called by LM. That function uses the parameters in 'pParam' to compute its output. The difference between the computed values and the measurement should be stored in 'pDiff'.
		 * \param vDiffsize : size of the pDiff array, it is equal to mYmeasu.size()
		 * \param vParamsize : size of the parameter array, ie the number of unknowns
		 * \param pParam : the parameters, to compute the model output
		 * \param pDiff : output array of double, will contain the differences between the model and the measure: used in the least sqr function
		 * \param vFlag: flag from the least sqr function.
		 * \return 0 if everything is ok
		 */
		int Exec(int vDiffsize, int vParamsize, const double* pParam, double* pDiff, int vFlag);
};


#endif
