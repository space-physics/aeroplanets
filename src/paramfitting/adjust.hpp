#ifndef ADJUST_HPP
#define ADJUST_HPP

#include <leastsq/least.hpp>
#include <atmo/atmo.hpp>
#include "runlsq.hpp"

class Adjustements
{
	private:
		/// The aeroplanet model
		Atmo* mpAtmo;
		/// The parameters for the adjustements xml file
		XmlParameters* mpParams;
		/// The report at the end of the fits
		std::string mReport;
		/// The report at the end of the fits, in xml form for studying the effects of MC inputs!
		std::string mXmlReport;
		/// The final data from the model, corrected
		std::deque< ublas::vector<double> > mFinalData;
		/// one adjustement
		void Adjust(TiXmlNode* vpNode);
		/**
		 * Fills the report for a node
		 * \param vRet : the output status of the fitting
		 * \param vNames : the name of the adjusted species
		 * \param vOrigParams : the first guess value parameters for these speecies
		 * \param vParams : the value of the parameters for these speecies
		 * \param vModels : the models for theses parameters
		 * \param vbIsElecPrecipAdjusted : if the electron precipitation is adjusted
		 * \param vElecPrecipModel : the electron precipitation model
		 * \param vElecPrecipAdditionalParameters : the additional parameters for the electron precipitation model
		 * \param vNbCalls : the number of calls to the function
		 * \param vChiSquarev : the value of Chi square divided by the number of free parameters (should be  close to 1)
		 * \param vComatrix  : covariance matrix, its elements are the variance and covariance of the fitted params
		 */
		void FillReport(int vRet, std::deque<std::string> vNames, std::deque<double> vOrigParams,  std::deque<double> vParams,  std::deque< std::deque<double> > vAdditionalParameters, std::deque<int> vModels, bool vbIsElecPrecipAdjusted, unsigned vElecPrecipModel, std::deque<double> vElecPrecipAdditionalParameters, unsigned vNbCalls, double vChiSquarev, ublas::matrix<double> vComatrix);
		/**
		 * Apply the final modifs to the atmosphere: pre-next computation
		 * \param vNames : the name of the adjusted species
		 * \param vParams : the value of the parameters for these speecies
		 * \param vModels : the models for theses parameters
		 * \param vAdditionalParameters : the non adjusted parameters
		 * \param vbIsElecPrecipAdjusted : if the electron precipitation is adjusted
		 * \param vElecPrecipModel : the electron precipitation model
		 * \param vElecPrecipAdditionalParameters : the additional parameters for the electron precipitation model
		 */
		void ApplyModifs(std::deque<std::string> vNames, std::deque<double> vParams, std::deque<int> vModels,std::deque< std::deque<double> > vAdditionalParameters, bool vbIsElecPrecipAdjusted, unsigned vElecPrecipModel, std::deque<double> vElecPrecipAdditionalParameters);
		/// If true, the calibration of the data can be adjusted
		bool mbUseCalib;

		/*
		/// The altitudes of the measurements
		ublas::vector<double> mYalts;
		/// The measurements
		ublas::vector<double> mYmeasu;
		/// The errors in the measurements
		ublas::vector<double> mYerrors;
		*/
		// deque with the measurements altitudes vectors
		std::deque< ublas::vector<double> > mMyalts;
		// deque with the measurements value vectors
		std::deque< ublas::vector<double> > mMymeasu;
		// deque with the measurements error vectors
		std::deque< ublas::vector<double> > mMyerrors;

	public:
			/** 
			 * Constructor
			 * \param vpAtmo : the aeroplanets model
			 * \param vpParams : the adjustements xml file parameters
			 */
			Adjustements(Atmo* vpAtmo,XmlParameters* vpParams):mpAtmo(vpAtmo),mpParams(vpParams)
	{
		mReport="";
		mbUseCalib=false;
	}

			/// Starts all the adjustements
			void StartAdjustements();

			/**
			 * Print the report of the fit
			 * \param vSuffix : suffix for the output file name
			 */
			void PrintReport(std::string vSuffix="");


};







#endif



