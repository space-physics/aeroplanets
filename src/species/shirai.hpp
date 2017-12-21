/**
 * \file shirai.hpp
 * \brief Defines the shirai class, used to translate the parameters of Shirai into actual cross sections 
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: shirai.hpp 881 2010-02-22 06:21:37Z gronoff $
 */


#ifndef SHIRAI_HPP
#define SHIRAI_HPP
#include <cppscixml/scixml.hpp>


/**
 * Reads the Shirai kind Cross sections.
 * Now, three papers allows us to work with these cross sections:
 * Shirai 2001 about CO2,CO and H20
 * Shirai 2002 about CH4, ...
 * Tabata 2006 about N2, N2+
 *
 */
class Shirai
{
	private:
		/// Cross section multiplicator
		double mSigmacm2;
		/// Rydberg energy 
		double mErkeV;


		/** The first kind elementary function
		 * \param vE : the energy vector
		 * \param vC1 : the first parameter
		 * \param vC2 : the second parameter
		 * \return the cross section in cm2
		 */
		ublas::vector<double> F1(ublas::vector<double> vE,double vC1,double vC2);
		/** The second kind elementary function
		 * \param vE : the energy vector
		 * \param vC1 : the first parameter
		 * \param vC2 : the second parameter
		 * \param vC3 : the third parameter
		 * \param vC4 : the fourth parameter
		 * \return the cross section in cm2
		 */
		ublas::vector<double> F2(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4);

		/** The third kind elementary function
		 * \param vE : the energy vector
		 * \param vC1 : the first parameter
		 * \param vC2 : the second parameter
		 * \param vC3 : the third parameter
		 * \param vC4 : the fourth parameter
		 * \param vC5 : the fifth parameter
		 * \param vC6 : the sixth parameter
		 * \return the cross section in cm2
		 */
		ublas::vector<double> F3(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4,double vC5,double vC6);

		/** The fourth kind elementary function (defined by taking into account the last S type of each paper)
		 * \param vE : the energy vector
		 * \param vC1 : the first parameter
		 * \param vC2 : the second parameter
		 * \param vC3 : the third parameter
		 * \param vC4 : the fourth parameter
		 * \return the cross section in cm2
		 */
		ublas::vector<double> F4(ublas::vector<double> vE,double vC1,double vC2,double vC3,double vC4);

		/** The first CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S1CH4(ublas::vector<double> vE);
		/** The 2 section CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S2CH4(ublas::vector<double> vE);
		/** The 3 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S3CH4(ublas::vector<double> vE);
		/** The 4 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S4CH4(ublas::vector<double> vE);
		/** The 5 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S5CH4(ublas::vector<double> vE);
		/** The 6 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S6CH4(ublas::vector<double> vE);
		/** The 7 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S7CH4(ublas::vector<double> vE);
		/** The 8 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S8CH4(ublas::vector<double> vE);
		/** The 9 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S9CH4(ublas::vector<double> vE);
		/** The 10 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S10CH4(ublas::vector<double> vE);
		/** The 11 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S11CH4(ublas::vector<double> vE);
		/** The 12 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S12CH4(ublas::vector<double> vE);
		/** The 13 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S13CH4(ublas::vector<double> vE);
		/** The 14 CH4 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S14CH4(ublas::vector<double> vE);


		/** The first CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S1CO2(ublas::vector<double> vE);
		/** The 2 section CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S2CO2(ublas::vector<double> vE);
		/** The 3 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S3CO2(ublas::vector<double> vE);
		/** The 4 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S4CO2(ublas::vector<double> vE);
		/** The 5 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S5CO2(ublas::vector<double> vE);
		/** The 6 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S6CO2(ublas::vector<double> vE);
		/** The 7 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S7CO2(ublas::vector<double> vE);
		/** The 8 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S8CO2(ublas::vector<double> vE);
		/** The 9 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S9CO2(ublas::vector<double> vE);
		/** The 10 CO2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S10CO2(ublas::vector<double> vE);

		/** The first N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S1N2(ublas::vector<double> vE);
		/** The 2 section N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S2N2(ublas::vector<double> vE);
		/** The 3 N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S3N2(ublas::vector<double> vE);
		/** The 4 N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S4N2(ublas::vector<double> vE);
		/** The 5 N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S5N2(ublas::vector<double> vE);
		/** The 6 N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S6N2(ublas::vector<double> vE);
		/** The 7 N2 main  function
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> S7N2(ublas::vector<double> vE);
		/** Return the interpolation of the equation in the case of a
		 * species in the CH4 paper. The threshold correction is not made here
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> ReturnCH4(ublas::vector<double> vE); 
		/** Return the interpolation of the equation in the case of a
		 * species in the CO2 paper. The threshold correction is not made here
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> ReturnCO2(ublas::vector<double> vE); 
		/** Return the interpolation of the equation in the case of a
		 * species in the N2 paper. The threshold correction is not made here
		 * \param vE : the energy vector
		 * \return the cross section in cm2
		 */
		ublas::vector<double> ReturnN2(ublas::vector<double> vE); 

		/** Loads the necessary parameters for the computation
		 * \return true if the cross section is correctly loaded
		 */ 
		bool LoadCrs();
		/// vector of the values necessary for the computation
		ublas::vector<double> mAvalues;
		/// The minimum defined energy
		double mMinEnerkeV;
		/// The maximum defined energy
		double mMaxEnerkeV;
		/// The threshold energy
		double mThreshkeV;

		/// The equation type
		unsigned mEquationType;
		/// The article Id, to select if we launch CH4, CO2, or N2 type of equation
		std::string mArticleId;

		/// The parameters file: to extract the system
		XmlParameters* mpParams;
		/// The node to extract
		TiXmlNode* mNode;
	public:

		/** Initialize the object.
		 * \param vpParams : the xml file parameter pointer
		 * \param vNode : the node where we extract the shirai cross section
		 * \param vThresholdeV : the threshold in eV
		 */
		Shirai(XmlParameters* vpParams,TiXmlNode* vNode,double vThresholdeV);
		/** Returns the result, with the energy range in parameter
		 * \param vEeV : the energy vector  in eV
		 * \return the cross section in cm2
		 */
		ublas::vector<double> GetCrs(ublas::vector<double> vEeV);






};




#endif
