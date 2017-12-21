/**
 * \file marsbinary.hpp
 * \brief Defines the martian atmosphere from binary lecture
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: marsbinary.hpp 927 2010-03-12 17:23:03Z gronoff $
 */
#ifndef MARS_BINARY_HPP
#define MARS_BINARY_HPP
#include "datamarinerviking.hpp"
#include "../planet.hpp"


class DumpData
{
	private:
		std::string mFilename;
		unsigned mSizeOf;
		void GetVal(std::ifstream&  rIfs,double& rTmp);
		void GetVect(std::ifstream& rIfs,unsigned vSize,ublas::vector<double>& vVect);
		void CopyVect(ublas::matrix<double>& rMat,unsigned line,const ublas::vector<double>& vVec,unsigned start,unsigned end);
	
	public:
		unsigned mNbAlt;
		unsigned mNbCol;
		unsigned mIAnnee;
		unsigned mIMois;
		unsigned mIJour;
		unsigned mIHeure;
		unsigned mIMinute;
		double mSeconde;
		unsigned mIPas;
		double mF107;
		double mChiDeg;
		ublas::vector<double> mAltitudeKm;
		ublas::matrix<double> mDensitecm_3;
		ublas::matrix<double> mVz;
		ublas::matrix<double> mVy;
		ublas::matrix<double> mTemperatureK;
		ublas::matrix<double> mFluxQ;
		ublas::matrix<double> mNeutreDensitycm_3;
		ublas::vector<double> mNeutreTemperatureK;
		ublas::matrix<double> mProdioncm_3s_1;
		ublas::matrix<double> mLossioncm_3s_1;
		ublas::vector<double> mHeat;
		ublas::vector<double> mChampB;
		ublas::matrix<double> mSupra;
		ublas::vector<double> mUy;
		ublas::vector<double> mUz;

		DumpData(std::string vFilename,unsigned vSize=8);
		virtual ~DumpData();
		void ReadTableau();

		/**
		 * Returns the neutral atmosphere
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return vector of vector : the neutral atmosphere
		 */
		std::deque< ublas::vector<double> > ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);
		/**
		 * Returns the  ionospheric composition and densities
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * The last value is the electron density!
		 * \return vector of vector : the ionosphere
		 * H+,O+,O2+,CO2+
		 */
		std::deque< ublas::vector<double> > ReturnIono(const ublas::vector<double>& vAltitudeGridKm);
		
		
		/**
		 * Returns the  neutral temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm);
		
		/**
		 * Returns the  electron temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnETemp(const ublas::vector<double>& vAltitudeGridKm);
			/**
		 * Returns the  ion temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnITemp(const ublas::vector<double>& vAltitudeGridKm);
			
};

#endif

/*

int main(int argc, char**argv)
{

	cout<<"Hello"<<endl;

	cout<<sizeof(int)<<" "<<sizeof(uint32_t)<<endl;
	cout<<sizeof(float)<<" "<<sizeof(double)<<endl;
	cout<<sizeof(char)<<endl;
	DumpData dump("filein.dat"); 
	dump.ReadTableau();
	return 0;
}
*/

