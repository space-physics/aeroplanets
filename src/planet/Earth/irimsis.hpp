/**
 * \file irimsis.hpp
 * \brief Wrapper around Iri and Msis models
 * Copyright G Gronoff Nov 2009
 * Last Modification $Id: irimsis.hpp 831 2010-01-08 17:58:49Z gronoff $ 
 */


#ifndef IRI_MSIS_HPP
#define IRI_MSIS_HPP
#include "../planet.hpp"




extern "C"
{
	/**
	 * \ingroup Planetes
	 * Load the msis model
	 * \param modatmos : - 1 if tinf computed by msis
	 *  		     - 2 if iterating on f107
	 *  		     - 3 if iterating on flux and Ap
	 * \param neutspe : number of species (max=8, N2,O2,O,H,He,N,A,NO)
	 * \param dayno : day number
	 * \param nan : year (yy)
	 * \param Apind : index Ap
	 * \param f107bar : the mean f107
	 * \param f107 : the f107 parameter
	 * \param alat: the geo latitude
	 * \param along : the geo longitude
	 * \param UT : the universal time
	 * \param z : array of altitudes
	 * \param nalt : number of altitudes
	 * \param sortie : bool : 0 (no output)
	 * \param dens : array(8,nalt) the output densities
	 * \param fcdens : array containing 8 multiplicative factors for the density (eg 1. because our program does the modification later)
	 * \param tinf : temperature at the bottom if not imposed by the atmosphere (modatmos!=1), put tinf below -300 to compute it by msis.
	 * \param tneutre : the output neutral temperature array(nalt) 
	 * \param fctemp : correction of the temperature : 1.
	 * \param xmasdens : array(nalt) containing  the total mass density (1E-3kg/cm3)
	 * \param impress : value for the output file... put 6 : it is the std output...
	 *
	 */

	extern void atmntr_(int* modatmos,int* neutspe,float* dayno,int* nan,float* Apind,float* f107bar,float* f107,float* alat,float* along,float* UT,float* z,int* nalt,int*sortie,float* dens,float* fcdens,float* tinf,float* tneutre,float* fctemp,float* xmasdens,int* impress );

	/**
	 * \ingroup Planetes
	 * Load the iri model
	 *
	 * \param jf: 1,0 values to compute the different options (see below)
	 *
	 *
	 *
	 *         JF(1:12)      TRUE/FALSE FLAGS FOR SEVERAL OPTIONS
          JF(1)=.TRUE.[.FALSE.]   ELECTRON DENSITY IS [NOT] CALCULATED
          JF(2)=T[F]    TEMPERATURES ARE [NOT] CALCULATED
          JF(3)=T[F]    ION COMPOSITION IS [NOT] CALCULATED
          JF(4)=T[F]    B0 FROM TABLE [FROM GULYEAVA 1987]
          JF(5)=T[F]    F2 PEAK FROM CCIR [FROM URSI]
          JF(6)=T[F]    ION COMP. STANDARD [DANILOV-YAICHNIKOV-1985]
          JF(7)=T[F]    STAND. IRI TOPSIDE [IRI-79]
          JF(8)=T[F]    NMF2 PEAK MODEL [INPUT VALUES]
          JF(9)=T[F]    HMF2 PEAK MODEL [INPUT VALUES]
          JF(10)=T[F]   TE MODEL [TE-NE MODEL WITH NE INPUT]
          JF(11)=T[F]   NE STANDARD [LAY-FUNCTIONS VERSION]
          JF(12)=T[F]   MESSAGE ARE WRITTEN TO UNIT=6 [=12]

  JF(1:11)=.TRUE. GENERATES THE STANDARD IRI-90 PARAMETERS.
  IF YOU SET JF(8)=.FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK
  NMF2/M-3 OR FOF2/MHZ IN OARR(1). SIMILARLY, IF YOU SET JF(9)=
  .FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK HEIGHT HMF2/KM IN
  OARR(2). IF YOU SET JF(10)=.FALSE., THAN YOU HAVE TO PROVIDE THE
  ELECTRON DENSITY IN M-3 AT 300KM AND/OR 400KM AND/OR 600KM IN
  OARR(3), OARR(4), AND OARR(5). IF YOU WANT TO USE THIS OPTION AT
  ONLY ONE OF THE THREE ALTITUDES, THAN SET THE DENSITIES AT THE
  OTHER TWO TO ZERO.
	 * \param jmag 0 if geographic coordinates, 1 if geomagnetic coordinates
	 * \param alati : the latitude in degree (with respect to jmag)
	 * \param along : the longitude in degree (with respect to jmag)
	 * \param minusf107 : f107 as a negative number
	 * \param minusdoy : day of year as a negative number
	 * \param UTp25 : universal time plus 25
	 * \param minh : minimal height in km
	 * \param maxh : maximal height in km
	 * \param stph : height step in km : there should be less than 50 steps!
	 * \param outf : array (11,50) containing the main outputs
	 * - (warning : fortran array, begin at 1; inverted in c++, where begins at 0) outf(1,*) : electron density m-3
	 * - outf(2,*) : T neutre K
	 * - outf(3,*) : T ion K
	 * - outf(4,*) : T electron K
	 * - outf(5,*) : % ion O+
	 * - outf(6,*) : % ion H+
	 * - outf(7,*) : % ion He+
	 * - outf(8,*) : % ion O2+
	 * - outf(9,*) : % ion NO+
	 * - outf(10,*) : % cluster ions
	 * - outf(11,*) : % N+ ions
	 *
	 * \param oarr : array(30), add output parameters
	 * 
            OARR(1:30)   ADDITIONAL OUTPUT PARAMETERS         
              OARR(1) = NMF2/M-3        OARR(2) = HMF2/KM
              OARR(3) = NMF1/M-3        OARR(4) = HMF1/KM
              OARR(5) = NME/M-3         OARR(6) = HME/KM
              OARR(7) = NMD/M-3         OARR(8) = HMD/KM
              OARR(9) = HHALF/KM        OARR(10) = B0/KM
              OARR(11) =VALLEY-BASE/M-3 OARR(12) = VALLEY-TOP/KM
              OARR(13) = TE-PEAK/K      OARR(14) = TE-PEAK HEIGHT/KM
              OARR(15) = TE-MOD(300KM)  OARR(16) = TE-MOD(400KM)/K
              OARR(17) = TE-MOD(600KM)  OARR(17) = TE-MOD(1400KM)/K
              OARR(18) = TE-MOD(3000KM) OARR(19) = TE(120KM)=TN=TI/K
              OARR(20) = TI-MOD(430KM)  OARR(21) = X/KM, WHERE TE=TI
              OARR(22) = SOLAR ZENITH ANGLE/DEG
              OARR(23) = SUN DECLINATION/DEG
              OARR(24) = DIP            OARR(25) = DIP LATITUDE
              OARR(26) = MODIFIED DIP LATITUDE
              OARR(27:30) FREE

 */
	extern void iris12_(int* jf, int* jmag,float* alati,float* along,float*minusf107,int* minusdoy,float* UTp25,float* minh,float* maxh,float*stph,float* outf,float* oarr);

};


/**
 * To have the terrestrial atmosphere-ionosphere
 * \todo correct the ion and electron temperature : same as neutral temperature at low altitude!!!
 *
 */
class EarthIriMsis
{
	private:


		/// Geographic latitude in degree
		double mLocalLatDeg;
		/// Geographic longitude in degree
		double mLocalLongDeg;
		/// Subsolar latitude in degree
		double mLatDeg;
		/// Subsolar longitude in degree
		double mLongDeg;
		/// Solar declination angle in degree
		double mSdaDeg;
		/// Solar zenith angle in degree
		double mSZADeg;
		/// Geomagnetic latitude in degree
		double mMagLatDeg;
		/// Geomagnetic longitude in degree
		double mMagLongDeg;

		/**
		 * Modification of the coordinate system. The passage from the two systems is given by the  vPsiDeg parameter
		 * if vPsiDeg= mSdaDeg this code changes the subsolar coordinate to local coordinate 
		 * if vPsiDeg= -mSdaDeg this code changes the local coordinates to subsolar coordinate
		 * if vPsiDeg = 90-latitude this codes changes the Hour coordinates (proportional to alpha-deg, but with a time dependance) into azimuth-dec
		 * if vPsiDeg = -(90-latitude) this codes changes the azimuth-deg into Hour coordinates 
		 * \param vLatDeg : the latitude in degree in the first frame
		 * \param vLonDeg : the longitude in degree in the first frame
		 * \param vPsiDeg : the Difference in degree between the two poles
		 * \param rLatDeg : the latitude in degree in the final frame -> result, reference
		 * \param rLongDeg : the longitude in degree in the final frame -> result, reference
		 */
		void ChgtCoord(double vLatDeg,double vLonDeg, double vPsiDeg, double& rLatDeg,double& rLongDeg);

		/**
		 * Compute the sda thanks to the day of the year
		 * \param vDayOfYear : day of the year
		 * The output is in the private variable mSdaDegree
		 */
		void DayToSda(double vDayOfYear);

		/// The type of the atmosphere if 1, tn computed by the atmosphere
		int mModAtmo;
		/// The f107 parameter
		double mF107;
		/// The f107  parameter : mean over 3 month
		double mF107bar;
		/// The year
		int mYear;
		/// The day of the year
		int mDayOfYear;
		/// The UT parameter
		double mUT;
		/// The Tinf (<-300 : computed by msis)
		double mTinf;

		/// The Ap parameter
		double mAp;



		/// True if Msis has already been launched
		bool mbIsMsisDone;
		/// True is Iri has already been launched
		bool mbIsIriDone;

		/// Launch Msis and store the result
		void LaunchMsis(const ublas::vector<double>& vAltGridKm);

		/// Launch Iri and store the result
		void LaunchIri(const ublas::vector<double>& vAltGridKm);

		/// Store the neutral temperature
		ublas::vector<double> mTempNeutre;

		/// Store the neutral densities
		std::deque< ublas::vector<double> > mDensities;

		/// Store the msis grid
		ublas::vector<double> mMsisAltGridKm;
		
		/// Store the iri grid
		ublas::vector<double> mIriAltGridKm;

		
		/// Store the iri densities
		std::deque< ublas::vector<double> > mIonDensities;

		/// Store the iri temperature (n,i,e)
		std::deque< ublas::vector<double> > mIriTemp;
		




	public:
		/**
		 * Init
		 * \param vModAtmo : model of atmosphere (see the atmntr_ parameters)
		 * \param vF107 : the f107 parameter
		 * \param vF107bar : the mean of f107 over 3 month
		 * \param vUT : the universal time in hour
		 * \param vTinf : the temperature at the bottom of the msis simu, if computed by msis : <-300
		 * \param vAp : the Ap index value
		 * \param vYear : the year
		 * \param vDayOfYear : the number of the day of the year
		 */
		EarthIriMsis(int vModAtmo,double vF107,double vF107bar,int vYear,int vDayOfYear,double vUT,double vTinf,double vAp);

		/**
		 * Init thanks to the local (geographical) coordinates
		 * \param vLocalLatDeg : the geographic latitude 
		 * \param vLocalLongDeg : the geographic longitude 
		 *
		 */
		void InitLocalCoords(double vLocalLatDeg,double vLocalLongDeg);

		/**
		 * Init thanks to the subsolar coordinates
		 * \param vLatDeg : the subsolar latitude 
		 * \param vLongDeg : the subsolar longitude 
		 *
		 */
		void InitSubsolarCoords(double vLatDeg, double vLongDeg);


		/**
		 * Returns the msis neutral atmosphere
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return vector of vector : the neutral atmosphere
		 */
		std::deque< ublas::vector<double> > ReturnNeutralAtmo(const ublas::vector<double>& vAltitudeGridKm);
		/**
		 * Returns the  Iri ionospheric composition and densities
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * The last value is the electron density!
		 * \return vector of vector : the ionosphere
		 */
		std::deque< ublas::vector<double> > ReturnIono(const ublas::vector<double>& vAltitudeGridKm);
		
		
		/**
		 * Returns the Msis  neutral temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnNTemp(const ublas::vector<double>& vAltitudeGridKm);
		/**
		 * Returns the IRI ion  temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnITemp(const ublas::vector<double>& vAltitudeGridKm);
		/**
		 * Returns the IRI electron temperatures
		 * \param vAltitudeGridKm : the final altitude grid. Log interpolation is done on that grid
		 * \return the neutral temperatures
		 */
		ublas::vector<double> ReturnETemp(const ublas::vector<double>& vAltitudeGridKm);

		double ReturnLatDeg()
		{
			return mLatDeg;
		}
		double ReturnLongDeg()
		{
			return mLongDeg;
		}
		double ReturnSZADeg()
		{
			return mSZADeg;
		}


};




#endif


