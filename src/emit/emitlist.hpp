/**
 * \file emitlist.hpp
 * \brief Defines the different emissions 
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emitlist.hpp 1241 2011-03-25 20:15:03Z gronoff $
 */


#ifndef EMIT_LIST_HPP
#define EMIT_LIST_HPP
#include "emit.hpp"


/**
 * Emission of the Vegard kaplan Bands
 * N2(A3S)
 *
 *
 */

class EmitN2A3S : public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the N2A3S density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitN2A3S(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();

};




/**
 *
 * Emission of the O(1S) species : green line 557.7 and UV 297.2
 *
 */
class EmitO1S:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the O1S density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitO1S(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};


/**
 * Emission of the Vegard kaplan Bands
 * N(2D)
 *
 *
 */

class EmitND : public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the N2D density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitND(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();

};


/**
 *
 * Emission of the O(1D) species : triplet at 630 nm
 *
 */
class EmitO1D:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the O1D density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitO1D(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};


/**
 *
 * Emission of the CO(a3Pi) species : Bands around 200 nm
 *
 */
class EmitCOa3Pi:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the  extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitCOa3Pi(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};


/**
 *
 * Emission of the CO2+(B) species : doublet at 289 nm
 *
 */
class EmitCO2pB:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the  extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitCO2pB(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};
/**
 *
 * Emission of the CH(A2D) species : doublet at 289 nm
 *
 */
class EmitCHA2D:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the  extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitCHA2D(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};
/**
 *
 * Emission of the O+(2P) species : 2 lines (doublet) in the red 7320 and 7330 and one line in the UV 2470 
 *
 */
class EmitOplus2P:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the O1S density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitOplus2P(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};

/**
 *
 * Emission of the O((3p)3P) species 844.6 nm (3P*)
 *
 */
class EmitO3p3P:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the O1S density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitO3p3P(ublas::vector<double> vAltGridKm);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};


/**
 *
 * Emission  template for the allowed transitions
 *
 */
class EmitAllowed:public Emit
{
	private:
		/// To write the extra information
		void WriteExtraInfo();
		/// Matrix containing the output of the O1S density extra information
		ublas::matrix<double> mExtraInfo;
		/// Deque containing the name of these extra informations
		std::deque< std::string > mExtraInfoName;
		/// String for the extra warnings
		std::string mWarnings;
		/// map containing  the emitting frequencies and their branching ratio
		std::map<double,double> mFreqBratio;
	public:
		/** The constructor
		 * \param vAltGridKm the altitude grid
		 */
		EmitAllowed(ublas::vector<double> vAltGridKm, unsigned vId,std::string vName,std::string vState, std::map<double,double> vFreqBratio);

		/*
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		*/
 		void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions,std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="");
		
	
		/**
		 * Function to give informations about the emission computed
		 */
		std::string Info();
};





#endif
