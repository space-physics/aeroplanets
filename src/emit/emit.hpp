/**
 * \file emit.hpp
 * \brief Define the emit class: it allows to compute and write the spectrum.
 * Copyright G Gronoff March 2010
 * Last Modification : $Id: emit.hpp 1191 2010-11-23 21:22:19Z gronoff $
 * \defgroup Emit Emissions
 */

#ifndef EMIT_HPP
#define EMIT_HPP
#include <observ/geoobs.hpp>
#include <chem/chem.hpp>

/**
 * The emission class: allows to compute the number of Rayleigh of an emission
 * Then the emission profile is used to compute the integrated value. And then,
 * the synthetic spectrum. This is an abstract class
 * \ingroup Emit
 */
class Emit
{
	protected:

		/// The filename for the column emission file
		std::string mColumnFileName;
		/// The filename for the profile emission file (limb observation)
		std::string mProfileFileName;
		/// Transforms the emission spectrum into limb view \param vPath : the geometry object
		void EmitSpectrumToLimb(boost::shared_ptr<GeoObs> vPath);
		/// Integrates the emit spectrum
		void EmitSpectrumIntegrate();

		/// Write the column file
		void WriteColumn();
		/// Write the profile file \param vPath : the geometry object
		void WriteProfile(boost::shared_ptr<GeoObs> vPath);

		/// The altitude grid in Km
		ublas::vector<double> mAltGridKm;


		/// If the profile system is used: modified by the vPath object
		bool mbPathUsed;

		/// Filename to write extra informations
		std::string mExtraInfoFilename;
		/// To write extra informations
		virtual void WriteExtraInfo()=0;

		/// To use the absorption for some lines.
		bool mbUseAbsorption;

		/// The id of the emitting specie
		SpecieId mEmittingSpecie;


		/**
		 * Compute the absorption in  cm from the absorption lists, and the chemistry
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 */
		void ComputeAbsorption(std::map< std::string, std::map<double,double> > vAbsorptions, boost::shared_ptr<Chem> vChem);




	public:
		/// The name of the reaction, initialized in the SetId
		std::string mName;

		/// The Id of the reaction
		unsigned mId;
		/// Returns the id \return id of the chemical reaction
		unsigned GetId()
		{
			return mId;
		}
	
		/// The wavelength of the emission in nm
		std::deque<double> mLambdanm;
		/// The  total integrated column emission in R
		double mIntegratedColEmitR;

		/// The integrated against the wavelength, in R
		std::map<double,double> mEmitSpectrumIntegratedR;
		

		/**
		 * The total emission in phcm_3s_1
		 */
		ublas::vector<double> mColEmitphcm_3s_1;
		/**
		 * The emission in phcm_3s_1 (with division in terms of wavelength)
		 */
		std::map<double, ublas::vector<double> > mEmitSpectrumcm_3s_1;


		/// The integrated profile, in R
		ublas::vector<double> mProfileIntegratedR;

		/// The integrated profile in R (with wavelength)
		std::map<double, ublas::vector<double> > mProfileSpR;

		/** The absorption in cm2/cm3, with respect to the altitude grid
		 * if the vector is not set (size != size of the column) then
		 * the absorption is not computed
		 */
		std::map<double, ublas::vector<double> > mColumnAbsorption_cm;



		/** The constructor
		 * \param vAltGridKm the altitude grid
		 * \param vId : the Id
		 * \param vName : the name of the emitting species
		 * \param vState : the state of the emitting species
		 */
		Emit(ublas::vector<double> vAltGridKm,unsigned vId,std::string vName,std::string vState):mAltGridKm(vAltGridKm),mEmittingSpecie(vName,vState),mId(vId)
		{
			mbPathUsed=false;// If the vPath object is correctly defined, it will automatically turn into true when emitspectrumtolimb is called
			mbUseAbsorption=false;
		}

		/**
		 * The destructor
		 */
		virtual ~Emit()
		{
		}


		/**
		 * Accessor to the Id
		 */
		SpecieId GetSpecieId()
		{
			return mEmittingSpecie;
		}

		/**
		 * Abstract class for the computation of the reaction
		 * \param vChem : shared pointer to the chemistry. It is used to compute the densities when needed
		 * \param vPath : shared pointer to the path: necessary to compute the limb emissions
		 * \param vbUseAbsorption : set if the absorption is used
		 * \param vAbsorptions : the map species, lambda (nm), absorption in cm_2
		 * \param vColumnFileName : if different than "", the name of the file for the profile of the emission
		 * \param vProfileFileName : the name of the file for the limb observation.
		 * \param vExtraInfoFilename : the name of the file for writing extra informations.  
		 * This function fills basically mEmitSpectrumcm_3s_1 and mColEmitphcm_3s_1
		 * then, it calls EmitSpectrumIntegrate to fill mIntegratedColEmitR, mEmitSpectrumIntegratedR
		 * and finally, it calls EmitSpectrumToLimb
		 */
		virtual void ComputeReaction(boost::shared_ptr<Chem> vChem,boost::shared_ptr<GeoObs> vPath,bool vbUseAbsorption,std::map< std::string, std::map<double,double> > vAbsorptions, std::string vColumnFileName="",std::string vProfileFileName="",std::string vExtraInfoFilename="")=0;

		/**
		 * Class to give informations about the emission computed
		 */
		virtual std::string Info()=0;


		/**
		 * Returns the volume emission rate (cm_3s_1) for that emission
		 * \param vPos : the position of the emission: if <0, it is the total emission, if > 0 it is determined by the position in the maps (mLamdanm)
		 */
		ublas::vector<double> GetVER(int vPos);
		/**
		 * Returns the limb emission rate (R) for that emission
		 * \param vPos : the position of the emission: if <0, it is the total emission, if > 0 it is determined by the position in the maps (mLamdanm)
		 */
		ublas::vector<double> GetLimb(int vPos);
};


#endif






