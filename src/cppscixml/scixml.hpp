/**
 * \defgroup file_utils File utilities
 * \file scixml.hpp
 * \brief Defines the scixml class : to retrieve xml data
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: scixml.hpp 1523 2012-06-07 20:43:29Z gronoff $
 */


#ifndef SCIXML_HPP
#define SCIXML_HPP

#include <math/mathfunction.hpp>
#include <tinyxpath/xpath_static.h>
#include <string>
#if BOOST_VERSION < 103500
	#include <boost/filesystem/convenience.hpp>
#else
	#include <boost/filesystem.hpp>
#endif
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/progress.hpp>
#include "error.hpp"


/**
 * \ingroup file utils
 */
namespace ggfileutils
{
	/**
	 * \ingroup file_utils
	 * To join two path, or a file and its path.
	 * Helps with the '/' problem
	 * \param vPath : the initial path
	 * \param vFilename : the filename or pathname to join
	 * \return final path
	 */

	std::string JoinPath(std::string vPath,std::string vFilename);


	/**
	 * \ingroup file_utils
	 * Helps to retrieve the folder name of a file.
	 * \param vPathfile : complete filename with its path. It can be relative or absolute
	 * \return : path relative or absolute, depending on the parameter
	 */
	std::string FilenameToPath(std::string vPathfile);

	/**
	 * \ingroup file_utils
	 * finds if the file (or the path) exists
	 * \param vPathfile : complete path to the file to check
	 * \return true if exists
	 */
	bool FileExists(std::string vPathfile);


	/**
	 * \ingroup file_utils
	 * Return the absolute path to the file
	 * \param vPathfile : complete path to the file to check
	 * \return absolute path
	 */
	std::string FullPath(std::string vPathfile);
};

// Ok, this is not the best thing wrote here...
using namespace ggfileutils;
/**
 * \ingroup file_utils
 * \~english Class to work with xml files. Allows to retrieve values...
 *
 * \~francais Classe pour travailler avec les fichiers xml: 
 *permet de récupérer des tableaux, et des valeurs
 */
class XmlParameters
{
	private:
		///  The name of the file to analyse \~francais Le nom du fichier a analyser
		std::string mFilename;
		/// The tixml pointer to our document
		TiXmlDocument* mpDoc;

		/// True if the error analysis through MonteCarlo is activated
		bool mMonteCarloActivated;

		/**
		 * Get the error parameters.
		 * Returns false if the error system is not active for that function
		 * (globally or locally not activated)
		 * Return true and fills the rPercent or rRange parameter with the actual sigma or range
		 * The fact parameter has no effect since we consider that the uncertainty
		 * is defined for the final value. The fact multiply the  value+MC factor... (but it is in the following).
		 * \param xPath : the path to the node
		 * \param rPercent : the percentage of uncertainty, returned as a reference if defined ("30%" ->30)
		 * \param rRange : the range of uncertainty, returned as a reference  ("10" -> equivalent to +/- 10)
		 * \param vNumber : the position, if necessary
		 * \param vUncertaintyKeyName : the name of the key where the uncertainty is searched. Defining it allows to choose a "factor_uncertainty", wich allows to put MC on the multiplicator of the value. (if necessary)
		 */
		bool GetMCParams(std::string xPath,double & rPercent,double&rRange,unsigned vNumber=0,std::string vUncertaintyKeyName="uncertainty");


		/**
		 * Get the error parameters.
		 * Returns false if the error system is not active for that function
		 * (globally or locally not activated)
		 * Return true and fills the rPercent or rRange parameter with the actual sigma or range
		 * The fact parameter has no effect since we consider that the uncertainty
		 * is defined for the final value. The fact multiply the  value+MC factor... (but it is in the following).
		 * \param vNode : the node where we search for MonteCarlo
		 * \param xPath : the path to the node
		 * \param rPercent : the percentage of uncertainty, returned as a reference if defined ("30%" ->30)
		 * \param rRange : the range of uncertainty, returned as a reference  ("10" -> equivalent to +/- 10)
		 * \param vNumber : the position, if necessary
		 * \param vUncertaintyKeyName : the name of the key where the uncertainty is searched. Defining it allows to choose a "factor_uncertainty", wich allows to put MC on the multiplicator of the value. (if necessary)
		 */
		bool GetMCParams(TiXmlNode* vNode,std::string xPath,double & rPercent,double&rRange,unsigned vNumber=0,std::string vUncertaintyKeyName="uncertainty");

	protected:
		/**
		 * Apply the monte carlo modification to rValue
		 * \param rValue : the value to be modified, is also the mean
		 * \param vPercent : the percent at 1sigma, negative if not active
		 * \param vRange : the range at 1 sigma, negative if not active
		 */

		template<class T>  void ApplyMC(T & rValue, double vPercent, double vRange)
		{
			double sigma=vRange;
			if(vPercent>0)
			{// We are in the percent case
				sigma=vPercent*rValue/100.;
			}
			double tmp=static_cast<double>(rValue);
			double ran=MathRandom::GetNormal(tmp,sigma);
			//tmp+=ran;
			rValue=static_cast<T>(ran);

		}

		/**
		 * Apply the monte carlo modification to the vector rvalue, based on the vUncert vector
		 * check if it is a Range (if not, it is a percent), and if it is a factor
		 * \param rValue : the vector of value to be modified, is also the mean
		 * \param vUncert : the vector for the uncertainties (default: range, else factor)
		 * \param vbIsRange : true if it is a range
		 * \param vbIsFactor : true if the uncertainty is computed through a factor
		 */
		template<class T> void VectorApplyMC(ublas::vector<T> & rValue, ublas::vector<double> vUncert, bool vbIsRange, bool vbIsFactor)
		{
			assert(rValue.size() == vUncert.size());
			if(vbIsFactor)
			{
				double factor = 0;
				ApplyMC(factor,-1,1.);// We use a normal law, centered at 0, with sigma = 1
			
				if(vbIsRange)
				{
					for(size_t i =0; i < rValue.size(); ++i)
					{
						rValue[i] += vUncert[i] * factor;
					}
				}else{
					for(size_t i =0; i < rValue.size(); ++i)
					{
						rValue[i] *=  (1 + vUncert[i] * factor / 100.);
					}
				}
			}else{
				if(vbIsRange)
				{
					for(size_t i =0; i < rValue.size(); ++i)
					{
						//rValue[i] += vUncert[i] * factor;
						ApplyMC(rValue[i], -1, vUncert[i]);
					}
				}else{
					for(size_t i =0; i < rValue.size(); ++i)
					{
						//rValue[i] *=  (1 + vUncert[i] * factor / 100.);
						ApplyMC(rValue[i], vUncert[i], -1);
					}
				}

			}

		}


	public:
		/**
		 * Apply the monte carlo modification to rValue, only if the monte carlo system is activated!
		 * \param rValue : the value to be modified, is also the mean
		 * \param vPercent : the percent at 1sigma, negative if not active
		 * \param vRange : the range at 1 sigma, negative if not active
		 */

		template<class T>  void ApplyMCP(T & rValue, double vPercent, double vRange)
		{
			if(mMonteCarloActivated)
			{
				ApplyMC(rValue, vPercent, vRange);
			}
		}

		/**
		 * Apply the monte carlo modification to the vector rvalue, based on the vUncert vector
		 * check if it is a Range (if not, it is a percent), and if it is a factor
		 * Only if the Monte Carlo system is activated!
		 * \param rValue : the vector of value to be modified, is also the mean
		 * \param vUncert : the vector for the uncertainties (default: range, else factor)
		 * \param vbIsRange : true if it is a range
		 * \param vbIsFactor : true if the uncertainty is computed through a factor
		 */
		template<class T> void VectorApplyMCP(ublas::vector<T> & rValue, ublas::vector<double> vUncert, bool vbIsRange, bool vbIsFactor)
		{
			if(mMonteCarloActivated)
			{
				VectorApplyMC(rValue, vUncert, vbIsRange, vbIsFactor);
			}
		}

		/**
		 * Applies the uncertainty in percent on a vector in argument
		 * \param vVal : the vector
		 * \param vPercent : the uncertainty in percent
		 * \return the vector corresponding to vVal with uncertaintites
		 */
		template<class T> ublas::vector<T> ReturnMC(ublas::vector<T> vVal,double vPercent)
		{
			if(mMonteCarloActivated)
			{
				for(unsigned i=0;i<vVal.size();i++)
				{
					ApplyMC(vVal[i],vPercent,0);
				}
			}
			return vVal;
		}
		/**
		 * The constructor : needs our file
		 * \param vFilename : the file to analyse with the xml parser
		 */ 
		XmlParameters(std::string vFilename);

		/**
		 * The destructor
		 */
		~XmlParameters();

		/// Set the Monte Carlo process as True -> The result are modified randomly
		void SetMonteCarloActive();
		/// Disable the Monte Carlo process.
		void SetMonteCarloInactive();

		/// Returns the status of the MonteCarlo process.
		bool GetMonteCarlo()
		{
			return mMonteCarloActivated;
		}


		/**
		 * Get the path to the file given in parameter with respect to the path for the present parameter file.
		 * Example: we have a main parameter file in
		 * ./superdossier/mainparameter.xml
		 * in this file, a seconde parameter file is defined:
		 * ./subdossier/secondparameter.xml
		 * if we open directly this defined path, we do not open
		 * ./superdossier/subdossier/secondparameter.xml but
		 * ./subdossier/secondparameter.xml which is wrong.
		 *
		 * But if we open get_subfilename("./subdossier/secondparameter.xml")
		 * then we have the desired file.
		 * \param vFilename : the filepath relative to the main parameter file.
		 * \return  string : the relative path, but relative to the execution folder.
		 *
		 */
		std::string SubFileName(std::string vFilename);



		/** Check if the xpath exists
		  \param  xPath : the path to check
		  \return  bool : true if the path exists
		  */
		bool Exists(std::string xPath);

		/** Check if the xpath, from the node exists
		  \param vNode : the initial node
		  \param  xPath : the path to check
		  \return  bool : true if the path exists
		  */
		bool Exists(TiXmlNode* vNode,std::string xPath);


		/** Check if the xpath exists
		 * if it does not exists -> print the non existing parameter
		 * 			 -> exit the program
		 \param  xPath : the path to check
		 \param vMsg : the message to print, more than just the name of the missing parameter
		 */
		void ExistsOrDie(std::string xPath,std::string vMsg);

		/** Check if the xpath, from the node exists
		 * if it does not exists -> print the non existing parameter
		 * 			 -> exit the program
		 \param vNode : the initial node
		 \param  xPath : the path to check
		 \param vMsg : the message to print, more than just the name of the missing parameter
		 */
		void ExistsOrDie(TiXmlNode* vNode,std::string xPath,std::string vMsg);





		/**
		  Directly processes a xpath, gives directly the result string.
		  Can be interesting when the xpath uses some functions.
		  As example "cout(/xml/numbers/)"
		  \param xPath : the path to process
		  */

		std::string Process(std::string xPath);

		/**
		  Directly processes a xpath,  which originates from a node,
		  gives directly the result string.
		  Can be interesting when the xpath uses some functions.
		  As example "cout(/xml/numbers/)"
		  \param vNode : the origine node
		  \param xPath : the path to process
		  */
		std::string Process(TiXmlNode* vNode, std::string xPath);



		/** Return the number of elements checking the path
		  \param  xPath : the path to check
		  \return  int : the number of elements checking the path
		  */
		unsigned Numbers(std::string xPath);

		/** Return the number of elements checking the path
		  \param vNode : the origine node
		  \param  xPath : the path to check
		  \return  int : the number of elements checking the path
		  */
		unsigned Numbers(TiXmlNode* vNode,std::string xPath);
		/**
		  Get the element, as a string, checking the path at the place number.
		  Put only the path, no / at the end
		  \param  xPath : the path  to the element
		  \param  vNumber : the place number for the element
		  \return std::string : the element
		  */
		std::string Elem(std::string xPath,unsigned vNumber=0);
		/**
		  Get the element, as a string, checking the path at the place number.
		  Put only the path, no / at the end
		  \param vNode : the origine node
		  \param  xPath : the path  to the element
		  \param  vNumber : the place number for the element
		  \return std::string : the element
		  */
		std::string Elem(TiXmlNode* vNode,std::string xPath,unsigned vNumber=0);


		/** Get the Node that satisfy the xpath
		 * \param xPath : the path to the node
		 * \param vNumber : the number of the element
		 * \return the node
		 */
		TiXmlNode* GetNode(std::string xPath,unsigned vNumber=0);

		/** Get the Nodes that satisfy the xpath
		 * \param xPath : the path to the node
		 * \return  vector with all the nodes
		 */
		std::vector<TiXmlNode*> GetNodes(std::string xPath);

		/** Get the Node that satisfy the xpath
		 * \param vNode : the origine node
		 * \param xPath : the path to the node
		 * \param vNumber : the number of the element
		 * \return the node
		 */
		TiXmlNode* GetNode(TiXmlNode* vNode,std::string xPath,unsigned vNumber=0);

		/** Get the Nodes that satisfy the xpath
		 * \param vNode : the origine node
		 * \param xPath : the path to the node
		 * \return  vector with all the nodes
		 */
		std::vector<TiXmlNode*> GetNodes(TiXmlNode* vNode,std::string xPath);


		/*
		   get_value -> template int double...)
		   get_array -> template int double...
		   get_nkey -> template int double
		   get_key -> 
		   */


		/**
		  Count the number of params
		  \param  xPath : the path where the params are
		  \param  vParam : the parameter
		  \param  vNumberatt : if numberatt>number of params, return 0
		  \return unsigned: the number of params
		  */
		unsigned NbParams(std::string xPath,std::string vParam,unsigned vNumberatt=0);

		/**
		  Count the number of params from a node
		  \param vNode : the origine node
		  \param  xPath : the path where the params are
		  \param  vParam : the parameter
		  \param  vNumberatt : if numberatt>number of params, return 0
		  \return unsigned: the number of params
		  */
		unsigned NbParams(TiXmlNode* vNode,std::string xPath,std::string vParam,unsigned vNumberatt=0);

		/**
		 * Get the value of the key param
		 * \param xPath : the path where the param is. No / at the end
		 * \param vParam : the parameter
		 * \param  vNumber : the number, 0 by default. 
		 */
		std::string GetKey(std::string xPath,std::string vParam,unsigned vNumber=0);
		/**
		 * Get the value of the key param
		 \param vNode : the origine node
		 * \param xPath : the path where the param is. No / at the end
		 * \param vParam : the parameter
		 * \param  vNumber : the number, 0 by default. 
		 */
		std::string GetKey(TiXmlNode* vNode,std::string xPath,std::string vParam,unsigned vNumber=0);

		/**
		 * Find if the key exists
		 * \param xPath : the path where the param is. No / at the end
		 * \param vParam : the parameter
		 */
		bool KeyExists(std::string xPath,std::string vParam);
		/**
		 * Find if the key exists, with a node
		 \param vNode : the origine node
		 * \param xPath : the path where the param is. No / at the end
		 * \param vParam : the parameter
		 */
		bool KeyExists(TiXmlNode* vNode,std::string xPath,std::string vParam);


		/**
		 * get the key into numerical value
		 * \param xPath : the path to the key
		 * \param vParam : the key name
		 * \param  rResu : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetNKey(std::string xPath,std::string vParam,T& rResu,unsigned vNumber=0)
		{
			std::string value=GetKey(xPath,vParam,vNumber);
			std::istringstream istr(value);
			istr>>rResu;
		}
		/**
		 * get the key into numerical value
		 * \param vNode : the origine node
		 * \param xPath : the path to the key
		 * \param vParam : the key name
		 * \param  rResu : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetNKey(TiXmlNode* vNode,std::string xPath,std::string vParam,T& rResu,unsigned vNumber=0)
		{
			std::string value=GetKey(vNode,xPath,vParam,vNumber);
			std::istringstream istr(value);
			istr>>rResu;
		}
		/**
		 * get a value into numerical value
		 * \param xPath : the path to the key
		 * \param  rMachin : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetValue(std::string xPath,T& rMachin,unsigned vNumber=0)
		{
			std::string value=Elem(xPath,vNumber);
			std::istringstream istr(value);
			istr>>rMachin;

			double percent,range;
			if(GetMCParams(xPath,percent,range,vNumber))
			{
				//Log::SetPriority(Log::CONFIG);
				Log::mL<<"Monte Carlo activated for "<<xPath<<std::endl;
				ApplyMC(rMachin,percent,range);
			}
			bool multiply=false;
			T facteur=static_cast<T>(1.);
			if(NbParams(xPath,"fact",vNumber) > 0)
			{
				//	T facteur;
				GetNKey(xPath.c_str(),"fact", facteur, vNumber);
				multiply=true;
			}

			if(GetMCParams(xPath,percent,range,vNumber,"fact_uncertainty"))
			{
				ApplyMC(facteur,percent,range);
				multiply=true;
			}
			if(multiply)
			{
				rMachin*=facteur;
			}
		}
		/**
		 * get a value into numerical value, if the path does not exists, let the default value
		 * \param xPath : the path to the key
		 * \param  rMachin : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetValueOrDefault(std::string xPath,T& rMachin,unsigned vNumber=0)
		{
			if(Exists(xPath))
				GetValue(xPath,rMachin,vNumber);
		}

		/**
		 * get a value into numerical value
		 * \param vNode : the origine node
		 * \param xPath : the path to the key
		 * \param  rMachin : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetValue(TiXmlNode* vNode,std::string xPath,T& rMachin,unsigned vNumber=0)
		{
			std::string value=Elem(vNode,xPath,vNumber);
			std::istringstream istr(value);
			istr>>rMachin;
			double percent,range;
			if(GetMCParams(vNode,xPath,percent,range,vNumber))
			{
				//Log::SetPriority(Log::CONFIG);
				Log::mL<<"Monte Carlo activated for "<<xPath<<std::endl;
				ApplyMC(rMachin,percent,range);
			}
			bool multiply=false;
			T facteur=static_cast<T>(1.);
			if(NbParams(vNode, xPath,"fact",vNumber)>0)
			{
				//	T facteur;
				GetNKey(vNode,xPath.c_str(),"fact",facteur,vNumber);
				multiply=true;
			}

			if(GetMCParams(vNode,xPath,percent,range,vNumber,"fact_uncertainty"))
			{
				ApplyMC(facteur,percent,range);
				multiply=true;
			}
			if(multiply)
			{
				rMachin*=facteur;
			}
		}

		/**
		 * get a value into numerical value, if the path does not exists, let the default value
		 * \param vNode : the origine node
		 * \param xPath : the path to the key
		 * \param  rMachin : The result, numerical (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void GetValueOrDefault(TiXmlNode* vNode,std::string xPath,T& rMachin,unsigned vNumber=0)
		{
			if(Exists(vNode,xPath))
				GetValue(vNode,xPath,rMachin,vNumber);
		}


		/**
		 * get a value into 1D array 
		 * \param xPath : the path to the key
		 * \param  rResult : The result, vector<numerical> (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void Get1DArray(std::string xPath,ublas::vector<T>& rResult,unsigned vNumber=0)
		{
			std::string value=Elem(xPath,vNumber);
			if(!MathString::LitToutString(value,rResult))
			{
			//	Log::SetPriority(Log::DEBUGG);// Because throw an error!
				Log::mE<<"Error while reading table"<<xPath<<std::endl;
				Error err("Get1DArray","lit string","Error while reading table "+xPath);
				throw err;
			}
			double percent,range;
			if(GetMCParams(xPath,percent,range,vNumber))
			{
				//Log::SetPriority(Log::CONFIG);
				Log::mL<<"Monte Carlo activated for "<<xPath<<std::endl;

				for(unsigned i=0;i<rResult.size();i++)
				{
					ApplyMC(rResult[i],percent,range);
				}
			}
			bool multiply=false;
			T facteur=static_cast<T>(1.);
			if(NbParams(xPath,"fact",vNumber)>0)
			{
				//	T facteur;
				GetNKey(xPath.c_str(),"fact",facteur,vNumber);
				multiply=true;
			}

			if(GetMCParams(xPath,percent,range,vNumber,"fact_uncertainty"))
			{
				ApplyMC(facteur,percent,range);
				multiply=true;
			}
			if(multiply)
			{
				for(unsigned i=0;i<rResult.size();i++)
				{
					rResult[i]*=facteur;
				}
			}

		}

		/**
		 * get a value into 1D array 
		 * \param vNode : the origine node
		 * \param xPath : the path to the key
		 * \param  rResult : The result, vector<numerical> (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */
		template<class T> void Get1DArray(TiXmlNode* vNode,std::string xPath,ublas::vector<T>& rResult,unsigned vNumber=0)
		{
			std::string value=Elem(vNode,xPath,vNumber);
			if(!MathString::LitToutString(value,rResult))
			{
				//Log::SetPriority(Log::DEBUGG);// Because throw an error!
				Log::mE<<"Error while reading table"<<xPath<<std::endl;
				Error err("Get1DArray (node)","lit string","Error while reading table "+xPath);
				throw err;
			}
			double percent,range;
			if(GetMCParams(vNode,xPath,percent,range,vNumber))
			{
				//Log::SetPriority(Log::CONFIG);
				Log::mL<<"Monte Carlo activated for "<<xPath<<std::endl;

				for(unsigned i=0;i<rResult.size();i++)
				{
					ApplyMC(rResult[i],percent,range);
				}
			}
			bool multiply=false;
			T facteur=static_cast<T>(1.);
			if(NbParams(vNode,xPath,"fact",vNumber)>0)
			{
				//	T facteur;
				GetNKey(vNode,xPath.c_str(),"fact",facteur,vNumber);
				multiply=true;
			}

			if(GetMCParams(vNode,xPath,percent,range,vNumber,"fact_uncertainty"))
			{
				ApplyMC(facteur,percent,range);
				multiply=true;
			}
			if(multiply)
			{
				for(unsigned i=0;i<rResult.size();i++)
				{
					rResult[i]*=facteur;
				}
			}
		}


		/**
		 * Get a string array (vector) of value into an 1D array
		 * \param xPath : the path
		 * \param vNumber : the number of the xpath
		 * \return the string vector
		 */
		std::vector<std::string> GetParamArray(std::string xPath,unsigned vNumber=0);

		/**
		 * Get a string deque of value into an 1D array
		 * \param xPath : the path
		 * \param vNumber : the number of the xpath
		 * \return the string vector
		 */
		std::deque<std::string> GetParamDeQue(std::string xPath,unsigned vNumber=0);


		/**
		 * Get a string array of value into an 1D array
		 * \param xPath : the path
		 * \param vNode : the origine node
		 * \param vNumber : the number of the xpath
		 * \return the string vector
		 */
		std::vector<std::string> GetParamArray(TiXmlNode* vNode,std::string xPath,unsigned vNumber=0);


		/**
		 * Get a string array of value into an 1D array
		 * \param xPath : the path
		 * \param vNode : the origine node
		 * \param vNumber : the number of the xpath
		 * \return the string vector
		 */
		std::deque<std::string> GetParamDeQue(TiXmlNode* vNode,std::string xPath,unsigned vNumber=0);



		/**
		 * Get a 2D array. The fact parameter is valid, but not the MC.
		 * The number of the column should be put after the fact, example fact0=14 fact1=12 
		 * \param xPath : the path to the key
		 * \param  rResu : The result, vector<vector<numerical>> (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */

		template<class T> void Get2DArray(std::string xPath,ublas::matrix<T> & rResu,unsigned vNumber=0)
		{
			std::string value=Elem(xPath,vNumber);
			//std::vector< std::vector<T> > result;
			ublas::matrix<T> result;
			if(!MathString::Lit2DString(value,result))
			{
				//Log::SetPriority(Log::DEBUGG);// Because throw an error!
				Log::mE<<"Error while reading table"<<xPath<<std::endl;
				Error err("Get2DArray","lit string","Error while reading table "+xPath);
				throw err;
			}

			unsigned nligne=result.size1();
			if(nligne==0)
			{// tableau vide
				rResu=result;
				return;
			}
			//unsigned ncol=result.size2();
			// On reorganise le tableau
			/*	for(unsigned c=0;c<ncol;++c)
				{
				std::vector<T> tmp;
				for(unsigned l=0;l<nligne;++l)
				{
				tmp.push_back(result[l][c]);
				}
				rResu.push_back(tmp);

				}*/
			rResu=trans(result);

			// On travaille maintenant sur le facteur numerique
			// On ne le met pas dans l'autre boucle en
			// cas de demande de tableau non numerique
			for(unsigned c=0;c<rResu.size1();++c)
			{
				std::string fac="fact"+ntostr(c);
				if(NbParams(xPath,fac,vNumber)>0)
				{
					T facteur;
					GetNKey((xPath.c_str()),fac,facteur,vNumber);
					for(unsigned i=0;i<rResu.size2();i++)
					{
						rResu(c,i)*=facteur;
					}

				}

			}

		}

		/**
		 * Get a 2D array. The fact parameter is valid, but not the MC.
		 * The number of the column should be put after the fact, example fact0=14 fact1=12 
		 * \param vNode : the origine node
		 * \param xPath : the path to the key
		 * \param  rResu : The result, vector<vector<numerical>> (it is a template function)
		 * \param vNumber : the number of the  xpath
		 */

		template<class T> void Get2DArray(TiXmlNode* vNode,std::string xPath,ublas::matrix<T> & rResu,unsigned vNumber=0)
		{
			std::string value=Elem(vNode,xPath,vNumber);
			ublas::matrix<T>  result;
			if(!MathString::Lit2DString(value,result))
			{
				//Log::SetPriority(Log::DEBUGG);// Because throw an error!
				Log::mE<<"Error while reading table"<<xPath<<std::endl;
				Error err("Get2DArray (node)","lit string","Error while reading table "+xPath);
				throw err;
			}

			unsigned nligne=result.size1();
			if(nligne==0)
			{// tableau vide| void array
				rResu=result;
				return;
			}
			/*
			   unsigned ncol=result[0].size();
			// On reorganise le tableau | reorganization of the array
			for(unsigned c=0;c<ncol;++c)
			{
			std::vector<T> tmp;
			for(unsigned l=0;l<nligne;++l)
			{
			tmp.push_back(result[l][c]);
			}
			rResu.push_back(tmp);

			}
			*/
			rResu=trans(result);

			// On travaille maintenant sur le facteur numerique
			// On ne le met pas dans l'autre boucle en
			// cas de demande de tableau non numerique
			//
			// Now, we are on a numerical array...
			for(unsigned c=0;c<rResu.size1();++c)
			{
				std::string fac="fact"+ntostr(c);
				if(NbParams(vNode,xPath,fac,vNumber)>0)
				{
					T facteur;
					GetNKey(vNode,(xPath.c_str()),fac,facteur,vNumber);
					for(unsigned i=0;i<rResu.size2();i++)
					{
						rResu(c,i)*=facteur;
					}

				}

			}

		}

		/**
		 * Return the value of the xpath as a valid filename
		 * (relative to the execution folder
		 * \param xPath : the path where the filename is
		 */
		std::string GetSubFileName(std::string xPath);

		/**
		 * Return the value of the xpath as a valid filename
		 * (relative to the execution folder
		 * \param vNode : the parent node
		 * \param xPath : the path where the filename is
		 */
		std::string GetSubFileName(TiXmlNode* vNode,std::string xPath);



		/*       _\|/_
			 (o o)
			 +----oOO-{_}-OOo----------------------------------------------------+
			 |                                                                   |
			 |Here, we will work on the automatic pointer to XmlParameter values.|
			 |                                                                   |
			 +------------------------------------------------------------------*/

	private:
		/// List of the open files
		static std::vector< std::string > mFileOpenList;
		/// Map giving the different opened files XmlParameter pointer
		static std::map< std::string, XmlParameters*> mParamMap;
		/// Map giving the number of files wich are using the pointer
		static std::map< std::string, unsigned > mOpenedMap;
	public:

		/**
		 * \page newprocess_xmlparameters Dynamic allocation of XmlParameters*
		 *
		 * If a XmlParameter object must be opened a lot a times, by different objects, and destroyed at the end, it is very difficult to  manage to pointer!!!
		 * The stupid solution is to open several times the object.
		 * This is not the supidest solution in fact!!! Because
		 * we are used of destroying at the end, and there is no
		 * problem to know if another process uses, or not, the file pointer!
		 * 
		 * The second solution is to use a pre-defined politic of pointer destruction... Great on the paper (whaaa la super traduction), but in reality, which politic I used last time? (Especially true if you have several projects with different politics in parallel).
		 *
		 * The third solution is to copy the pointer without deleting it, and to use a garbage collector to do the deleting work (I do not think, the GC does it...). Welcome to the java world ;-). Modern garbage collectors does not suffer from their ancestors problems, and sometimes, they are more efficient than hand freeing memory!
		 * I write this page the 2009 oct 2-> C++Oxx (or whatever its name) is still a project, and I want to run my code on existing computers, and not only on an experimental machine. So I can't use the C++Oxx garbage collector!
		 *(But in a near (??? C++Oxx) future, the GC will be implemented by default in  c++... So I do not want to break the compatibility with it)
		 * The GC is a huge think to implement in a code, so, considering the future of C++, and my use of memory, I do not implement one (difficult to choose). Moreover, I have a much better solution for our problem!!!!
		 *
		 * The solution used here looks like a mini garbage collector, and allows to retrieve the opened pointers without need for copy!!!!!!!
		 * The pointers are stored in a container, when we have the filename, we can retrieve the pointer!
		 * Each time we no more use the pointer, we detach it!
		 * when nobody is using the pointer (attach counter == 0)
		 * then the pointer is deleted!!!
		 * Thanks to that system, when my file was opened before and is not closed, I can retrieve the pointer without knowing where it was opened!!!
		 *
		 * Here is an example of the process.
		 * \code
		 *
		 XmlParameters* file1=XmlParameters::AttachFileParameter("test4.xml");
		 XmlParameters* file2=XmlParameters::AttachFileParameter("test4.xml");
		 XmlParameters* file3=XmlParameters::AttachFileParameter("test4.xml");

		 file3->DetachFileParameter();
		 file1->DetachFileParameter();
		 file2->DetachFileParameter();
		 XmlParameters* file4=XmlParameters::AttachFileParameter("test4.xml");
		 file4->DetachFileParameter();

		 *\endcode
		 * In the real life file1...n can be in different objects!
		 *
		 */

		/**
		 * Allows to attach a pointer to the parameter file
		 * Create the pointer if it does not exists
		 * \warning no thread safe
		 * \param vFilename : the name of the file
		 *
		 * Please note that this is a static function
		 */
		static XmlParameters* AttachFileParameter(std::string vFilename);
		/** Allows to detach the pointer to the parameter file. 
		 * if no more class is attached, the pointer is destroyed
		 * Please note that this is NOT a static function...
		 */
		void DetachFileParameter();

		/// Check if there is attached parameters! return true if so.
		static bool DoINeedAGarbageCollector();

};







#endif
