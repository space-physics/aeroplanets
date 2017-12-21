/**
 * \ingroup file_utils 
 * \file logging.hpp
 * \brief defines  a static logging class
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: logging.hpp 1111 2010-08-12 19:43:33Z gronoff $
 *
 */
#ifndef LOGGING_HPP
#define LOGGING_HPP
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <climits>
#include <config.h> 
#include <map>
#include <algorithm>
#include <deque>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp> // Test multi-array
#include <boost/progress.hpp>

// usage of smart pointers now
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


namespace ublas = boost::numeric::ublas;

class Log;


/**
 * Performs the basic logging operations
 * It allows to define the level of the message
 * (Warning, Debug, Info...)
 * And to display it on the screen and/on on the Log file!
 *
 * On the programmer side, if you consider the default options
 * you can just consider Log as the equivalent of cout.
 *
 * You can then choose the level of warning for 
 * the file and the screen outputs....
 *
 *
 */
class Log
{
	private:

		/// The name of the priorities: they are used to display the importance... Priorities are enumerated below, here, it is juste a string array
		static const std::string msPriorityNames[];
		static bool mbIsOfstream;
	public:
		/**
		 * The constructor
		 */
		Log();

		virtual ~Log()
		{
			Log::Close();
		}

		static void Close();

		/// My log instance, very important
		static Log mL;

		/// Log for debugging
		static Log mD;

		/// Log for errors
		static Log mE;
	
		/// Log for Warnings
		static Log mW;

		/// Log for infos
		static Log mI;

		/// Log for the sections
		static Log mS;


		/** The file output stream
		 * When the program starts, msStream is defined
		 * as cerr...
		 * But when the initialization is finished, it
		 * is likely to be a file stream.
		 *
		 */
		static std::ostream* msStream;
		
		/**
		 * Priority enumeration
		 * Allows to check the priority of
		 * the message....
		 * Be careful, it is DEBUGG with 2 G because DEBUG is
		 * already defined, or not, for the Makefile...
		 */
		enum Priority
		{
			DEBUGG,
			CONFIG,
			INFO,
			SECTION,
			WARNING,
			ERROR
		};
	
		/**
		 * Initialization of the Log system
		 * \param filename : the name of the file where we write the log
		 * \param vCoutPr : the cout priority, defaulted
		 * \param vFilePr : the priority for the file
		 * \param vCerrPr : the priority for cerr
		 *
		 * Rem: the priorities are min priorities. If the 
		 * priority of the message is more important than the minimum
		 * (for example ERROR is the most important)
		 * then the message is used in respectively cout, the file, cerr...
		 */
		static void Init(std::string filename,Log::Priority vCoutPr=Log::DEBUGG,Log::Priority vFilePr=Log::CONFIG,Log::Priority vCerrPr=Log::ERROR);

		/**
		 * Returns the date with the local format
		 */
		static std::string ReturnDate();
	
		/// The priority of the following messages, see SetPriority. By default ERROR
		Priority mNextMessagePriority;

		/// The priority for the cout. By default DEBUG when debug is activated, CONFIG else
		static Priority msCoutPriority;

		/// The priority for the file log. By default CONFIG
		static Priority msFilePriority;

		/// The priority for the Cerr. By default ERROR
		static Priority msCerrPriority;



		/** The Log message to add at the end of the different files
		 * It stores the messages add by the user in the file
		 * It DOES NOT stores the errors and warning messages
		 * (Program has to put it if you want that feature)
		 *
		 * It is mainly used to remember the version of the program 
		 * in the output files
		 *
		 */
		static std::string msMessageLog;

		/**
		 *  Add a part to the message log
		 * Please consider using # at the beggining
		 *
		 * It is mainly intended to add informations about authors
		 * and versions of the different configurations files
		 * It can be used to remember important informations
		 * such as f107, bend line...
		 *
		 */
		static void AddVersionInfo(std::string vMessage);



		/**
		 * Allows to modify the priority of the messages following the call
		 * for example 
		 * \code
		 * Log::SetPriority(Log::ERROR,"Error in function ...");
		 * Log<<"Your array size is -1, Houston, we've got a problem"<<endl;
		 *
		 * Log::SetPriority(Log::DEBUG);
		 * Log::mL<<"He, je ne pensais pas arriver jusque là"<<endl;
		 * Log::mL<<"He, where am I? Hello???"<<endl;
		 * \endcode
		 * \param vPriority : the priority to consider
		 * \param vMessage : the optional message. Generally the name of the function...
		 */
		void SetPriority(Log::Priority vPriority,std::string vMessage="");


		/**
		 * Here we use the possibilities of osteam
		 * but the possibility to use std::endl is 
		 * created in the other operator overload.
		 * \param vParam : the parameter to pass
		 */
		template<class T > Log& operator<<(const T& vParam)
		{
			if(mNextMessagePriority>=Log::msCoutPriority)
			{
				std::cout<<vParam;
			}

			if(mNextMessagePriority>=Log::msFilePriority)
			{
				*Log::msStream<<vParam;
			}

			
			if(mNextMessagePriority>=Log::msCerrPriority)
			{
				std::cerr<<vParam;
			}
			return *this;
		}
		/**
		 * This one allows to use std::endl;
		 * \param pF the std::cout....
		 */
		Log& operator<<(std::ostream& (*pF)( std::ostream& ) )
		{
			if(mNextMessagePriority>=Log::msCoutPriority)
			{
				std::cout<<pF;
			}

			if(mNextMessagePriority>=Log::msFilePriority)
			{
				*Log::msStream<<pF;
			}


			if(mNextMessagePriority>=Log::msCerrPriority)
			{
				std::cerr<<pF;
			}
			return *this;

		}

};





#endif
