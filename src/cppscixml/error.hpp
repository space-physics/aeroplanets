/** 
 * \file error.hpp
 * \brief Defined and implements the Error class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: error.hpp 1111 2010-08-12 19:43:33Z gronoff $
 *
 */


#ifndef ERROR_HPP
#define ERROR_HPP

#include <iostream>
#include <string>
/**
 * Class to throw exceptions
 * Example :
 * \code
 *   // there is an error in the function InitError
 *   Error err("InitError","beginning","Just to test throw")
 *   throw err
 * \endcode
 *
 * Then you can catch the error by doing
 * \code
 * catch(Error&err)
 * {
 * 	err.Affiche();
 * }
 * \endcode
 *
 */
class Error
{
	public:
		// The function where the error is discovered
		std::string mFunction;
		// The position of the error (it can be function specific, like missing pathà
		std::string mPosition;
		// The message to the end-user : how to define parameters to avoid the error.
		std::string mMessage;


		/**
		 * The error is constructed-defined
		 * \param vFunction : function where discovered
		 * \param vPosition : for the mPosition
		 * \param vMessage : for the mMessage
		 */
		Error(std::string vFunction,std::string vPosition,std::string vMessage): mFunction(vFunction),mPosition(vPosition),mMessage(vMessage){}

		void Affiche()
		{
			Log::mE<<"================================================"<<std::endl;
			Log::mE<<" ------- Error in function "<<mFunction<<std::endl;
			Log::mE<<" ------- This error is caused at the position : "<<std::endl<<mPosition<<std::endl;
			Log::mE<<" -----------------------------------------------"<<std::endl;
			Log::mE<<" ------- And  the cause is :"<<std::endl<<mMessage<<std::endl;
			Log::mE<<"================================================"<<std::endl;
		}

};

#endif
