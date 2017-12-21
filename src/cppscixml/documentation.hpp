/** 
 * \file cppscixml/documentation.hpp
 * \brief Documentation of the scixml class
 * Copyright G Gronoff Oct 2009
 * Last Modification $Id: documentation.hpp 872 2010-02-04 04:16:57Z gronoff $
 */


/**
 * \page scixml_doc Documentation of SciXml
 *
 * \section scixml_brief The XmlParameters Object
 *
 * This object allows to retrieve the different options 
 * and parameters in xml files.
 *
 * \subsection xmlparam_init How to initialize?
 *
 *
 * \subsection xmlparam_element How to get a string option?
 *
 *
 * \subsection xmlparam_key How to get a key (numerical or not)?
 *
 *  The keys are the elements included inside brackets, for example
 *  the type key \<value type="float" /\>.
 *
 *  GetKey allows to get the value of the key as strings:
 *  string maclef=xmlparams.GetKey("/value","type");
 *
 *  If you want a numerical value, you should use GetNKey (Get Numerical Key)
 *  for example \<machin fact="212.2"/\>
 *  double factor=0;
 *  xmlparams.GetNKey("/machin","fact",factor);
 *  The factor is now 212.2
 *  Please not that it is a template class (if you do not understand this
 *  it is not a problem, except if you want to modify the code).
 *  Therefore, you have to define the type of the value.
 *
 *
 *
 * \subsection xmlparam_vals How to get numerical values
 *
 * To get a string element, we used the Elem method, and sometimes the GetSubFile method.
 * For the array, we used the GetParamArray method.
 *
 * For numerical values, the thing is a little bit different:
 * - It is interesting to have a correction factor for values: example, we want to multiply by 1E-18 
 *   because it is a cross section...
 * - We want to have different types of array, (with common multiplicator)
 * - We can want to do a parameter error propagation!!! \ref xmlparam_erroranalysis
 *
 * A standard key for the numerical options will be << fact >> this will be used to multiply
 * the result.
 *
 * To get a single value: \<val fact="22.2"\>1\</val\>
 *
 * double maval;
 * xmlparams.GetValue("/val",maval); // maval=22.2
 * Note that an inexisting fact is, by default, equal to 1 (seems logical, it is!)
 *
 * To get a list of values as an 1D array 
 * \<val fact="12.2"\>1 2 3 4 5\</val\> or 
 * \<val fact="12.2"\>1
 * 2
 * 3
 * 4 5\</val\>
 *
 * vector<double> myarray;
 * xmlparams.Get1DArray("/val",myarray);
 * 
 * 
 * If you want to work with a 2D array
 * \<val fact0="12" fact1="42" \>
 * 1	2	3
 * 2	4	9
 * 3	6	27
 * 4	8	81
 * \</val\>
 * \warning For the 2D, you have to precise the column for the factor -> fact0, fact1 ...
 * 	    The columns starts at 0
 * 	    The MonteCarlo process is not implemented.
 * 	    fact has no effet
 *
 * 
 *
 *
 *
 * \subsection xmlparam_erroranalysis I want to perform an error analysis with a Monte Carlo technique!!!
 * 		
 *	 The process of getting the values from XmlParameters allows you to define variable
 *	 parameter, witout modifying the xml file by yourself!!! (Ideal for cluster launch!)	
 *	 For the GetValue and Get1DArray, you can add an uncertainty parameter which can
 *	 be used for the monte carlo simulations.
 *
 * 	 \warning The MonteCarlo is not implemented for the 2D array. Please contact me if you want it.
 *
 * 	 \warning The MonteCarlo model use only normal law, the mean is the value, and the sigma is you parameter.
 * 	 	 Please contact me if you want more laws! It can be done really easily (especially if it is integrated
 * 	 	 in boost random).
 *
 * 	To use the Monte Carlo system:
 *	- 1) The Monte Carlo should be activated in your XmlParameter class, by xmlparam.SetMonteCarloActive();
 *	- 2) The Monte Carlo should be activated in your element by the setMC="active" key
 *	- 3) You must have defined uncertainty="-12.5" or uncertainty="30%" (if it is not defined -> only WARNING)
 *	- 4) You can also define fact_uncertainty, as uncertainty. In that case, there is also an uncertainty on 
 *	    the factor. Useful if you want to test the uncertainty of the global value of a cross section instead
 *	    of testing the uncertainty of each parameters.
 *	\warning The uncertainty is at 1 sigma!!!
 *
 *	This double checking is very important, it allows you to disable the MonteCarlo in the global scale
 *	or in the local scale. DO NOT MAIL ME SAYING THIS IS TOO MUCH WORK, AND TOO VERBOSE IN THE LOG!
 *	=> You will thank me later.
 *
 *
 *
 *	SUMMARY:
 *	\code
 *	XmlParameters xmlparameter("myfile.xml")
 *	xmlparameter.SetMonteCarloActive();
 *	....
 *	*In myfile.xml :
 *	<value setMC="active" uncertainty="30%"> 12 322 13 </value> <!-- each value has a 30% uncertainty (1 sigma) -->
 *	<value setMC="active" uncertainty="1.4"> 12 322 13 </value> <!-- each value has a 1.4 uncertainty (1 sigma) -->
 *	<value setMC="active" fact_uncertainty="30%"> 12 322 13 </value> <!-- the multiplicative value is 1. with an uncertainty of 30% (1 sigma)-->
 *	<value setMC="active" fact="3" fact_uncertainty="30%"> 12 322 13 </value> <!-- the multiplicative value is 3. with an uncertainty of 30% (1 sigma) -->
 *	\endcode
 * 		
 *
 * 
 *
 *
 *
 *
 */


