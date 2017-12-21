/**
 * \ingroup Chem
 * \file reaclist.hpp
 * \brief Defines the implementation of the different chemical reactions. Inherit from ChemReact.
 * These reactions will be used in the Chem class to compute the densities, and
 * in the emission class to compute the emissions
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: reaclist.hpp 1270 2011-06-06 14:05:50Z gronoff $
 * \defgroup Reacts Reactions Group for all the chemical reactions
 */


#ifndef REACTION_CHEM_LIST_HPP
#define REACTION_CHEM_LIST_HPP
#include "reaction.hpp"




/**
 * \ingroup Chem
 * \ingroup Reacts
 * The  Einstein coefficient for the emission of O(1D) 630nm
 *
 */
class Reac0 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac0(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 * The  Einstein coefficient for the emission of O(1D) 636nm
 *
 */
class Reac1 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac1(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 * The  Einstein coefficient for the emission of O(1D) 639nm
 *
 */
class Reac2 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac2(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1D) + O2 -> O2 + O
 *
 */
class Reac3 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac3(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1D) + O -> O + O
 *
 */
class Reac4 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac4(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1D) + N2 -> O + N2
 *
 */
class Reac5 : public ChemReact
{
	private:
		void SetId()
		{
			mId=5;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac5(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1D) + e -> O + e
 *
 */
class Reac6 : public ChemReact
{
	private:
		void SetId()
		{
			mId=6;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac6(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Reacts
 * \ingroup Chem
 *  O(1D) + CO2 -> O + CO2
 *
 */
class Reac7 : public ChemReact
{
	private:
		void SetId()
		{
			mId=7;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac7(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1D) + CO -> O + CO
 *
 */
class Reac8 : public ChemReact
{
	private:
		void SetId()
		{
			mId=8;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac8(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};




/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O2+ + e -> O(1D) + O*
 *  (Thermic electrons)
 *  The coefficient is the coefficient of O(1D) production. It is the sum of the reactions
 *  O2+ + e -> O(1D) + O*
 *  O2+ + e -> 2 O(1D) 
 */
class Reac9 : public ChemReact
{
	private:
		void SetId()
		{
			mId=9;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac9(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};




/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O + e -> O(1D)
 *  (Thermic electrons)
 */
class Reac10 : public ChemReact
{
	private:
		void SetId()
		{
			mId=10;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac10(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N+ + O2 -> O(1D)
 */
class Reac11 : public ChemReact
{
	private:
		void SetId()
		{
			mId=11;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac11(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO+ + e -> O(1D) + C
 *  (thermal ion)
 */
class Reac12 : public ChemReact
{
	private:
		void SetId()
		{
			mId=12;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac12(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};





/**
 * \ingroup Chem
 * \ingroup Reacts
 * The  Einstein coefficient for the emission of O(1S) 557.7nm
 * \todo find uncertainty
 */
class Reac13 : public ChemReact
{
	private:
		void SetId()
		{
			mId=13;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac13(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 * The  Einstein coefficient for the emission of O(1S) 297.2nm
 * \todo find uncertainty
 */
class Reac14 : public ChemReact
{
	private:
		void SetId()
		{
			mId=14;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac14(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};




/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + O2 -> O2 +  O(1D) | O
 */
class Reac15 : public ChemReact
{
	private:
		void SetId()
		{
			mId=15;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac15(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + O -> 2 O 
 */
class Reac16 : public ChemReact
{
	private:
		void SetId()
		{
			mId=16;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac16(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + CO -> CO + O 
 */
class Reac17 : public ChemReact
{
	private:
		void SetId()
		{
			mId=17;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac17(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + CO2 -> CO2 + O | O(1D) 
 */
class Reac18 : public ChemReact
{
	private:
		void SetId()
		{
			mId=18;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac18(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + N2 -> N2 + O  
 */
class Reac19 : public ChemReact
{
	private:
		void SetId()
		{
			mId=19;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac19(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + e -> e + O  
 */
class Reac20 : public ChemReact
{
	private:
		void SetId()
		{
			mId=20;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac20(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O2+ + e -> O1S  
 */
class Reac21 : public ChemReact
{
	private:
		void SetId()
		{
			mId=21;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac21(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O2+ + N -> O1S  + NO + 
 *  Frederick/Kopp needs a verification, see G Gronoff PhD... Gronoff et al 2008 
 */
class Reac22 : public ChemReact
{
	private:
		void SetId()
		{
			mId=22;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac22(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};




/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O(1S) + e -> e + O(1D)  
 */
class Reac23 : public ChemReact
{
	private:
		void SetId()
		{
			mId=23;
		}
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac23(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2+ + O -> N(2D)  
 */
class Reac24 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac24(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2+ + e -> N(2D)  
 */
class Reac25 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac25(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  NO+ + e -> N(2D)  
 */
class Reac26 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac26(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N+ + O2 -> N(2D)  + O2+
 */
class Reac27 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac27(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + O2 ->  NO + O | O 1D
 */
class Reac28 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac28(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + O -> ?
 */
class Reac29 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac29(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + e -> N
 */
class Reac30 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac30(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) Emission
 */
class Reac31 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac31(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};





/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + O2+ -> 
 */
class Reac32 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac32(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + O+ -> 
 */
class Reac33 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac33(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N(2D) + CO2 -> NO + CO 
 */
class Reac34 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac34(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2(A3S)  -> N2 + VK
 */
class Reac35 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac35(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2(A3S) + O  ->  N2  + O | O1S
 */
class Reac36 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac36(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2(A3S) + CO  ->  N2  + CO 
 */
class Reac37 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac37(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2(A3S) + CO2  ->  N2  + CO2 
 */
class Reac38 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac38(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2(A3S) + O2 ->  N2  + O2 
 */
class Reac39 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 * \todo find the uncertainty
		 */
		Reac39(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO2+ + e- ->  CO(A1PI | a3Pi)  + O 
 */
class Reac40 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac40(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  -> hv cameron bands
 *   
 */
class Reac41 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac41(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  - CO2 quenching
 *   
 */
class Reac42 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac42(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  - CO quenching
 *   
 */
class Reac43 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac43(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  - NO quenching
 *   
 */
class Reac44 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac44(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  - N2 quenching
 *   
 */
class Reac45 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac45(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO(a3Pi)  - O2 quenching
 *   
 */
class Reac46 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac46(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) - O quenching
 *   
 */
class Reac47 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac47(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) - N2 quenching
 *   
 */
class Reac48 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac48(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) - e quenching
 *   
 */
class Reac49 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac49(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) -  7320 emission
 *   
 */
class Reac50 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac50(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) -  7330 emission
 *   
 */
class Reac51 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac51(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) -  247 emission
 *   
 */
class Reac52 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac52(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O+(2P) -  CO2 quenching
 *   
 */
class Reac53 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac53(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  N2
 *   
 */
class Reac54 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac54(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};
/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  CO2
 *   
 */
class Reac55 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac55(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  O
 *   
 */
class Reac56 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac56(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  CO
 *   
 */
class Reac57 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac57(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  He
 *   
 */
class Reac58 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac58(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  He
 *   
 */
class Reac59 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac59(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  e-
 *   
 */
class Reac60 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac60(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  O++ -  H
 *   
 */
class Reac61 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac61(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2++ -  N2
 *   
 */
class Reac62 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac62(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2++ -  O
 *   
 */
class Reac63 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac63(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2++ -  CO2
 *   
 */
class Reac64 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac64(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};

/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2++ -  e dissociative recombination
 */
class Reac65 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac65(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  N2++ autodissociation
 */
class Reac66 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac66(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO2++ autodissociation
 */
class Reac67 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac67(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO2++  + e-  dissociative recombination
 */
class Reac68 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac68(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};


/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO2++  + CO2
 */
class Reac69 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac69(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};



/**
 * \ingroup Chem
 * \ingroup Reacts
 *  CO2++  + O
 */
class Reac70 : public ChemReact
{
	public:
		/** The constructor \param vpParam Pointer to the xml file
		 */
		Reac70(XmlParameters* vpParam);
		/** Function to get the reaction rate value, for all the altitudes at the same time!
		 * \param vTe : the electron temperature
		 * \param vTi : the ion temperature
		 * \param vTn : the neutral temperature
		 */
		ublas::vector<double> GetReactionRate(ublas::vector<double> vTe,ublas::vector<double> vTi,ublas::vector<double> vTn);
};




#endif
