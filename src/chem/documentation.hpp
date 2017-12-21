/**
 * \ingroup Chem
 *
 * \file chem/documentation.hpp
 * \brief Documentation for the chemical reaction processes
 * Copyright G Gronoff Feb 2010
 * Last Modification : $Id: documentation.hpp 888 2010-02-22 23:21:04Z gronoff $
 */

/**
 * \page chem_doc Documentation of the Chemistry
 *
 * The chemistry is divided in three main parts: computing the chemical reactions, computing the densities, and organising the whole thing (neutral densities, call,...)
 *
 * \section chem_react The chemical reaction
 * The chemical reactions are based  (derived from) on the class ChemReact. The idea of this class is to have a logic system to define the reaction, modify its 
 * value from the xml file, and be able to allow the uncertainty for all the chemical reactions.
 *
 * The chemical reactions have an Id, and this Id is related to their position in the list of chemical reactions. (this is checked in the test, so if your modifications
 * goes through that, it works!)
 *
 * The different chemical reactions are defined in the reaclist. It takes some time to fill this database, but it ensure an improved flexibility
 *
 * \section chem_densi The photochemical  equilibrium density computation
 *
 * The CalcDensPhEq class is the basis for computing photochemical equilibriums. This class also computes the matrix of the importance of each processes, interesting
 * for sensitivity studies
 *
 *
 * \section chem_chem The Chem class
 * The Chem class allows to read the xml file to check the error and define the densities that are missing
 * It stores the different chemical reactions, and density computation functions
 *
 *
 *
 *
 *
 */


