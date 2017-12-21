 /** 
 * \file species/documentation.hpp
 * \brief Documentation for the species
 * Copyright G Gronoff Sept 2009
 * Last Modification : $Id: documentation.hpp 1336 2011-11-08 22:11:28Z gronoff $
 *
 */



/**
 * \page ee_ionization Electron impact ionization processes
 *
 *   \section ee_crs_ionization Electron impact cross sections
 *                The electron impact cross sections are divided into 3 parts
 *                (4 if we take into account the total cross section).
 *
 *                - 1 - The elastic cross section. Which is needed to compute the fluxes
 *                - 2 - The Ionization cross section. The total cross section for
 *                           the different ionization processes.
 *                - 3 - The Excitation cross section. The total cross section for 
 *                           the different excitation processes.
 *                - (4- the total cross section is the sum of ionization and excitation)
 *
 *
 *            Total cross section and elastic cross sections are needed to compute the fluxes of electrons. Nothing more.
 *            But ionization and excitations are needed to compute the redistribution factor. Because the ionization
 *            creates an electron, it must be separated from the excitation, which only redistributes electrons.
 *            \warning Excitation and Ionization cross sections are used to compute the corrected cross section, mCrsCincm2 !!!
 *                     as in trans*, this cross section is used for later computations.
 *                     An issue now (23/9/09) is the computation of the cinex cross section. It is only used once in trans
 *                     and maybe the use of relative crs will be better....
 *            THUS, FOR COMPUTING ELECTRON DEGRADATION AND SECONDARY ELECTRONS, IONIZATION AND EXCITATION ARE NECESSARY.
 *
 *            \warning The  secondary electron function as an energy limit at 5600eV. To check
 *
 *
 *
 *            In the cross section file, the elastic cross section must be defined. 
 *            For ionization and excitation, it is different. 
 *            An unique, and total cross section can be defined, with one characteristic energy.
 *            But it is better to define the different processes, because the different threshold, for these processes
 *            means different degradation.
 *            To do that, we added two properties in the processes: \<Ionization/\> and \<Excitation/\>.
 *            If one is defined, the cross section is counted in the degradation process. (If two: error, because
 *            when the ionization is defined, the degradation of the impacting electron is computed).
 *            If nothing is defined, the process is not counted in the degradation process. This is useful when
 *            we want to compute a subproduct of another (counted) process!
 *
 *            In the original Trans* code, the ionization was only ONE process. Therefore all other processes 
 *            were counted as excitation. Unfortunately, in martian and venusian cases, some ionization processes
 *            were added as excitation(!) leading to a big mistake. The right way to do that in Trans* was to use
 *            branching ratios. Unfortunately, branching ratio are not accurate for these processes.
 *            To correct the problem of the different branching ratio thresholds, (giving different redistributions)
 *            these one were added in the process, allowing the redistribution to be more accurate.
 *
 *            In the present code, there is a kind of branching ratio, given by the \<species\>\</species\> parameters.
 *            The main goal of this parameter is not to be a branching ratio, but to give the different created species.
 *            This is why there is no threshold for the different productions => If you want a new threshold, it means that
 *            you want a new production => you create a new process!!! This is the only way to do the same in the new code.
 *            The logic behind is that the branching ratio is only a ration, so only an approximation: new cross sections
 *            could give something very different from a branching ratio.  
 *		
 *
 *
 *            \warning Please note that nothing was done to correct the fact that double ionization creates two electrons.
 *
 *
 * 		" In the good old days physicists repeated each other’s experiments, just to be sure. Today they stick to FORTRAN, so that they can share each other’s programs, bugs included. " Edsger W. Dijkstra
 * 		http://blogs.zdnet.com/Murphy/?p=568
 *
 *
 *
 *
 * \page ph_ionization Proton/Hydrogen impact ionization processes
 *   \section ph_crs_ionization Proton/Hydrogen impact cross sections
 *               The impact of proton has several possible effects on a species:
 *               	- Elastic collision: just a modification of the trajectory
 *               	- Ionization or excitation of the species: which do not modify the proton itself
 *               	- Charge exchange: the proton takes the electron, and becomes an hydrogen atom
 *               The impact of the Hydrogen atom has the same effect: the charge exchange works the other way; 
 *               or the electron can just be arrached from the hydrogen.
 *               Therefore, the cross section system for the proton impact will also take into account the hydrogen impact.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */


