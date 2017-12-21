/**
 * \file constphys.hpp
 * \ingroup Math
 * \brief To define physical constants
 * Copyright G Gronoff Sept 2009
 *  Last modification : $Id: constphys.hpp 1344 2011-11-10 21:06:50Z gronoff $
 */


#ifndef CONST_PHYS_HPP
#define CONST_PHYS_HPP


// Atomic mass unit in g
#define AMU 1.66054E-24
// Boltzman constant in g*m2/(K*s2)
#define BOLTZMANN 1.3806505E-20
// Boltzmann constant in SI kg*m2/(K*s2)==Joule/K
#define BOLTZMANNJ_K 1.3806505E-23


// Electron mass in kg
#define ELECTRON_MASS_KG 9.10938291E-31

// PROTON mass in kg
#define PROTON_MASS_KG 1.672621777E-27


// Electron / Proton mass ration
#define ELECTRON_PROTON_MASS_RATIO 5.446170219E-4


// k/amu in m2/(K.s2)
#define BOLTZ_DIV_AMU 8.3144727E3

// k/amu*1K/(1.m/s2) in cm (when g expressed in m/s2)

#define SCALEH_CONST 8.3144727E5


// eV in Joule
#define eV_IN_JOULE 1.60217653E-19

// BOLTMAN div eV
#define BOLTZ_DIV_eV 8.61734318E-05

// Avogadro constant
#define AVOGADRO 6.02214179E23


// Conversion of eV to erg
#define EV_TO_ERG 1.602177E-12

// Conversion of eV to nm lambda (nm) = EV_TO_NM / E(eV) = 1.24E3 / E
#define EV_TO_NM 1.239842E3

// Conversion of nm to eV E (eV) = NM_TO_EV / Lambda (nm) = 1.24E3 / Lambda (nb: it is the same than EV_TO_NM)
#define NM_TO_EV 1.239842E3

#define EV_TO_CMM1 8.065547E3

#define CMM1_TO_EV 0.1239842E-3



// c'est beaucoup plus joli de faire PI au lieu de M_PI
#define PI M_PI
#define DEG_TO_RAD ( 180. / PI )
#define RAD_TO_DEG ( PI / 180. )
// Sinon l'autre solution c'est :
//
// Que j'aime à faire apprendre un nombre utile aux sages
// 3 . 1  4   1   5           9  2      6     5  3      5
// Immortel Archimède, artiste, ingénieur,
//        8         9        7  9
// Qui de ton jugement peut priser la valeur ?
// 3   2   3         8    4      6  2      6
// Pour moi ton problème eut de pareils avantages 
// 4    3    3         8  3   2       7         9
// 3.1415926535897932384626433832795
// avec ce truc, on peut presque calculer le diamètre de l'univers  observable (R=15My) a l'atome près!




#endif
