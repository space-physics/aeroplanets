/**
 * \file main/documentation.hpp
 * \brief The main documentation
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: documentation.hpp 667 2009-10-19 13:05:45Z gronoff $
 *
 */

#ifndef DOCUMENTATION_HPP
#define DOCUMENTATION_HPP

#include "timeline.hpp"
/**
 * \mainpage
 * \section mainpage_project The project.
 * Hello world, you are on the main page of the Trans++ project.
 * This \ref project is intended to understand the impact of energy
 * precipitation on the atmospheres. It is based on the old Transsolo 
 * code, written in fortran, which objective was to understand the 
 * impact of photoionization and electron impact ionization on the
 * earth upper atmosphere. It was the improve for other planets
 * (Trans*), and included, for the earth, proton precipitation (TransCube).
 * It was also coupled with a fluid model (TransCar, or TransMars).
 *
 * Unfortunately, fortran code grows in complexity, and in un-readability
 * very easily. And the Trans* codes were buggy!!!
 * Considering the amount of time needed to implement very small features
 * (for example, monthes to implement a new specie, because size of 
 * array had to be modified, and specie-dependant parameters were to be
 * added...) it seemed preferable to implement the same features of the
 * Trans* into a new kind a OO model. D'où Trans++!
 * 
 * \section mainpage_modifications I want to modify the code!!!
 *
 * Good idea!!!
 *
 * A living code needs always improvements: new features, bug corrections.
 * For the Trans++ code, it is easy to look at these new features:
 * - more particle impact models (protons, cosmics... it was done in Trans*!)
 *                 but you can look at IR processes, bent geometry, radiative transfer...
 * - more planets! It is easy to add a new atmosphere model! But it also need:
 * - more species. Ionization, excitation... All these things are done through
 *                 cross sections! But are the cross sections updated
 * - more physical processes : fluid model can use Trans++, as atmosphere retrieval
 *                 models.
 * - Gimme ideas please!!!! Trans++ can also be used for things I have no ideas!
 *   		  
 *
 * But, before putting your hands into dust, please consider some parameters
 *  - \ref versioning_system VERY IMPORTANT, FOR YOUR OWN SAFETY. Use the version control system.
 *  		Firstly, it allows you to see the activity on the project, and to participate!
 *  		Secondly, it avoids you to lose all your improvements in an unsafe rm -fr
 *  		And MAINLY, it allows you to retrieve bug corrections!
 *  - \warning a bug report system is not implemented now, but I  wish to have one soon
 *  - \ref coding_style Please help me to have a standardized code. It has it advantages.
 *  		I don't want this new code growing like it ancestors!!!
 *  - Please CODE IN C++ : the C++ has so many advantage among the Fortran.
 *              Of course, some call of Fortran routines are made in the code (because
 *              I am lasy to re-write all!).  BUT, I may modify it by running 
 *              f2c (fortran to c).
 *              So, If you have to re-use some very very very specific blackbox, ok.
 *              But if you start a new project -> C++, if the you are integrating
 *              a code that manage some species properties -> C++ (I worked hard
 *              to have something consistent in species properties management!).
 *              If you want some advices, curses... In C++, it is there \ref cpp_curses .
 *
 *
 * 
 * Then, you can consider the \ref architecture of the code.
 *  
 *
 *
 *
 *
 */




/**
 * \page project Trans++
 * Trans++ is a project to implement the Trans* code in C++
 * The main objective is to have a versatile programme, easy to modify
 * with less bugs coming from memory problems.
 *
 * The difficult point is, of course, that you have to know C++.
 * \section Install
 * The installation is very easy: it is based on the autotools.
 * If you got the program through the svn server, you have
 * to launch ./prepare_conf.sh to set up the configure program.
 * Then ./configure ; make ; make check ;  make install
 * In the majority of cases, you will modify the program
 * therefore, the make install process is not necessary (you will
 * execute in local). But you have to consider the developing 
 * options : 
 * ./configure --enable-maintainer-mode --enable-debug
 * the first options makes an auto-configuration when the Makefile.am
 * and configure.ac are modified (when you add a new library).
 * The second one set up the debug options, and the DEBUG macro.
 * It is very interesting to put an "#ifdef DEBUG" when you are debugging
 * moreover, you can use asserts.
 *
 * \section Architecture
 * A more detailed view of the architecture can be foud at \ref architecture
 * The basis of the code is the Atmo object.
 * Atmo embedded a planet object, which set up the local conditions.
 * It has also the photoionization model.
 *
 * One important thing, is that all these objects need the XmlParameters one.
 * This parameter reads the XML page where all the options are defined .
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


/**
 * \page architecture Architecture
 * The architecture of the project is based on Atmo. This object has the charge to set up the other models, and
 * to call the computation functions (embedded in other objects, of course).
 * 
 * \section architecture_class Architecture of the different classes and namespaces
 *
 *
 * \subsection architecture_math The Math
 *		The math functions are given under the Math module.
 *		Basically, it can be splitted in several parts:
 *		- MathString : reads the string values, and gives numerical outputs
 *		             (integer, double... or arrays -in fact vectors-)
 *		- MathFunction : gives interpolation, integration and mean functions
 *		- MathGrid : allows to create specific grids, necessary to the 
 *				numerical discretisation
 *		- MathFlux : allows to fill the grids with typical functions.
 *			     It is mainly used to fill electron precipitation
 *			     input fluxes. Hence the name. 
 *
 *
 *
 *
 * \subsection architecture_xml_basis The Xml Basis
 *             The Trans++ code uses the tinyxpath library to read the xml files.
 *             To avoid the difficult task of learning that library and the xPath
 *             fundamentals, all these points are hidden behind the XmlParameters 
 *             class.
 *             To have a better understanding of the processes, reading the 
 *             cppscixml/test.cpp  test file is a good idea.
 *
 *
 *             The File Utilities namespace has also some functions to work
 *             with the filenames. JoinPath allows to join to path, without
 *             problems with the number of /. FilenameToPath allows to 
 *             retrieve the path name of a file (please read the source 
 *             code if you want to do a lot of work with that function,
 *             it only extracts the path from the filename!).
 *             FileExists allows to check the presence of a file.
 *             And FullPath allows to retrieve the absolute path of a
 *             file if the parameter name is relative to the execution dir.
 *
 * \subsection architecture_xml_errors The Exceptions (Errors)
 * 		Trans++ uses the standard try throw catch model.
 * 		It throw two types of errors, Error and MathError. The last
 * 		one is reserved for the math library.
 * 		Both error classes have an Affiche() member function that
 * 		allows the user to understand why the exception were raised.
 *
 *
 *
 *
 *
 *
 */

/**
 *
 * \page coding_style Coding style
 *
 * \section why_a_coding_style  Why setting up  a coding style ?
 * Coding style is very important for several reasons.
 * 	- A good looking code is easy to read (look at python! Indentation implies that no bracket are necessary)
 * 	- A good coding style avoid confusions (between members and functions for example...)
 * 	- A smart coding style avoid basis mistakes, for example, to avoid units confusions : units should be embedded inside  the variable.
 * 	- And also an automatic documentation avoid problems during re-reading of the code. For example, put a \f$\vec{F}=m\vec{a} \f$ is
 *         much more useful
 *
 * \section coding_rules_for_transpp Coding rules for Trans++
 * The coding rules are basically the one found in http://www.possibility.com/Cpp/CppCodingStandard.html
 * I put here a resume of the more important choices. (Sometimes, there are too
 * many rules in the precedent document)
 * 
 * \subsection coding_rules_name Make Names Fit
 *
Names are the heart of programming. In the past people believed knowing someone's true name gave them magical power over that person. If you can think up the true name for something, you give yourself and the people coming after power over the code. Don't laugh!

A name is the result of a long deep thought process about the ecology it lives in. Only a programmer who understands the system as a whole can create a name that "fits" with the system. If the name is appropriate everything fits together naturally, relationships are clear, meaning is derivable, and reasoning from common human expectations works as expected.

If you find all your names could be Thing and DoIt then you should probably revisit your design.

\subsection coding_rules_class_name Class Names

    * Name the class after what it is. If you can't think of what it is that is a clue you have not thought through the design well enough.
    * Compound names of over three words are a clue your design may be confusing various entities in your system. Revisit your design. Try a CRC card session to see if your objects have more responsibilities than they should.
    * Avoid the temptation of bringing the name of the class a class derives from into the derived class's name. A class should stand on its own. It doesn't matter what it derives from.
    * Suffixes are sometimes helpful. For example, if your system uses agents then naming something DownloadAgent conveys real information. 

\subsection  coding_rules_function_name Method and Function Names

    * Usually every method and function performs an action, so the name should make clear what it does: CheckForErrors() instead of ErrorCheck(), DumpDataToFile() instead of DataFile(). This will also make functions and data objects more distinguishable.

      Classes are often nouns. By making function names verbs and following other naming conventions programs can be read more naturally.

    * Suffixes are sometimes useful:
          - Max - to mean the maximum value something can have.
          - Cnt - the current count of a running count variable.
          - Key - key value. 

      For example: RetryMax to mean the maximum number of retries, RetryCnt to mean the current retry count.

    * Prefixes are sometimes useful:
          - Is - to ask a question about something. Whenever someone sees Is they will know it's a question.
          - Get - get a value.
          - Set - set a value. 

      For example: IsHitRetryLimit.

\subsection  coding_rules_units Include Units in Names
If a variable represents time, weight, or some other unit then include the unit in the name so developers can more easily spot problems. For example:
\code
uint32 mTimeoutMsecs;
uint32 mMyWeightLbs;
\endcode 
Better yet is to make a variable into a class so bad conversions can be caught.
GG -> If negative exponents, please use _ :
\code
vector<double> mResuFluxPhcm_2s_1
\endcode


\subsection coding_rules_uppercaseabbr No All Upper Case Abbreviations

    * When confronted with a situation where you could use an all upper case abbreviation instead use an initial upper case letter followed by all lower case letters. No matter what. 

\subsubsection cr_upperjustif Justification

    * People seem to have very different intuitions when making names containing abbreviations. It's best to settle on one strategy so the names are absolutely predictable.

      Take for example NetworkABCKey. Notice how the C from ABC and K from key are confused. Some people don't mind this and others just hate it so you'll find different policies in different code so you never know what to call something. 

\subsubsection cr_upperexamp Example
\code
   class FluidOz             // NOT FluidOZ
   class NetworkAbcKey       // NOT NetworkABCKey
\endcode
\subsection cr_classname Class Names

    * Use upper case letters as word separators, lower case for the rest of a word
    * First character in a name is upper case
    * No underbars ('_') 

\subsubsection cr_cn_justif Justification

    * Of all the different naming strategies many people found this one the best compromise. 

\subsubsection cr_cn_exa Example
\code
   class NameOneTwo
  
   class Name
\endcode 
\subsection cr_cln  Class Library Names

    * Now that name spaces are becoming more widely implemented, name spaces should be used to prevent class name conflicts among libraries from different vendors and groups.
    * When not using name spaces, it's common to prevent class name clashes by prefixing class names with a unique string. Two characters is sufficient, but a longer length is fine. 

\subsubsection cr_cln_example Example
John Johnson's complete data structure library could use JJ as a prefix, so classes would be:

\code
   class JjLinkList
   {
   }
\endcode 

\subsection cr_mn Method Names

    * Use the same rule as for class names. 

\subsubsection cr_mn_justif Justification

    * Of all the different naming strategies many people found this one the best compromise. 

\subsubsection cr_mn_exa Example
\code
   class NameOneTwo
   {
   public:
      int                   DoIt();
      void                  HandleError();
   }
\endcode 
\subsection cr_can  Class Attribute Names

    * Attribute names should be prepended with the character 'm'.
    * After the 'm' use the same rules as for class names.
    * 'm' always precedes other name modifiers like 'p' for pointer. 

\subsubsection cr_can_just Justification

    * Prepending 'm' prevents any conflict with method names. Often your methods and attribute names will be similar, especially for accessors. 

\subsubsection cr_can_exa Example
\code
   class NameOneTwo
   {
   public:
      int                   VarAbc();
      int                   ErrorNumber();
   private:
      int                   mVarAbc;
      int                   mErrorNumber;
      String*               mpName;
   }
\endcode 
\subsection cr_man Method Argument Names

    * The first character should be lower case.
    * All word beginnings after the first letter should be upper case as with class names. 

\subsubsection cr_man_justif Justification

    * You can always tell which variables are passed in variables.
    * You can use names similar to class names without conflicting with class names. 

\subsubsection cr_man_exa Example
\code 
   class NameOneTwo
   {
   public:
      int                   StartYourEngines(
                               Engine& rSomeEngine, 
                               Engine& rAnotherEngine,
			       double vAmountOfFuel,
			       int * pIHaveNoIdeaWhatAPointerCouldDoHereButIHaveToShowAnExampleEvenVeryBadOfATooLongNameAndAPointerSoHereAreBoth);
   }
\endcode 
(
Be careful here: we have an r because these are references. 
That's why I added vAmountOfFuel, and also the, really too long, pointer variable
GG 13 sept 09.
)

\subsection cr_vns Variable Names on the Stack

    * use all lower case letters
    * use '_' as the word separator. 

\subsubsection cr_vns_justif Justification

    * With this approach the scope of the variable is clear in the code.
    * Now all variables look different and are identifiable in the code. 

\subsubsection cr_vns_exa Example
\code 
   int
   NameOneTwo::HandleError(int errorNumber)
   {
      int            error= OsErr();
      Time           time_of_error;
      ErrorProcessor error_processor;
      Time*          p_out_of_time= 0;
   }
\endcode 

The standard pointer notation is not entirely satisfactory because it doesn't look quite right, but it is consistent.

How do you handle statics? There's never a reason to have a static local to a function so there's no reason to invent a syntax for it. But like for most absolute rules, there is an exception, that is when making singletons. Use a "s_" prefix in this case. Take a look at Singleton Pattern for more details.

\subsection cr_pvar Pointer Variables

    * pointers should be prepended by a 'p' in most cases
    * place the * close to the pointer type not the variable name 

\subsubsection cr_pvar_justif Justification

    * The idea is that the difference between a pointer, object, and a reference to an object is important for understanding the code, especially in C++ where -> can be overloaded, and casting and copy semantics are important.
    * Pointers really are a change of type so the * belongs near the type. One reservation with this policy relates to declaring multiple variables with the same type on the same line. In C++ the pointer modifier only applies to the closest variable, not all of them, which can be very confusing, especially for newbies. You want to have one declaration per line anyway so you can document each variable. 

\subsubsection cr_pvar_exa Example
\code 
  String* pName= new String;

  String* pName, name, address; // note, only pName is a pointer.
\endcode 
\subsection cr_rvfrr Reference Variables and Functions Returning References

    * References should be prepended with 'r'. 

\subsubsection cr_rvfrr_justif Justification

    * The difference between variable types is clarified.
    * It establishes the difference between a method returning a modifiable object and the same method name returning a non-modifiable object. 

\subsubsection cr_rvfrr_exa Example
\code
   class Test
   {
   public:
      void               DoSomething(StatusInfo& rStatus);

      StatusInfo&        rStatus();
      const StatusInfo&  Status() const;

   private:
      StatusInfo&        mrStatus;
   }
\endcode 
\subsection cr_gv Global Variables

    * Global variables should be prepended with a 'g'. 

\subsubsection cr_gv_justif Justification

    * It's important to know the scope of a variable. 

\subsubsection cr_gv_exa Example
\code
    Logger  gLog;
    Logger* gpLog;
\endcode 
\subsection cr_gc Global Constants

    * Global constants should be all caps with '_' separators. 

\subsubsection cr_gc_justif Justification
It's tradition for global constants to named this way. You must be careful to not conflict with other global \#defines and enum labels.
\subsubsection cr_gc_exa Example
\code 
    const int A_GLOBAL_CONSTANT= 5;
\endcode

\subsection cr_sv Static Variables

    * Static variables may be prepended with 's'. 

\subsubsection cr_sv_justi Justification

    * It's important to know the scope of a variable. 

\subsubsection cr_sv_exa Example
\code
   class Test
   {
   public:
   private:
      static StatusInfo msStatus;
   }
\endcode
\subsection cr_tn Type Names

    * When possible for types based on native types make a typedef.
    * Typedef names should use the same naming policy as for a class with the word Type appended. 

\subsubsection cr_tn_ju Justification

    * Of all the different naming strategies many people found this one the best compromise.
    * Types are things so should use upper case letters. Type is appended to make it clear this is not a class. 

\subsubsection cr_tn_exa Example
\code
   typedef uint16  ModuleType;
   typedef uint32  SystemType;
\endcode 
\subsection cr_enumn Enum Names
Labels All Upper Case with '_' Word Separators
This is the standard rule for enum labels.
\subsubsection cr_enumn_exa Example
\code 
   enum PinStateType
   {
      PIN_OFF,
      PIN_ON
   };
\endcode

\subsection cr_enums Enums as Constants without Class Scoping
Sometimes people use enums as constants. When an enum is not embedded in a class make sure you use some sort of differentiating name before the label so as to prevent name clashes.
\subsection cr_enums_exa Example
\code
   enum PinStateType            If PIN was not prepended a conflict 
   {                            would occur as OFF and ON are probably
      PIN_OFF,                  already defined.
      PIN_ON
   };
\endcode
\subsection cr_enumsc Enums with Class Scoping
Just name the enum items what you wish and always qualify with the class name: Aclass::PIN_OFF.
Make a Label for an Error State
It's often useful to be able to say an enum is not in any of its valid states. Make a label for an uninitialized or error state. Make it the first label if possible.
\subsubsection cr_enumsc_exa Example
\code 
enum { STATE_ERR,  STATE_OPEN, STATE_RUNNING, STATE_DYING};
\endcode 

\subsection cr_def_macro Define and Macro Names

    * Put defines and macros in all upper using '_' separators. 

 \subsubsection cr_def_macro_exa Justification
This makes it very clear that the value is not alterable and in the case of macros, makes it clear that you are using a construct that requires care.

Some subtle errors can occur when macro names and enum labels use the same name.
 \subsubsection cr_def_macro_exa Example
\code
#define MAX(a,b) blah
#define IS_ERR(err) blah
\endcode

\subsection cr_cfna C Function Names (and Fortran, but for Fortran, you have to. GG)

    * In a C++ project there should be very few C functions.
    * For C functions use the GNU convention of all lower case letters with '_' as the word delimiter. 
   
  \subsubsection cr_cfna_just Justification

    * It makes C functions very different from any C++ related names. 

  \subsubsection cr_cfna_exa Example
 \code
   int
   some_bloody_function()
   {
   }
 \endcode

 *
 *
 *
 *
 *
 *
 *
 *
 *
 * \section coding_rules_documentation Documentation rules
 *
 * It is also interesting to have some documentation rules.
 * The CppCodingStandard (\ref coding_rules_for_transpp) has a lot of advices
 * (sometimes too strong!). The most important point are
 * - Documentation of all functions, classes, namespaces, members...
 *   	It is very important. Unfortunately, it is difficult to have
 *   	an updated view of this documentation if we say that it is called
 *   	by ... to do ... Indeed, the aim of these functions is to be
 *   	reused!
 * - To avoid the problem of the too important documentation of the functions,
 *      the good idea is to update the \ref architecture page: the idea is to 
 *      explain what is the aim of the class in the scope of the project,
 *      and to give an example of the call through the timeline.
 *      This timeline is a MAJOR feature of this documentation.
 * - Timeline : every step of the code has its functionnalities, and solve
 *   	a given issue. Unfortunately, the major idea, the major issue, can
 *   	be difficult to discover when we re-read the code. The timeline
 *   	is here to make this major documentation, and the explain the
 *   	scope of the function in a global point of view.(\ref architecture_timeline)
 * - Modules : The possibility to organise the code into library is shown in
 *   	the documentation through the modules. To create a new module, you have to use the defgroup - ingroup commands.
 * - Class-related page : some classes are virtual, and need to be derived!
 *   	WHY? Because these classes are  templates: a general view that need
 *   	to be precised. The reason why we use templates is that many kind
 *   	of specialization are possible, while not compatible!
 *   	For example, the Planet class is derived into Venus or Mars. 
 *   	Venus and Mars have different sizes, distances to the sun... That
 *   	are defined in the herited object. Therefore, if we want to work on
 *   	Titan, we have to create that planet, inherited from Planet!
 *   	But how to do that??? Mainly by copying the other planets, but it 
 *   	can be tricky sometimes (we don't want tricky things, but...).
 *   	Therefore, we create, thanks to the page function, the documentation
 *   	for the heritage (and more).
 *
 *
 *
 */


/**
 * \page versioning_system Versioning System
 * \section vers_sys_why Why using a version control system?
 * Versioning is very important!!!
 *
 * Explanation : (http://www.softpanorama.org/Lang/Cpp_rama/humor.shtml)
 * 
 * how to shot yourself in the foot in unix :
 *
 * Unix:
 *\code
 * % ls
 * foot.c foot.h foot.o toe.c toe.o
 * % rm * .o
 * rm:.o no such file or directory
 * % ls
 * %
 * \endcode
 *
 *
 *
 */

/**
 * \page cpp_curses Learn the C++
 *
 * C++ is not just another language, it is another paradigm!
 * Fortran and C just use functions, C++ use classes.
 * The basic idea is that some calculs can be very close.
 * For example, computing the orbit of two planets is basically
 * the same calcul, but with different initializations.
 *
 * In C or Fortran, you just do a function to compute the orbit
 * and you do a switch (a test) to see if it is Mars or Venus.
 * Unfortunately, you will have to add these tests all along your 
 * code, if you need, for example, the radius of the planet, its 
 * gravity.
 *
 * In C++, you just create a Planet class, which gives you the 
 * different functions, (needing some parameters, for example 
 * private Radius, Position, Orbit initalization). And you
 * inherit this class to have a Mars or Venus class (basically
 * these classes initializes the private parameters).
 * You just have to initialize your planet variable with Mars
 * or Venus to run the whole process, without need for tests, 
 * switches...
 * C++ is not easier, creating good classes is hard, and
 * need a long reflexion: but it is really MORE FLEXIBLE.
 * -> to add a new planet, you just have to create another 
 *  heritage, when it is complete, your code works! In Fortran
 *  you have to update all the switches: generally, you shot 
 *  yourself in the foot!!!
 *
 *  But don't forget:
 *  "In C++ it's harder to shoot yourself in the foot, but when you do, you blow off your whole leg."    — Bjarne Stroustrup.
 *
 *
 *
 *
 * See http://cpp.developpez.com/cours/cpp/ for a french curse.
 *
 * Or, in english, http://www.mindview.net/Books/DownloadSites
 * thinking in C++
 *
 * Then you will be able to understand:
 * "C++ : Where friends have access to your private members."    — Gavin Russell Baker. 
 *
 *
 * Some tips can be seen in the section \ref coding_style . This part 
 * is very important for the integration of your work.
 *
 * And don't forget http://www.softpanorama.org/Lang/Cpp_rama/humor.shtml
 * http://www.pbm.com/~lindahl/real.programmers.html
 */

/**
 * \page debugging
 *
 * \section how_to_debug How to debug this program correctly?
 * Firstly the code should be configured with the option --enable-debug
 * this option set up the debug option for the main code and the libraries.
 *
 * Since the code uses the autotools to create libraries, and then libtool
 * to  make the execution of the code (when it is not installed), the
 * link is not perfectly finished (the real installation finish the link,
 * please refer to libtool documentation to understand).
 *
 * This particularity prevent us to directly launch gdb or ddd on the 
 * executable.
 *
 * To do so, you have to call libtool!
 * It is better to give an example: to debug the mstest with ddd, you
 * simply call
 * <code> libtool --mode=execute ddd mstest </code> 
 * To be able to debug the library, you have to load it. (at this point,
 * you can only see the header!).
 * Two solutions: the first one is to load the library (I forgot the 
 * command) inside gdb or ddd. The second is to run the code one time,
 * that load the library, and then you can load the sources files of
 * that library, and debugg it (as you can see, I'm partial to that 
 * choice).
 *
 * For more example please see the curse of Emmanuel Fleury.
 *
 *
 *
 */



#endif
