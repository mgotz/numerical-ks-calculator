#ifndef SETTINGS
#define SETTINGS

//use to define some global preprocessor switches

//turn off the recalculation of the electric field, i.e. deactivate shield by charges
//#define CONSTFIELD

// logging levels are 0, 1 and 2
// at 0 only the inital settings and the results are outputted
// at 1 a status every 1000 steps is added
// at 2 the calculation of the stepsize is logged at every iteration
// in order to allow a definition of DEBUGLEVEL at the compilation command line this is only defined if it has not happened already
#ifndef DEBUGLEVEL
    #define DEBUGLEVEL 2
#endif

//output files containing E-Field, electron and ion distributions, every 1000 steps
//#define CHARGETIMELINE

//undefine saves computations compared to just setting electron recombination coefficient=0
//#define ELECTRONCOMBO


// define if using constant ion mobilities, constant ion recombination rate or constant electron recombination rate (no E-field dependence)
// significantly speeds up the calculation
#define CONSTIONMOB
#define CONSTRECOMB
#define CONSTELECTRONCOMBO



#endif // SETTINGS
