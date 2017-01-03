# IC_Simulation
numerical solution of charge transport, creatiion and recombination inside a 1D ionization chamber model

Building:
run cmake from a directory you want to build in and give a path to this directory as argument
cmake will create a makefile, which can be used to build by simply runn make

Requirements: standard library and program_options from the boost library
Tested on gcc 4.9.2


Alternatively compile manually. Source files are in src, headers are in include directories.
Link all files and don't forget -lboost_program_options

Usage:
The program's parameters can be given either on the command line or in a configuration file (or any combination of the two)
Run IK_Simluation -h to see a list of the parameters.

Beyond the straightfoward settings, the program needs a chamber file and an active medium file.
The chamber file gives the dimensions and properties (mainly calibration factor to calculate liberated charge from a given dose) 
of the chamber.
See AdvMarkus_kammer.txt for an example.

The active_medium file specifies the properties of the filling gas.
Those can be simple constants, function definitions or point to additional files, if they are given in a spline format
See active_medium_air_with_direct.txt for an example.