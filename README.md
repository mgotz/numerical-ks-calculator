# numerical-ks-calculator
Numerical solution of charge transport, creation and recombination inside a 1D ionization chamber model.

Building:
Run cmake from a directory you want to build in and give a path to this directory as argument.
Cmake will create a makefile, which can be used to build by simply running make.

Requirements: standard library and program_options from the boost library.
Tested on gcc 4.9.2

Alternatively compile manually. Source files are in src, headers are in include directories.
Link all files and don't forget -lboost_program_options

Usage:
The program's parameters can be given either on the command line or in a configuration file (or any combination of the two)
Run IK_Simluation -h to see a list of the parameters.

Beyond the straightfoward settings, the program needs a chamber file and an active medium file.
The chamber file gives the dimensions and properties (mainly calibration factor to calculate liberated charge from a given dose) 
of the chamber.
See AdvMarkus_chamber.txt for an example.
The chamber file must maintain this order of the properties and no blank lines maybe added.
The programm expects each property in a certain position. In each line the content up to the first whitespace is discarded,
then the data is read an the following discription is also discarded.

The active_medium file specifies the properties of the filling gas.
Those can be simple constants, function definitions or point to additional files, if they are given in a spline format
See active_medium_air_with_direct.txt for an example and a detailed discription of its format.
