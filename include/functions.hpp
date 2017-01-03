#ifndef FUNCTIONS_H
#define FUNCTIONS_H
/*
 assorted helper functions, so far only the reading of the input paramters is here
*/


#include "active_medium.hpp"
#include "beam_model.hpp"
#include <stdexcept>

//struct to return all the user input variables from parse_cmd_parameters to main
struct variable_package{
    double dose;
    double voltage;
    double pulseduration;
    int numberOfBins;
    double binSize;
    double rateLimit;
    double speedLimit;
    double electrodeArea;
    double Q;
    double E0;
    streambuf* logBuffer;
    string summaryFileName;
    int summaryIndepVar;
    string mediumName;
    active_medium* fillgas;
    beam_model* this_beam;
    string settingsLogFileName;
};

//exceptions calls to signal an abort
class aborted: public std::exception
{};

//return the parsed command line arguments as a struct of the variable_package form
variable_package parse_cmd_parameters(int, char**, ofstream &pLogStream, string);

#endif // FUNCTIONS_H
