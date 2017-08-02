#ifndef FUNCTIONS_H
#define FUNCTIONS_H
/*
 * helper functions to parse user input and log those parameters
 * also defines structs used to pass user input to and from the helper functions and
 * return the results from the loop_template
*/


#include "active_medium.hpp"
#include "beam_model.hpp"
#include "chamber.hpp"

#include <stdexcept>

//struct to return all the user input variables from parse_cmd_parameters to main
struct variable_package{
    double dose;
    double voltage;
    int polarity;
    double pulseduration;
    unsigned int numberOfBins;
    unsigned int margin;
//    double binSize;
    double rateLimit;
    double speedLimit;
//    double electrodeArea;
    double Q;
//    double E0;
    streambuf* logBuffer;
    string summaryFileName;
    int summaryIndepVar;
    string mediumName;
    active_medium* fillgas;
    beam_model* this_beam;
    chamber* this_chamber;
    string settingsLogFileName;
};

//exceptions calls to signal an abort
class aborted: public std::exception
{};

//return the parsed command line arguments as a struct of the variable_package form
variable_package parse_cmd_parameters(int, char**, ofstream &pLogStream, string);

//write the parameters to the log
void output_parameters(variable_package &userInput, ostream &log, string startTime);

//struct to return the calculation results from the loop to the main
struct loop_return{
    double collectedCharges[3];
    double remainingCharges[3];
    double liberatedCharges;
    double ionRecombinedCharges;
    double directlyRecombinedCharges;
    double EmaxAll;
    double EminAll;
    double elapsedTime;
    double electronCollectionTime;
    double attachedCharge;
    double beamOnTime;
    bool positivityWarning;
    size_t steps;

    //define constructor do make initialization simple
    loop_return(){
        for  (size_t i=0;i<3;i++)
        {
            collectedCharges[i] = 0;
            remainingCharges[i] = numeric_limits<double>::max();
        }

        liberatedCharges = 0;
        attachedCharge = 0;
        directlyRecombinedCharges=0;
        ionRecombinedCharges=0;

        EmaxAll = 0;
        EminAll = numeric_limits<double>::max();

        elapsedTime=0;
        electronCollectionTime=-1;
        beamOnTime = -1;

        positivityWarning=false;
        steps = 0;
    }
};
#endif // FUNCTIONS_H
