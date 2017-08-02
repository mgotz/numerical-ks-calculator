/*
Copyright (c) 2016 Malte Gotz

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


//standard headers
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

//my headers
#include "active_medium.hpp"
#include "chamber.hpp"
#include "functions.hpp"
#include "config.hpp"
#include "loop_template.hpp"

using namespace std;

int main(int argc, char **argv){
// charges are calculated considering the entire chamber cross-section.
// to obtain a charge density the values must be divided by the area

variable_package userInput;
ofstream logStream;

//define the start time

time_t now = time(0);
struct tm* tm_now;
tm_now = localtime(&now);
char startTime[20];
strftime(startTime, 20, "%Y-%m-%d_%H-%M-%S",tm_now);

//easier formatting of the time, that needs a c++11 function not implemented on the gcc I used
//ostringstream startTimeStream;
//startTimeStream<<std::put_time(localtime(&now),"%Y-%m-%d_%H-%M-%S");
//string startTime = startTimeStream.str();


//get the user input
try
    {
    userInput = parse_cmd_parameters(argc, argv, logStream,startTime);
    }
catch (aborted)
    {
        return 0;
    }
catch (invalid_argument e)
    {
        cerr<<"Error parsing the input parameters"<<endl;
        cerr<<e.what()<<endl;
        return 1;
    }

//build log stream, it recieves extra logging from during the simulation, while starting info and results are written to cout (can be used to declutter the output if logging is not disabled)
ostream log(userInput.logBuffer);

//show the parameters in use
output_parameters(userInput, log, startTime);


loop_return * results;



/*
 * different versions of the loop template are constructed for the different combinations of chamber and polarity
 * this essentially creates duplicate code: one loop version for each call to the template, but has the advantage
 * that branching inside the inner loops of the calculation is avoided.
*/
if (userInput.this_chamber->get_type() == 1)
{
    results = time_loop<1,1>(userInput, log, startTime);
}
else if (userInput.this_chamber->get_type() == 2)
{
    if (userInput.polarity == 1)
    {
        results = time_loop<1,2>(userInput, log, startTime);
    }
    else
    {
        results =time_loop<-1,2>(userInput, log, startTime);
    }
}
else
{
    cerr<<"unkown chamber type: "<<userInput.this_chamber->get_type()<<endl;
    return 1;
}


if (!results){
    cerr<<"main calculation aborted"<<endl;
    return 1;
}

//final output
cout<<"***Results***"<<endl;
cout<<"remaining charge (I+/Q;I-/Q,e-/Q): "<<results->remainingCharges[0]/userInput.Q<<" "<<results->remainingCharges[1]/userInput.Q<<" "<<results->remainingCharges[2]/userInput.Q<<endl;

//if there are hardly any negative charges left, assume that all the remaining positive charges were collected without issue
if (results->remainingCharges[1]<userInput.Q*1.0e-06 && results->remainingCharges[2]<userInput.Q*1.0e-06)
    {
    results->collectedCharges[0] += results->remainingCharges[0];
    results->remainingCharges[0]=0;
    }

if (results->remainingCharges[0]>userInput.Q*1.0e-04 || results->remainingCharges[1]>userInput.Q*1.0e-04 || results->remainingCharges[2]>userInput.Q*1.0e-04)
    {
    cerr<<"too many charges at the end of the loop --> increase calculation time!! (increase reactionRateLimit or NumberOfSteps) "<<std::endl;
    }


if (fabs(results->liberatedCharges/userInput.Q- 1) > 0.001)
{
    cerr<<"actually liberated charge deviates from target values, check calculation";
    results->collectedCharges[0] = 0;
    results->collectedCharges[1] = 0;
    results->collectedCharges[2] = 0;
}

cout<<endl;
cout<<"liberated charge (Q): "<<results->liberatedCharges<<" pC, attached: "<<results->attachedCharge<<" pC"<<endl;
cout<<"collected charge (I+; I-; e-; I-+e-): "
    <<results->collectedCharges[0]<<" "
    <<results->collectedCharges[1]<<" "
    <<results->collectedCharges[2]<<" "
    <<results->collectedCharges[1]+results->collectedCharges[2]<<endl;
cout<<"recombined charge (I+&I-; I+&e-; Q-recombined): "
    <<results->ionRecombinedCharges<<" "
    <<results->directlyRecombinedCharges<<" "
    <<results->liberatedCharges-results->ionRecombinedCharges-results->directlyRecombinedCharges<<endl;
cout<<endl;

cout<<"collection times: "<<endl;
cout<<"total: "<<results->elapsedTime/1000.<<" µs; " <<results->steps<<" steps"<<endl;
cout<<"electrons: "<<results->electronCollectionTime/1000.<<" µs, Beam: "<<results->beamOnTime/1000.<<" µs"<<endl;


cout<<"fraction of free electrons: ";
if (results->collectedCharges[1]+results->collectedCharges[2]>0)
    cout<<results->collectedCharges[2]/results->liberatedCharges<<endl;
else
    cout<<0<<endl;


cout.precision(8);


cout<<"ksat (Q/collected +; Q/collected -): ";

if (results->collectedCharges[0]>0)
    cout<<results->liberatedCharges/(results->collectedCharges[0])<<"; ";
else
    cout<<0<<"; ";

if ((results->collectedCharges[1]+results->collectedCharges[2]) >0)
    cout<<results->liberatedCharges/(results->collectedCharges[1]+results->collectedCharges[2]);
else
    cout<<0;

cout<<endl;

if (results->positivityWarning)
{
    cerr << "Waring: Negative charge density was cut off, conservation of total charge not guaranteed" << endl;
    }


cout<<"Emax: "<<results->EmaxAll<<", Emin: "<<results->EminAll
   <<", total charge conservation (1?): "<<(results->ionRecombinedCharges+results->directlyRecombinedCharges)/(results->liberatedCharges-results->collectedCharges[0])<<endl;
cout<<endl;


//compact output (summary)
double indep_variable;
if (userInput.summaryIndepVar)
{
    ofstream compact_out(userInput.summaryFileName.c_str(), ios::out|ios::app);
    switch (userInput.summaryIndepVar)
    {
    case 'D':
        indep_variable = userInput.dose;
        log<<"result summary with dose printed to "<<userInput.summaryFileName<<endl;
        break;
    case 'Q':
        indep_variable = userInput.Q/1000.;
        log<<"result summary with total charge printed to "<<userInput.summaryFileName<<endl;
        break;
    case 't':
        indep_variable = userInput.pulseduration;
        log<<"result summary with pulse duration printed to "<<userInput.summaryFileName<<endl;
        break;
    case 'N':
        indep_variable = ceil(userInput.pulseduration/77.0);
        log<<"result summary with number of pulses printed to "<<userInput.summaryFileName<<endl;
        break;
    case 'V':
        indep_variable = userInput.voltage;
        log<<"result summary with voltage printed to "<<userInput.summaryFileName<<endl;
        break;
    }

    compact_out.precision(16);
    compact_out<<scientific;
    compact_out<<indep_variable<<", "<<results->liberatedCharges/(results->collectedCharges[0])<<endl;

}

tm_now = localtime(&now);
char endTime[20];
strftime(endTime, 20, "%Y-%m-%d_%H-%M-%S",tm_now);

log<<"done at "<<endTime<<endl;
log<<endl;

return 0;
}
