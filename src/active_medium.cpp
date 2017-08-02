#include "active_medium.hpp"
#include "active_medium_functiods.hpp"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>

#include <vector>
#include <stdexcept>
#include <string>

using namespace std;
namespace pt = boost::property_tree;

//the number of required and optional sections and thus functiods in the active medium file
const size_t numberOfRequired = 6;
const size_t numberOfOptional = 1;


//construct active medium properties from an INI-file with boost::property_tree
active_medium::active_medium(std::string pMediumFile)
{
    string tempString;

    //vector detailing the sections of the INI-File, this allows to simply iterate the vector
    //instead of repeating code for each section, all sections here are required
    //when adding more section: increase vector size, add alements and also define appropriate methods
    vector<string> names(numberOfRequired);
    names[0] = "volume_recombination";
    names[1] = "positive_ion_mobility";
    names[2] = "negative_ion_mobility";
    names[3] = "electron_attachment_rate";
    names[4] = "electron_mobility";
    names[5] = "direct_recombination";
    //only iterating up to 6 currently, additional ones could be defined here

    //define additional, optional parameters (so far only simple constants)
    //their default values are defined in a the second vector (super inelegant, but everything else is so complicated)
    vector<string> optionalNames(numberOfOptional);
    optionalNames[0] = "relative_permittivity";

    vector<double> optionalValues(numberOfOptional);
    optionalValues[0] = 1.0;




    pt::ptree properties;
    pt::read_ini(pMediumFile,properties);

    string modelType;

    for (size_t i=0; i < names.size();i++){
        // read the type of model and construct the corresponding functiod
        modelType = properties.get<string>(names[i]+".function");
        if (modelType == "constant")
        {
            mFunctiods.push_back(new constant (properties.get<double>(names[i]+".param0")));
        }
        else if (modelType == "line")
        {
            mFunctiods.push_back(new linear (properties.get<double>(names[i]+".param0"),
                                             (properties.get<double>(names[i]+".param1")-properties.get<double>(names[i]+".param0"))/4000.));
        }
        else if (modelType == "exp")
        {
            mFunctiods.push_back(new exponential (properties.get<double>(names[i]+".param0"),
                                                  properties.get<double>(names[i]+".param1"),
                                                  properties.get<double>(names[i]+".param2")));
        }
        else if (modelType == "spline")
        {
            tempString = properties.get<string>(names[i]+".file_name");
            ifstream parameterFile(tempString.c_str() );
            try{
                mFunctiods.push_back(new spline_eval_only (parameterFile));
                }
            catch(ifstream::failure e){
                tempString = "failure parsing parameter file for "+names[i]+": "+e.what();
                throw ifstream::failure(tempString);
                }
        }
        else if (modelType == "poly4")
        {
            mFunctiods.push_back(new pol4 (properties.get<double>(names[i]+".param0"),
                                           properties.get<double>(names[i]+".param1"),
                                           properties.get<double>(names[i]+".param2"),
                                           properties.get<double>(names[i]+".param3"),
                                           properties.get<double>(names[i]+".param4")));

        }
        else if (modelType == "poly5")
        {
            mFunctiods.push_back(new pol5 (properties.get<double>(names[i]+".param0"),
                                           properties.get<double>(names[i]+".param1"),
                                           properties.get<double>(names[i]+".param2"),
                                           properties.get<double>(names[i]+".param3"),
                                           properties.get<double>(names[i]+".param4"),
                                           properties.get<double>(names[i]+".param5")));

        }
        else
        {
            tempString = "Model type: "+modelType +", for "+names[i];
            throw unknown_arg_error(tempString);
        }

        //check for the smalles range of the functions
        mMaxE = std::numeric_limits<double>::max();
        for (size_t i=0; i < mFunctiods.size();i++){
            if (mFunctiods[i]->mMaxE < mMaxE){
                mMaxE = mFunctiods[i]->mMaxE;
            }
        }

    }
    //iterate over the optional values
    for (size_t i=0; i < optionalNames.size();i++){
        mFunctiods.push_back(new constant (properties.get(optionalNames[i]+".param0",optionalValues[i])));
    }


}


//define functions for the different properties, by pointing to the created array pointers
double active_medium::volume_recombination(double pEField){
    return (mFunctiods[0]->evaluate(pEField));
}
double active_medium::mobility_pos(double pEField){
    return (mFunctiods[1]->evaluate(pEField));
}

double active_medium::mobility_neg(double pEField){
    return (mFunctiods[2]->evaluate(pEField));
}

double active_medium::attachment(double pEField){
    return (mFunctiods[3]->evaluate(pEField));
}

double active_medium::mobility_e(double pEField){
    return (mFunctiods[4]->evaluate(pEField));
}

double active_medium::direct_recombination(double pEField){
    return (mFunctiods[5]->evaluate(pEField));
}

//the optional ones
double active_medium::relative_permittivity(double pEField){
    return (mFunctiods[numberOfRequired+0]->evaluate(pEField));
}


void active_medium::write_log(ofstream &pLogFile){
    pLogFile << "E in V/cm\t"
             << "µ_I+ in m^2/Vs\t"
             << "µ_I- in m^2/Vs\t"
             << "µ_e in m^2/Vs\t"
             << "electron attachment in 1/ns\t"
             << "Ion-ion-recombination in m^3/s\t"
             << "Ion-electron-recombination in m^3/s\t"
             << endl;
    for(size_t i=0; i < 8001; i++)
    {
        pLogFile << i << "\t"
                 << this->mobility_pos(i) << "\t"
                 << this->mobility_neg(i) << "\t"
                 << this->mobility_e(i) << "\t"
                 << this->attachment(i) << "\t"
                 << this->volume_recombination(i) << "\t"
                 << this->direct_recombination(i) << "\t"
                 << endl;
     }
    pLogFile << endl;

}


