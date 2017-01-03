#include "functions.hpp"
#include "chamber.hpp"
#include "beam_model.hpp"

#include <boost/program_options.hpp>
#include <boost/property_tree/exceptions.hpp>

#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <sstream>

using namespace std;
namespace po = boost::program_options;

variable_package parse_cmd_parameters(int argc,char **argv,ofstream &pLogFile, string startTime)
{
    string config_file;
    variable_package user_input;

    //declare cmd only options
    po::options_description cmdOnly("command line options");
    cmdOnly.add_options()
            ("help,h","produce help message")
            ("config,c",po::value<string>(&config_file)->default_value("config.cfg"),
                "name of the configuration file")
            ;

    //declare options that can be in file as well as command line
    po::options_description config("configuration");
    config.add_options()
            ("dose,D",po::value<double>(),"applied pulse dose in mGy")
            ("charge,Q",po::value<double>(),"liberated charge in nC (alternative to and overwrites dose value)")
            ("number-of-pulses,N",po::value<int>(),"number of pulses in the pulse train")
            ("micro-pulsed,P",po::bool_switch()->default_value(false),"set switch to use a micro pulse structure like on ELBE")
            ("voltage,U",po::value<double>(),"chamber voltage in V")
            ("chamber-file",po::value<string>(),"file giving chamber specifications")
            ("active-medium",po::value<string>(),"file specifying the active medium")
            ("rate-limit,r",po::value<double>(),"maximum reaction in one timestep")
            ("speed-limit,v",po::value<double>(),"maximum movement in bins per timestep")
            ("number-of-bins,n",po::value<int>(),"number of bins for spatial discretization")
            ("log-file",po::value<string>(),"file for non critical output, defaults to cou")
            ("summary-file",po::value<string>(),"file to write a condensed result to, optional")
            ("summary-variable",po::value<char>(),"independent variable for the summary output, optional")
            ("settings-log",po::value<string>(),"file to write the active medium settings to, empty string to write no log")
            ;

    //declare positional alternative, to ensure backwards compatibility
    po::positional_options_description posOptions;
    posOptions.add("chamber-file",1);
    posOptions.add("active-medium",1);
    posOptions.add("voltage",1);
    posOptions.add("dose",1);
    posOptions.add("number-of-pulses",1);
    posOptions.add("rate-limit",1);
    posOptions.add("number-of-bins",1);
    posOptions.add("summary-file",1);
    posOptions.add("summary-variable",1);
    posOptions.add("settings-log",1);
    posOptions.add("speed-limit",1);

    po::options_description cmdOptions;
    cmdOptions.add(cmdOnly).add(config);

    po::options_description configFileOptions;
    configFileOptions.add(config);

    po::variables_map vm;
    store(po::command_line_parser(argc,argv).options(cmdOptions).positional(posOptions).run(),vm);
    notify(vm);

    int numberOfPulses = 0; //local variable

    //read and parse config file
    if (vm.count("help"))
    {
        cout<<cmdOptions<<endl;
        throw aborted();
    }

    ifstream ifs(config_file.c_str());
    if (!ifs)
    {
        cout << "can not open config file: " << config_file << "\n";
    }
    else
    {
        store(parse_config_file(ifs, configFileOptions), vm);
        notify(vm);
    }

    //check all the options and put appropriate values in user_input struct

    if (vm.count("charge"))
    {
        user_input.Q=vm["charge"].as<double>()*1000.;
        user_input.dose=0.0;
    }
    else if (vm.count("dose"))
    {
        user_input.dose=vm["dose"].as<double>();
    }
    else
    {
        throw invalid_argument("neither dose nor charge specified");
    }

    if (vm.count("number-of-pulses"))
    {
        numberOfPulses = vm["number-of-pulses"].as<int>();
        user_input.pulseduration = (numberOfPulses-1)*77.0+5.0e-3; //in ns

    }
    else
    {
        throw invalid_argument("number of pulses not specified");
    }

    if (vm.count("voltage"))
    {
        user_input.voltage=vm["voltage"].as<double>();
    }
    else
    {
        throw invalid_argument("voltage not specified");
    }

    if (vm.count("rate-limit"))
    {
        user_input.rateLimit=vm["rate-limit"].as<double>();
    }
    else
    {
        user_input.rateLimit=0.1; //default value
    }

    if (vm.count("speed-limit"))
    {
        user_input.speedLimit=vm["speed-limit"].as<double>();
        if (user_input.speedLimit >= 1.0)
        {
            throw invalid_argument("specified speedlimit too large, must be < 1");
        }
    }
    else
    {
        user_input.speedLimit=0.5; //default value
    }

    if(vm.count("number-of-bins"))
    {
        user_input.numberOfBins=vm["number-of-bins"].as<int>();
    }
    else
    {
        user_input.numberOfBins=1000; //default value
    }

    if (vm.count("chamber-file"))
    {
        ifstream chamberFile((vm["chamber-file"].as<string>()).c_str());
        if (!chamberFile)
            {
            throw invalid_argument("failed to open chamber file");
            }
        chamber thisChamber(chamberFile);
        user_input.electrodeArea=thisChamber.get_area(); //in mm2
        user_input.E0=thisChamber.calc_E0(user_input.voltage); //in V/cm
        user_input.binSize=thisChamber.calc_bin_size(user_input.numberOfBins);//in nm
        //if Q is give, dose is set 0 and Q takes precedence
        if (user_input.dose > 0)
        {
            user_input.Q=thisChamber.calc_charge_liberated(user_input.dose);//total liberated charge in pC
        }
    }
    else
    {
        throw invalid_argument("chamber file not specified");
    }

    if (vm["micro-pulsed"].as<bool>())
    {
        user_input.this_beam = new elbe_micropulses(user_input.Q/(double)(numberOfPulses)/(double)(user_input.numberOfBins),user_input.pulseduration);
    }
    else
    {
        user_input.this_beam = new continuous_pulse(user_input.Q/user_input.pulseduration/(double)(user_input.numberOfBins),user_input.pulseduration);
    }

    string errorString;
    if (vm.count("active-medium"))
    {
        //The properties of the filling-gas are dealt with in the active_medium_class
        //it is created by reading an ini-style file
        user_input.mediumName=vm["active-medium"].as<string>();
        try{
            user_input.fillgas = new active_medium(user_input.mediumName);
        }
        catch (boost::property_tree::ptree_bad_path e)
        {
            errorString = "Missing a required parameter in active medium file: ";
            errorString.append(e.what());
            throw invalid_argument(errorString.c_str());
        }
        catch (boost::property_tree::ptree_bad_data e)
        {
            errorString = "Malformed data in active medium file: ";
            errorString.append(e.what());
            throw invalid_argument(errorString.c_str());
        }

        catch (ifstream::failure e){
            errorString = "Error reading or parsing file for active medium:\n";
            errorString.append(e.what());
            throw invalid_argument(errorString.c_str());
        }
        catch (unknown_arg_error e){
            errorString = "can not interpret parameter in active mediume file: \n";
            errorString.append(e.what());
            throw invalid_argument(errorString.c_str());
            }
        catch (file_open_error e){
            errorString = "failed to oper file: ";
            errorString.append(e.what());
            throw invalid_argument(errorString.c_str());
            }
        //catch (...){
        //    errorString = "unknown exception during active medium construction";
        //    throw invalid_argument(errorString.c_str());
        //    }
    }
    else
    {
        throw invalid_argument("active medium file not specified");
    }

    if (vm.count("log-file"))
    {
        pLogFile.open(vm["log-file"].as<string>().c_str(),ios::out|ios::app);
        if (!pLogFile)
        {
            cerr<<"failed to open log file: "<<vm["log-file"].as<string>()<<endl;
            user_input.logBuffer=cout.rdbuf();
        }
        else
        {
            user_input.logBuffer=pLogFile.rdbuf();
        }

    }
    else
    {
        user_input.logBuffer=cout.rdbuf();
    }

    if (vm.count("summary-file") || vm.count("summary-variable"))
    {
        if (!(vm.count("summary-file")))
        {
            user_input.summaryIndepVar=0;
            cerr<<"Inconsistent options: a variable for summary was given but no file"<<endl;
        }
        else if (!(vm.count("summary-variable")))
        {
            user_input.summaryIndepVar=0;
            cerr<<"Inconsistent options: a summary file was given but no independent variable"<<endl;
        }
        else
        {
            user_input.summaryIndepVar=vm["summary-variable"].as<char>();
            user_input.summaryFileName=vm["summary-file"].as<string>();
        }
    }

    if (vm.count("settings-log"))
    {
        user_input.settingsLogFileName=vm["settings-log"].as<string>();
    }
    else
    {
        ostringstream fileName;
        fileName << "settings_"
                 << startTime
                 <<".log";
        user_input.settingsLogFileName=fileName.str();
    }

    return user_input;
}

