#ifndef ACTIVE_MEDIUM
#define ACTIVE_MEDIUM
/*
 * This class defines an active medium, with several methods to return its properties,
 * such as ion mobility or electron attachment, depending on the electric field strength.
 * It is constructed by reading an ini-style file with boost::property_tree
*/


#include "active_medium_functiods.hpp"

#include <fstream>
#include <limits>
#include <vector>
#include <string>
#include <stdexcept>


using namespace std;


class active_medium
{
	public:
        active_medium(std::string); //constructor takes file path to the active medium ini file

        //methods to get the different properties:
        double mobility_e (double);
        double mobility_pos (double);
        double mobility_neg (double);
        double attachment (double);
        double volume_recombination (double);
        double direct_recombination (double);

        double relative_permittivity(double);

        void write_log(ofstream&); //write table of function values in stream

        double mMaxE;
	private:
        //functiods for the different properties
        vector<evaluation_functiod*> mFunctiods;
};

class unknown_arg_error : public runtime_error
{
	public:
		unknown_arg_error(const string &message) : runtime_error(message) {}
};
class file_open_error : public runtime_error
{
	public:
		file_open_error(const string &message) : runtime_error(message) {}
};




#endif
