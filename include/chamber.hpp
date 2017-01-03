#ifndef CHAMBER
#define CHAMBER
/*
 Takes care of loading the chamber file and returning the chamber parameters pertinent to the calculation
*/


#include <fstream>

const static double pi = 3.14159265358979323846;

class chamber
{
	private:
        double plateDistance;
        double diameter;
        double kQ;
        double Nw;
	public:
		chamber(std::ifstream& );
        double calc_bin_size(int numberOfBins) {return (plateDistance*1.0e6/numberOfBins);} //in nm
        double calc_E0(double voltage) { return (voltage*10./plateDistance);} //in [Voltage]/cm
        double calc_charge_liberated(double dose) {return (1.0e9*dose/(Nw*kQ));} // in pC for a dose given in mGy
        double get_area() { return (pi/4.*diameter*diameter);} //in mm^2

};

#endif
