#ifndef CHAMBER
#define CHAMBER
/*
 * Takes care of loading the chamber file and returning the chamber parameters pertinent to the calculation.
 * One abstract base class: chamber, provides the interface to call the methods relevant for the calculation,
 * while two implementations for cylinder and plane-parallel provide the specific implementations of those methods.
*/


#include <fstream>

const static double pi = 3.14159265358979323846;

class chamber
{
protected:
    double binLength; //nm
    double* binVolume; //mm^3
    double* boundaryArea; //mm^2
    double* EField; //V/cm
    double* r; //nm
    double E0;
    double kQ;
    double Nw;
    double U; //V
    int type; //1 for plane-parallel, 2 for clyinder, maybe more later
    int polarity;
    size_t margin;
    size_t numberOfBins;
public:
    double volume; //mm^3
    virtual ~chamber() = 0;
    virtual double return_binLength() const {return (binLength);} //in nm
    virtual double* return_radius() const {return (r);} // in nm
    virtual double* return_binVolume() const {return (binVolume);} //in mm^3
    virtual int get_type() const {return type;}
    virtual double calc_liberated_charge(double pDose) const {return (1.0e9*pDose/(Nw*kQ));} //in pC for dose in mGy
    virtual double return_average_E() const {return (E0);} // in V/cm
    virtual double* calc_E(const double * const  rhoPosIon, const double * const rhoNegIon, const double * const  rhoE, const double permittivity) = 0;
};

inline chamber::~chamber(){}

class cylinder : public chamber{
private:
    double r1;
    double r2;
    double length;
    double QOffset;
    double UFromCharges;
public:
    cylinder(double pR1, double pR2, double pLength, double pNw, double pKq, double pU, int pPolarity,size_t pNumberOfBins, size_t pMargin);
    virtual double* calc_E(const double * const  rhoPosIon, const double * const  rhoNegIon, const double * const  rhoE, const double permittivity) ;
    ~cylinder();

};

class plane_parallel : public chamber{
private:
    double plateDistance;
    double diameter;
    double UFromCharges;
    double EOffset;
public:
    plane_parallel(double pPlateDistance, double pDiameter, double pNw, double pKq, double pU, int pPolarity, size_t pNumberOfBins, size_t pMargin);
    virtual double* calc_E(const double * const  rhoPosIon, const double * const  rhoNegIon, const double * const  rhoE,const double permittivity);
    ~plane_parallel();
};

chamber* read_chamber_file(std::ifstream& chamberFile,double voltage, int polarity, size_t pNumberOfBins,size_t pMargin);

#endif
