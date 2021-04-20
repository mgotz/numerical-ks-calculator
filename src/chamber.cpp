#include "chamber.hpp"
#include "config.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

chamber* read_chamber_file(std::ifstream& chamberFile, double voltage, int polarity, size_t pNumberOfBins, size_t pMargin)
{
string temp;
string chamberType;
chamber* thisChamber;
double value0, value1, value2, value3, value4;
chamberFile>>temp>>chamberType;
if (chamberType == "cylindrical")
    {
    #if LOGLEVEL > 1
        cout<<"cylindrical type chamber"<<endl;
    #endif
    chamberFile>>temp>>value0>>temp;
    chamberFile>>temp>>value1>>temp;
    chamberFile>>temp>>value2>>temp;
    chamberFile>>temp>>value3>>temp;
    chamberFile>>temp>>value4>>temp;
    thisChamber = new cylinder(value0,value1, value2, value3, value4, voltage, polarity, pNumberOfBins, pMargin);
}
else if (chamberType == "plane-parallel")
{
    #if LOGLEVEL > 1
        cout<<"plane-parallel type chamber"<<endl;
    #endif
    chamberFile>>temp>>value0>>temp;
    chamberFile>>temp>>value1>>temp;
    chamberFile>>temp>>value2>>temp;
    chamberFile>>temp>>value3>>temp;
    thisChamber = new plane_parallel(value0, value1, value2, value3, voltage, polarity, pNumberOfBins, pMargin);
}
else
{
    #if LOGLEVEL > 1
        cout<<"defaulting to plane-parallel chamber"<<endl;
    #endif
    value0 = std::strtod(chamberType.c_str(),NULL);
    //value0 = std::stod(chamberType); //c++11
    chamberFile>>temp>>value1>>temp;
    chamberFile>>temp>>value2>>temp;
    chamberFile>>temp>>value3>>temp;
    thisChamber = new plane_parallel(value0, value1, value2, value3, voltage, polarity, pNumberOfBins, pMargin);
}

return (thisChamber);

}

cylinder::cylinder(double pR1, double pR2, double pLength, double pNw, double pKq, double pU, int pPolarity, size_t pNumberOfBins, size_t pMargin)
{
    type = 2;
    r1 = pR1;
    r2 = pR2;
    length = pLength;
    volume = length*(r2*r2-r1*r1)*pi;
    Nw = pNw;
    kQ = pKq;
    numberOfBins = pNumberOfBins;
    margin = pMargin;
    U = pU;
    polarity = pPolarity;

    E0 = 10.*U/(r2-r1);

    //create the bin arrays
    binVolume = new double[numberOfBins+2*margin];
    boundaryArea = new double[numberOfBins+2*margin+1];

    //r are the radii of the cell boundaries
    r = new double[numberOfBins+2*margin+1];

    //uniform bins in r
    binLength = (r2-r1)/numberOfBins;
    for (size_t i = 0; i<(numberOfBins+2*margin+1); i++)
    {
        r[i] = binLength*int(i-margin)+r1; // in mm
    }

    binLength = binLength *1.0e6; //convert to nm
    //calculate the size of the area at the boundary
    for (size_t i = 0; i<(numberOfBins+2*margin+1); i++)
    {
        boundaryArea[i] = r[i]*2.0*pi*length; //in mm^2
    }
    //calculate the bin volume
    for (size_t i = 0; i<(numberOfBins+2*margin); i++)
    {
        binVolume[i] = pi*length*(r[i+1]*r[i+1]-r[i]*r[i]); //in mm^3
    }

    EField = new double[numberOfBins+2*margin+1]();

}

cylinder::~cylinder()
{
    delete[] binVolume;
    delete[] boundaryArea;
    delete[] EField;
}

double* cylinder::calc_E(const double * const rhoPosIon,const double * const rhoNegIon,const double * const rhoE, double const permittivity)
{
    UFromCharges = 0;
    QOffset = 0;

    for (size_t i = 0; i<(margin+1); i++)
    {
        EField[i] = 0;
    }

    for (size_t i = margin; i < (numberOfBins+margin);i++)
    {
        EField[i+1] = (EField[i]*boundaryArea[i]*permittivity+polarity*(rhoPosIon[i]-rhoNegIon[i]-rhoE[i])*binVolume[i]*10000.)/(boundaryArea[i+1]*permittivity); //permittivity in pF/m, rho in pC, area in mm^3 -> E in V/cm
        UFromCharges += (EField[i+1]+EField[i])/2.0*binLength/1.0e7; //in V, calculate the potential difference, i.e. the E-field in the center of the bin times the bin width
    }

    for (size_t i = margin+numberOfBins+1; i < (numberOfBins+2*margin+1);i++)
    {
        EField[i] = EField[numberOfBins+margin];
    }

    //caluclate the difference between potential from charges and the applied voltage
    double UDiff = U - UFromCharges;
    //calucalte the charge on central wire resulting from that difference
    QOffset = UDiff*2.0*pi*length*1.0e-3*permittivity/log(r2/r1); //in pC
    //normalize the E-Field using that charge
    for (size_t i = 0; i <(numberOfBins+2*margin+1); i++)
    {
        EField[i] += QOffset*10000./(boundaryArea[i]*permittivity);
    }


    return EField;
}


plane_parallel::plane_parallel(double pPlateDistance, double pDiameter, double pNw, double pKq, double pU, int pPolarity, size_t pNumberOfBins, size_t pMargin)
{
    //initalize all the parameters of the chamber
    type = 1;
    plateDistance = pPlateDistance;
    diameter = pDiameter;
    volume = pi/4.*diameter*diameter*plateDistance;
    Nw = pNw;
    kQ = pKq;
    numberOfBins = pNumberOfBins;
    margin = pMargin;
    U = pU;
    polarity = pPolarity;
    E0 = U*10./plateDistance;

    //create the bin arrays
    binVolume = new double[numberOfBins+2*margin];

    binLength  = plateDistance*1.0e6/pNumberOfBins; //binLength in nm

    for(size_t i = 0; i < (numberOfBins + 2*margin);i++)
    {
        binVolume[i] = plateDistance/pNumberOfBins*pi/4.*diameter*diameter;
    }

    EField = new double[numberOfBins+2*margin+1]();

}

double* plane_parallel::calc_E(const double * const  rhoPosIon, const double * const  rhoNegIon, const double * const  rhoE,const double permittivity)
{
    EOffset = 0;
    UFromCharges = 0;
    EField[0] = 0;

    for (size_t i = 0; i<(margin+1); i++)
    {
        EField[i] = 0;
    }

    for(size_t i = margin; i < (numberOfBins+margin);i++)
    {
        EField[i+1] = EField[i]+(rhoPosIon[i]-rhoNegIon[i]-rhoE[i])*binLength/(permittivity*100.);
        // permittivity in pF/m, rho in pC/mm^3, binLength in nm -> E in V/cm
        // this is the same calculation as for the cylinder chamber, but all boundaryArea[i] and binVolumes[i] are the same, which allows a lot of simplification.
        UFromCharges += (EField[i+1]+EField[i])/2.0*binLength/1.0e7; //in V, calculate the potential difference, i.e. the E-field in the center of the bin times the bin width
    }


    for (size_t i = margin+numberOfBins+1; i < (numberOfBins+2*margin+1);i++)
    {
        EField[i] = EField[numberOfBins+margin];
    }

    //determine the offset, so that the potential difference of the chamber stays constant (boundary condition)
    EOffset=10.*(U-UFromCharges)/plateDistance;

    for(size_t i = 0;i< (numberOfBins+2*margin+1);i++)
    {
        EField[i] += EOffset;
    }

    return (EField);
}

plane_parallel::~plane_parallel(){
    delete[] binVolume;
    delete[] EField;
}
