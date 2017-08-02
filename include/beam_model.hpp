#ifndef BEAM_MODEL
#define BEAM_MODEL
/*
describe the beam, for now either a conitnuous rate over the entire pulse duration or a micropulse structure like at ELBE
*/
#include "config.hpp"

#include <string>

class beam_model{
public:
    virtual double ionization(double, double) = 0; //=0 makes this an abstract class, i.e. method is not defined and needs to be defined in the children
    virtual void print_characeristics() = 0;
    virtual ~beam_model() = 0;
protected:
    double mPulseDuration;
};

inline beam_model::~beam_model() {}

class continuous_pulse : public beam_model{
public:
    virtual double ionization(double, double);
    virtual void print_characeristics();
    continuous_pulse (double, double); //give the charge created per ns in pC/mm^3 and the pulse duration
private:
    double mIonizationRate;
};


class elbe_micropulses : public beam_model{
public:
    virtual double ionization(double , double );
    virtual void print_characeristics();
    elbe_micropulses (double, double); //give the charge created per micro pulse in pC/mm^3 and the pulse duration
private:
    double mChargePerPulse;


};

#endif // BEAM_MODEL
