#include "beam_model.hpp"

#include <cmath>
#include <iostream>

#include "config.hpp"

continuous_pulse::continuous_pulse(double pIonizationRate, double pPulseDuration) //give the charge created per ns and bin
{
    mIonizationRate = pIonizationRate;
    mPulseDuration = pPulseDuration;
}

//return ionizations in the time interval from now-timestep to now (now = totalTime)
double continuous_pulse::ionization(double pTotalTime, double pTimeStep)
{
    if (pTotalTime <= mPulseDuration){
        return pTimeStep*mIonizationRate;
    }
    else if (pTotalTime-pTimeStep <= mPulseDuration) {
        return (mPulseDuration-pTotalTime+pTimeStep)*mIonizationRate;
    }
    else{
        return 0.0;
    }
}

void continuous_pulse::print_characeristics()
{
    std::cout << "continous pulse with charge creation rate in pC/mm^3: " << mIonizationRate<<std::endl;
}



elbe_micropulses::elbe_micropulses (double pChargePerPulse, double pPulseDuration) //give the charge created per micro pulse in pC/mm^3
{
    mChargePerPulse = pChargePerPulse;
    mPulseDuration = pPulseDuration;
}

//return the amount of created charge, based on the number of micro pulses in the interval (now-timestep, now]
double elbe_micropulses::ionization(double pTotalTime, double pTimeStep){
    double endTime;
    double startTime = pTotalTime - pTimeStep;
    if (pTotalTime > mPulseDuration){
        endTime = mPulseDuration;
    }
    else{
        endTime = pTotalTime;
    }
    int pulsesInInterval = (int)std::ceil((endTime - 2.5e-3) / 77) - (int)std::ceil((startTime -2.5e-3) / 77);
    #if LOGLEVEL > 0
        if (pulsesInInterval > 0){
            std::streamsize oldPrecision = std::cout.precision(16);
            std::cout<<pulsesInInterval<<" micropulses at "<<pTotalTime<<" ns"<<std::endl<<std::endl;
            std::cout.precision(oldPrecision);
        }
    #endif
    return mChargePerPulse*pulsesInInterval;

}

void elbe_micropulses::print_characeristics()
{
    std::cout << "micro-pulsed pulse with micro-pulse-charge in pC/mm^3: " << mChargePerPulse<<std::endl;
}
