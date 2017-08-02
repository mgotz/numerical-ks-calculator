#ifndef LOOP_TEMPLATE_H
#define LOOP_TEMPLATE_H

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "functions.hpp"
#include "config.hpp"

/*
 * The core loop of the numeric calculation as a template for different chamber and polarity types
*/

using namespace std;

template<int Polarity, int ChamberType>
loop_return* time_loop(variable_package &userInput, ostream &log, char* startTime){


    //transfer user input to local variables
    const double reactionRateLimit=userInput.rateLimit;         //maximum allowed reaction rate per time step and space bin
    const double speedLimit = userInput.speedLimit;             //maximum allowed speed in bins/timestep -> determines the timestep size
    const double pulseDuration=userInput.pulseduration;         //in ns
    double E0=userInput.this_chamber->return_average_E();       //in V/cm
    const double Q=userInput.Q;                                 //total liberated charge in pC
    active_medium* fillgas = userInput.fillgas;
    const double epsilon=8.8542*fillgas->relative_permittivity(E0);// PikoFarad/Meter
    const int margin=userInput.margin;
    const int NumberOfPoints=userInput.numberOfBins+2*margin;
    double binLength = userInput.this_chamber->return_binLength(); //length in direction of the movement of each bin in nm
    double* binVolume = userInput.this_chamber->return_binVolume(); //volume of the bin in mm^3
    double* r = userInput.this_chamber->return_radius(); //radii of cell bin boundaries

    //create the return struct and make reference to it, for convinience and to minimize the changes from an older version
    loop_return * results = new loop_return ();

    double& electronCollectionTime = results->electronCollectionTime;
    double& beamOnTime = results->beamOnTime;
    double& elapsedTime = results->elapsedTime;

    double& Eminall = results->EminAll;
    double& Emaxall = results->EmaxAll;

    bool& PositivityWarning = results->positivityWarning;

    double& liberatedCharge = results->liberatedCharges;
    double& collectedCharge1 = results->collectedCharges[0];
    double& collectedCharge2 = results->collectedCharges[1];
    double& collectedCharge3 = results->collectedCharges[2];
    double& recombinedCharge = results->ionRecombinedCharges;
    double& directlyRecombinedCharge = results->directlyRecombinedCharges;
    double& attachedCharge=results->attachedCharge;

    double& sum1 = results->remainingCharges[0];
    double& sum2 = results->remainingCharges[1];
    double& sum3 = results->remainingCharges[2];

    //initialize all the required variables
    double Emin=E0,Emax=E0,maxEMedium=fillgas->mMaxE;
    double chargeMax=0;
    double timestep=0.0;
    double alpha=0;   //recombination: I+ I-
    double beta = 0;  //recombination: e- I+
    double gamma=0;   //attachment: e- --> I-

    const size_t maxNumberOfSteps=1000000000;
    //const size_t maxNumberOfSteps=100;

    // electric field and velocities shall be defined at the cell boundaries of the space discretization | x_i-1 | x_i | x_i+1 |,
    // therefore field and velocities need one more point than x (charge densities)
    // () initalizes all values to 0

    double* v1 = new double[NumberOfPoints+1](); //velocity pos. ions in nm/ns, i.e. m/s
    double* v2 = new double[NumberOfPoints+1](); //velocity neg. ions
    double* v3 = new double[NumberOfPoints+1](); //velocity e-

    double* rho1 = new double[NumberOfPoints](); //in pC/mm^3
    double* rho2 = new double[NumberOfPoints]();
    double* rho3 = new double[NumberOfPoints]();
    double* rho1temp = new double[NumberOfPoints]();
    double* rho2temp = new double[NumberOfPoints]();
    double* rho3temp = new double[NumberOfPoints]();

    double* E = userInput.this_chamber->calc_E(rho1,rho2,rho3,epsilon);   //E-field in V/cm


    double ionizationPerVolume = 0;

    double v1min=numeric_limits<double>::max(),v1max=0;
    double v2min=numeric_limits<double>::max(),v2max=0;
    double v3min=numeric_limits<double>::max(),v3max=0;
    double vIonMax=0;
    double GammaMax=1.0,AlphaMax=0; //save the maximum reaction rates to calculate the timestep

    //helper index variables
    int j=0;
    int k=0;

    ostringstream errorStrStream;


    #if CONSTIONMOB == 1
        double mobPos = fillgas->mobility_pos(E0);
        double mobNeg = fillgas->mobility_neg(E0);
    #endif

    #if LOGLEVEL > 1
        cout<<"begin time-loop"<<endl;
    #endif

    #if CHARGETIMELINE == 1
        ostringstream fileName;
        fileName << "Pos_Ion_TimeEvol_"<< startTime<<".txt";
        std::ofstream PosIonTimeline(fileName.str().c_str());

        fileName.str("");
        fileName << "Neg_Ion_TimeEvol_"<< startTime<<".txt";
        std::ofstream NegIonTimeline(fileName.str().c_str());

        fileName.str("");
        fileName << "Electron_TimeEvol_"<< startTime<<".txt";
        std::ofstream ElectronTimeline(fileName.str().c_str());

        fileName.str("");
        fileName << "EField_TimeEvol_"<< startTime<<".txt";
        std::ofstream EField(fileName.str().c_str());
    #endif

    size_t step=0; //counts the number of timesteps taken, declared here to be available after the loop has exited

    try
    {
    for (step=0;step<maxNumberOfSteps && (sum1>Q*1.0e-06 || sum2>Q*1.0e-06 || sum3>Q*1.0e-06 || beamOnTime < 0) ;step++)
        //abort when either maxNumberOfSteps exceeded (avoid infinite run), or when there are less than 10^-6 of the created charge left, but not before irradiation is complete

        {
        //do some regular logging every 1000 time steps
        #if LOGLEVEL > 0
        if  ((step-1)%1000 == 0)
            {
            log<<"Now at step "<<step<<" time "<<elapsedTime<<" ns; timestep " <<timestep<<" ns"<<"; timestep from electron speed limit: "<<speedLimit/fabs(v3max)<<endl;
            log<< "extreme velocities (I+, I-, e, I) "<<v1min*timestep<<" "<<v1max*timestep<<" "<<v2min*timestep<<" "<<v2max*timestep<<" "<<v3min*timestep<<" "<<v3max*timestep<<" "<<vIonMax*timestep<<endl;
            log<<"E-field: min: " << Emin << ", max: " << Emax<<" V/cm" << endl;
            log<<"in chamber: I+: "<<sum1<<", I-: "<<sum2<<", e: "<<sum3<<" pC"<<endl;
            log<<"collected: I+: "<<collectedCharge1<<", I-: "<<collectedCharge2<<", e:"<<collectedCharge3<<", attached: "<<attachedCharge<<" pC"<<endl;
            log<<"positivity warning: "<< bool(PositivityWarning) <<endl;
            #if ELECTRONCOMBO == 1
                log<<"beta value: "<< beta << "; ";
            #endif
            log<<"alpha value: "<<alpha << "; directLosses: "<<directlyRecombinedCharge<<endl;
            log<<endl;
            }
        #endif


        //primary ionizations
        if (beamOnTime < 0){
            ionizationPerVolume = userInput.this_beam->ionization(elapsedTime,timestep);
            for (int i=margin; i<NumberOfPoints-margin;i++)
            {
            rho1temp[i]+= ionizationPerVolume;
            rho3temp[i]+= ionizationPerVolume;
            liberatedCharge+= ionizationPerVolume*binVolume[i];
            }
        }

        //use one alpha for every spot in the chamber (independent of E)
        #if CONSTRECOMB == 1
            alpha = fillgas->volume_recombination(E0) * timestep / (1.6e-07);
        #endif

        // recombination and drift of ions
        for (int i=0; i<(NumberOfPoints-1);i++)
            {
            //recombination, calculate alpha at each bin from E-field if CONSTRECOMB = 0
            #if CONSTRECOMB == 0
                alpha=fillgas->volume_recombination(0.5*(E[i]+E[i+1])) * timestep /(1.6e-07);
            #endif

            rho1temp[i]-=alpha*rho1[i]*rho2[i];
            rho2temp[i]-=alpha*rho1[i]*rho2[i];
            recombinedCharge+=alpha*rho1[i]*rho2[i]*binVolume[i];

            #if LOGLEVEL > 1
            if (rho2temp[i]<0||rho1temp[i]<0)
            {
                log<<"positivity violation at time "<< elapsedTime << " during recombination at bin "<<i<<endl;
                log<<" E: "<<(0.5*(E[i]+E[i+1]))<<", alpha: "<<alpha<<", rho1: "<<rho1[i]<<", rho2: "<<rho2[i]<<endl;
                log<<" recombinedCharge: "<<alpha*rho1[i]*rho2[i]<<", temp1: "<<rho1temp[i]<<", temp2: "<<rho2temp[i]<<endl;
            }
            #endif

            //drift is discretized in 1st order upwind: -rho[i]*v[i]+rho[i-1]*v[i-1] direction may need to be adapted depending on sign of v


            //new version which assumes the directionality of v through Polarity template parameter and consequently takes care of the directionality
            if (ChamberType == 1) //no need for polarity, because in plane-parallel geometry it is perfectly symmetric
            {
                rho1temp[i+1] += timestep*(v1[i+1]*rho1[i] - v1[i+2]*rho1[i+1]);
                rho2temp[i] += timestep*(v2[i+1]*rho2[i+1] - v2[i]*rho2[i]);
            }
            else if (ChamberType == 2 && Polarity == -1)
            {
                rho1temp[i] += timestep*(v1[i+1]*r[i+1]*rho1[i+1] - v1[i]*r[i]*rho1[i]) * 2.0/(r[i]+r[i+1]);
                rho2temp[i+1] += timestep*(v2[i+1]*r[i+1]*rho2[i] - v2[i+2]*r[i+2]*rho2[i+1]) * 2.0/(r[i+1]+r[i+2]);
            }
            else if (ChamberType == 2 && Polarity == 1)
            {
                rho1temp[i+1] += timestep*(v1[i+1]*r[i+1]*rho1[i] - v1[i+2]*r[i+2]*rho1[i+1]) * 2.0/(r[i+1]+r[i+2]);
                rho2temp[i] += timestep*(v2[i+1]*r[i+1]*rho2[i+1] - v2[i]*r[i]*rho2[i]) * 2.0/(r[i]+r[i+1]);
            }
            else
            {
                errorStrStream.str("unkown chamber type: ");
                errorStrStream<<ChamberType;
                throw invalid_argument(errorStrStream.str());
            }


            #if LOGLEVEL > 1
                if (Polarity == -1)
                {
                    j = i; //index for rho1
                    k = i+1; //index for rho2
                }
                else
                {
                    j = i+1;
                    k = i;
                }

                if (rho1temp[j]<0)
                {
                    log<<"positivity violation at step "<< step <<"("<<elapsedTime<<" ns)"<<" during rho1 move at bin "<<j<<endl;
                    log<<"v1[i]: "<<v1[j]*timestep<<" v1[i+1]: "<<v1[j+1]*timestep<<endl;
                    log<<"rho1[i]: "<<rho1[j]<<showpos<<" rho1[i"<<Polarity<<"]: "<<noshowpos<<rho1[k]<<endl;
                }
                if (rho2temp[k]<0)
                {
                    log<<"positivity violation at step "<< step <<"("<<elapsedTime<<" ns)"<<" during rho2 move at bin "<<k<<endl;
                    log<<"v2[i]: "<<v2[k]*timestep<<" v2[i+1]: "<<v2[k+1]*timestep<<endl;
                    log<<"rho2[i]: "<<rho2[k]<<showpos<<" rho2[i"<<Polarity<<"]: "<<noshowpos<<rho2[j]<<endl;
                }
            #endif
            }


        //electron attachment, drift and direct recombination (e- + I+ -> I) (only executed if there are still electrons left)
        if (electronCollectionTime < 0){
            #if (CONSTELECTRONCOMBO == 1) && (ELECTRONCOMBO == 1)
                beta = fillgas->direct_recombination(E0) * timestep / (1.6e-07);
            #endif
            for(int i=0; i<(NumberOfPoints-1); i++)
            {
                //attachment
                gamma=timestep*fillgas->attachment(0.5*(E[i]+E[i+1]));	//Umw. e- --> I-
                rho2temp[i]+=gamma*rho3[i];
                rho3temp[i]-=gamma*rho3[i];
                attachedCharge += gamma*rho3[i]*binVolume[i];

                #if LOGLEVEL > 1
                    if (rho3temp[i]<0){log<<"positivity violation at time "<< elapsedTime << " during reaction at bin "<<i<<std::endl;}
                    if (rho2temp[i]<0){
                        log<<"positivity violation at time "<< elapsedTime << " during reaction at bin "<<i<<std::endl;
                        log<<"E: "<<(0.5*(E[i]+E[i+1]))<<"V/cm, gamma: "<<gamma<<" rho3: "<<rho3[i]<<std::endl;}
                #endif

                //do the direct recombination if desired
                #if ELECTRONCOMBO == 1
                    #if CONSTELECTRONCOMBO == 0
                        beta = fillgas->direct_recombination(0.5*(E[i]+E[i+1])) * timestep / (1.6e-07);
                    #endif

                    rho1temp[i] -= beta*rho3[i]*rho1[i];
                    rho3temp[i] -= beta*rho3[i]*rho1[i];
                    directlyRecombinedCharge += beta*rho3[i]*rho1[i]*binVolume[i];
                #endif

                //new version which assumes all v are positive (for both ion species and electrons) and takes care of directionalty throught the expression her

                if (ChamberType == 1)
                {
                    rho3temp[i] += timestep*(v3[i+1]*rho3[i+1] - v3[i]*rho3[i]);
                }
                else if (ChamberType == 2 && Polarity == -1)
                {
                    rho3temp[i+1] += timestep*(v3[i+1]*r[i+1]*rho3[i] - v3[i+2]*r[i+2]*rho3[i+1]) * 2.0/(r[i+1]+r[i+2]);
                }
                else if(ChamberType == 2 && Polarity == 1)
                {
                    rho3temp[i] += timestep*(v3[i+1]*r[i+1]*rho3[i+1] - v3[i]*r[i]*rho3[i]) * 2.0/(r[i]+r[i+1]);
                }
                else
                {
                    errorStrStream.str("unkown chamber type: ");
                    errorStrStream<<ChamberType;
                    throw invalid_argument(errorStrStream.str());
                }


                #if LOGLEVEL > 1
                    if (rho3temp[i]<0)
                    {
                        log<<"positivity violation at step "<< step <<"("<<elapsedTime<<" ns)"<<" during rho3 move at bin "<<i<<endl;
                        log<<"v3[i]: "<<v3[i]*timestep<<" v3[i+1]: "<<v3[i+1]*timestep<<" v3[i-1]: "<<v3[i-1]*timestep<<endl;
                        log<<"rho3[i]: "<<rho3[i]<<showpos<<" rho3[i"<<Polarity<<"]: "<<noshowpos<<rho3[i+Polarity]<<endl;
                    }
                #endif
            }
        }

        //sum all the charges
        sum1=sum2=sum3=chargeMax=0;

        for (int i=0; i<NumberOfPoints;i++)
            {
            if (rho1temp[i]>chargeMax) {chargeMax=rho1temp[i];}
            if (rho2temp[i]>chargeMax) {chargeMax=rho2temp[i];}
            #if LOGLEVEL > 0
            if (rho1temp[i]<0) {PositivityWarning=true;rho1temp[i]=0;}
            if (rho2temp[i]<0) {PositivityWarning=true;rho2temp[i]=0;}
            if (rho3temp[i]<0) {PositivityWarning=true;rho3temp[i]=0;}
            #endif
            sum1+=rho1temp[i]*binVolume[i];
            sum2+=rho2temp[i]*binVolume[i];
            sum3+=rho3temp[i]*binVolume[i];

            }

        //zero the margin = collected charge
        for (int i=0; i<margin;i++)
            {
            collectedCharge1 = collectedCharge1 + Polarity*rho1temp[NumberOfPoints-i-1]*binVolume[NumberOfPoints-i-1];
            rho1temp[NumberOfPoints-i-1]=0;
            collectedCharge1 = collectedCharge1 - Polarity*rho1temp[i]*binVolume[i];
            rho1temp[i]=0;
            collectedCharge2 = collectedCharge2 + Polarity*rho2temp[i]*binVolume[i];
            rho2temp[i]=0;
            collectedCharge2 = collectedCharge2 - Polarity*rho2temp[NumberOfPoints-i-1]*binVolume[NumberOfPoints-i-1];
            rho2temp[NumberOfPoints-i-1]=0;
            }


        if (electronCollectionTime<0)
            {
            for (int i=0; i<margin;i++)
                {
                collectedCharge3 = collectedCharge3 + Polarity*rho3temp[i]*binVolume[i];
                rho3temp[i]=0;
                collectedCharge3 = collectedCharge3 - Polarity*rho3temp[NumberOfPoints-i-1]*binVolume[NumberOfPoints-i-1];
                rho3temp[NumberOfPoints-i-1]=0;
                }
            }

        //transfer the temp data
        for(int i=0;i<NumberOfPoints;i++)
            {
            rho1[i] =rho1temp[i];
            rho2[i] =rho2temp[i];
            rho3[i] =rho3temp[i];
    //        cout<<i<<", 1: "<<rho1[i]<<" ,2: "<<rho2[i]<<" ,3: "<<rho3[i]<<endl;
            }




    #if CHARGETIMELINE == 1
        if (elapsedTime<pulseDuration){
            if  (step/1000 == step/1000.){
                PosIonTimeline << elapsedTime << "\t";
                NegIonTimeline << elapsedTime << "\t";
                ElectronTimeline << elapsedTime << "\t";
                EField << elapsedTime << "\t";
                for(int i=0; i<NumberOfPoints;i++){
                    PosIonTimeline << rho1[i] << "\t";
                    NegIonTimeline << rho2[i] << "\t";
                    ElectronTimeline << rho3[i] << "\t";
                    EField << E[i] << "\t";
                }
                EField << E[NumberOfPoints] << "\t";
                PosIonTimeline << std::endl;
                NegIonTimeline << std::endl;
                ElectronTimeline << std::endl;
                EField << std::endl;
            }
        }
        else if (step%10 == 0){
            PosIonTimeline << elapsedTime << "\t";
            NegIonTimeline << elapsedTime << "\t";
            ElectronTimeline << elapsedTime << "\t";
            EField << elapsedTime << "\t";
            for(int i=0; i<NumberOfPoints;i++)
            {
                PosIonTimeline << rho1[i] << "\t";
                NegIonTimeline << rho2[i] << "\t";
                ElectronTimeline << rho3[i] << "\t";
                EField << E[i] << "\t";
            }
            EField << E[NumberOfPoints] << "\t";
            PosIonTimeline << std::endl;
            NegIonTimeline << std::endl;
            ElectronTimeline << std::endl;
            EField << std::endl;
        }

    #endif

        //calculate the E-field
        #if CONSTFIELD != 1
            E = userInput.this_chamber->calc_E(rho1,rho2,rho3,epsilon);
        #endif
        Emin=numeric_limits<double>::max();
        Emax=0;
        v3max = 0;
        v2max = 0;
        v1max = 0;
        v3min = numeric_limits<double>::max();
        v2min = numeric_limits<double>::max();
        v1min = numeric_limits<double>::max();
        //calculate Emax, Emin and calcualte the velocities
        for (int i=0; i<NumberOfPoints+1;i++)
            {
            if (i>margin && i<(NumberOfPoints-margin+1) && E[i]<Emin) {Emin=E[i];}
            if (i>margin && i<(NumberOfPoints-margin+1) && E[i]>Emax) {Emax=E[i];}

            //ion velocities, if constant mobilities are used, the function calls can be reduced
            //mobility in m^2/Vs, E in V/cm, binLength in nm -> v in bin/ns


            #if CONSTIONMOB == 1
                v1[i]=(mobPos*E[i])/binLength*100.;
                v2[i]=(mobNeg*E[i])/binLength*100.;
            #else
                v1[i]=(fillgas->mobility_pos(E[i])*E[i])/binLength*100.;
                v2[i]=(fillgas->mobility_neg(E[i])*E[i])/binLength*100.;
            #endif
            v3[i]=(fillgas->mobility_e(E[i])*100.*E[i])/binLength;

            if (v1[i]<v1min){v1min = v1[i];}
            if (v1[i]>v1max){v1max = v1[i];}
            if (v2[i]<v2min){v2min = v2[i];}
            if (v2[i]>v2max){v2max = v2[i];}
            if (v3[i]<v3min){v3min = v3[i];}
            if (v3[i]>v3max){v3max = v3[i];}

    //        cout<<i<<": "<<E[i]<<endl;
            }
        if (Emin<Eminall) { Eminall = Emin; }
        if (Emax>Emaxall) { Emaxall = Emax; }
        if (Emax > maxEMedium)
        {
            errorStrStream.str("The E-field exeeds the range of one of the active medium functions. Emax: ");
            errorStrStream<<Emax;
            throw domain_error(errorStrStream.str().c_str());
        }
        //determine extreme velocities
        vIonMax = (v2max > v1max)? v2max : v1max;

        //determine extreme reaction rates
        GammaMax = fillgas->attachment(Emin);
        AlphaMax = chargeMax*fillgas->volume_recombination(Emin)/(1.6e-07);


        //realize the end of electron collection time and signal via electronCollectionTime
        //but dont stop before irradiation is finished (for micro-pulsed important, where electrons might disappear between micro-pulses
        if (sum3<Q*1.0e-07 && beamOnTime>0 && electronCollectionTime<0)
            {
            for (int i=0; i<NumberOfPoints;i++)
                {
                rho3[i]=0;
                rho3temp[i]=0;
                }

            gamma = 0.0;
            electronCollectionTime=elapsedTime;
            log<<"electrons collected after "<<elapsedTime/1000.<<" µs (step:"<<step<<")"<<endl;
            log<<"Alpha (Limit): " << AlphaMax<<" ("<<reactionRateLimit<<"), chargeMax: "<<chargeMax<<endl;
            log<<"total charges (I+.I-,e): " <<sum1<<" "<<sum2<<" "<<sum3<<endl;

            timestep = (reactionRateLimit/AlphaMax > speedLimit/vIonMax)? speedLimit/vIonMax : reactionRateLimit/AlphaMax;
            log<<"new timestep "<<timestep<<" ns"<<endl<<endl;
            }

        //determine timestep
        if (electronCollectionTime < 0)
        {
            timestep = (reactionRateLimit / GammaMax > speedLimit / v3max) ? speedLimit / v3max : reactionRateLimit / GammaMax;
            //limit the timestep during irradiation (relevant for very short pulses or very slow electrons (LIC))
            if (beamOnTime < 0)
            {
                if (timestep > pulseDuration/10.)
                {
                    timestep = pulseDuration/10.;
                }
            }
            #if LOGLEVEL > 1
                log << "step: " << step << " timestep: " << timestep << " time: " << elapsedTime << " electron v: " << timestep*v3max << " gamma: " << GammaMax*timestep<< std::endl;
                log << "   v3min: "<<v3min<< " v3max: "<<v3max<<"; Emax: "<<Emax<<" Emin: "<<Emin<<endl;
                log << "   sum1: "<<sum1<<", sum2: "<<sum2<<", sum3: "<<sum3<<endl;
            #endif
        }
        else
            {
            timestep = (reactionRateLimit/AlphaMax > speedLimit/vIonMax)? speedLimit/vIonMax : reactionRateLimit/AlphaMax;
            #if LOGLEVEL > 1
                log << "step: "<< step << " timestep: " << timestep << " time: " << elapsedTime << " ion v: "<< timestep*vIonMax << " alpha: " << AlphaMax*timestep << std::endl;
                log << v1max << " " << v1min << " "<< v2max << " "<< v2min << std::endl;
                log << Emin << " "<< Emax << std::endl;
            #endif
            }


        //realize irradiation end
        if(elapsedTime>=pulseDuration && beamOnTime<0)
            {
            beamOnTime=elapsedTime;
            log<<"beam off after "<<elapsedTime/1000.<<" µs (step: "<<step<<")"<<endl<<endl;
            }

        elapsedTime += timestep;



        #if LOGLEVEL > 1
            log<<endl;
        #endif


        }

    //time loop complete
    }
    /*handle aborts from inside the loop,
     *free the results memory and set the returned pointer to NULL to inform caller of failure
     *using the exception instead of returning directly at the point of failure, ensures that the cleanup is run
    */
    catch (domain_error e)
    {
        delete results;
        results = NULL;
        cerr<<e.what()<<endl;
    }
    catch (invalid_argument e)
    {
        delete results;
        results = NULL;
        cerr<<"invalid argument: "<<e.what()<<endl;
    }


    //cleanup
    delete[] rho1;
    delete[] rho2;
    delete[] rho3;
    delete[] rho1temp;
    delete[] rho2temp;
    delete[] rho3temp;
    delete[] v1;
    delete[] v2;
    delete[] v3;

    return results;


}


#endif
