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
#include "settings.hpp"

using namespace std;

int main(int argc, char **argv){
// charges are calculated considering the entire chamber cross-section.
// to obtain a charge density the this values must be divided by the area

variable_package userInput;
ofstream logStream;

//define the start time
time_t now = time(0);
struct tm *tm_now = localtime(&now);
char startTime[20];
strftime(startTime, 20, "%Y-%m-%d_%H-%M-%S",tm_now);

//easier formatting of the time, that needs a c++11 function not implemented on the gcc I used
//ostringstream startTime;
//startTime<<std::put_time(localtime(&t),"%Y-%m-%d_%H-%M-%S");


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

//transfer user input to local variables
const double reactionRateLimit=userInput.rateLimit;         //maximum allowed reaction rate per time step and space bin
const double speedLimit = userInput.speedLimit;             //maximum allowed speed in bins/timestep -> determines the timestep size
const double pulseDuration=userInput.pulseduration;         //in ns
const double area=userInput.electrodeArea;     //in mm2
double E0=userInput.E0;                                     //in V/cm
const double binSize=userInput.binSize;                   //in nm
const double Q=userInput.Q;                                 //total liberated charge in pC
active_medium* fillgas = userInput.fillgas;
const double epsilon=8.8542*fillgas->relative_permittivity(E0);// PikoFarad/Meter

//build log stream
ostream log(userInput.logBuffer);


//output start time and also write time to the log, if the log is something else than cout

cout<<"start time: "<<startTime<<endl<<endl;
if (userInput.logBuffer != cout.rdbuf())
{
    log<<"start time: "<<startTime<<endl<<endl;
}


// output the settings to stdout
cout<<"beam settings: "<<endl;
cout<<"dose: "<<userInput.dose<<" mGy"
    <<", pulse duration: "<<userInput.pulseduration<<" ns"
    <<", total charge: "<<Q<<" pC"
    <<", max rection rate: "<<reactionRateLimit
    <<", max velocity: "<<speedLimit
    <<", bins: "<<userInput.numberOfBins<<endl;
userInput.this_beam->print_characeristics();
cout<<"chamber properties: "<<endl;
cout<<"voltage: "<<userInput.voltage<<" V"
    <<", E0: "<<E0<<" V/cm"
    <<", area: "<<area<<" mm2"
    <<", electrode distance: "<<binSize*userInput.numberOfBins/1.0e6<<" mm"
    <<", binsize: "<<binSize<<" nm"<<endl;
cout<<"active medium file: " << userInput.mediumName<<endl;
cout<<"values at E0="<<E0<<" V/cm:"<<endl;
cout<<"volume recombination: "<<fillgas->volume_recombination(E0)
    <<", attachment rate: "<<fillgas->attachment(E0);
cout<<", mu von e: "<<fillgas->mobility_e(E0)
    <<", mu von neg: "<<fillgas->mobility_neg(E0)
    <<", mu von pos: "<<fillgas->mobility_pos(E0);
cout<<endl;
#ifdef ELECTRONCOMBO
    cout<<"using direct recombination: "<<fillgas->direct_recombination(E0)<<" m^3/s"<<endl;
    double beta=0;    //Rek. I+ e-
#else
    cout<<"no direct recombination"<<endl;
#endif
cout<<endl;
cout<<"permittivity: "<<epsilon<<endl;


//log the settings of the active medium
if (userInput.settingsLogFileName != "")
{
    ofstream altOut(userInput.settingsLogFileName.c_str(), ios::out|ios::app);
    if (!altOut)
    {
        cerr<<"failed to open settigns log: "<<userInput.settingsLogFileName.c_str()<<endl;
    }
    else{
        fillgas->write_log(altOut);
    }
}
else
{
    log<<"no settings written"<<endl<<endl;
}

//initialize all the required variables
double E_offset=0;
double Emin=E0,Eminall=numeric_limits<double>::max(),Emax=E0,Emaxall=0,maxEMedium=fillgas->mMaxE;
double chargeMax=0;
double timestep=0.0;
double beamOnTime=-1;
double electronCollectionTime=-1;
double elapsedTime=0;
double alpha=0;   //recombination: I+ I-
double gamma=0;   //attachment: e- --> I-


const int margin=2;
const int NumberOfPoints=userInput.numberOfBins+2*margin;
const size_t maxNumberOfSteps=1000000000;

// electric field and velocities shall be defined at the cell boundaries of the space discretization | x_i-1 | x_i | x_i+1 |,
// therefore field and velocities need one more point than x (charge densities)
// () initalizes all values to 0

double* v1 = new double[NumberOfPoints+1](); //velocity pos. ions in bins/timestep
double* v2 = new double[NumberOfPoints+1](); //velocity neg. ions
double* v3 = new double[NumberOfPoints+1]();  //velocity e-
double* E = new double[NumberOfPoints+1](); //E-field in V/cm
double* E_int = new double[NumberOfPoints+1](); //help variable to calculate E-field

double* rho1 = new double[NumberOfPoints]();
double* rho2 = new double[NumberOfPoints]();
double* rho3 = new double[NumberOfPoints]();
double* rho1temp = new double[NumberOfPoints]();
double* rho2temp = new double[NumberOfPoints]();
double* rho3temp = new double[NumberOfPoints]();

bool PositivityWarning=false;
double ionizationPerBin = 0;
double liberatedCharge=0;
double collectedCharge1=0;
double collectedCharge2=0;
double collectedCharge3=0;
double recombinedCharge=0;
double directlyRecombinedCharge=0;

double help;
double v1min=numeric_limits<double>::max(),v1max=0;
double v2min=numeric_limits<double>::max(),v2max=0;
double v3min=numeric_limits<double>::max(),v3max=0;
double vIonMax=0;
double sum1=numeric_limits<double>::max(),sum2=numeric_limits<double>::max(),sum3=numeric_limits<double>::max();
double GammaMax=1.0,AlphaMax=0; //save the maximum reaction rates to calculate the timestep

#ifdef CONSTIONMOB
    double mobPos = fillgas->mobility_pos(E0);
    double mobNeg = fillgas->mobility_neg(E0);
#endif

#if DEBUGLEVEL > 1
    cout<<"begin time-loop"<<endl;
#endif

#ifdef CHARGETIMELINE
    ostringstream fileName;
    fileName << "Pos_Ion_TimeEvol_"<< startTime<<".txt";
    std::ofstream PosIonTimeline(fileName.str());

    fileName.str("");
    fileName << "Neg_Ion_TimeEvol_"<< startTime<<".txt";
    std::ofstream NegIonTimeline(fileName.str());

    fileName.str("");
    fileName << "Electron_TimeEvol_"<< startTime<<".txt";
    std::ofstream ElectronTimeline(fileName.str());

    fileName.str("");
    fileName << "EField_TimeEvol_"<< startTime<<".txt";
    std::ofstream EField(fileName.str());
#endif

size_t step=0; //counts the number of timesteps taken, declared here to be available after the loop has exited

for (step=0;step<maxNumberOfSteps && (sum1>Q*1.0e-06 || sum2>Q*1.0e-06 || sum3>Q*1.0e-06 || beamOnTime < 0) ;step++)
    //abort when either maxNumberOfSteps exceeded (avoid infinite run), or when there are less than 10^-6 of the created charge left, but not before irradiation is complete

	{
    //do some regular logging every 1000 time steps
    #if DEBUGLEVEL > 0
    if  ((step-1)%1000 == 0)
		{
        log<<"Now at step "<<step<<" time "<<elapsedTime<<" ns; timestep " <<timestep<<" ns"<<"; timestep from electron speed limit: "<<speedLimit/fabs(v3max)<<endl;
        log<< "extreme velocities (I+, I-, e, I) "<<v1min*timestep<<" "<<v1max*timestep<<" "<<v2min*timestep<<" "<<v2max*timestep<<" "<<v3min*timestep<<" "<<v3max*timestep<<" "<<vIonMax*timestep<<endl;
        log<<"E-field: min: " << Emin << ", max: " << Emax<<" V/cm" << endl;
        log<<"charge: I+: "<<sum1<<", I-: "<<sum2<<", e: "<<sum3<<", collected: "<<collectedCharge3<<" pC"<<endl;
        log<<"positivity warning: "<< bool(PositivityWarning) <<endl;
        #ifdef ELECTRONCOMBO
            log<<"beta value: "<< beta << "; ";
        #endif
        log<<"alpha value: "<<alpha << "; directLosses: "<<directlyRecombinedCharge<<endl;
        log<<endl;
		}
    #endif
	
	
    //primary ionizations
    if (beamOnTime < 0){
        ionizationPerBin = userInput.this_beam->ionization(elapsedTime,timestep);
        for (int i=margin; i<NumberOfPoints-margin;i++)
		{
        rho1temp[i]+= ionizationPerBin;
        rho3temp[i]+= ionizationPerBin;
        liberatedCharge+= ionizationPerBin;
		}
	}
	
    //use one alpha for every spot in the chamber (independent of E)
    #ifdef CONSTRECOMB
        alpha = fillgas->volume_recombination(E0) / area*timestep / binSize / (1.6e-07)*(1.0e+06);
    #endif

    // recombination and drift of ions
    for (int i=1; i<(NumberOfPoints-1);i++)
		{
        //recombination, calculate alpha at each bin from E-field unless CONSTRECOMB is defined
        #ifndef CONSTRECOMB
            alpha=fillgas->volume_recombination(0.5*(E[i]+E[i+1]))/area*timestep/binSize/(1.6e-07)*(1.0e+06);
        #endif
		
		rho1temp[i]-=alpha*rho1[i]*rho2[i];	
		rho2temp[i]-=alpha*rho1[i]*rho2[i];	
        recombinedCharge+=alpha*rho1[i]*rho2[i];

        #if DEBUGLEVEL > 1
        if (rho2temp[i]<0||rho1temp[i]<0)
        {
            std::cerr<<"positivity violation at time "<< elapsedTime << " during recombination at bin "<<i<<" E: "<<(0.5*(E[i]+E[i+1]))<<std::endl;
        }
        #endif

        //drift is discretized in 1st order upwind: -rho[i]*v[i]+rho[i-1]*v[i-1] direction may need to be adapted depending on sign of v

        //positive ions

        //an old version, which evaluates a conditional for each v, which is very time consuming (branching) but safer
        //rho1temp[i] += (v1[i]>0) ? (timestep*(v1[i] * rho1[i-1] - v1[i+1] * rho1[i])) : (timestep*(v1[i]*rho1[i] - v1[i+1]*rho1[i+1]));

        //new version which assumes all v are positive (for both ion species and electrons) and takes care of directionalty throught the expression her
        rho1temp[i] += timestep*(v1[i] * rho1[i-1] - v1[i+1] * rho1[i]);

        #if DEBUGLEVEL > 1
		if (rho1temp[i]<0){
            log<<"positivity violation at time "<< elapsedTime<<" during rho1 move at bin "<<i<<std::endl;
            log<<"v1[i]: "<<v1[i]*timestep<<" v1[i+1]: "<<v1[i+1]*timestep<<" v1[i-1]: "<<v1[i-1]*timestep<<std::endl;
            log<<"rho1[i]: "<<rho1[i]<<" rho1[i+1]: "<<rho1[i+1]<<" rho1[i-1]: "<<rho1[i-1]<<std::endl;
			}
		#endif

        //negative ions

        //an old version, which evaluates a conditional for each v, which is very time consuming (branching) but safer
        //rho2temp[i] += (v2[i]>0) ? (timestep*(v2[i] * rho2[i-1] - v2[i+1] * rho2[i])) : (timestep*(v2[i]*rho2[i] - v2[i+1]*rho2[i+1]));

        //new version which assumes all v are positive (for both ion species and electrons) and takes care of directionalty throught the expression her
        rho2temp[i] += timestep*(v2[i]*rho2[i] - v2[i+1]*rho2[i+1]);

        #if DEBUGLEVEL > 1
		if (rho2temp[i]<0){
            log<<"positivity violation at time "<< elapsedTime<<" during rho2 move at bin "<<i<<std::endl;
            log<<"v2[i]: "<<v2[i]*timestep<<" v2[i+1]: "<<v2[i+1]*timestep<<" v2[i-1]: "<<v2[i-1]*timestep<<std::endl;
            log<<"rho2[i]: "<<rho2[i]<<" rho2[i+1]: "<<rho2[i+1]<<" rho2[i-1]: "<<rho2[i-1]<<std::endl;
			}
		#endif
		}
	
	

    //electron attachment, drift and direct recombination (e- + I+ -> I) (only executed if there are still electrons left)
    if (electronCollectionTime < 0){
        #if defined(CONSTELECTRONCOMBO) && defined (ELECTRONCOMBO)
            beta = fillgas->direct_recombination(E0)/area*timestep/binSize/(1.6e-07)*(1.0e+06);
        #endif
        for(int i=1; i<(NumberOfPoints-1); i++){
            //attachment
            gamma=timestep*fillgas->attachment(0.5*(E[i]+E[i+1]));	//Umw. e- --> I-
			rho2temp[i]+=gamma*rho3[i];	
			rho3temp[i]-=gamma*rho3[i];

            #if DEBUGLEVEL > 1
                if (rho3temp[i]<0){log<<"positivity violation at time "<< elapsedTime << " during reaction at bin "<<i<<std::endl;}
                if (rho2temp[i]<0){
                    log<<"positivity violation at time "<< elapsedTime << " during reaction at bin "<<i<<std::endl;
                    log<<"E: "<<(0.5*(E[i]+E[i+1]))<<"V/cm, gamma: "<<gamma<<" rho3: "<<rho3[i]<<std::endl;}
			#endif

            //do the direct recombination if desired
            #ifdef ELECTRONCOMBO
                #ifndef CONSTELECTRONCOMBO
                    beta = fillgas->direct_recombination(0.5*(E[i]+E[i+1]))/area*timestep/binSize/(1.6e-07)*(1.0e+06);
                #endif

                rho1temp[i] -= beta*rho3[i]*rho1[i];
                rho3temp[i] -= beta*rho3[i]*rho1[i];
                directlyRecombinedCharge += beta*rho3[i]*rho1[i];
            #endif

            //an old version, which evaluates a conditional for each v, which is very time consuming (branching) but safer
            //rho3temp[i] += (v3[i]>0) ? (timestep*(v3[i] * rho3[i-1] - v3[i+1] * rho3[i])) : (timestep*(v3[i]*rho3[i] - v3[i+1]*rho3[i+1]));

            //new version which assumes all v are positive (for both ion species and electrons) and takes care of directionalty throught the expression her
            rho3temp[i] += timestep*(v3[i]*rho3[i] - v3[i+1]*rho3[i+1]);

            #if DEBUGLEVEL > 1
			if (rho3temp[i]<0){
                log<<"positivity violation at step "<< step <<"("<<elapsedTime<<" ns)"<<" during rho3 move at bin "<<i<<std::endl;
                log<<"v3[i]: "<<v3[i]*timestep<<" v3[i+1]: "<<v3[i+1]*timestep<<" v3[i-1]: "<<v3[i-1]*timestep<<std::endl;
                log<<"rho3[i]: "<<rho3[i]<<" rho3[i+1]: "<<rho3[i+1]<<" rho3[i-1]: "<<rho3[i-1]<<std::endl;
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
        #if DEBUGLEVEL > 0
		if (rho1temp[i]<0) {PositivityWarning=true;rho1temp[i]=0;}
		if (rho2temp[i]<0) {PositivityWarning=true;rho2temp[i]=0;}
		if (rho3temp[i]<0) {PositivityWarning=true;rho3temp[i]=0;}
        #endif
		sum1+=rho1temp[i];
		sum2+=rho2temp[i];
		sum3+=rho3temp[i];
		
		}

    //zero the margin = collected charge
    for (int i=0; i<margin;i++)
		{
        collectedCharge1+=rho1temp[NumberOfPoints-i-1];
		rho1temp[NumberOfPoints-i-1]=0;		
        collectedCharge1-=rho1temp[i];
		rho1temp[i]=0;	
        collectedCharge2+=rho2temp[i];
		rho2temp[i]=0;		
        collectedCharge2-=rho2temp[NumberOfPoints-i-1];
		rho2temp[NumberOfPoints-i-1]=0;		
		}


    if (electronCollectionTime<0)
		{		
        for (int i=0; i<margin;i++)
			{
            collectedCharge3+=rho3temp[i];
			rho3temp[i]=0;		
            collectedCharge3-=rho3temp[NumberOfPoints-i-1];
			rho3temp[NumberOfPoints-i-1]=0;
			}
		}
	
    //transfer the temp data
    for(int i=0;i<NumberOfPoints;i++)
		{
		rho1[i] =rho1temp[i];	 
		rho2[i] =rho2temp[i];	 
		rho3[i] =rho3temp[i];
		}	




#ifdef CHARGETIMELINE
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

#endif

    //calculate the E-field
	E[0]=0;
	E_int[0]=0;
    for (int i=0; i<NumberOfPoints;i++)
		{
        E[i+1]=E[i]+(rho1[i]-rho2[i]-rho3[i])*10000./(epsilon*area);
        E_int[i+1]=E_int[i]+E[i+1]; //the integrated E-field from 0 to i
        }

    //determine the offset, so that the potenial difference of the chamber stays constant (boundary condition)
    //the electric field in the center of the cell is approximated by the average of its two boundaries (linear interpolation)
    E_offset=E0-(E_int[NumberOfPoints-margin-1]-E_int[margin+1]+0.5*(E[margin]+E[NumberOfPoints-margin]))/(NumberOfPoints-2*margin);
	Emin=100000;
	Emax=0;


    //apply boundary condition, calculate Emax, Emin and calcualte the velocities
    for (int i=0; i<NumberOfPoints+1;i++)
		{
        #ifndef CONSTFIELD
            E[i]+=E_offset;
        #else
            E[i] = E0;
        #endif
//		Eabs = fabs(E[i]); not needed, because E is defined as positiv and flipping due to charges should not be possible
        if (i>margin && i<(NumberOfPoints-margin+1) && E[i]<Emin) {Emin=E[i];}
        if (i>margin && i<(NumberOfPoints-margin+1) && E[i]>Emax) {Emax=E[i];}

        //ion velocities, if constant mobilities are used, the function calls can be reduced
        #ifdef CONSTIONMOB
            v1[i]=(mobPos*E[i])*100./binSize;
            v2[i]=-(mobNeg*E[i])*100./binSize;
        #else
            v1[i]=(fillgas->mobility_pos(E[i])*E[i]*100.)/binSize;
            v2[i]=-(fillgas->mobility_neg(E[i])*E[i]*100.)/binSize;
        #endif

        v3[i]=-(fillgas->mobility_e(E[i])*100.*E[i])/binSize;		//mobility in m^2/Vs, E in V/cm, binsize in nm -> v in bin/ns
		}
	if (Emin<Eminall) { Eminall = Emin; }
    if (Emax>Emaxall) { Emaxall = Emax; }
    if (Emax > maxEMedium){
        cerr<<"abort, because Emax exeeds range of one of the active medium functions"<<endl;
        return 1;
    }
    //determine extreme velocities
    v3max = -(fillgas->mobility_e(Emax)*100.*Emax)/binSize;
    v2max = -(fillgas->mobility_neg(Emax)*Emax*100.)/binSize;
    v1max = (fillgas->mobility_pos(Emax)*Emax*100.)/binSize;
    vIonMax = ((-1.0*v2max) > v1max)? (-1.0*v2max) : v1max;

    v3min = -(fillgas->mobility_e(Emin)*100.*Emin)/binSize;
    v2min = -(fillgas->mobility_neg(Emin)*Emin*100.)/binSize;
    v1min = (fillgas->mobility_pos(Emin)*Emin*100.)/binSize;

    //determine extreme reaction rates
    GammaMax = fillgas->attachment(Emin);
    AlphaMax = chargeMax*fillgas->volume_recombination(Emin)/area/binSize/(1.6e-07)*(1.0e+06);


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
        timestep = (reactionRateLimit / GammaMax > speedLimit / (-1.0*v3max)) ? speedLimit / (-1.0*v3max) : reactionRateLimit / GammaMax;
        //limit the timestep during irradiation (relevant for very short pulses or very slow electrons (LIC))
        if (beamOnTime < 0)
        {
            if (timestep > pulseDuration/10.)
            {
                timestep = pulseDuration/10.;
            }
        }
        #if DEBUGLEVEL > 1
            log << "step: " << step << " timestep: " << timestep << " time: " << elapsedTime << " electron v: " << timestep*v3max << " gamma: " << GammaMax*timestep<< std::endl;
            log << "   v3min: "<<v3min<< " v3max: "<<v3max<<"; Emax: "<<Emax<<" Emin: "<<Emin<<endl;
        #endif
    }
	else
		{
        timestep = (reactionRateLimit/AlphaMax > speedLimit/vIonMax)? speedLimit/vIonMax : reactionRateLimit/AlphaMax;
        #if DEBUGLEVEL > 1
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

    #if DEBUGLEVEL > 1
        log<<endl;
    #endif


	}

//time loop complete


//final output
cout<<"***Results***"<<endl;
cout<<"remaining charge (I+/Q;I-/Q,e-/Q): "<<sum1/Q<<" "<<sum2/Q<<" "<<sum3/Q<<endl;
if (sum2<Q*1.0e-06 && sum3<Q*1.0e-06)
	{
    collectedCharge1+=sum1;
	sum1=0;
	}

if (sum1>Q*1.0e-04 || sum2>Q*1.0e-04 || sum3>Q*1.0e-04)
	{
    std::cerr<<"too many charges at the end of the loop --> increase calculation time!! (increase reactionRateLimit or NumberOfSteps) "<<std::endl;
	}
if (fabs(liberatedCharge/Q- 1) > 0.001)
{
    cerr<<"actually liberated charge deviates from target values, check calculation";
    collectedCharge1 = 0;
    collectedCharge2 = 0;
    collectedCharge3 = 0;
}
if (((fabs(v1min*timestep)>margin/2. || fabs(v1max*timestep)>margin/2.
    || fabs(v2min*timestep)>margin/2. || fabs(v2max*timestep)>margin/2.)
    && (electronCollectionTime > 0)) ||
    ((fabs(v3min*timestep)>margin/2. || fabs(v3max*timestep)>margin/2.)
    && electronCollectionTime < 0)) //check max v, if electrons are collected only v1 and v2 are relevant, otherwise v3
	{
    cerr<<"max velocity exceeds margin "<<margin<<" "
		<<v1min<<" "<<v2min<<" "<<v3min<<" "<<v1max<<" "<<v2max<<" "<<v3max<<std::endl;
	}


cout<<endl;
cout<<"liberated charge (Q): "<<liberatedCharge<<" pC"<<endl;
cout<<"collected charge (I+,I-,e-,I-+e-): "
    <<collectedCharge1<<" "
    <<collectedCharge2<<" "
    <<collectedCharge3<<" "
    <<collectedCharge2+collectedCharge3<<endl;
cout<<"recombined charge (I+&I-;I+&e-;Q-recombined): "
    <<recombinedCharge<<" "
    <<directlyRecombinedCharge<<" "
    <<liberatedCharge-recombinedCharge-directlyRecombinedCharge<<endl;
cout<<endl;

cout<<"collection times: "<<endl;
cout<<"total: "<<elapsedTime/1000.<<" µs; " <<step<<" steps"<<endl;
cout<<"electrons: "<<electronCollectionTime/1000.<<" µs, Beam: "<<beamOnTime/1000.<<" µs"<<endl;

if (collectedCharge2+collectedCharge3>0)
    help=collectedCharge3/liberatedCharge;
else
    help=0;
cout<<"fraction of free electrons: "<<help<<endl;

cout.precision(8);
if (collectedCharge1>0)
    help=liberatedCharge/(collectedCharge1);
else
    help=0;
cout<<"ksat (Q/collected +;Q/collected -): "<<help;

if (collectedCharge2+collectedCharge3>0)
    help=liberatedCharge/(collectedCharge2+collectedCharge3);
else
    help=0;
cout<<" "<<help<<endl;

if (PositivityWarning){
    cerr << "Waring: Negative charge density was cut off, conservation of total charge not guaranteed" << endl;
	}

cout<<"Emax: "<<Emaxall<<", Emin: "<<Eminall<<", total charge conservation (1?): "<<(recombinedCharge+directlyRecombinedCharge)/(liberatedCharge-collectedCharge1)<<endl;
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
		indep_variable = Q/1000.;
        log<<"result summary with total charge printed to "<<userInput.summaryFileName<<endl;
		break;
	case 't':
        indep_variable = pulseDuration;
        log<<"result summary with pulse duration printed to "<<userInput.summaryFileName<<endl;
		break;
    case 'N':
        indep_variable = ceil(pulseDuration/77.0);
        log<<"result summary with number of pulses printed to "<<userInput.summaryFileName<<endl;
        break;
	case 'V':
        indep_variable = userInput.voltage;
        log<<"result summary with voltage printed to "<<userInput.summaryFileName<<endl;
		break;
	}
	
    compact_out.precision(16);
    compact_out<<scientific;
    compact_out<<indep_variable<<", "<<liberatedCharge/(collectedCharge1)<<endl;
	
}
log<<"done"<<endl;
log<<endl;


return 0;
}
