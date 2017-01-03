#include "active_medium_functiods.hpp"
#include <iostream>
#include <string>

using namespace std;

void solve_tridiag_inplace(size_t , const double*, const double*, double*, double*); //function to solve linear equation system

//contructor for splines with a parameter file, only evaluates the splines.
//the actual splines must be calculated and tabulated from outside.

spline_eval_only::spline_eval_only(ifstream& pParameterFile){

    //trash variable
    string tempString;
    //store some array sizes
    size_t NofSplinesGiven;
    size_t NofSplinesEval;
    //an additional counter (i already defined in header, because it's used in eval)
    size_t j = 0;
    //set file stream to throw exceptions
    pParameterFile.exceptions (fstream::failbit | fstream::badbit);
    //read first line, get number of splines (number of knots - 1) in the file
    // and mGridspacing (minimal distance between the splines)
    pParameterFile>>tempString>>NofSplinesGiven>>tempString>>mGridspacing;
    pParameterFile>>tempString;

    // read the knots
    double* knots = new double[NofSplinesGiven+1];
    for (i=0;i<NofSplinesGiven+1;i++){
        pParameterFile>>knots[i];
        }
    mX0 = knots[0];

    mMaxE = knots[NofSplinesGiven];
    //define number of splines for evaluation (create new equidistant grid)
    NofSplinesEval = (size_t)((knots[NofSplinesGiven]-knots[0])/mGridspacing);

    mA=new double[NofSplinesEval];
    mB=new double[NofSplinesEval];
    mC=new double[NofSplinesEval];
    mD=new double[NofSplinesEval];

    //read the spline parameters from the file and fill with duplicates where the uneven intervals require it

    for(i=0;i<NofSplinesEval;i++){
        if (knots[j] <= mGridspacing*i+mX0){
            pParameterFile>>mD[i]>>mC[i]>>mB[i]>>mA[i];
            j++;
        }
        else{
            mA[i]=mA[i-1];
            mB[i]=mB[i-1];
            mC[i]=mC[i-1];
            mD[i]=mD[i-1];
        }

    }

    //cleanup
    delete[] knots;

}

spline_eval_only::~spline_eval_only(){
delete[] mA;
delete[] mB;
delete[] mC;
delete[] mD;
}







void solve_tridiag_inplace(size_t N, const double a[], const double b[], double c[], double d[])
{

size_t i; //index of same type as size of array N
double temp; //help variable

/* solves a tri-diagonal linear equation system of the form:
b c 0   x   d
a b c . x = d
0 a c   x   d

the 3 diagnoals a b and c should be given in an array, d ist the right hand side. 
c and d are overwritten in the process and the solution is returned in d.
*/

c[0]=c[0]/b[0];
d[0]=d[0]/b[0];


for(i=1;i<N;i++){
	temp = b[i]-c[i-1]*a[i];
	c[i]=c[i]/temp;
	d[i]=(d[i]-d[i-1]*a[i])/temp;
	}


for (i=N-1;i-->0;){ //checks if i>0 then substracts 1 and then executes, thus first iteration i=N-2, last i=0
	d[i]=d[i]-c[i]*d[i+1];
	}

}


