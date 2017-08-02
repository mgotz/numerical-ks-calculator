#ifndef ACTIVE_MEDIUM_FUNCTIODS
#define ACTIVE_MEDIUM_FUNCTIODS
/*
 * These classes provide simple containers for a few parameters defining a function (e.g., a simple line) and
 * an evaluation method, implementing that function.
 * After construction one can therefore call the function repeatedly using the same parameters.
 * There is one abstract base class: the evaluation_functiod, which specifies an interface to call the evaluation method,
 * independent of which specific implementation is used.
*/


#include <fstream>
#include <cmath>
#include <limits>

class evaluation_functiod{
	public:
        virtual double evaluate(double) = 0; //=0 makes this an abstract class, i.e. method is not defined and needs to be defined in the children
		virtual ~evaluation_functiod() = 0;
        double mMaxE;
};

inline evaluation_functiod::~evaluation_functiod() {}



class constant : public evaluation_functiod{
	public:
        virtual double evaluate(double ){ return(mIntercept);}
        constant(double pIntercept)
        {
            mIntercept=pIntercept;
            mMaxE = std::numeric_limits<double>::max();
        }
	private:
        double mIntercept;

	};

class linear : public evaluation_functiod{
	public:
        virtual double evaluate(double pEField){ return(mIntercept+mSlope*pEField);}
        linear(double pIntercept, double pSlope)
        {
            mIntercept=pIntercept;
            mSlope = pSlope;
            mMaxE = std::numeric_limits<double>::max();
        }
	private:
        double mIntercept;
        double mSlope;

	};

class exponential : public evaluation_functiod{
	public:
        virtual double evaluate(double pEField){ return(m0+m1*exp(m2*pEField));}
        exponential (double p0, double p1, double p2)
        {
            m0 = p0;
            m1 = p1;
            m2 = p2;
            mMaxE = std::numeric_limits<double>::max();
        }
	private:
        double m0;
        double m1;
        double m2;

	};

class pol4 : public evaluation_functiod{
	public:
        virtual double evaluate(double pEField){
        return(m0+m1*pEField+m2*pEField*pEField+m3*pEField*pEField*pEField+m4*pEField*pEField*pEField*pEField);
		}
        pol4 (double p0, double p1, double p2, double p3, double p4)
        {
            m0 = p0;
            m1 = p1;
            m2 = p2;
            m3 = p3;
            m4 = p4;
            mMaxE = std::numeric_limits<double>::max();
        }
	private:
        double m0;
        double m1;
        double m2;
        double m3;
        double m4;
};

class pol5 : public evaluation_functiod{
	public:
        virtual double evaluate(double pEField){
        return(m0+m1*pEField+m2*pEField*pEField+m3*pEField*pEField*pEField+m4*pEField*pEField*pEField*pEField+m5*pEField*pEField*pEField*pEField*pEField);
		}
        pol5 (double p0, double p1, double p2, double p3, double p4, double p5)
        {
            m0 = p0;
            m1 = p1;
            m2 = p2;
            m3 = p3;
            m4 = p4;
            m5 = p5;
            mMaxE = std::numeric_limits<double>::max();
        }
	private:
        double m0;
        double m1;
        double m2;
        double m3;
        double m4;
        double m5;

};

class spline_eval_only : public evaluation_functiod{
	public:
        virtual double evaluate(double pEField){
            i = (size_t)((pEField-mX0)/mGridspacing);
            return (mA[i]*pEField*pEField*pEField+mB[i]*pEField*pEField+mC[i]*pEField+mD[i]);
			}
		spline_eval_only(std::ifstream&);
		~spline_eval_only();

	private:
        double* mA;
        double* mB;
        double* mC;
        double* mD;
        double mGridspacing;
        double mX0;
        size_t i;
};


#endif
