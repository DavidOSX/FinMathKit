#include "BSM.hpp"


//using namespace SiriusFM;

double BSMPxCall(double a_S0, 
                double a_K, 
                double a_TTE, 
                double rateA, 
                double rateB, 
                double a_sigma) 
    {
        assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
        if (a_TTE <= 0)
        // Return PayOff:
        return std::max<double>(a_S0 - a_K, 0);

        double xd = a_sigma * sqrt(a_TTE);
        double x1 = (log(a_S0 / a_K) + (rateB - rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
        double x2 = x1 - xd;
        double px = a_S0 * exp(-rateA * a_TTE) * Phi(x1) - a_K  * exp(-rateB * a_TTE) * Phi(x2);
        return px;
    }
    
double BSMPxPut(double a_S0, 
                double a_K, 
                double a_TTE, 
                double rateA, 
                double rateB, 
                double a_sigma)
    {
        double px = BSMPxCall(a_S0, a_K, a_TTE, rateA, rateB, a_sigma) - a_S0 + exp(-rateB * a_TTE) * a_K;
        assert(px > 0.);
        return px;
    }
