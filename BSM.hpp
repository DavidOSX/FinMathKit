#pragma once
#include <cmath>
#include <cassert>
#include <algorithm>

//namespace  {
    //BSM Pricer:
    inline double Phi(double x) {
        
        return 0.5 * (1. + erf(x / M_SQRT2));
        
    }
    
    double BSMPxCall(double a_S0, 
                     double a_K, 
                     double a_TTE, 
                     double rateA, 
                     double rateB, 
                     double a_sigma); 

    
    double BSMPxPut(double a_S0, 
                    double a_K, 
                    double a_TTE, 
                    double rateA, 
                    double rateB,
                    double a_sigma);
    
    inline double BSMDeltaCall(double a_S0, 
                               double a_K, 
                               double a_TTE, 
                               double rateA, 
                               double rateB, 
                               double a_sigma) 
    {
        assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
        if (a_TTE <= 0)
        return (a_S0 < a_K) ? 0 : (a_S0 > a_K) ? 1 : 0.5;

        double xd = a_sigma * sqrt(a_TTE);
        double x1 = (log(a_S0 / a_K) + (rateB - rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
        return Phi(x1);
    }
    
    inline double BSMDeltaPut(double a_S0, 
                              double a_K, 
                              double a_TTE, 
                              double rateA, 
                              double rateB, 
                              double a_sigma) 
    {
        return BSMDeltaCall(a_S0, a_K, a_TTE, rateA, rateB, a_sigma) - 1.;
    }
    
//};



