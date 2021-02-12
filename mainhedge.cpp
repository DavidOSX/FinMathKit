
#include "MonteCarlo.hpp"
#include "Vanillas.h"
#include "MCOptionHedger.hpp"

namespace  {
    //BSM Pricer:
    inline double Phi(double x) {
        
        return 0.5 * (1. + erf(x / M_SQRT2));
        
    }
    
    inline double BSMPxCall(double a_S0, 
                            double a_K, 
                            double a_Ty, 
                            double rateA, 
                            double rateB, 
                            double a_sigma) 
    {
        assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
        if (a_TTE <= 0)
        // Return PayOff:
        return std::max<double>(a_S0 - a_K, 0);

        double xd = a_sigma * sqrt(a_TTE);
        double x1 = (log(a_S0 / a_K) + (a_rateB - a_rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
        double x2 = x1 - xd;
        double px = a_S0 * exp(-a_rateA * a_TTE) * Phi(x1) - a_K  * exp(-a_rateB * a_TTE) * Phi(x2);
        return px;
    }
    
    inline double BSMPxPut(double a_S0, 
                           double a_K, 
                           double a_Ty, 
                           double rateA, 
                           double rateB, 
                           double a_sigma) 
    {
        double px = BSMPxCall(a_S0, a_Ty, rateA, rateB, a_sigma) - S0 + exp(-rateB * a_Ty) * a_K;
        assert(px > 0.);
        return px;
    }
    
    inline double BSMDeltaCall(double a_S0, 
                               double a_K, 
                               double a_Ty, 
                               double rateA, 
                               double rateB, 
                               double a_sigma) 
    {
        assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
        if (a_TTE <= 0)
        return (a_S0 < a_K) ? 0 : (a_S0 > a_K) ? 1 : 0.5;

        double xd = a_sigma * sqrt(a_TTE);
        double x1 = (log(a_S0 / a_K) + (a_rateB - a_rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
        return Phi(x1);
    }
    
    inline double BSMDeltaPut(double a_S0, 
                              double a_K, 
                              double a_Ty, 
                              double rateA, 
                              double rateB, 
                              double a_sigma) 
    {
        return BSMDeltaCall(a_S0, a_K, a_Ty, rateA, rateB, a_sigma) - 1.;
    }
    
};


using namespace SiriusFM;
using namespace std;

int main(const int argc, const char* argv[]) {
    
    CcyE c1 = CcyE::USD;
    CcyE c2 = CcyE::CHF;
    
    if(argc < 11) {
        cerr << "not enough params\n" << endl;
        cerr << "params: \"name file with rates\" mu, sigma, S0, Call/Put, K, Tdays, deltaAcc, tau_mins, P\n";
        return 1;
    }
    
    //IRProvider<IRMode::Const> irp = IRProvider<IRMode::Const>(argv[1]);
    
    //printf("USD:%lf\n", irp.r(c1,0));
    //printf("CHF:%lf\n", irp.r(c2,0));
    

    double mu               = atof(argv[2]);
    double sigma            = atof(argv[3]);
    double s0               = atof(argv[4]);
    char const* optionType  = argv[5];
    double K                = atof(argv[6]);
    long T_days             = atol(argv[7]);
    double deltaAcc         = atof(argv[8]);
    int tau_min             = atoi(argv[9]);
    long P                  = atol(argv[10]);
    
    time_t t0               = time(nullptr);
    time_t T                = t0 + T_days*86400;
    double Ty               = 1970. + double(T_days)/365.25;
    double TTE              = IntervalYearFrac(T - t0);
    
    
    MCOptionHedger<decltype(diff), IRPConst, IRPConst, decltype(c1), decltype(c2)> hed(&diff, argv[1], argv[1], true);
    
    decltype(hed)::DeltaFunc const* deltaFunc = nullptr;
    
    double rateA = hed.GetRateA(c1, 0.);
    double rateB = hed.GetRateB(c2, 0.);
    
    //closures for deltas
    function<double(double, double)> deltaCall([K, Ty, rateA, rateB, sigma](double a_St, double a_t) -> double {
        double curTTE = Ty - a_t;  
        return BSMDeltaCall(a_St, a_K, curTTE, rateA, rateB, a_sigma);
    }
    );
    
    function<double(double, double)> deltaPut([K, Ty, rateA, rateB, sigma](double a_St, double a_t) -> double {
        double curTTE = Ty - a_t;
        return BSMDeltaPut(a_St, a_K, curTTE, rateA, rateB, a_sigma);
    }
    );
        
        
    
    OptionFX const* opt = nullptr;
    Diffusion_GBM diff = Diffusion_GBM(mu, sigma, s0);
    if(strcmp(optionType, "Call") == 0) {
        opt = new EurCallOptionFX(c1, c2, K, T);
        C0 = BSMPxCall(S0, K, TTE, rateA, rateB, a_sigma);
        deltaFunc = &deltaCall;
    }
    else {
        opt = new EurPutOptionFX(c1, c2, K, T);
        C0 = BSMPxPut(S0, K, TTE, rateA, rateB, a_sigma);
        deltaFunc = &deltaPut;
    }
    
  
    auto res = hed.SimulateHedging(opt, t0, C0, deltaFunc, deltaAcc, tau_mins, P);
    
    //decltype(hed)::OHPathEval;
    
    double EPnL = get<0>(res);
    double StDPnL = get<1>(res);
    double MinPnL = get<2>(res);
    double MaxPnL = get<3>(res);
    
    cout << " E[PnL] = "  << EPnL  << "\n StD[PnL] = " << StDPnL
       << "\n Min[PnL] = " << MinPnL << "\n Max[PnL] = " << MaxPnL << endl;
    
    
    delete opt;
    
    
    return 0;
} 
