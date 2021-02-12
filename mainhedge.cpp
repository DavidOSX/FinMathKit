
#include "MonteCarlo.hpp"
#include "Vanillas.h"
#include "MCOptionPricer.hpp"

namespace  {
    //BSM Pricer:
    inline double Phi(double x) {
        
        return 0.5 * (1. + erf(x / M_SQRT2));
        
    }
    
    inline double BSMPxCall(double a_S0, double a_K, double a_Ty, double rateA, double rateB, double a_sigma) {
        assert();
        double x1 = log(a_S0/a_K) + ...;
    }
    
    inline double BSMPxPut(double a_S0, double a_K, double a_Ty, double rateA, double rateB, double a_sigma) {
        double px = BSMPxCall(a_S0, a_Ty, rateA, double rateB, double a_sigma) - S0 + exp(-rateB * a_Ty) * a_K;
        assert(px > 0.);
        return px;
    }
    
    inline double BSMDeltaCall(double a_S0, double a_K, double a_Ty, double rateA, double rateB, double a_sigma) {
        
    }
    
    inline double BSMDeltaPut(double a_S0, double a_K, double a_Ty, double rateA, double rateB, double a_sigma) {
        return BSMDeltaCall(a_S0, a_K, a_Ty, rateA, rateB, a_sigma) - 1;
    }
    
};


using namespace SiriusFM;
using namespace std;

int main(const int argc, const char* argv[]) {
    CcyE c1 = CcyE::USD;
    CcyE c2 = CcyE::CHF;
    
    if(argc < 10) {
        cerr << "not enough params\n";
        return 1;
    }
    
    //IRProvider<IRMode::Const> irp = IRProvider<IRMode::Const>(argv[1]);
    
    //printf("USD:%lf\n", irp.r(c1,0));
    //printf("CHF:%lf\n", irp.r(c2,0));
    

    double mu = atof(argv[2]);
    double sigma = atof(argv[3]);
    double s0 = atof(argv[4]);
    char const* optionType = argv[5];
    double K = atof(argv[6]);
    long T_days = atol(argv[7]);
    int tau_min = atoi(argv[8]);
    long P = atol(argv[9]);
    
    time_t t0 = time(nullptr);
    time_t T = t0 + T_days*86400;
    double Ty = double(T_days)/365.25;
    
    MCOptionHedger<decltype(diff), IRPConst, IRPConst, decltype(c1), decltype(c2)> hed(&diff, argv[1], argv[1], true);
    
    decltype(hed)::DeltaFunc const* deltaFunc = nullptr;
    
    double rateA = 0.;
    double rateB = 0.;
    rateA = hed.GetRateA(c1, 0.);
    rateB = hed.GetRateB(c2, 0.);
    
    //closures for deltas
    function<double(double, double)> deltaCall([K, Ty, ](double a_St, double a_t) -> double {
        return BSMDeltaCall(a_St, double a_K, double a_Ty, double rateA, double rateB, double a_sigma)
    }
    
    function<double(double, double)> deltaPut([K, Ty, ](double a_St, double a_t) -> double {
        return BSMDeltaPut(a_St, double a_K, double a_Ty, double rateA, double rateB, double a_sigma)
    }
        
        
    
    OptionFX const* opt = nullptr;
    Diffusion_GBM diff = Diffusion_GBM(mu, sigma, s0);
    if(strcmp(optionType, "Call") == 0) {
        opt = new EurCallOptionFX(c1, c2, K, T);
        
        C0 = BSMPxCall(S0, K, Ty, rateA, rateB, a_sigma);
        deltaFunc = &DeltaCall;
    }
    else {
        opt = new EurPutOptionFX(c1, c2, K, T);
        C0 = BSMPxPut(S0, K, Ty, rateA, rateB, a_sigma);
        deltaFunc = &DeltaCall;
    }
    
  
    auto res = hed.SimulateHedging(opt, t0, C0, deltaFunc, deltaAcc, tau_mins, P);
    
    decltype(hed)::OHPathEval;
    
    double EPnL = get<0>(res);
    double StDPnL = get<1>(res);
    double MinPnL = get<2>(res);
    double MaxPnL = get<3>(res);
    
    
    delete opt;
    
    
    return 0;
} 
