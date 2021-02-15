#include "Diffusion.h"
#include "IRProviderConst.h"
#include "Vanillas.h"
#include "MCOptionHedger.hpp"
#include "BSM.hpp"


using namespace SiriusFM;
using namespace std;

int main(const int argc, const char* argv[]) {
    
    
    if(argc < 12) {
        cerr << "not enough params\n" << endl;
        cerr << "params: \"name file with rates\" mu, sigma, S0, Call/Put, K, Tdays, deltaAcc, tau_mins, P, isAmerican\n";
        return 1;
    }
    
    CcyE c2 = CcyE::USD;
    CcyE c1 = CcyE::CHF;
    
    
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
    int tau_mins            = atoi(argv[9]);
    long P                  = atol(argv[10]);
    bool isAmerican         = bool(atoi(argv[11]));
    
    time_t t0               = time(nullptr);
    time_t T                = t0 + T_days * 86400;
    double Ty               = 1970. + double(T_days)/365.25;
    double TTE              = IntervalYearFrac(T - t0);
    
    OptionFX const* opt = nullptr;
    Diffusion_GBM diff = Diffusion_GBM(mu, sigma, s0);
    
    
    MCOptionHedger<decltype(diff), IRPConst, IRPConst, decltype(c1), decltype(c2)> hed(&diff, argv[1], argv[1], true);
    
    decltype(hed)::DeltaFunc const* deltaFunc = nullptr;
    
    double rateA = hed.GetRateA(c1, 0.);
    double rateB = hed.GetRateB(c2, 0.);
    
    //closures for deltas
    function<double(double, double)> deltaCall([K, Ty, rateA, rateB, sigma](double a_St, double a_tt) -> double {
        double curTTE = Ty - a_tt;  
        return BSMDeltaCall(a_St, K, curTTE, rateA, rateB, sigma);
    }
    );
    
    function<double(double, double)> deltaPut([K, Ty, rateA, rateB, sigma](double a_St, double a_tt) -> double {
        double curTTE = Ty - a_tt;
        return BSMDeltaPut(a_St, K, curTTE, rateA, rateB, sigma);
    }
    );
        
    double C0 = 0.;
    
    
    if(strcmp(optionType, "Call") == 0) {
        opt = new CallOptionFX(c1, c2, K, T, isAmerican);
        C0 = BSMPxCall(s0, K, TTE, rateA, rateB, sigma); 
        deltaFunc = &deltaCall;
    }
    else {
        opt = new PutOptionFX(c1, c2, K, T, isAmerican);
        C0 = BSMPxPut(s0, K, TTE, rateA, rateB, sigma); 
        deltaFunc = &deltaPut;
    } 
    
  
    auto res = hed.SimulateHedging(opt, t0, C0, deltaFunc, deltaAcc, tau_mins, P);
     
    //decltype(hed)::OHPathEval;
    
    double EPnL     = get<0>(res);
    double StDPnL   = get<1>(res);
    double MinPnL   = get<2>(res);
    double MaxPnL   = get<3>(res);
    
    cout << " E[PnL]   = "  << EPnL  
       << "\n StD[PnL] = " << StDPnL
       << "\n Min[PnL] = " << MinPnL 
       << "\n Max[PnL] = " << MaxPnL << endl;
    
    
    delete opt;
    
    
    return 0;
} 
