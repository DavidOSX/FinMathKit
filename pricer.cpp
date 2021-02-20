#include "Diffusion.h"
#include "IRProviderConst.h"

#include "Vanillas.h"
#include "MCOptionPricer.hpp"
#include "BSM.hpp"


using namespace SiriusFM;
using namespace std;

constexpr double beta = 0.5; 

int main(const int argc, const char* argv[]) {
    
    if(argc < 12) {
        cerr << "not enough params\n";
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
    int tau_min             = atoi(argv[8]);
    long P                  = atol(argv[9]);
    char const* diffType    = argv[10];
    bool isAmerican         = bool(atoi(argv[11]));
    
    
    time_t t0 = time(nullptr);
    time_t T = t0 + T_days*86400;
    double Ty = double(T_days)/365.25;
    double TTE   = IntervalYearFrac(T - t0);
    
    Diffusion *diff = nullptr;
    if (strcmp(diffType, "GBM") == 0)  
         diff = new Diffusion_GBM(mu, sigma, s0);
    else if (strcmp(diffType, "CEV") == 0)  diff = new Diffusion_CEV(mu, sigma, beta, s0);
        else throw invalid_argument("bad diffusion type");//if (strcmp(diffType, "Lipton") == 0) Lipton(mu, s);
    
    
    OptionFX const* opt;
    if(strcmp(optionType, "Call") == 0)     opt = new CallOptionFX(c1, c2, K, T, isAmerican);
    else if(strcmp(optionType, "Put") == 0) opt = new PutOptionFX(c1, c2, K, T, isAmerican);
        else throw invalid_argument("bad option type");
    
    
    MCOptionPricer<Diffusion, IRPConst, IRPConst, decltype(c1), decltype(c2)> pr(diff, argv[1], argv[1], true);
    double px = pr.PX(opt, t0, tau_min, P);
    
    double bsmpx = 0.;
    if(strcmp(optionType, "Call") == 0)  bsmpx = BSMPxCall(s0, K, TTE, pr.GetRateA(c1, 0), pr.GetRateB(c2, 0), sigma);
    else bsmpx = BSMPxPut(s0, K, TTE, pr.GetRateA(c1, 0), pr.GetRateB(c2, 0), sigma);
    
    

    //double VarST = (est2 - (double(nvp))*est*est)/double(nvp-1);
    //double sigma2E = VarST/Ty;
    //double muE = (est+VarST/2.0)/Ty;
    
    cout << "MonteCarlo Price = " << px << endl;
    if (strcmp(diffType, "GBM") == 0) cout << "BSM Price        = " << bsmpx << endl;
    
    delete diff;
    delete opt;
    
    
    return 0;
} 
