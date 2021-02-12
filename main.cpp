
#include "MonteCarlo.hpp"
#include "Vanillas.h"
#include "MCOptionPricer.hpp"

//namespace SiriusFM {
    //path evaluator for option pricing:
    
//};

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
    
    
    OptionFX const* opt;
    Diffusion_GBM diff = Diffusion_GBM(mu, sigma, s0);
    if(strcmp(optionType, "Call") == 0) opt = new EurCallOptionFX(c1, c2, K, T);
    else opt = new EurPutOptionFX(c1, c2, K, T);
    
    
    MCOptionPricer<decltype(diff), IRPConst, IRPConst, decltype(c1), decltype(c2)> pr(&diff, argv[1], argv[1], true);
    double px = pr.PX(opt, t0, tau_min, P);
    //double err = res.second;
    
    //OPPathEval pathEval(opt);
    
    
    
    //priceOp *= exp((-1)*irp.r(c1,0)*Ty);
    //double VarST = (est2 - (double(nvp))*est*est)/double(nvp-1);
    //double sigma2E = VarST/Ty;
    //double muE = (est+VarST/2.0)/Ty;
    
    //cout << "mu = " << mu << ", muE = " << muE << endl;
    //cout << "sigma2 = " << sigma*sigma << ", sigma2E = " << sigma2E << endl;
    
    //if(strcmp(optionType, "Call") == 0) cout << "price(Call) = " << priceOp << endl;
    //else cout << "price(Put) = " << priceOp << endl; 
    //auto res = pathEval.GetPxStats();
    
    
    cout << "Price = " << px << endl;
    //cout << "Err = " << err << endl;
    
    delete opt;
    
    
    return 0;
} 
