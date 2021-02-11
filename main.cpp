
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
    
    
    Option const* opt;
    Diffusion_GBM diff = Diffusion_GBM(mu, sigma, s0);
    if(strcmp(optionType, "Call") == 0) opt = new EurCallOption(K, T);
    else opt = new EurPutOption(K, T);
    
    
    MCOptionPricer<decltype(diff), IRPConst, IRPConst, decltype(c1), decltype(c2), OPPathEval> pr(&diff, argv[1], argv[1], true);
    double px = pr.PX(opt, c1, c2, t0);
    //double err = res.second;
    
    //OPPathEval pathEval(opt);
    
    
    //MCEngine<decltype(diff), decltype(irp), decltype(irp), decltype(c1), decltype(c2), decltype(pathEval)> mce(20'000, 20'000);
    //mce.Simulate<true>(t0, T, tau_min, P, true, &diff, &irp, &irp, c1, c2, &pathEval);
    //mce.printPaths();
    //auto res  = mce.GetPaths();
    //long L1 = get<0>(res);
    //long P1 = get<1>(res);
    //double const* paths = get<2>(res);
    //double est = 0., est2 = 0.;
    //int nvp = 0;
    /*for(int p = 0; p < P1; ++p) {
        double const* path = paths + p*L1;
        double st = path[L1-1];
        if(st <= 0.) continue;
        ++nvp;
        double rt = log(st/s0);
        est += rt; est2 += rt*rt;
        //priceOp += opt -> payoff(L1, nullptr, path);
       
    }*/
    //assert(nvp > 1);
    //est /= double(nvp);
    //priceOp /= double(nvp);
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
