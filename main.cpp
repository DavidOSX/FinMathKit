//#include "Diffusion.h"
//#include <time.h>
#include "MonteCarlo.hpp"

using namespace SiriusFM;
using namespace std;

int main(const int argc, const char* argv[]) {
    CcyE c1 = CcyE::USD;
    CcyE c2 = CcyE::CHF;
    //Diffusion_GBM diff = Diffusion_GBM(0, 0.1);
    //IRMode IRM = IRM::Const;
    if(argc < 8) {
        cerr << "not enough params\n";
        return 1;
    }
    
    IRProvider<IRMode::Const> irp = IRProvider<IRMode::Const>(argv[1]);
    printf("USD:%lf\n", irp.r(c1,0));
    printf("CHF:%lf\n", irp.r(c2,0));
    
    if(argc != 8) {
        cerr << "not enough params\n";
        return 1;
    }
    double mu = atof(argv[2]);
    double sigma = atof(argv[3]);
    double s0 = atof(argv[4]);
    long T_days = atol(argv[5]);
    int tau_min = atoi(argv[6]);
    long P = atol(argv[7]);
     Diffusion_GBM diff = Diffusion_GBM(mu, sigma);
    time_t t0 = time(nullptr);
    time_t T = t0 + T_days*86400;
    double Ty = double(T_days)/365.25;
    
    
    MCEngine<decltype(diff),decltype(irp), decltype(irp), decltype(c1), decltype(c2)> mce(20'000, 20'000);
    mce.Simulate<false>(t0, T, tau_min, s0, P, &diff, &irp, &irp, c1, c2/* false*/);
    //mce.printPaths();
    auto res  = mce.GetPaths();
    long L1 = get<0>(res);
    long P1 = get<1>(res);
    double const* paths = get<2>(res);
    double est =0, est2 = 0;
    int nvp = 0;
    for(int p = 0; p < P1; ++p) {
        double const* path = paths + p*L1;
        double st = path[L1-1];
        if(st <= 0.) continue;
        ++nvp;
        double rt = log(st/s0);
        est += rt; est2 += rt*rt;
    }
    assert(nvp > 1);
    est /= double(nvp);
    double VarST = (est2 - double(nvp)*est*est)/double(nvp-1);
    double sigma2E = VarST/Ty;
    double muE =(est+VarST/2.0)/Ty;
    cout << "mu="<<mu<<",muE="<<muE<<endl;
    cout << "sigma2="<<sigma*sigma<<",sigma2E="<<sigma2E<<endl;
    
    
    
    
    
    return 0;
}
