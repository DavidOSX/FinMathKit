#include "Diffusion.h"

using namespace SiriusFM;

int main(const int argc, const char* argv[]) {
    CcyE c1 = CcyE::USD;
    CcyE c2 = CcyE::CHF;
    //IRMode IRM = IRM::Const;
    
    IRProvider<IRMode::Const> irp = IRProvider<IRMode::Const>(argv[1]);
    printf("USD:%lf\n", irp.r(c1,0));
    printf("CHF:%lf\n", irp.r(c2,0));
    
    return 0;
}
