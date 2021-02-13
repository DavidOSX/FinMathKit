#include "IRProviderConst.h"
#include "Diffusion.h"
#include "Vanillas.h"
#include "GridOP_S3.hpp"
#include "BSM.hpp"


using namespace SiriusFM;
using namespace std;

int main(const int argc, const char* argv[])
{
	if(argc != 9)
	{
		cerr << "params: \"file with rates\" sigma, S0, {Call/Put}, K, Tdays, NS, tauMins\n";
		return 1;
	}
	
	CcyE c1 = CcyE::USD;
    CcyE c2 = CcyE::CHF;
    
	double sigma        = atof(argv[2]);
	double S0           = atof(argv[3]);
	const char* OptType = argv[4];
	double K            = atof(argv[5]);
	long   Tdays        = atol(argv[6]);
    int    NS           = atol(argv[7]);
	int    tauMins      = atoi(argv[8]);

	assert(sigma > 0 && S0 > 0 && K > 0 && Tdays > 0 && NS > 0 && tauMins > 0);


	Diffusion_GBM diff(0.0, sigma, S0);     

  
	time_t t0  = time(nullptr);           
	time_t T   = t0 + 86400 * Tdays;  

	OptionFX const* opt = nullptr;

    if (strcmp(OptType, "Call") == 0)     opt = new EurCallOptionFX(c1, c2, K, T);
    else
        if (strcmp(OptType, "Put")  == 0) opt = new EurPutOptionFX (c1, c2, K, T);
        
            else throw invalid_argument("Bad option type");

  
  GridNOP<decltype(diff), IRPConst, IRPConst, CcyE, CcyE> grid(argv[1], argv[1]);

  // Presto! 
  grid.RunBI(opt, &diff, S0, t0, NS, tauMins);

  delete opt;
	return 0;
}
