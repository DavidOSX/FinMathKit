#include "IRProviderConst.h"
#include "Diffusion.h"
#include "Vanillas.h"
#include "GridOP_S3.hpp"
#include "BSM.hpp"


using namespace SiriusFM; 
using namespace std;

int main(const int argc, const char* argv[])
{ 
	if(argc != 10)
	{
		cerr << "params: \"file with rates\" sigma, S0, {Call/Put}, K, Tdays, NS, tauMins, isAmerican\n";
		return 1;
	}
	
	CcyE c2 = CcyE::USD; 
    CcyE c1 = CcyE::CHF;
    
	double sigma        = atof(argv[2]);
	double S0           = atof(argv[3]);
	const char* OptType = argv[4];
	double K            = atof(argv[5]);
	long   Tdays        = atol(argv[6]);
    int    NS           = atol(argv[7]);
	int    tauMins      = atoi(argv[8]);
    bool isAmerican     = bool(atoi(argv[9]));

	assert(sigma > 0 && S0 > 0 && K > 0 && Tdays > 0 && NS > 0 && tauMins > 0);


	Diffusion_GBM diff(0.0, sigma, S0);     

  
	time_t t0  = time(nullptr);           
	time_t T   = t0 + 86400 * Tdays;  

	OptionFX const* opt = nullptr;

    if (strcmp(OptType, "Call") == 0)     opt = new CallOptionFX(c1, c2, K, T, isAmerican);
    else
        if (strcmp(OptType, "Put")  == 0) opt = new PutOptionFX (c1, c2, K, T, isAmerican);
        
            else throw invalid_argument("Bad option type");

  
  GridNOP<decltype(diff), IRPConst, IRPConst, CcyE, CcyE> grid(argv[1], argv[1]);

  // Presto! 
  grid.RunBI<false>(opt, &diff, S0, t0, NS, tauMins);
  
  auto   res   = grid.GetPriceDeltaGamma();
  
  double px    = get<0>(res);
  double delta = get<1>(res);
  double gamma = get<2>(res);
  
  cout << " Px    = " << px 
     << "\n Delta = " << delta 
     << "\n Gamma = " << gamma << endl;


  delete opt;
	return 0;
}
