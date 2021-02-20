#include "IRProviderConst.h"
#include "Diffusion.h"
#include "Vanillas.h"
#include "GridOP_S3.hpp"
#include "BSM.hpp"


using namespace SiriusFM; 
using namespace std;


constexpr double beta = 0.5; 

int main(const int argc, const char* argv[])
{ 
	if(argc != 11)
	{
		cerr << "params: \"file with rates\" sigma, S0, {Call/Put}, K, Tdays, NS, tauMins, isAmerican\n";
		return 1;
	}
	
	CcyE c2 = CcyE::USD; 
    CcyE c1 = CcyE::CHF;
    
	double sigma           = atof(argv[2]);
	double S0              = atof(argv[3]);
	const char* OptType    = argv[4];
	double K               = atof(argv[5]);
	long   Tdays           = atol(argv[6]);
    int    NS              = atol(argv[7]);
	int    tauMins         = atoi(argv[8]);
    char const* diffType   = argv[9];
    bool isAmerican        = bool(atoi(argv[10]));
    

	assert(sigma > 0 && S0 > 0 && K > 0 && Tdays > 0 && NS > 0 && tauMins > 0);

    Diffusion *diff = nullptr;
	if (strcmp(diffType, "GBM") == 0)  
        diff = new Diffusion_GBM(0., sigma, S0);
    else if (strcmp(diffType, "CEV") == 0) diff = new Diffusion_CEV(0., sigma, beta, S0);
        else throw invalid_argument("bad diffusion type");     

  
	time_t t0  = time(nullptr);           
	time_t T   = t0 + 86400 * Tdays;  

	OptionFX const* opt = nullptr;

    if (strcmp(OptType, "Call") == 0)     opt = new CallOptionFX(c1, c2, K, T, isAmerican);
    else
        if (strcmp(OptType, "Put")  == 0) opt = new PutOptionFX (c1, c2, K, T, isAmerican);
            else throw invalid_argument("bad option type");

  
  GridNOP<Diffusion, IRPConst, IRPConst, CcyE, CcyE> grid(argv[1], argv[1]);

  // Presto! 
  grid.RunBI<false>(opt, diff, S0, t0, NS, tauMins);
  
  auto   res   = grid.GetPriceDeltaGamma(opt, t0);
  
  double px    = get<0>(res);
  double delta = get<1>(res);
  double gamma = get<2>(res);
  
  cout << " Feynman-Kac:\n"
       << " Price = " << px 
     << "\n Delta = " << delta 
     << "\n Gamma = " << gamma << endl;
     
  grid.RunBI<true>(opt, diff, S0, t0, NS, tauMins);
  res = grid.GetPriceDeltaGamma(opt, t0);
  px = get<0>(res);
  cout << " Fokker-Planck:\n"
       << " Price = " << px << endl;
    
  delete diff;
  delete opt;
	return 0;
}
