#include "Time.h"
#include "GridOp_S3.h"

// ================================
//      ""
// Grid Pricer for Non-IR Options using 3-Point Stencils and 1-st order Runge-Kutta-Chebyshev time marshalling

namespace SiriusFM {



    template <typename Diffusion, 
            typename AProvider, 
            typename BProvider, 
            typename AssetClassA, 
            typename AssetClassB>
    inline void class GridNOP::RunBI(Option<AssetClassA, AssetClassB> const* a_option,
                                    Diffusion const*    a_diff,
                                    char const*        a_irsfileA,
                                    char const*        a_irsfileB,
                                    double              a_S0,
                                    time_t              a_t0,
                                    long                a_N,
                                    int                 a_tauMins,
                                    double              a_BFactor
                                    ) 
                {
                    double TTE = IntervalYearFrac(a_option->ExpirTime() - a_t0);
                    
                    if(TTE <= 0) throw std::invalid_argument("option has already expired");
                    
                    //fill in the timeline
                    long Mints = (a_option -> ExpirTime() - a_t0) / (a_tauMins * 60);
                    if(Mints <= 0 || TTE <= 0) throw std::invalid_argument("Option has already expired");
                    
                    double M = Mints + 1;
                    double tau = TTE / double(Mints);
                    
                    for (int j = 0; j < M; ++j) {
                        m_ts[j] = YearFrac(a_t0 + j * tauSecs);
        
                        double rA = m_airp.r(a-option -> assetA(), m_ts[j]); 
                        double rB = m_birp.r(a-option -> assetB(), m_ts[j]);
                        
                        double rateDiff = 
                    
                    
                    
                    
                }
    };
};
                    
          
