
#pragma once
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
    inline void GridNOP<Diffusion, 
                AProvider, 
                BProvider,
                AssetClassA, 
                AssetClassB>::RunBI(Option<AssetClassA, AssetClassB> const* a_option,
                                    Diffusion const*    a_diff,
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
                    
                    long M = Mints + 1;
                    double tau = TTE / double(Mints);
                    double integAB = 0.;
                    m_ES[0] = a_S0;
                    m_VarS[0] = 0;
                    
                    for (int j = 0; j < M; ++j) {
                        m_ts[j] = YearFrac(a_t0 + j * (a_tauMins * 60));
        
                        double rA = m_airp.r(a_option -> assetA(), m_ts[j]); 
                        double rB = m_birp.r(a_option -> assetB(), m_ts[j]);
                        
                        double rateDiff = std::max<double>(rB-rA, 0);
                    
                    
                        if(j < M - 1) {
                            integAB += rateDiff * tau;
                            m_ES[j+1] = a_S0 * exp(integAB);
                            double sigma = a_diff->sigma(m_ES[j], m_ts[j]);
                            m_VarS[j+1]  = m_VarS[j] + sigma * sigma * tau;
                        }
                    }
                    // Upper Bound for S (the Lower Bound is 0):
                    double B = m_ES[M-1] + a_BFactor * sqrt(m_VarS[M-1]);

                    // Generate the S-Line:
                    double h = B / double(a_N-1);
                    
                    
                }
};

                    
          
