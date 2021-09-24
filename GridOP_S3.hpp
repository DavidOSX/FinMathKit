
#pragma once
#include "GridOp_S3.h"
#include <tuple>
//#include <omp.h>


// ================================
//      ""
// Grid Pricer for Non-IR Options using 3-Point Stencils and 1-st order Runge-Kutta-Chebyshev time marshalling

namespace SiriusFM {



    template<typename Diffusion, 
             typename AProvider, 
             typename BProvider, 
             typename AssetClassA, 
             typename AssetClassB>
    template<bool isFwd> 
    inline void GridNOP<Diffusion, 
                AProvider, 
                BProvider,
                AssetClassA, 
                AssetClassB>::RunBI(Option<AssetClassA, AssetClassB> const* a_option, 
                                    Diffusion const*    a_diff,
                                    double              a_S0,
                                    time_t              a_t0,
                                    long                a_Nints,
                                    int                 a_tauMins,
                                    double              a_BFactor
                                    ) 
                {  
                    double TTE = IntervalYearFrac(a_option->ExpirTime() - a_t0);
                    
                    if(TTE <= 0) throw std::invalid_argument("option has already expired");
                    
                    //fill in the timeline
                    long Mints = (a_option -> ExpirTime() - a_t0) / (a_tauMins * 60);
                    if(Mints <= 0 || TTE <= 0) throw std::invalid_argument("Option has already expired");
                    
                    m_isFwd = isFwd;
                    m_M = Mints + 1;
                    double tau = TTE / double(Mints);
                    double integAB = 0.;
                    m_ES[0] = a_S0;
                    m_VarS[0] = 0;
                    
                    for (int j = 0; j < m_M; ++j) {
                        m_ts[j] = YearFrac(a_t0 + j * (a_tauMins * 60));
        
                        double rA = m_airp.r(a_option -> assetA(), m_ts[j]); 
                        double rB = m_birp.r(a_option -> assetB(), m_ts[j]);
                        
                        double rateDiff = std::max<double>(rB - rA, 0);
                    
                    
                        if(j < m_M - 1) {
                            integAB += rateDiff * tau;
                            m_ES[j + 1] = a_S0 * exp(integAB);
                            double sigma = a_diff->sigma(m_ES[j], m_ts[j]);
                            m_VarS[j + 1]  = m_VarS[j] + sigma * sigma * tau;
                        }
                    }
                    // Upper Bound for S (the Lower Bound is 0):
                    double B = m_ES[m_M - 1] + a_BFactor * sqrt(m_VarS[m_M - 1]);

                    // Generate the S-Line:
                    
                    double h = B / double(a_Nints-1);
                    
                    m_i0 = int(round(a_S0 / h)); 
                    h = a_S0 / double(m_i0);
                    if(!std::isfinite(h)) 
                        throw std::invalid_argument("S0 too small");
                    B = h * double(a_Nints);
                    m_N = a_Nints + 1;
                    if(m_N > m_maxN) 
                        throw std::invalid_argument("N too large");
                    
                    double* payOff = (!isFwd) ? m_grid + (m_M - 1) * m_N : nullptr;
                    
                    for(int i = 0; i < m_N; ++i) {
                        m_S[i] = double(i) * h; // again , low bounds is 0
                    
                    if (!isFwd) payOff[i] = a_option -> payoff(1, m_ts + (m_M - 1), m_S + i); 
                    if(isFwd) m_grid[i] = 0.; 
                    }
                    if(isFwd) m_grid[m_i0] = 1./h;
                    
                    bool isNeumann = false;
                    double UBC = 0.;
                    double fa = isFwd ? 0 : payOff[0];
                    if(!isFwd) {
                        isNeumann = (payOff[m_N - 1] != 0);
                        UBC = isNeumann ? (payOff[m_N - 1] - payOff[m_N - 2])  : 0;
                    }
                    
                    for (int j = 0; j < m_M - 1; ++j)  m_grid[j * m_N] = fa;
                        
                    
                    double D2 = 2 * h * h; //Denom in the diffusive term
                    
                    
                    
                    for (int j = isFwd ? 0 : m_M-1; isFwd ? (j <= m_M-2) : (j >= 1); j += (isFwd ? 1 : -1)) {
                        
                        double const* fj = m_grid + j * m_N;
                        double* fj1 = const_cast<double*>(isFwd ? (fj + m_N) : (fj - m_N));
                        
                        double tj = m_ts[j];
                        double rateAj = m_airp.r(a_option -> assetA(), tj);
                        double rateBj = m_birp.r(a_option -> assetB(), tj);
                        double C1 = (rateBj - rateAj) / (2 * h);
                        fj1[0] = fa;
                        //#pragma acc loop copyin(fj[0:N]) copyout()
                        //#pragma omp parallel for
                        for (int i = 1; i <= m_N - 2; ++i) {
                            double Si       = m_S[i];
                            double fjim1    = fj[i - 1];
                            double fji      = fj[i];
                            
                            double fjip1    = fj[i + 1];
                            double sigma    = a_diff -> sigma(Si, tj);
                            
                            double DfDt = 0.;
                            
                            if (isFwd) {
                                double SiM    = m_S[i - 1];
                                double SiP    = m_S[i + 1];
                                double sigmaM = a_diff -> sigma(SiM, tj);
                                double sigmaP = a_diff -> sigma(SiP, tj);

                                DfDt = - C1 * (SiP * fjip1 - SiM * fjim1)
                                + (sigmaP * sigmaP * fjip1 - 2 * sigma * sigma * fji 
                                +  sigmaM * sigmaM * fjim1) / D2;
                            }
                            else DfDt = rateBj * fji - C1 * Si * (fjip1 - fjim1) - sigma * sigma/D2 *(fjip1 - 2 * fji + fjim1);
                            
                            if(!isFwd) fj1[i] = fji - tau * DfDt;
                            else fj1[i] = fji + tau * DfDt;
                        }
                        
                        fj1[m_N - 1] = (!isFwd && isNeumann) ? (fj1[m_N - 2] + UBC) : UBC;
                        
                        if (a_option -> isAmerican() && !isFwd) {
                            for (int i = 0; i < m_N; ++i) {
                                
                                double  intrVal = a_option -> payoff(1, &tj, m_S + i);

                                fj1[i] = std::max<double>(fj1[i],  intrVal);
                            }
                        }
                    }
                    
                }
                
                
    template <typename Diffusion, 
             typename AProvider, 
             typename BProvider, 
             typename AssetClassA, 
             typename AssetClassB>
             inline std::tuple<double , double, double> 
             GridNOP<Diffusion, 
                    AProvider, 
                    BProvider,
                    AssetClassA, 
                    AssetClassB>::GetPriceDeltaGamma(Option<AssetClassA, AssetClassB> const* a_option, time_t a_t0) const 
                    {
                        if (m_M == 0 || m_N == 0)
                        throw std::runtime_error("RunBI first!");

                        assert(0 <= m_i0 && m_i0 < m_N);
                        
                        double h     = m_S[1] - m_S[0];
                        double px    = 0.;
                        double delta = 0.;
                        double gamma = 0.;
                        if(!m_isFwd) {
                        
                        px  = m_grid[m_i0];    // j=0
                        
                        if (0 < m_i0 && m_i0 <= m_N-2) {
                            delta = (m_grid[m_i0 + 1] - m_grid[m_i0-1]) / (2*h);
                            gamma = (m_grid[m_i0 + 1] - 2 * m_grid[m_i0] + m_grid[m_i0 - 1]) / (h * h);
                        } else if (m_i0 == 0) delta = (m_grid[1]   - m_grid[0])   / h;    // gamma remains 0
                            else {
                                assert(m_i0  == m_N-1);
                                delta = (m_grid[m_N-1] - m_grid[m_N-2]) / h; // gamma remains 0
                            }
                            
                        }
                        else {
                            for (int i = 0; i < m_N; ++i) {
                                px += (a_option -> payoff(1, m_ts, m_S + i)) * m_grid[(m_M - 1) * m_N + i] * h;
                            }
                            px *= m_birp.DF(a_option->assetB(), a_t0, a_option -> ExpirTime());
                        }
                        return std::make_tuple(px, delta, gamma);
                    }
                        
};

                    
          
