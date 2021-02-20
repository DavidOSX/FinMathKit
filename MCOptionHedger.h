
#pragma once
#include <functional>
#include "Options.h"
#include "MonteCarlo.hpp"

namespace SiriusFM {


          template <typename Diffusion, 
                    typename AProvider, 
                    typename BProvider, 
                    typename AssetClassA, 
                    typename AssetClassB>
          class MCOptionHedger {
          public:
                    using DeltaFunc = std::function<double(double, double)>;
    
          private:
          class OHPathEval {
        
          
          private :
                    Option<AssetClassA, AssetClassB> const* const m_option;
                    AProvider const*                              m_airp;
                    BProvider const*                              m_birp;
                    double*                                       m_rsA;
                    double*                                       m_rsB;
                    double const                                  m_C0; //init option premimum
                    DeltaFunc const* const                        m_deltaFunc;     
                    double const                                  m_deltaAcc;
                    long                                          m_P; // total paths evaluated
                    double                                        m_sumPnL; //sum of residual P&Ls
                    double                                        m_sumPnL2; // sum of P&Ls^2
                    double                                        m_minPnL; // min PayOff
                    double                                        m_maxPnL; // m_maxPO
          public:
                    OHPathEval(Option<AssetClassA, AssetClassB> const * a_option, 
                               AProvider const* a_airp,
                               BProvider const* a_birp,
                               DeltaFunc const* a_deltaFunc,
                               double           a_C0,
                               double           a_deltaAcc) 
                            : m_option    (a_option),
                              m_airp      (a_airp),
                              m_birp      (a_birp),
                              m_rsA       (nullptr),
                              m_rsB       (nullptr),
                              m_C0        (a_C0),
                              m_deltaFunc (a_deltaFunc),
                              m_deltaAcc  (a_deltaAcc),
                              m_P         (0),
                              m_sumPnL    (0),
                              m_sumPnL2   (0),
                              m_minPnL    (INFINITY),
                              m_maxPnL    (-INFINITY)
                    { 
                              assert(m_option != nullptr && 
                              m_deltaFunc != nullptr && 
                              m_deltaAcc >= 0.0 && 
                              m_airp != nullptr &&
                              m_birp != nullptr
                              ); 
            
                    }
                    void operator() (long a_L, 
                                     long a_PM, 
                                     double const* a_paths, 
                                     double const* a_ts)
                            {
                              if (m_rsA == nullptr) {
                                        m_rsA = new double[a_L];
                                        for (long l = 0; l < a_L; ++l)   m_rsA[l] = m_airp -> r(m_option->assetA(), a_ts[l]);
                              }
                              if (m_rsB == nullptr) {
                                        m_rsB = new double[a_L];
                                        for (long l = 0; l < a_L; ++l)   m_rsB[l] = m_birp -> r(m_option->assetB(), a_ts[l]);
                              }
                              for (long p = 0; p < a_PM; ++p) {
                
                                        double const* path = a_paths + p * a_L;
                                        
                                        // delta-hedging along this path:
                                        double M     = - m_C0;
                                        double delta    = 0.0; //cur delta
                                        for (long l = 0; l < a_L; ++l) {
                    
                                                  double St   = path[l];
                                                  double t    = a_ts[l];
                    
                                                  //manage money account
                                             if(l > 0) {
                                                double tau = t - a_ts[l - 1];
                                                double Sp = path[l - 1];
                                                M += M * tau * m_rsB[l - 1];
                        
                                                M += Sp * tau * m_rsA[l - 1];
                        
                                             }
                                            if(l < a_L - 1) {
                                            double deltaN = (*m_deltaFunc)(St, t);
                                            // also deltaN changes sign as we long the option
                                            deltaN = - round(deltaN / m_deltaAcc) * m_deltaAcc;
                                                if (deltaN != delta) {
                                                  // re-hedge
                                                  M -= (deltaN - delta) * St;
                                                  delta = deltaN;
                                                }
                                            }
                                      }
                
                              double PnL = M + delta * path[a_L - 1] + m_option -> payoff(a_L, a_ts, path);
            
                
                              m_sumPnL += PnL;
                              m_sumPnL2 += PnL * PnL;
                              m_minPnL = std::min<double>(m_minPnL, PnL);
                              m_maxPnL = std::max<double>(m_maxPnL, PnL);
                              }  m_P += a_PM;
                            }
        
                    std::tuple<double, double, double, double> GetStats() {
                              if (m_P < 2)   throw std::runtime_error("empty OHPathEval");
                              double mean =  m_sumPnL  / double(m_P);
                              double var  = (m_sumPnL2 - double(m_P) * mean * mean) / double(m_P - 1);
                              assert(var >= 0);
                              return std::make_tuple(mean, sqrt(var), m_minPnL, m_maxPnL);
                    }
        
                    ~OHPathEval() { 
                              delete[] (m_rsA); 
                              delete[] (m_rsB);
                              m_rsA = nullptr;
                              m_rsB = nullptr;
                    }
          };
          Diffusion const* const      m_diff;
          AProvider                   m_airp;
          BProvider                   m_birp;
          MCEngine<Diffusion, 
                    AProvider, 
                    BProvider, 
                    AssetClassA, 
                    AssetClassB, 
                    OHPathEval>       m_mce;
          bool                        m_useTimerSeed;
    
      public:  
    
          MCOptionHedger(Diffusion const* a_diff,
                   const char*      a_fileA,
                   const char*      a_fileB,
                   bool             a_useTimerSeed):
                   m_diff(a_diff),
                   m_airp(a_fileA),
                   m_birp(a_fileB),
                   m_mce(102'271, 4'096),
                   m_useTimerSeed(a_useTimerSeed)
                   {}
                   
          std::tuple<double, double, double, double> 
          SimulateHedging(Option<AssetClassA, AssetClassB> const* a_option,
                          time_t                                  a_t0,
                          double                                  a_C0,
                          DeltaFunc const*                        a_deltaFunc,
                          double                                  a_deltaAcc,
                          int                                     a_tauMins,
                          long                                    a_P);
                         
           inline double GetRateA(AssetClassA a_A, double a_ty) const { 
                    return m_airp.r(a_A, a_ty); 
           }

           inline double GetRateB(AssetClassB a_B, double a_ty) const { 
                    return m_birp.r(a_B, a_ty); 
          }
      };

};   
        
