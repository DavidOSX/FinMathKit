
#include <cassert>
#include <random>
#include "MonteCarlo.h"


namespace SiriusFM
{
    
inline double IntervalYearFrac(time_t t) {
    constexpr double SecY = 365.25*86400.;
    return (double) t/SecY;
}

inline double YearFrac(time_t t) {
    return 1970.+ IntervalYearFrac(t);
}



template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB
          typename PathEvaluator>  
template<bool a_isRN>
inline void MCEngine <Diffusion,
                      AProvider, 
                      BProvider, 
                      AssetClassA, 
                      AssetClassB
                      PathEvaluator> :: Simulate(time_t                           a_t0, 
                                               time_t                           a_T, 
                                               int                              a_tau_min,
                                               /*double                           a_S0,*/
                                               long                             a_P,
                                               bool                             a_useTimerSeed,
                                               Diffusion_GBM const*             a_diff, 
                                               IRProvider<IRMode::Const> const* a_ap,
                                               IRProvider<IRMode::Const> const* a_bp,
                                               CcyE                             a_A,
                                               CcyE                             a_B
                                               )
                      {
                          assert(a_diff != nullptr &&
                                 a_ap != nullptr && 
                                 a_bp != nullptr && 
                                 a_t0 <= a_T && 
                                 a_tau_min > 0 &&
                                 a_pathEval != nullptr);
                          
                          time_t Tsec = a_T - a_t0;
                          int tau_sec = a_tau_min*60;
                          
                          double tau = IntervalYearFrac(tau_sec);
                          long L = (Tsec % tau_sec == 0) ? Tsec/tau_sec : Tsec/tau_sec + 1;
                          
                          double y0 = YearFrac(a_t0);
                        
                          long P = 2*a_P;
                          double y = y0; 
                          
                          if(L > m_MaxL) std::invalid_argument("too many steps");
                          
                          double stau = sqrt(tau);
                          
                          std::normal_distribution<> nd(0.0, 1.0);
                          std::mt19937_64 u(a_useTimerSeed ? time(nullptr) : 0);
                          
                          long PM = (m_MaxL * m_MaxPM) / L;//how many paths we can store in memory 
                          
                          if(P % 2 != 0) --P;
                          assert(PM > 0 && P % 2 == 0);
                          long PMh = PM / 2;
                          // number of outer P iters:
                          long PI = (P % PM == 0) ? P/PM : (P/PM + 1);
                          
                          //adlust P
                          P = PI * PM;
                          
                          m_ts[L - 1] = m_ts[L-2]+tlast;
                          
                          for(long l = 0; l < L - 1;  ++l) m_ts[l] = y0 + double(l) * tau;
                          
                          
                          
                          
                          double tlast = (Tsec%tau_sec==0)? tau :IntervalYearFrac(Tsec-(L-1)*tau_sec);
                          //yT - y0 - double(L - 2)*tau;
                          assert(tlast <= tau && 0 < tlast);
                          double slast = sqrt(tlast);
                          ++L;
                          for(long i = 0; i < PI; ++i) {
                          
                            for(long p = 0; p < PMh; ++p) {
                              double* path0 = m_paths + 2 * p * L;
                              double* path1 = path0 + L;
                              path0[0] = a_diff -> GetStartPoint();
                              path1[0] = a_diff -> GetStartPoint();
                              //double y = y0;
                              double Sp0 = a_diff -> GetStartPoint();
                              double Sp1 = a_diff -> GetStartPoint();
                              for(long l = 1; l < L; ++l) {
                                  double mu0 = 0.;
                                  double mu1 = 0.;
                                  double y = m_ts[l - 1];
                                  if(a_isRN) {
                                      double delta_r = a_bp->r(a_B, y) - a_ap->r(a_A, y);
                                      mu0 = delta_r * Sp0;//+
                                      mu1 = delta_r * Sp1;
                                      
                                    } else {
                                        mu0 = a_diff -> mu(Sp0, y);
                                        mu1 = a_diff -> mu(Sp1, y);
                                    }
                                double sigma0 = a_diff -> sigma(Sp0,y);
                                double sigma1 = a_diff -> sigma(Sp1,y);
                                double Z = nd(u);
                                double Sn0, Sn1;
                                if(l == L - 1) {
                                    Sn0 = Sp0 + mu0*tlast + sigma0*Z*slast;
                                    Sn1 = Sp1 + mu1*tlast - sigma1*Z*slast;
                                } else {
                                    Sn0 = Sp0 + mu0*tau + sigma0*Z*stau;
                                    Sn1 = Sp1 + mu1*tau - sigma1*Z*stau;
                                }
                                path0[l] = Sn0;
                                path1[l] = Sn1;
                                Sp0 = Sn0;
                                Sp1 = Sn1;
                            }
                          }  (*a_pathEval)(L, PM, m_paths, m_ts) //evaluate in-memory paths
                        } //m_L = L; m_P = P;
                      }
};
