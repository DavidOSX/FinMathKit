#include "MonteCarlo.h"
#include <cassert>
#include <random>


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
          typename AssetClassB>  
template<bool a_isRN>
inline void MCEngine <Diffusion,
                      AProvider, 
                      BProvider, 
                      AssetClassA, 
                      AssetClassB> ::
                                       Simulate(time_t           a_t0, 
                                               time_t           a_T, 
                                               int              a_tau_min,
                                               double           a_S0,
                                               long             a_P,
                                               Diffusion_GBM const* a_diff, 
                                               IRProvider<IRMode::Const> const* a_ap,
                                               IRProvider<IRMode::Const> const* a_bp,
                                               CcyE             a_A,
                                               CcyE             a_B
                                               /*bool             a_isRN*/
                                               )
                      {
                          assert(a_diff != nullptr && a_ap != nullptr && a_bp != nullptr && a_t0 <= a_T&& a_tau_min > 0);
                          
                          time_t Tsec = a_T - a_t0;
                          int tau_sec = a_tau_min*60;
                          
                          double tau = IntervalYearFrac(tau_sec);
                          long L = (Tsec % tau_sec == 0) ? Tsec/tau_sec : Tsec/tau_sec + 1;
                          
                          double y0 = YearFrac(a_t0);
                          //double yT = YearFrac(a_T);
                          
                          //dou
                          //(long) ceil(yT-y0)/tau + 1;
                          long P = 2*a_P;
                          double y = y0; 
                          //double Sp0 = a_S0, Sp1 = a_S0;
                          
                          if(L * P > m_MaxL * m_MaxP) std::invalid_argument("...");
                          
                          double stau = sqrt(tau);
                          
                          std::normal_distribution<> nd(0.0, 1.0);
                          std::mt19937_64 u;
                          
                          
                          //printf("%lf\n", tau)
                          
                          double tlast = (Tsec%tau_sec==0)? tau :IntervalYearFrac(Tsec-(L-1)*tau_sec);
                          //yT - y0 - double(L - 2)*tau;
                          assert(tlast <= tau && 0 < tlast);
                          double slast = sqrt(tlast);
                          ++L;
                          
                          //assert(slast <= stau && 0 < slast);
                          
                          for(long p = 0; p < a_P; ++p) {
                              double* path0 = m_paths + 2*p*L;
                              double* path1 = path0 + L;
                              path0[0] = a_S0;
                              path1[0] = a_S0;
                              double y = y0;
                              double Sp0 = a_S0;
                              double Sp1 = a_S0;
                              for(long l = 1; l < L; ++l) {
                                  double mu0 = 0.;
                                  double mu1 = 0.;
                                  if(a_isRN) {
                                      double delta_r = a_bp->r(a_B, y) - a_ap->r(a_A, y);
                                      mu0 = delta_r * Sp0;//+
                                      mu1 = delta_r * Sp1;
                                      
                                }
                              else {
                                  mu0 = a_diff -> mu(Sp0, y);
                                  mu1 = a_diff -> mu(Sp1, y);
                                  
                              }
                              double sigma0 = a_diff -> sigma(Sp0,y);
                              double sigma1 = a_diff -> sigma(Sp1,y);
                              double Z = nd(u);
                              double Sn0, Sn1;
                               if(l == L - 1) {
                                 Sn0 = Sp0 + mu0*tlast + sigma0*Z*slast;
                                 Sn1 = Sp1 + mu1*tlast + sigma1*Z*slast;
                                }
                                else {
                                    Sn0 = Sp0 + mu0*tau + sigma0*Z*stau;
                                    Sn1 = Sp1 + mu1*tau + sigma1*Z*stau;
                                }
                                path0[l] = Sn0;
                                path1[l] = Sn1;
                                Sp0 = Sn0;
                                Sp1 = Sn1;
                            }
                          } m_L = L; m_P = P;
                         // return std::make_pair(L, p);
                      }
}
