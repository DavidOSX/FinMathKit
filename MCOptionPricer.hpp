#include "MCOptionPricer.h"

namespace SiriusFM {
    
template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB,
          typename PathEvaluator>  
inline double MCOptionPricer <Diffusion, 
                              AProvider, 
                              BProvider, 
                              AssetClassA, 
                              AssetClassB,
                              PathEvaluator>::PX(Option const* a_option,
                                                 AssetClassA a_A,
                                                 AssetClassB a_B,
                                                 time_t a_t0,
                                                 int a_tauMins,
                                                 long a_P
                                                )  {
              OPPathEval pathEval(a_option);
              m_mce.template Simulate<true>(a_t0, a_option -> ExpirTime(), a_tauMins, a_P, m_useTimerSeed, m_diff, &m_airp, &m_birp, a_A, a_B, &pathEval);
              
              auto res = pathEval.GetPxStats();
              double px = res.first; 
              
              px *= m_birp.DF(a_A, a_t0, a_option -> ExpirTime());
              
              return px;
              
              
          }

};
