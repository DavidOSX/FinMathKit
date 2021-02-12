#include "MCOptionHedger.h"

namespace SiriusFM {
    
template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB>  
inline std::tuple<double, double, double, double> 
              MCOptionHedger <Diffusion, 
                              AProvider, 
                              BProvider, 
                              AssetClassA, 
                              AssetClassB>::SimulateHedging(Option<AssetClassA, AssetClassB> const*      a_option,
                                                 time_t                                     a_t0,
                                                 double                                     a_C0,
                                                 DeltaFunc const*                           a_deltaFunc,
                                                 double                                     a_deltaAcc,
                                                 int                                        a_tauMins,
                                                 long                                       a_P
                                              ) 
                              {
                                OHPathEval pathEval(a_option, &m_airp, &m_birp, a_deltaFunc, a_C0, a_deltaAcc);
                                m_mce.template Simulate<true>(a_t0, a_option -> ExpirTime(), a_tauMins, a_P, m_useTimerSeed, m_diff, &m_airp, &m_birp, a_option->assetA(), a_option->assetB(), &pathEval);
                                
                                //m_mce.printPaths(m_diff -> GetStartPoint());
                                
              
                                //auto res = pathEval.GetStats();
                                //double px = res.first; 
              
                                //px *= m_birp.DF(a_option->assetB(), a_t0, a_option -> ExpirTime());
                    
                                return pathEval.GetStats();
              
              
                              }
          
                            

};
