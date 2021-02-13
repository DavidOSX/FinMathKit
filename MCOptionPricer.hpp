#pragma once
#include "MCOptionPricer.h"

namespace SiriusFM {
    
template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB>  
inline double MCOptionPricer <Diffusion, 
                              AProvider, 
                              BProvider, 
                              AssetClassA, 
                              AssetClassB>::PX(Option<AssetClassA, AssetClassB> const*      a_option,
                                                 time_t                                     a_t0,
                                                 int                                        a_tauMins,
                                                 long                                       a_P
                                              ) 
                              {
                                OPPathEval pathEval(a_option);
                                m_mce.template Simulate<true>(a_t0, a_option -> ExpirTime(), a_tauMins, a_P, m_useTimerSeed, m_diff, &m_airp, &m_birp, a_option->assetA(), a_option->assetB(), &pathEval);
                                
                                m_mce.printPaths(m_diff -> GetStartPoint());
                                
              
                                auto res = pathEval.GetStats();
                                double px = std::get<0>(res); 
                                //std::cout << 
              
                                px *= m_birp.DF(a_option->assetB(), a_t0, a_option -> ExpirTime());
                    
                                return px;
              
              
                              }
          
                            

};
