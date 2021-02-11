#pragma once
#include "Options.h"

namespace SiriusFM {

class OPPathEval {
    private :
        Option const* const m_option;
        long m_P; // total paths
        double m_sum; //sum of payoff
        double m_sum2; // sum of  payoff^2
    public:
        OPPathEval(Option const * a_option):
        m_option(a_option),
        m_P     (0),
        m_sum   (0),
        m_sum2  (0)
        //m_minPO (infty),
        //m_maxPO (-infty)
        { 
            assert(m_option != nullptr); 
            
        }
        void operator() (long a_L, 
                         long a_PM, 
                         double const* a_paths, 
                         double const* a_ts)
        {
            for (long p = 0; p < a_PM; ++p) {
                double const* path = a_paths + p * a_L;
                double payOff      = m_option -> payoff(a_L, a_ts, path);
                m_sum += payOff;
                m_sum2 += payOff * payOff;
                //m_minPO = std::min<double>(m_minPO, payOff);
                //m_maxPO = std::max<double>(m_maxPO, payOff);
            }
            m_P += a_PM;
        }
        
        std::pair<double, double> GetPxStats() const { //double GetPx
            if(m_P < 2) throw std::runtime_error("empty OPPathEval");
            double px = m_sum / double(m_P);
            double var = m_sum2 - double(m_P) * px * px;
            double err = (px != 0) ?  (fabs(px) / sqrt(var)) : sqrt(var);
            return std::make_pair(px, err);
        }
        //std::tuple GetStats
};

template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB,
          typename PathEvaluator>
class MCOptionPricer {
    
private:
    Diffusion const* const      m_diff;
    AProvider                   m_airp;
    BProvider                   m_birp;
    MCEngine<Diffusion, 
             AProvider, 
             BProvider, 
             AssetClassA, 
             AssetClassB, 
             PathEvaluator>     m_mce;
    bool                        m_useTimerSeed;
    
public:
    
    MCOptionPricer(Diffusion const* a_diff,
                   const char*      a_fileA,
                   const char*      a_fileB,
                   bool             a_useTimerSeed):
                   m_diff(a_diff),
                   m_airp(a_fileA),
                   m_birp(a_fileB),
                   m_mce(102'271, 4'096),
                   m_useTimerSeed(a_useTimerSeed)
                   { /* m_mce(102271,4096)*/ 
                       
                   }
                   
    double PX(Option const* a_option,
              AssetClassA a_A,
              AssetClassB a_B,
              time_t a_t0,
              int a_tauMins = 60,
              long a_P = 100'000
           );
};

};   
        
