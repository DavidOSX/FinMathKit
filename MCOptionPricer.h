
#pragma once
#include "Options.h"
#include "MonteCarlo.hpp"

namespace SiriusFM {



template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB>
class MCOptionPricer {
    
private:
    class OPPathEval {
    private :
        Option<AssetClassA, AssetClassB> const* const m_option;
        long m_P; // total paths
        double m_sum; //sum of payoff
        double m_sum2; // sum of  payoff^2
        double m_minPO; // min PayOff
        double m_maxPO; // m_maxPO
    public:
        OPPathEval(Option<AssetClassA, AssetClassB> const * a_option):
        m_option(a_option),
        m_P     (0),
        m_sum   (0),
        m_sum2  (0),
        m_minPO (INFINITY),
        m_maxPO (-INFINITY)
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
                m_minPO = std::min<double>(m_minPO, payOff);
                m_maxPO = std::max<double>(m_maxPO, payOff);
            }   m_P += a_PM;
        }
         
        std::tuple<double, double, double, double> GetStats() const { //double GetPx
            if(m_P < 2) throw std::runtime_error("empty OPPathEval");
            double px = m_sum / double(m_P);
            double err = m_sum2 - double(m_P) * px * px;
            //double err = (px != 0) ?  (fabs(px) / sqrt(var)) : sqrt(var);
            return std::make_tuple(px, sqrt(err), m_minPO, m_maxPO); 
        }
        //std::tuple GetStats
};
    Diffusion const* const      m_diff;
    AProvider                   m_airp;
    BProvider                   m_birp;
    MCEngine<Diffusion, 
             AProvider, 
             BProvider, 
             AssetClassA, 
             AssetClassB, 
             OPPathEval>        m_mce;
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
                   {}
                   
    double PX(Option<AssetClassA, AssetClassB> const* a_option,
              time_t a_t0,
              int a_tauMins,
              long a_P
             );
    inline double GetRateA(AssetClassA a_A, double a_ty) const { 
                    return m_airp.r(a_A, a_ty); 
    }

    inline double GetRateB(AssetClassB a_B, double a_ty) const { 
                    return m_birp.r(a_B, a_ty); 
    } 
};

};   
        
