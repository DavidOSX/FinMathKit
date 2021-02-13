#pragma once 
#include "Time.h"
#include "Options.h"


// ================================
//      ""
// Grid Pricer for Non-IR Options using 3-Point Stencils and 1-st order Runge-Kutta-Chebyshev time marshalling

namespace SiriusFM {



    template<typename Diffusion, 
            typename AProvider, 
            typename BProvider, 
            typename AssetClassA, 
            typename AssetClassB>
    class GridNOP {
    private:
        AProvider  m_airp;
        BProvider  m_birp;
        long        m_maxM;
        long        m_maxN;
        double* const m_grid;
        double* const m_ts;
        double* const m_S;
        double* const m_ES;
        double* const m_VarS;
        
    public:
        GridNOP(char const * a_fileA, char const * a_fileB, long a_maxM = 210384, long a_maxN = 2048): 
        m_airp(a_fileA),
        m_birp(a_fileB),
        m_maxM(a_maxM),
        m_maxN(a_maxN),
        m_grid(new double[m_maxN * m_maxM]),
        m_ts(new double[m_maxM]),
        m_S(new double[m_maxN]),
        m_ES(new double[m_maxM]),
        m_VarS(new double[m_maxM])
        {}
        
        ~GridNOP() {
            delete[] (m_grid);
            delete[] (m_ts);
            delete[] (m_S);
            delete[] (m_ES);
            delete[] (m_VarS);
            const_cast<double*&>(m_grid) = nullptr;
            const_cast<double*&>(m_S)    = nullptr;
            const_cast<double*&>(m_ts)   = nullptr;
            const_cast<double*&>(m_ES)   = nullptr;
            const_cast<double*&>(m_VarS) = nullptr;
            
        }
        
        void RunBI (Option<AssetClassA, AssetClassB> const* a_option,
                    Diffusion const* a_diff,
                    double      a_S0,
                    time_t      a_t0,
                    long        a_N = 500,
                    int         a_tauMins = 30,
                    double      a_BFactor = 4.5
                   );
    };
};
                    
          
