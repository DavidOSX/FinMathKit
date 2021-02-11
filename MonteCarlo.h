#pragma once
#include <iostream>
#include <tuple>
#include <time.h>
#include "Diffusion.h"
#include "IRProviderConst.h"

namespace SiriusFM {

template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB
          typename PathEvaluation> class MCEngine {
private:
    long const         m_MaxLP; //long const         m_MaxP;
    double*            m_paths;
    long               m_L;
    long               m_P;/*
    double             m_dt;// timestep
    double             m_t0;
    Diffusion_GBM const*   m_diff;
    IRProvider const*   m_ap; // AProvider
    IRProvider const*   m_bp; // BProvider
    CcyE                 m_b;//AE
    CcyE                 m_a;//BE*/
    
public:
    MCEngine(long a_MaxL, long a_MaxP):
    m_MaxLP(a_MaxLP),
    //m_MaxP(a_MaxP),
    m_paths(new double[m_MaxLP]),
    m_L(0),
    m_P(0)
    { 
              if(m_MaxLP <= 0) throw std::invalid_argument("..."); 
              for(int lp = 0; lp < m_MaxLP; ++lp) m_paths[lp] = 0; 
    }
    
    template<bool> void Simulate(time_t                       a_t0, 
                                 time_t                       a_T, 
                                 int                          a_tau_min,
                                 double                       a_S0,
                                 long                         a_P,
                                 Diffusion const*             a_diff, 
                                 AProvider const*             a_ap,
                                 BProvider const*             a_bp,
                                 AssetClassA                  a_A,
                                 AssetClassB                  a_B
                                 PathEvaluator*               a_pathEval
                                );
                    
    void printPaths() { 
        FILE *f = fopen("PATHS.txt","w");
        fprintf(f,"%ld %ld\n", m_L, m_P);
        for(int i = 0; i < m_L; ++i) fprintf(f,"%lf\n", m_paths[i]);
        fclose(f);
    }
                    
    std::tuple<long, long, double const*> GetPaths() const {
        return (m_L <= 0 || m_P < 0) ? std::make_tuple(0,0,nullptr) : std::make_tuple(m_L,m_P,m_paths);
    }
                    
    ~MCEngine() { delete[] m_paths; }
                    
    MCEngine(MCEngine&) = delete;
                    
    MCEngine& operator= (MCEngine const&) = delete; 
};
    
};    
        
