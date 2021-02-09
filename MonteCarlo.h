#pragma once
#include <iostream>
#include <tuple>
#include <time.h>
#include "Diffusion.h"
#include "IRProviderConst.h"

using namespace std;

namespace SiriusFM {

template <typename Diffusion, 
          typename AProvider, 
          typename BProvider, 
          typename AssetClassA, 
          typename AssetClassB> class MCEngine {// 
private:
    long const         m_MaxL;
    long const         m_MaxP;
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
    m_MaxL(a_MaxL),
    m_MaxP(a_MaxP),
    m_paths(new double[m_MaxL*m_MaxP]),
    m_L(0),
    m_P(0)/*
    m_dt(nan),
    m_t0(nan),
    m_diff(nullptr),
    m_ap(nullptr),
    m_bp(nullptr),
    m_a(),//AssetClassA::...
    m_b()//AssetClassB::...*/
    { if(m_MaxL <= 0|| m_MaxP <= 0) throw std::invalid_argument("..."); }
    
    template<bool> void Simulate(time_t           a_t0, 
                  time_t           a_T, 
                  int              a_tau_min,
                  double           a_S0,
                  long             a_P,
                  Diffusion_GBM const* a_diff, 
                  IRProvider<IRMode::Const> const* a_ap,
                  IRProvider<IRMode::Const> const* a_bp,
                  CcyE      a_A,
                  CcyE      a_B
                  /*bool             a_isRN*/
                  );
    void printPaths() { 
        FILE *f = fopen("PATHS.txt","w");
        //fprintf(f,"%ld %ld\n", m_L, m_P);
        for(int i = 0; i < m_L; ++i) fprintf(f,"%lf\n", m_paths[i]);
        fclose(f);
    }
    std::tuple<long,long,double const*> GetPaths() const {
        return (m_L <= 0||m_P<=0)?std::make_tuple(0,0,nullptr):std::make_tuple(m_L,m_P,m_paths);
    }
    ~MCEngine() { delete[] m_paths; }
    MCEngine(MCEngine&) = delete;
    MCEngine& operator= (MCEngine const&) = delete; 
};
    
};    
        
