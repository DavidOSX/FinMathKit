#pragma once
#include "IRProvider.h"
 

namespace SiriusFM {

template<> class IRProvider<IRMode::Const> {
private:
    double m_IRs[(int)CcyE::N];
    
public:
    IRProvider(char const* a_file); 
    
    double r(CcyE a_ccy, double a_t) const { 
        return m_IRs[(int)a_ccy]; 
    }
    
    double DF(CcyE a_ccy,time_t a_t0, time_t a_t1) const 
    {
        double y = (a_t1 - a_t0)/(365.25*86400.);
        return exp(-m_IRs[(int) a_ccy]*y);
    } 
};

using IRPConst = IRProvider<IRMode::Const>;
    
};
    
