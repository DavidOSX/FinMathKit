#pragma once
#include "IRProvider.h"


using namespace SiriusFM;

template<> class IRProvider<IRMode::Const> {
private:
    double m_IRs[(int)CcyE::N];
    
public:
    IRProvider(char const* a_file);
    double r(CcyE a_ccy, double a_t) const { return m_IRs[(int)a_ccy]; }
};
    
