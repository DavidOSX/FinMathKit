#pragma once
#include <algorithm>
#include <cassert>
#include "Options.h"

namespace SiriusFM {

class EurCallOption: public Option {
private:
    double const m_K;
    
public:
    EurCallOption(double a_K, time_t a_Tdays):
    Option(a_Tdays, false),
    m_K(a_K)
    {
        if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
    }
    
    double payoff(long a_L, double const* a_t, double const* a_path) const override {
        assert(a_L > 1 && a_path != nullptr);
        return std::max<double>(a_path[a_L - 1] - m_K, 0);
    }
    ~EurCallOption() override {}
};

class EurPutOption: public Option {
private:
    double const m_K;
    
public:
    EurPutOption(double a_K, time_t a_Tdays):
    Option(a_Tdays, false),
    m_K(a_K)
    {
        if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
    }
    
    double payoff(long a_L, double const* a_t, double const* a_path) const override {
        assert(a_L > 1 && a_path != nullptr);
        return std::max<double>(-a_path[a_L - 1] + m_K, 0);
    }
    ~EurPutOption() override {}
    
    
};

//using EurCallOptionFX = EurCallOption<CcyE, CcyE>;
//using EurPutOptionFX = EurPutOption<CcyE, CcyE>;
};

