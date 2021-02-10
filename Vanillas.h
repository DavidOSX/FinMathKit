
#include <algorithm>
#include "Options.h"

namespace SiriusFM {

class EurCallOption: public Option {
private:
    double const m_K;
    
public:
    EurCallOption(double a_K, long a_Tdays):
    Option(false, a_Tdays),
    m_K(a_K)
    {
        if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
    }
    
    double payoff(long a_L, double const* a_t, double const* a_S) const override {
        return std::max<double>(a_S[a_L - 1] - m_K, 0);
    }
    ~EurCallOption() override {}
};

class EurPutOption: public Option {
private:
    double const m_K;
    
public:
    EurPutOption(double a_K, long a_Tdays):
    Option(false, a_Tdays),
    m_K(a_K)
    {
        if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
    }
    
    double payoff(long a_L, double const* a_t, double const* a_S) const override {
        return std::max<double>(-a_S[a_L - 1] + m_K, 0);
    }
    ~EurPutOption() override {}
};
};

