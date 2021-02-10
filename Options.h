#pragma once

namespace SiriusFM {

class Option {
protected:
    bool const m_isAmerican;
    long const m_Tdays;
    Option(bool a_isAmerican, long a_Tdays):
    m_isAmerican(a_isAmerican),
    m_Tdays(a_Tdays)
    {}
public:
    virtual double payoff(long a_L, double const* a_t, double const* a_S) const = 0;
    
    virtual ~Option() {}
};
};
