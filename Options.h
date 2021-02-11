#pragma once

namespace SiriusFM {
//template <AssetClassA, AssetClassB>
class Option {
protected:
    bool const m_isAmerican;
    //long const m_Tdays;
    time_t m_expirTime;
    Option(time_t a_Tdays, bool a_isAmerican):
    m_expirTime(a_Tdays),
    m_isAmerican(a_isAmerican)
    {}
public:
    virtual double payoff(long a_L, double const* a_t, double const* a_S) const = 0;
    time_t ExpirTime() const {
        return m_expirTime;
    }
    virtual ~Option() {}
    
};

 //using OptionFX = Option<CcyE, CcyE>;
};
