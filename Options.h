#pragma once

 
namespace SiriusFM { 
    
template <typename AssetClassA, typename AssetClassB> 
class Option {
protected:
    AssetClassA const m_assetA;
    AssetClassB const m_assetB;
    time_t            m_expirTime;
    bool const        m_isAmerican;
    bool const        m_isAsian;
    
    
    Option(AssetClassA a_assetA, 
           AssetClassB a_assetB, 
           time_t a_Tdays, 
           bool a_isAmerican,
           bool a_isAsian):
    m_assetA        (a_assetA),
    m_assetB        (a_assetB),
    m_expirTime     (a_Tdays),
    m_isAmerican    (a_isAmerican),
    m_isAsian       (a_isAsian)
    {}
public: 
    
    virtual double payoff(long a_L, double const* a_t, double const* a_S) const = 0;
    
    inline time_t ExpirTime() const {
        return m_expirTime;
    }
    inline AssetClassA assetA() const {
        return m_assetA;
    }
     
    inline AssetClassB assetB() const {
        return m_assetB;
    }
    inline bool isAmerican() const {
        return m_isAmerican;
    }
    
    virtual ~Option() {} 
    
};

 using OptionFX = Option<CcyE, CcyE>;
};
