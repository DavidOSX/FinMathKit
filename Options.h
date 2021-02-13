#pragma once

 
namespace SiriusFM {
    
template <typename AssetClassA, typename AssetClassB>
class Option {
protected:
    AssetClassA const m_assetA;
    AssetClassB const m_assetB;
    bool const        m_isAmerican;
    time_t            m_expirTime;
    
    Option(AssetClassA a_assetA, 
           AssetClassB a_assetB, 
           time_t a_Tdays, 
           bool a_isAmerican):
    m_assetA        (a_assetA),
    m_assetB        (a_assetB),
    m_expirTime     (a_Tdays),
    m_isAmerican    (a_isAmerican)
    {}
public: 
    
    virtual double payoff(long a_L, double const* a_t, double const* a_S) const = 0;
    time_t ExpirTime() const {
        return m_expirTime;
    }
    AssetClassA assetA() const {
        return m_assetA;
    }
    
    AssetClassB assetB() const {
        return m_assetB;
    }
    
    virtual ~Option() {} 
    
};

 using OptionFX = Option<CcyE, CcyE>;
};
