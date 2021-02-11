#pragma once
#include <algorithm>
#include <cassert>
#include "Options.h"

namespace SiriusFM {

    template <typename AssetClassA, typename AssetClassB>
    class EurCallOption: public Option<AssetClassA, AssetClassB> {
    private:
        double const m_K;
    
    public:
        EurCallOption(AssetClassA   a_A, 
                      AssetClassB   a_B, 
                      double        a_K, 
                      time_t        a_Tdays):
        Option<AssetClassA, AssetClassB>(a_A, a_B, a_Tdays, false),
        m_K(a_K)
        {
            std::cout << "Call:" << std::endl;
            if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
        }
    
        double payoff(long a_L, double const* a_t, double const* a_path) const override {
            assert(a_L > 1 && a_path != nullptr);
            return std::max<double>(a_path[a_L - 1] - m_K, 0);
        }
        ~EurCallOption() override {}
    };
    
    template <typename AssetClassA, typename AssetClassB>
    class EurPutOption: public Option<AssetClassA, AssetClassB> {
    private:
        double const m_K;
    
    public:
        EurPutOption(AssetClassA    a_A,
                     AssetClassB    a_B,
                     double         a_K, 
                     time_t         a_Tdays):
        Option<AssetClassA, AssetClassB>(a_A, a_B, a_Tdays, false),
        m_K(a_K)
        {
            std::cout << "Put:" << std::endl;
            if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
        }
    
        double payoff(long a_L, double const* a_t, double const* a_path) const override {
            assert(a_L > 1 && a_path != nullptr);
            return std::max<double>(-a_path[a_L - 1] + m_K, 0);
        }
        ~EurPutOption() override {}
    
    
    };
    // Aliases ==================================
    using EurCallOptionFX = EurCallOption<CcyE, CcyE>;
    using EurPutOptionFX = EurPutOption<CcyE, CcyE>;
    ///==========================================
};

