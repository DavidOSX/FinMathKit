#pragma once
#include <algorithm>
#include <cassert>
#include "Options.h"

namespace SiriusFM {

    template <typename AssetClassA, typename AssetClassB>
    class CallOption: public Option<AssetClassA, AssetClassB> {
    private:
        double const m_K;
    
    public:
        CallOption(AssetClassA   a_A, 
                   AssetClassB   a_B, 
                   double        a_K, 
                   time_t        a_Tdays,
                   bool          a_isAmerican ):
        Option<AssetClassA, AssetClassB>(a_A, a_B, a_Tdays, a_isAmerican, false),
        m_K(a_K)
        {
            std::cout << "Call:" << std::endl;
            if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
        }
    
        double payoff(long a_L, double const* a_t, double const* a_path) const override {
            assert(a_L > 0 && a_path != nullptr);
            return std::max<double>(a_path[a_L - 1] - m_K, 0);
        }
        ~CallOption() override {}
    };
    
    template <typename AssetClassA, typename AssetClassB>
    class PutOption: public Option<AssetClassA, AssetClassB> {
    private:
        double const m_K;
    
    public:
        PutOption(AssetClassA    a_A,
                  AssetClassB    a_B,
                  double         a_K, 
                  time_t         a_Tdays,
                  bool           a_isAmerican):
        Option<AssetClassA, AssetClassB>(a_A, a_B, a_Tdays, a_isAmerican, false),
        m_K(a_K)
        {
            std::cout << "Put:" << std::endl;
            if(a_K <= 0 || a_Tdays <= 0) std::invalid_argument("bad option");
        }
    
        double payoff(long a_L, double const* a_t, double const* a_path) const override {
            assert(a_L > 0 && a_path != nullptr);
            return std::max<double>(-a_path[a_L - 1] + m_K, 0);
        }
        ~PutOption() override {}
    
    
    };
    // Aliases ==================================
    using CallOptionFX = CallOption<CcyE, CcyE>;
    using PutOptionFX = PutOption<CcyE, CcyE>;
    ///==========================================
};

