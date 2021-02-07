#pragma once
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
namespace SiriusFM {

enum class CcyE : int {
    USD = 0,
    EUR = 1,
    GBP = 2,
    CHF = 3,
    N = 4

};

enum class IRMode:int {
    Const = 0,
    ForwCurve = 1,
    Stock = 2
};

inline char const* CcyE2Str(CcyE a_ccy) {
    switch(a_ccy) {
        case CcyE::USD:
            return "USD";
            break;
        case CcyE::EUR:
            return "EUR";
            break;
        case CcyE::GBP:
            return "GBP";
            break;
        case CcyE::CHF:
            return "CHF";
            break;
        default:
            throw std::invalid_argument("...");
            break;
    }
};

inline int aux(const char* a_ccy) {
    if(strcmp(a_ccy, "USD") == 0) return 1;
    if(strcmp(a_ccy, "EUR") == 0) return 2;
    if(strcmp(a_ccy, "GBP") == 0) return 3;
    if(strcmp(a_ccy, "CHF") == 0) return 4;
    return 0;
};

inline  CcyE Str2CcyE(const char* a_ccy) {
    switch(aux(a_ccy)) {
        case 1:
            return CcyE::USD;
            break;
        case 2:
            return CcyE::EUR;
            break;
        case 3:
            return CcyE::GBP;
            break;
        case 4:
            return CcyE::CHF;
            break;
        default:
            throw std::invalid_argument("...");
            break;
    }
};


template<IRMode IRM> class IRProvider;


};
