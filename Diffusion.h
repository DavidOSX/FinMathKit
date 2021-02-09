#pragma once
#include <iostream>
#include <cmath>
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

//using namespace std;
namespace SiriusFM {

class Diffusion_GBM
{
private:
    double const m_mu;
    double const m_sigma;
public:
    double mu(double a_s, double t) const { return m_mu*a_s; }
    double sigma(double a_s, double t) const { return (a_s<0) ? 0 : m_sigma*a_s; }
    Diffusion_GBM(double a_m, double a_sigma):
    m_mu(a_m),
    m_sigma(a_sigma)
    { if(m_sigma <= 0) throw std:: invalid_argument("..."); }

};

class Diffusion_CEV {
private:
    double const m_mu;
    double const m_sigma;
    double const m_beta;
public:
    double mu(double a_s, double t) const { return m_mu*a_s; }
    double sigma(double a_s, double t) const { return (a_s<0) ? 0 : m_sigma*a_s*m_beta; }
    Diffusion_CEV(double a_m, double a_sigma, double a_beta):
    m_mu(a_m),
    m_sigma(a_sigma),
    m_beta(a_beta)
    { if(m_sigma <= 0) throw std:: invalid_argument("..."); }

};

class Uhlenbek {
    private:
        double const m_cappa;
        double const m_sigma;
        double const m_theta;
    public:
        double mu(double a_s, double t) const { return m_cappa*(m_theta-a_s); }
        double sigma(double a_s, double t) const { return (a_s<0) ? 0 : m_sigma; }
        Uhlenbek(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { if(m_sigma <= 0) throw std:: invalid_argument("..."); }

};

class Lipton {
    private:
        double const m_sigma_1;
        double const m_sigma_2;
        double const m_sigma_3;
    public:
        //double mu(double a_s, double t)  { return m_cappa*(m_theta-a_s); }
        double sigma(double a_s, double t) const {
            return (a_s<0) ? 0 : m_sigma_1 + m_sigma_2*a_s + m_sigma_3*a_s*a_s;
        }
        Lipton(double a_sigma1, double a_sigma2, double a_sigma3):
        m_sigma_1(a_sigma1),
        m_sigma_2(a_sigma2),
        m_sigma_3(a_sigma3)
        { if(m_sigma_1 <= 0) throw std:: invalid_argument("...");
          if(m_sigma_3 <= 0) throw std:: invalid_argument("...");
          if(m_sigma_2*m_sigma_2-4*m_sigma_1*m_sigma_3 > 0) throw std:: invalid_argument("...");
        }
};

class CIR {
    private:
        double const m_cappa;
        double const m_sigma;
        double const m_theta;
    public:
        double mu(double a_s, double t) const { return m_cappa*(m_theta-a_s); }
        double sigma(double a_s, double t) const { return (a_s<0) ? 0 : m_sigma*sqrt(a_s); }
        CIR(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { if(m_sigma <= 0) throw std:: invalid_argument("..."); }

};

/*enum class CcyE : int {
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

template<> class IRProvider<IRMode::Const> {
private:
    double m_IRs[(int)CcyE::N];
    
public:
    IRProvider(char const* a_file);
    double r(CcyE a_ccy, double a_t) { return m_IRs[(int)a_ccy]; }
};*/
    
};

//class
                       
//Uhlenbek, Lipton
