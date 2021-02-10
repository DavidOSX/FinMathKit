#pragma once
#include <iostream>
#include <cmath>
#include <stdio.h>

namespace SiriusFM {

class Diffusion_GBM
{
private:
    double const m_mu;
    double const m_sigma;
public:
    double mu(double a_s, double t) const { 
        return m_mu*a_s; 
    }
    double sigma(double a_s, double t) const { 
        return (a_s<0) ? 0 : m_sigma*a_s; 
    }
    Diffusion_GBM(double a_m, double a_sigma):
    m_mu(a_m),
    m_sigma(a_sigma)
    { 
        if(m_sigma <= 0) throw std::invalid_argument("..."); 
    }

};

class Diffusion_CEV {
private:
    double const m_mu;
    double const m_sigma;
    double const m_beta;
public:
    double mu(double a_s, double t) const { 
        return m_mu*a_s; 
    }
    double sigma(double a_s, double t) const { 
        return (a_s<0) ? 0 : m_sigma*a_s*m_beta; 
    }
    Diffusion_CEV(double a_m, double a_sigma, double a_beta):
    m_mu(a_m),
    m_sigma(a_sigma),
    m_beta(a_beta)
    { 
        if(m_sigma <= 0) throw std:: invalid_argument("bad sigma"); 
    }

};

class Uhlenbek {
    private:
        double const m_cappa;
        double const m_sigma;
        double const m_theta;
    public:
        double mu(double a_s, double t) const { 
            return m_cappa*(m_theta-a_s); 
        }
        double sigma(double a_s, double t) const { 
            return (a_s<0) ? 0 : m_sigma; 
        }
        Uhlenbek(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { 
            if(m_sigma <= 0) throw std:: invalid_argument("bad sigma"); 
        }

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
        { 
          if(m_sigma_1 <= 0) throw std:: invalid_argument("...");
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
        double mu(double a_s, double t) const { 
            return m_cappa*(m_theta-a_s); 
        }
        double sigma(double a_s, double t) const { 
            return (a_s<0) ? 0 : m_sigma*sqrt(a_s); 
        }
        CIR(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { 
            if(m_sigma <= 0) throw std:: invalid_argument("..."); 
        }

};

