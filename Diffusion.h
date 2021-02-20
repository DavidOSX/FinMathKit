#pragma once
#include <iostream>
#include <cmath>


namespace SiriusFM {
    
class Diffusion {
protected:
    double const m_mu;
    double const m_sigma;
    double const m_S0;
public:
    
    Diffusion(double a_m, double a_sigma, double a_S0):
    m_mu    (a_m),
    m_sigma (a_sigma),
    m_S0    (a_S0)
    {
        if(m_sigma <= 0) throw std::invalid_argument("bad sigma");
    }
    
    virtual double mu(double a_s, double t) const = 0;
    
    virtual double sigma(double a_s, double t) const = 0;
    
    inline double GetStartPoint() const {
        return m_S0;
    }   
    
    virtual ~Diffusion() {}
    
};    

class Diffusion_GBM: public Diffusion {
//private:
  //  double const m_mu;
  //  double const m_sigma;
  //  double const m_S0;
public:
    
    Diffusion_GBM(double a_m, double a_sigma, double a_S0):
    Diffusion(a_m, a_sigma, a_S0)
    {}
    
    double mu(double a_s, double t) const override { 
        return m_mu*a_s; 
    }
    
    double sigma(double a_s, double t) const override { 
        return (a_s<0) ? 0 : m_sigma*a_s; 
    }
    
     ~Diffusion_GBM() override {}

};

class Diffusion_CEV: public Diffusion {
//private:
//    double const m_mu;
 //   double const m_sigma;
    double const m_beta;
    //double const m_S0;
    
public:
    
    Diffusion_CEV(double a_m, double a_sigma, double a_beta, double a_S0):
    Diffusion(a_m, a_sigma, a_S0),
    m_beta  (a_beta)
    { 
        //if(m_sigma <= 0) throw std:: invalid_argument("bad sigma"); 
    }
    
    double mu(double a_s, double t) const override { 
        return m_mu*a_s; 
    }
    
    double sigma(double a_s, double t) const override { 
        return (a_s<0) ? 0 : m_sigma * pow(a_s, m_beta); 
    }
    
    ~Diffusion_CEV() override {}

};

//================

class Lipton {
    private:
        double const m_mu;
        double const m_sigma_1;
        double const m_sigma_2;
        double const m_sigma_3;
    public:
        
        Lipton(double a_mu, double a_sigma1, double a_sigma2, double a_sigma3):
        m_mu        (a_mu),
        m_sigma_1   (a_sigma1),
        m_sigma_2   (a_sigma2),
        m_sigma_3   (a_sigma3)
        { 
          if(m_sigma_1 <= 0) throw std:: invalid_argument("bad sigma 1");
          if(m_sigma_3 <= 0) throw std:: invalid_argument("bad sigma 2");
          if(m_sigma_2*m_sigma_2-4*m_sigma_1*m_sigma_3 > 0) throw std:: invalid_argument("bad discriminant");
        }
        
        double mu(double a_s, double t) const { 
            return m_mu * a_s; 
        }
        double sigma(double a_s, double t) const {
            return (a_s<0) ? 0 : m_sigma_1 + m_sigma_2*a_s + m_sigma_3*a_s*a_s;
        }
        
};

//models for IR:

class CIR {
    private:
        double const m_cappa;
        double const m_sigma;
        double const m_theta;
        
    public:
        
        CIR(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { 
            if(m_sigma <= 0) throw std:: invalid_argument("bad sigma"); 
        }
        
        double mu(double a_s, double t) const { 
            return m_cappa*(m_theta-a_s); 
        }
        
        double sigma(double a_s, double t) const { 
            return (a_s<0) ? 0 : m_sigma*sqrt(a_s); 
        }
        

};


class Uhlenbek {
    private:
        double const m_cappa;
        double const m_sigma;
        double const m_theta;
        
    public:
        
        Uhlenbek(double a_c, double a_sigma, double a_theta):
        m_cappa(a_c),
        m_sigma(a_sigma),
        m_theta(a_theta)
        { 
            if(m_sigma <= 0) throw std:: invalid_argument("bad sigma"); 
        }
        
        double mu(double a_s, double t) const { 
            return m_cappa*(m_theta - a_s); 
        }
        
        double sigma(double a_s, double t) const { 
            return (a_s<0) ? 0 : m_sigma; 
        }
        

};

};

