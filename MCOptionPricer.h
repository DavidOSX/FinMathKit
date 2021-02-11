template <Diffusion, 
          AProvider, 
          BProvider, 
          AssetClassA, 
          AssetClassB>
class MCOptionPricer {
    
private:
    AProvider                   m_airp;
    BProvider                   m_birp;
    MCEngine<Diffusion, 
             AProvider, 
             BProvider, 
             AssetClassA, 
             AssetClassB, 
             PathEvaluator>     m_mce;
    
public:
    
    MCOptionPricer(const char* a_fileA,
                   const char* a_fileB);
    
        
