#include "Diffusion.h"


using namespace SiriusFM; 
    
IRProvider<IRMode::Const>::IRProvider(char const* a_file) {
        char str[16];
        int len = 4, i = 0;
        CcyE cur;
        FILE *f = fopen(a_file,"r");
        for(int j = 0; j < (int)CcyE::N; j++) m_IRs[j] = 0;
        while(fgets(str, len, f)) {
            if(i%2 == 0) {
                len = 10;
                cur = Str2CcyE(str);
            }
            else {
                len = 4;
                m_IRs[(int) cur] = atof(str);
            }
            i++;
        }
              
}
    

