#include "data/reader.hpp"
#include "fusion/f_fusion.hpp"

using namespace Spike;

FFusion::Stats FFusion::analyze(const std::string &file, const Options &options)
{
    Reader r(file);
    std::vector<std::string> tokens;
    
    while (r.nextTokens(tokens, "@"))
    {
        assert(tokens.size() > 1);
        
        /*
         * chrT-chrT  2082667  4107441  fr  223  37  86  0  47545  81  0.054435
         */
        
        std::cout << tokens[0] << std::endl;
        
        
    }
    
    //            bool nextTokens(std::vector<std::string> &, const std::string &c) const;
    
    //            bool nextLine(std::string &) const;
    
    
    
    
    
    return FFusion::Stats();
}