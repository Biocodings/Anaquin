#include <assert.h>
#include "data/tokens.hpp"
#include "parsers/parser_fusion.hpp"

using namespace Spike;

void ParserFusion::parse(const Reader &r, Callback x)
{
    Fusion f;
    ParserProgress p;
    
    std::vector<std::string> temp, tokens;
    
    while (r.nextTokens(tokens, "@"))
    {
        assert(tokens.size() > 1);

        /*
         * chrT-chrT  2082667  4107441  fr  223  37  86  0  47545  81  0.054435
         */

        Tokens::split(tokens[0], "\t", tokens);

        // "chrT" and "chrT"
        Tokens::split(tokens[0], "\t", temp);
        
        assert(temp.size() == 2);
        
        f.chr_1 = temp[0];
        f.chr_2 = temp[1];
        f.pos_1 = stoi(tokens[1]);
        f.pos_2 = stoi(tokens[2]);
        
        const auto orient = tokens[3];
        
        f.ori_1 = orient[0] == 'f' ? Forward : Backward;
        f.ori_2 = orient[0] == 'f' ? Forward : Backward;;

        p.i++;
        x(f, p);
    }
}