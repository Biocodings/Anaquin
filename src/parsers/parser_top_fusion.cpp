#include <assert.h>
#include "data/tokens.hpp"
#include "parsers/parser_top_fusion.hpp"

using namespace Anaquin;

void ParserTopFusion::parse(const Reader &r, Callback x)
{
    Data data;
    ParserProgress p;
    
    std::vector<std::string> temp, tokens;
    
    while (r.nextTokens(tokens, "@"))
    {
        assert(tokens.size() > 1);

        // chrT-chrT  2082667  4107441  fr  223  37  86  0  47545  81  0.054435
        Tokens::split(std::string(tokens[0]), "\t", tokens);
        
        // "chrT" and "chrT"
        Tokens::split(std::string(tokens[0]), "-", temp);
        
        assert(temp.size() == 2);

        data.cID_1 = temp[0];
        data.cID_2 = temp[1];
        
        // Starting position of the first chromosome
        data.l1 = stoi(tokens[1]) + 1;
        
        // Starting position of the secodn chromosome
        data.l2 = stoi(tokens[2]) + 1;

        data.s1    = tokens[3][0] == 'f' ? Forward : Backward;
        data.s2    = tokens[3][1] == 'f' ? Forward : Backward;;
        data.reads = stoi(tokens[4]);

        p.i++;
        x(data, p);
    }
}