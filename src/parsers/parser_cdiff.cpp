#include "data/standard.hpp"
#include "parsers/parser_cdiff.hpp"

using namespace Anaquin;

bool ParserCDiff::isCDiff(const Reader &r, unsigned &nTrans, unsigned &nGenes)
{
    std::string line;
    std::vector<Token> toks;
    
    // Read the header
    if (r.nextLine(line))
    {
        Tokens::split(line, "\t", toks);
        
        if (toks.size() == 14               &&
            toks[0]  == "test_id"           &&
            toks[1]  == "gene_id"           &&
            toks[2]  == "gene"              &&
            toks[3]  == "locus"             &&
            toks[4]  == "sample_1"          &&
            toks[5]  == "sample_2"          &&
            toks[6]  == "status"            &&
            toks[7]  == "value_1"           &&
            toks[8]  == "value_2"           &&
            toks[9]  == "log2(fold_change)" &&
            toks[10] == "test_stat"         &&
            toks[11] == "p_value"           &&
            toks[12] == "q_value"           &&
            toks[13] == "significant")
        {
            const auto &rs = Standard::instance().r_rna;
            
            // Sequin transcripts
            const auto &l1 = rs.seqsL1();
            
            // Sequin genes
            const auto &l2 = rs.seqsL2();

            nTrans = nGenes = 0;
            
            ParserCDiff::parse(Reader(r), [&](const ParserCDiff::Data &x, const ParserProgress &)
            {
                if      (l1.count(x.iID)) { nTrans++; }
                else if (l2.count(x.iID)) { nGenes++; }
            });

            return true;
        }
    }
    
    return false;
}

