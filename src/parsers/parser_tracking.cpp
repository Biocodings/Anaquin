#include <map>
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser_tracking.hpp"

using namespace Anaquin;

enum TrackingField
{
    T_TrackID = 0,
    T_GeneID  = 3,
    T_Locus   = 6,
    T_FPKM    = 9,
    T_FPKM_LO = 10,
    T_FPKM_HI = 11,
    T_Status  = 12,
};

void ParserTracking::parse(const std::string &file, std::function<void (const Tracking &, const ParserProgress &)> f)
{
    Reader i(file);

    static const std::map<std::string, TrackingStatus> mapper =
    {
        { "OK", OK },
        { "HIDATA", HIData }
    };

    Tracking t;
    ParserProgress p;
    
    std::string line;
    std::vector<std::string> temp, tokens;
    
    while (i.nextLine(line))
    {
        p.i++;
        Tokens::split(line, "\t", tokens);
       
        /*
         * tracking_id  code  nearest_ref  gene_id  gene_short  tss_id  locus  length  coverage  FPKM  FPKM_conf_lo  FPKM_conf_hi  FPKM_status
         */

        if (!mapper.count(tokens[T_Status]))
        {
            continue;
        }
        
        assert(!tokens[T_GeneID].empty());
        assert(!tokens[T_FPKM].empty());
        assert(!tokens[T_FPKM_LO].empty());
        assert(!tokens[T_FPKM_HI].empty());
        assert(!tokens[T_Status].empty());

        t.trackID = tokens[T_TrackID];
        t.geneID  = tokens[T_GeneID];

        // Eg: chrT:1082119-1190836
        Tokens::split(tokens[T_Locus], ":", temp);
        
        // Eg: chrT
        const auto c = temp[0];
        
        // Eg: 1082119-1190836
        Tokens::split(temp[1], "-", temp);
        
        // Eg: 1082119, 1190836
        t.l = Locus(stod(temp[0]), stod(temp[1]));
        
        t.fpkm    = stod(tokens[T_FPKM]);
        t.lFPKM   = stod(tokens[T_FPKM_LO]);
        t.uFPKM   = stod(tokens[T_FPKM_HI]);
        t.status  = mapper.at(tokens[T_Status]);

        if (t.status != TrackingStatus::HIData)
        {
            f(t, p);
        }
    }
}