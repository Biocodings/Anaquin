#include <map>
#include "file.hpp"
#include "tokens.hpp"
#include "parsers/parser_tracking.hpp"

using namespace Spike;

enum TrackingField
{
    T_TrackID = 0,
    T_GeneID  = 3,
    T_FPKM    = 9,
    T_FPKM_LO = 10,
    T_FPKM_HI = 11,
    T_Status  = 12,
};

void ParserTracking::parse(const std::string &file, std::function<void (const Tracking &)> f)
{
    File i(file);

    static const std::map<std::string, TrackingStatus> mapper =
    {
        { "OK", OK },
        { "HIDATA", HIData }
    };

    Tracking t;
    
    std::string line;
    std::vector<std::string> tokens;
    
    while (i.nextLine(line))
    {
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
        t.fpkm    = stof(tokens[T_FPKM]);
        t.lFPKM   = stof(tokens[T_FPKM_LO]);
        t.uFPKM   = stof(tokens[T_FPKM_HI]);
        t.status  = mapper.at(tokens[T_Status]);

        if (t.status != TrackingStatus::HIData)
        {
            f(t);
        }
    }
}


