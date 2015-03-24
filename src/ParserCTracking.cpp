#include <map>
#include <fstream>
#include <vector>
#include <assert.h>
#include "ParserCTracking.hpp"
#include <boost/algorithm/string.hpp>

bool ParserCTracking::parse(const std::string &file, std::function<void (const CTracking &)> f)
{
    std::string line;
    std::ifstream i(file);
    
    if (!i)
    {
        return false;
    }
    
    std::map<std::string, CTrackingStatus> mapper =
    {
        { "OK", OK },
        { "HIDATA", HIData }
    };

    CTracking t;
    std::vector<std::string> tokens;
    
    while (std::getline(i, line))
    {
        boost::split(tokens, line, boost::is_any_of("\t"));
       
        /*
         * tracking_id  code  nearest_ref  gene_id  gene_short  tss_id  locus  length  coverage  FPKM  FPKM_conf_lo  FPKM_conf_hi  FPKM_status
         *   R_9_1_V     -        -        R_9_1_V      -         -       -       -        -     1188     1131          1244           OK
         */

        if (!mapper.count(tokens[12]))
        {
            continue;
        }
        
        assert(!tokens[3].empty());
        assert(!tokens[9].empty());
        assert(!tokens[10].empty());
        assert(!tokens[11].empty());
        assert(!tokens[12].empty());

        t.geneID  = tokens[3];
        t.fpkm    = stof(tokens[9]);
        t.lFPKM   = stof(tokens[10]);
        t.uFPKM   = stof(tokens[11]);
        t.trackID = tokens[0];
        t.status  = mapper[tokens[12]];

        if (t.status != CTrackingStatus::HIData)
        {
            f(t);
        }
    }
    
    return true;
}


