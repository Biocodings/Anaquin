#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "data/convert.hpp"
#include "parsers/parser_cufflink.hpp"

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

void ParserCufflink::parse(const FileName &file, std::function<void (const ParserCufflink::Data &, const ParserProgress &)> f)
{
    Reader i(file);

    static const std::map<std::string, TrackingStatus> mapper =
    {
        { "OK", TrackingStatus::OK },
        { "HIDATA", TrackingStatus::HIData }
    };

    ParserCufflink::Data t;
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

        t.id  = tokens[T_GeneID];
        t.tID = tokens[T_TrackID];

        // Eg: chrT:1082119-1190836
        Tokens::split(tokens[T_Locus], ":", temp);
        
        t.cID = temp[0];
        
        // Eg: 1082119-1190836
        Tokens::split(std::string(temp[1]), "-", temp);

        // Eg: 1082119, 1190836
        t.l = Locus(stod(temp[0]), stod(temp[1]));
        
        t.status = mapper.at(tokens[T_Status]);
        
        try
        {
            t.abund = s2d(tokens[T_FPKM]);
            t.lFPKM = s2d(tokens[T_FPKM_LO]);
            t.uFPKM = s2d(tokens[T_FPKM_HI]);
        }
        catch (const std::out_of_range &)
        {
            t.abund = t.lFPKM = t.uFPKM = 0;
        }
        
        if (t.status != TrackingStatus::HIData)
        {
            f(t, p);
        }
    }
}