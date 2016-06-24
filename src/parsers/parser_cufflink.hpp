#ifndef PARSER_CUFFLINK_HPP
#define PARSER_CUFFLINK_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserCufflink
    {
        enum class TrackingStatus
        {
            OK,
            HIData,
            NoTest
        };
        
        struct Data
        {
            ChrID cID;
            
            GeneID id;
            
            // It's also known as trackID
            IsoformID tID;
            
            Locus l;
            
            // The expression for the abundance
            FPKM abund;

            // Lower 95% confidence
            FPKM lFPKM;
            
            // Upper 95% confidence
            FPKM uFPKM;
            
            TrackingStatus status;
        };

        static void parse(const FileName &, std::function<void (const Data &, const ParserProgress &)>);
    };
}

#endif