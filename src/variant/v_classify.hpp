#ifndef GI_V_CLASSIFY_HPP
#define GI_V_CLASSIFY_HPP

#include <boost/format.hpp>
#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    template <typename Stats, typename Options, typename F> void classify
                (Stats &stats, const std::string &file, const Options &o, F f)
    {
        const auto &r = Standard::instance().r_var;
        
        o.info("Parsing VCF file");
        o.writer->open("VarAllele_false.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->write((boost::format(format) % "start"
                                               % "matched"
                                               % "type"
                                               % "alt"
                                               % "ref").str());

        ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
        {
            if (v.id != Standard::instance().id)
            {
                stats.n_hg38++;
                return;
            }
            
            stats.n_chrT++;
            const Variation *match;

            if (classify(stats.m, v, [&](const VCFVariant &)
            {
                const auto found = (match = r.findVar(v.l)) != nullptr;
                const auto type  = (match && match->type == v.type);
                const auto alt   = (match && match->alt  == v.alt);
                const auto ref   = (match && match->ref  == v.ref);
                
                if (!found || !type || !alt || !ref)
                {
                    o.writer->write((boost::format(format) % v.l.start
                                                           % found
                                                           % type
                                                           % alt
                                                           % ref).str());
                    return Negative;
                }
                
                stats.detected++;
                
                f(v, match);

                return Positive;
            }))
            {
                stats.h.at(match->id)++;
            }
        });
        
        o.writer->close();
    }
}

#endif