#ifndef V_CLASSIFY_HPP
#define V_CLASSIFY_HPP

#include <boost/format.hpp>
#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    struct VClassify
    {
        struct Stats
        {
            struct FalsePositive
            {
                inline bool operator<(const Locus &l) const  { return this->l < l;    }
                inline bool operator!=(const Locus &l) const { return !operator==(l); }
                inline bool operator==(const Locus &l) const { return this->l == l;   }
                
                Locus l;
                
                // Matching for position?
                bool pos;
                
                // Matching for variant type?
                bool type;
                
                // Matchinf for allele?
                bool alt;
                
                // Matching for reference?
                bool ref;
            };
            
            std::set<FalsePositive> fps;
        };
    };
    
    template <typename Stats, typename Options, typename F> void classify
                (Stats &stats, const FileName &file, const Options &o, F f)
    {
        const auto &r = Standard::instance().r_var;
        
        o.info("Parsing VCF file");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->write((boost::format(format) % "Start"
                                               % "Match"
                                               % "Type"
                                               % "Alt"
                                               % "Rqef").str());

        ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
        {
            if (v.id != Standard::instance().id)
            {
                stats.n_expT++;
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
                                                           % (found ? "t" : "f")
                                                           % (type ? "t" : "f")
                                                           % (alt ? "t" : "f")
                                                           % (ref ? "t" : "f")).str());
                    return Negative;
                }
                
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