#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/variant.hpp"
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

    struct VariantMatch
    {
        // Matched by position?
        const Variant *match;
        
        // Match by type (SNP, Indel etc)?
        bool type;
        
        // Matched by variant allele?
        bool alt;
        
        // Matched by reference allele?
        bool ref;
    };
    
    /*
     * Common framework for parsing and matching a VCF variant file
     */

    template <typename F> void parseVCF(const FileName &file, F f)
    {
        const auto &r = Standard::instance().r_var;

        ParserVCF::parse(file, [&](const ParserVCF::VCFVariant &v, const ParserProgress &)
        {
            VariantMatch m;

            if (v.chrID != ChrT)
            {
                f(v, nullptr);
            }
            else
            {
                m.match = r.findVar(v.l);
                m.alt   = m.match && m.match->alt  == v.alt;
                m.ref   = m.match && m.match->alt  == v.ref;
                m.type  = m.match && m.match->type == v.type;
                
                f(v, m.match ? &m : nullptr);
                
                /*
                 const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
                 
                 o.writer->write((boost::format(format) % "Start"
                 % "Match"
                 % "Type"
                 % "Alt"
                 % "Rqef").str());
                 */
                

                //if (!found || !type || !alt || !ref)
            //    {
                    // o.writer->write((boost::format(format) % v.l.start
                    //                                      % (found ? "t" : "f")
                    //                                    % (type ? "t" : "f")
                    //                                  % (alt ? "t" : "f")
                    //                                % (ref ? "t" : "f")).str());
                    //return Negative;
             //   }
                
//                if (classify(stats.data.m, v, [&](const VCFVariant &)
                    //stats.chrT->h.at(match->id)++;
            }
        });
    }
}

#endif