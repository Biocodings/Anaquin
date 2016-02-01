#ifndef VARIANT_HPP
#define VARIANT_HPP

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

    /*
     * Common framework for parsing a VCF variant file
     */

    template <typename Stats, typename F> void parseVCF(const FileName &file, Stats &stats, F f)
    {
        const auto &r = Standard::instance().r_var;

/*
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->write((boost::format(format) % "Start"
                                               % "Match"
                                               % "Type"
                                               % "Alt"
                                               % "Rqef").str());
*/

        ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
        {
            
            
            
            
          //  if (v.id != ChrT)
           // {
                //stats.chrT->n_endo++;
            //    return;
            //}

            //stats.chrT->n_chrT++;

            /*
             * The following information is required:
             *
             *     - 
             *
             *
             */
             
            
            
            
            
            
            const Variation *match;

            if (classify(stats.data.m, v, [&](const VCFVariant &)
            {
                const auto found = (match = r.findVar(v.l)) != nullptr;
                const auto type  = (match && match->type == v.type);
                const auto alt   = (match && match->alt  == v.alt);
                const auto ref   = (match && match->ref  == v.ref);
                
                if (!found || !type || !alt || !ref)
                {
                   // o.writer->write((boost::format(format) % v.l.start
                     //                                      % (found ? "t" : "f")
                       //                                    % (type ? "t" : "f")
                         //                                  % (alt ? "t" : "f")
                           //                                % (ref ? "t" : "f")).str());
                    return Negative;
                }
                
                f(v, match);

                return Positive;
            }))
            {
                //stats.chrT->h.at(match->id)++;
            }
        });
        
        //o.writer->close();
    }
}

#endif