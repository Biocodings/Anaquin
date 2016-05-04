#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include <boost/format.hpp>
#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_varscan.hpp"
#include <boost/algorithm/string/predicate.hpp>

// Defined in main.cpp
extern std::string mixture();

namespace Anaquin
{
    struct VariantStats
    {
        // Number of SNPs detected
        Counts n_snp;

        // Number of indels detected
        Counts n_ind;
    };
    
    /*
     * This class represents variant matching to the synthetic chromosome.
     */
    
    struct VariantMatch
    {
        CalledVariant query;

        // Matched by position?
        const Variant *match = nullptr;
        
        // Matched by sequin region?
        const Variant *seq = nullptr;

        /*
         * Defined only if seq is defined.
         */
        
        Proportion eFold;
        Proportion eAllFreq;
        
        /*
         * Defined only if there's a match
         */

        // Matched by variant allele?
        bool alt;
        
        // Matched by reference allele?
        bool ref;
    };

    inline long var2hash(const SequinID &id, Mutation type, const Locus &l)
    {
        const auto str = (boost::format("%1%_%2%_%3%_%4%") % id
                                                           % type
                                                           % l.start
                                                           % l.end).str();
        return std::hash<std::string>{}(str);
    }
    
    inline bool isRefID(const SequinID &id)
    {
        if (boost::algorithm::ends_with(id, "_R"))
        {
            return true;
        }
        else if (boost::algorithm::ends_with(id, "_V"))
        {
            return false;
        }
        
        throw std::runtime_error("Unknown sequin: " + id);
    }
    
    // Eg: D1_1_R to D1_1
    inline SequinID baseID(const SequinID &id)
    {
        auto tmp = id;

        boost::replace_all(tmp, "_R", "");
        boost::replace_all(tmp, "_V", "");
        
        return tmp;
    }
    
    inline SequinID refID(const SequinID &id)
    {
        return baseID(id) + "_R";
    }
    
    inline SequinID varID(const SequinID &id)
    {
        return baseID(id) + "_V";
    }

    /*
     * Common framework for parsing and matching a called variant
     */

    template <typename F, typename Soft> void parseVariant(const FileName &file, Soft soft, F f)
    {
        const auto &r = Standard::instance().r_var;

        VariantMatch m;

        auto match = [&](const CalledVariant &query)
        {
            m.query = query;
            m.seq   = nullptr;
            m.match = nullptr;

            if (query.cID == ChrT)
            {
                // Can we match by position?
                m.match = r.findVar(query.l, Exact);

                if (m.match)
                {
                    m.seq = m.match;
                    m.ref = m.match->ref == query.ref;
                    m.alt = m.match->alt == query.alt;
                }
                else
                {
                    m.seq = r.findVar(query.l, Contains);
                }
                
                if (m.seq && !mixture().empty())
                {
                    m.eFold    = r.matchFold(baseID(m.seq->id));
                    m.eAllFreq = r.matchAlleleFreq(baseID(m.seq->id));
                }
            }

            return m;
        };

        switch (soft)
        {
            case Soft::GATK:
            {
                ParserVCF::parse(file, [&](const ParserVCF::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }

            case Soft::VarScan:
            {
                ParserVarScan::parse(file, [&](const ParserVarScan::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }
                
            default : { break; }
        }
    }
}

#endif