#ifndef VARQUIN_HPP
#define VARQUIN_HPP

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

    enum class Caller
    {
        GATK,
        VarScan,
    };

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
    
    inline SequinID baseID(const SequinID &id)
    {
        auto tmp = id;

        boost::replace_all(tmp, "_R", "");
        boost::replace_all(tmp, "_V", "");
        
        return tmp;
    }
    
    /*
     * Common framework for parsing and matching a variant output
     */

    template <typename F> void parseVariant(const FileName &file, Caller caller, F f)
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
                    m.eFold    = r.fold(baseID(m.seq->id));
                    m.eAllFreq = r.alleleFreq(baseID(m.seq->id));
                }
            }

            return m;
        };

        switch (caller)
        {
            case Caller::GATK:
            {
                ParserVCF::parse(file, [&](const ParserVCF::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }

            case Caller::VarScan:
            {
                ParserVarScan::parse(file, [&](const ParserVarScan::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }
        }
    }
}

#endif