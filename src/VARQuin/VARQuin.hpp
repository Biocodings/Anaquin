#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/variant.hpp"
#include "data/standard.hpp"

#include "parsers/parser_vcf.hpp"
#include "parsers/parser_varscan.hpp"

namespace Anaquin
{
    /*
     * This class represents variant matching to the synthetic chromosome.
     */
    
    struct VariantMatch
    {
        const CalledVariant * query;

        // Matched by position?
        const Variant *match = nullptr;
        
        // Matched by sequin region?
        const Variant *seq = nullptr;

        /*
         * Expected allele frequency. Defined only if seq is defined.
         */
        
        double eAllFreq;
        
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

    /*
     * Common framework for parsing and matching a variant output
     */

    template <typename F> void parseVariant(const FileName &file, Caller caller, F f)
    {
        const auto &r = Standard::instance().r_var;

        VariantMatch m;

        auto match = [&](const CalledVariant &query)
        {
            m.query = &query;
            m.seq   = nullptr;
            m.match = nullptr;

            if (query.chrID == ChrT)
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
                
                if (m.seq)
                {
                    m.eAllFreq = r.alleleFreq(m.seq->id);
                }
            }

            return m;
        };

        switch (caller)
        {
            case Caller::GATK:
            {
                ParserVCF::parse(file, [&](const ParserVCF::VCFVariant &v, const ParserProgress &)
                {
                    //f(match(v));
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