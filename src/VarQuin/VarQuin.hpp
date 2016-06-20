#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_variants.hpp"
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
     * This class represents variant matching to the synthetic chromosome
     */
    
    struct VariantMatch
    {
        CalledVariant query;

        // Matched by position?
        const Variant *match = nullptr;
        
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
     * Common framework for parsing and matching called variants
     */

    template <typename F, typename Input> void parseVariants(const FileName &file, Input input, F f)
    {
        const auto &r = Standard::instance().r_var;

        VariantMatch m;

        auto match = [&](const CalledVariant &query)
        {
            m.query = query;
            m.match = nullptr;

            if (Standard::isSynthetic(query.cID))
            {
                // Can we match by position?
                m.match = r.findVar(query.cID, query.l);

                if (m.match)
                {
                    m.ref = m.match->ref == query.ref;
                    m.alt = m.match->alt == query.alt;
                }
                
                if (m.match && !mixture().empty())
                {
                    m.eFold    = r.findAFold(baseID(m.match->id));
                    m.eAllFreq = r.findAFreq(baseID(m.match->id));
                }
            }

            return m;
        };

        switch (input)
        {
            case Input::VCF:
            {
                ParserVCF::parse(file, [&](const ParserVCF::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }

            case Input::Text:
            {
                ParserVariant::parse(file, [&](const ParserVariant::Data &d, const ParserProgress &)
                {
                    f(match(d));
                });

                break;
            }
        }
    }
}

#endif