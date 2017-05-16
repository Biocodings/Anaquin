#ifndef V_CANCER_HPP
#define V_CANCER_HPP

#include <set>
#include <vector>
#include "data/vData.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    // Hask key for mapping a variant
    typedef long VarHashKey;
    
    typedef std::map<VarHashKey, Counts> VarHashTable;

    struct VCancer
    {
        struct VariantMatch
        {
            // The called variant
            Variant query;
            
            // Sequin matched by position?
            const Variant *seqByPos = nullptr;
            
            // Matched by variant allele? Only if position is matched.
            bool alt;
            
            // Matched by reference allele? Only if position is matched.
            bool ref;
            
            // Does the variant fall into one of the reference regions?
            SequinID rReg;
        };

        typedef AnalyzerOptions Options;
        
        struct Stats
        {
            // True positives
            std::vector<VariantMatch> tps;
            
            // False postivies
            std::vector<VariantMatch> fps;

            /*
             * Performance statistics
             */
            
            // Performance for different context (only sequins)
            std::map<SeqVariant::Context, Confusion> g2c;

            // Performance for different variation
            std::map<Variation, Confusion> m2c;

            // Overall performance
            Confusion oc;

            /*
             * Statistics for allele frequency
             */
            
            struct AlleleStats : public SequinStats, public LimitStats {};

            // Performance for each variation
            std::map<Variation, AlleleStats> m2a;

            // Overall performance
            AlleleStats oa;
            
            inline const VariantMatch * findTP(const SequinID &id) const
            {
                for (auto &i : tps)
                {
                    if (i.seqByPos->name == id)
                    {
                        return &i;
                    }
                }

                return nullptr;
            }
        };

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
