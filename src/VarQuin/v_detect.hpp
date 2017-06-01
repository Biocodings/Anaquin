#ifndef V_DETECT_HPP
#define V_DETECT_HPP

#include "data/vData.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    // Hask key for mapping a variant
    typedef long VarHashKey;
    
    typedef std::map<VarHashKey, Counts> VarHashTable;

    struct VDetect
    {
        struct Match
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

        enum class Method
        {
            NotFiltered,
            Passed,
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            Method meth = Method::NotFiltered;
        };

        struct EStats
        {
            unsigned found = 0;
        };
        
        struct SStats
        {
            // True positives
            std::vector<Match> tps;
            
            // False postivies
            std::vector<Match> fps;

            /*
             * Performance statistics
             */
            
            // Performance by context (only sequins)
            std::map<SequinVariant::Context, Confusion> g2c;

            // Performance by variation
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
            
            inline const Match * findTP(const SequinID &id) const
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

        static EStats analyzeE(const FileName &, const Options &o);
        static SStats analyzeS(const FileName &, const Options &o);

        // Report for both endogenous and sequins
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
