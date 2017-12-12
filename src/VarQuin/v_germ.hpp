#ifndef V_GERM_HPP
#define V_GERM_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct EStats
    {
        std::set<Variant> vs;
        
        // Counts by context
        std::map<SequinVariant::Context, Counts> c2c;
        
        // Counts by variation
        std::map<Variation, Counts> v2c;
        
        // Counts by genotype
        std::map<Genotype, Counts> g2c;

        Counts lGC  = 0;
        Counts hGC  = 0;
        Counts vlGC = 0;
        Counts vhGC = 0;
        Counts sRep = 0;
        Counts sSom = 0;
        Counts lRep = 0;
    };

    struct BaseCallerStats
    {
        EStats es; // Statistics for endogenous
    };
    
    struct VGerm
    {
        struct Match
        {
            // The called variant
            Variant qry;
            
            // Sequin matched by position?
            const Variant *var = nullptr;
            
            // Matched by variant allele? Only if position is matched.
            bool alt;
            
            // Matched by reference allele? Only if position is matched.
            bool ref;
            
            // Does the variant fall into one of the reference regions?
            SequinID rID;
        };

        struct Options : public AnalyzerOptions
        {
            Options() : filter(VCFFilter::NotFiltered) {}
            VCFFilter filter;
        };

        struct SStats
        {
            std::vector<Match> tps, fns, fps;

            // Performance by context
            std::map<SequinVariant::Context, Confusion> c2c;

            // Performance by variation
            std::map<Variation, Confusion> v2c;

            // Performance by genotype
            std::map<Genotype, Confusion> g2c;
            
            // Performance for low/high GC contents
            Confusion lGC2c, hGC2c;

            // Performance for very low/high GC contents
            Confusion vlGC2c, vhGC2c;

            // Performance for simple repeats
            Confusion sr2c, lr2c;
            
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
                    if (i.var->name == id)
                    {
                        return &i;
                    }
                }

                return nullptr;
            }
        };
        
        struct Stats : public BaseCallerStats
        {
            SStats ss;
        };

        static EStats analyzeE(const FileName &, const Options &o);
        static SStats analyzeS(const FileName &, const Options &o);

        static Stats report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
