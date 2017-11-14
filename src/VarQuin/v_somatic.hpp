#ifndef V_SOMATIC_HPP
#define V_SOMATIC_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/v_germ.hpp"

namespace Anaquin
{
    struct VSomatic
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

            // Performance by allele frequency group
            std::map<float, Confusion> f2c;
            
            // Performance by context
            std::map<SequinVariant::Context, Confusion> c2c;
            
            // Performance by variation
            std::map<Variation, Confusion> v2c;
            
            // Performance by genotype
            std::map<Genotype, Confusion> g2c;
            
            // Overall performance
            Confusion oc;
            
            /*
             * Caller specific fields
             */
            
            std::map<std::string, std::map<long, int>>   si;
            std::map<std::string, std::map<long, float>> sf;
            
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
