#ifndef V_STRUCTURE_HPP
#define V_STRUCTURE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VStructure
    {
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

        struct Match
        {
            // Called variant
            Variant qry;
            
            // Matched sequin
            const Variant *var = nullptr;

            // Does the variant fall into reference regions?
            SequinID rID;
        };
        
        struct EStats
        {
            std::map<Variation, Counts> v2c;
        };
        
        struct SStats
        {
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
           
            // Overall performance
            Confusion oc;
            
            std::vector<Match> tps, fps;
            
            // Performance by variation
            std::map<Variation, Confusion> v2c;
        };

        static EStats analyzeE(const FileName &, const Options &o);
        static SStats analyzeS(const FileName &, const Options &o);

        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
