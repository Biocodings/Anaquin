#ifndef F_DISCOVER_HPP
#define F_DISCOVER_HPP

#include "stats/analyzer.hpp"
#include "FusQuin/FUSQuin.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller caller;
        };

        struct Stats : public FusionStats
        {
            typedef FUSQuin::Match ChrTData;

            #define TP(x) data.at(x).tps
            #define FP(x) data.at(x).fps
            #define FN(x) data.at(x).fns
            
            inline Counts countTP(const ChrID &id) const
            {
                return TP(id).size();
            }

            inline Counts countFP(const ChrID &id) const
            {
                return FP(id).size();
            }
            
            inline Counts countKnown(const ChrID &id) const
            {
                return TP(id).size() + FN(id).size();
            }

            inline Counts countDetect(const ChrID &id) const
            {
                return TP(id).size() + FP(id).size();
            }

            inline Proportion sn(const ChrID &id) const
            {
                return static_cast<Proportion>(TP(id).size()) / (TP(id).size() + FN(id).size());
            }

            inline Proportion pc(const ChrID &id) const
            {
                return static_cast<Proportion>(TP(id).size()) / countDetect(id);
            }

            struct Data : public SequinStats
            {
                // List of missing fusions
                std::vector<FusionRef::KnownFusion> fns;

                std::vector<FUSQuin::Match> fps, tps;
            };

            std::map<ChrID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif