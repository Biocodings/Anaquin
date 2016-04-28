#ifndef F_CLASSIFY_HPP
#define F_CLASSIFY_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_top_fusion.hpp"
#include "parsers/parser_star_fusion.hpp"

namespace Anaquin
{
    enum FusionCaller
    {
        StarFusion,
        TopHatFusion,
    };

    struct FUSQuin
    {
        enum Label
        {
            Geno     = -2,
            GenoChrT = -1,
            Negative = 0,
            Positive = 1,
        };

        struct Match
        {
            // State of the query
            Label label;

            CalledFusion query;

            // Only defined if the label is Positive
            const FusionRef::KnownFusion *known;
        };
        
        typedef CalledFusion FalsePositive;

//        struct Results
//        {
//            Results(Code code, const FusionRef::FusionPoint *match = nullptr) : code(code), match(match) {}
//
//            Code code;
//
//            // Where the fusion matches
//            const FusionRef::FusionPoint *match;
//        };

//        template <typename Options, typename T> static Results classifyFusion(const T &f, Options &o)
//        {
//            if (f.chr_1 != ChrT || f.chr_2 != ChrT)
//            {
//                if (f.chr_1 != ChrT && f.chr_2 != ChrT)
//                {
//                    return Results(Endo);
//                }
//                else
//                {
//                    return Results(EndoChrT);
//                }
//            }
//
//            const auto min = std::min(f.l1, f.l2);
//            const auto max = std::max(f.l1, f.l2);
//            const auto m   = Standard::instance().r_fus.find(min, max, f.s1, f.s2, o.fuzzy);
//
//            return Results(m ? Code::Positive : Code::Negative, m);
//        }

        template <typename T> static Match classifyFusion(const T &f, double fuzzy)
        {
            Match match;
            match.query = f;
            
            if (f.cID_1 != ChrT || f.cID_2 != ChrT)
            {
                if (f.cID_1 != ChrT && f.cID_2 != ChrT)
                {
                    match.label = Geno;
                }
                else
                {
                    match.label = GenoChrT;
                }
            }
            else
            {
                match.label = Negative;

                const auto min = std::min(f.l1, f.l2);
                const auto max = std::max(f.l1, f.l2);
                
                const auto &r = Standard::instance().r_fus;

                // Any match by position with fuzzy?
                if ((match.known = r.findFusion(min, max, f.s1, f.s2, fuzzy)))
                {
                    match.label = Positive;
                }
            }

            return match;
        }

        template <typename Options, typename F> static void analyze(const FileName &file, const Options &o, F f)
        {
            o.info("Fuzzy level: " + std::to_string(o.fuzzy));
            o.info("Parsing alignment file");

            switch (o.caller)
            {
                case StarFusion:
                {
                    ParserStarFusion::parse(Reader(file), [&](const CalledFusion &t, const ParserProgress &)
                    {
                        f(classifyFusion(t, o.fuzzy));
                    });

                    break;
                }
                    
                case TopHatFusion:
                {
                    ParserTopFusion::parse(Reader(file), [&](const CalledFusion &t, const ParserProgress &)
                    {
                        f(classifyFusion(t, o.fuzzy));
                    });

                    break;
                }
            }
        }
    };
}

#endif