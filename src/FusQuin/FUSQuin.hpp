#ifndef F_CLASSIFY_HPP
#define F_CLASSIFY_HPP

#include "data/standard.hpp"
#include "stats/analyzer.hpp"
#include "data/reference.hpp"
#include "parsers/parser_top_fusion.hpp"
#include "parsers/parser_star_fusion.hpp"
#include <boost/algorithm/string/replace.hpp>

namespace Anaquin
{
    enum FusionCaller
    {
        StarFusion,
        TopHatFusion,
    };

    struct FUSQuin
    {
        static SequinID normToFusion(const SequinID &id)
        {
            auto x = id;
            
            // Eg: NG1_1_P1 to FG1_1_P1
            boost::replace_all(x, "NG", "FG");
            
            return x;
        }

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

            switch (o.soft)
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
//                    ParserTopFusion::parse(Reader(file), [&](const CalledFusion &t, const ParserProgress &)
//                    {
//                        f(classifyFusion(t, o.fuzzy));
//                    });

                    break;
                }
            }
        }
    };
}

#endif