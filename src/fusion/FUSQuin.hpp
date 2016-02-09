#ifndef F_CLASSIFY_HPP
#define F_CLASSIFY_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_top_fusion.hpp"
#include "parsers/parser_star_fusion.hpp"

namespace Anaquin
{
    enum FusionCaller
    {
        Star,
        TopHat,
    };
    
    struct FUSQuin
    {
        enum State
        {
            EndoChrT = -1,
            Negative = 0,
            Positive = 1, // Either chrT or endogenous
        };

        struct Match
        {
            // State of the query
            State state;
            
            CalledFusion query;
            
            /*
             * Only defined if the code is Positive
             */

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

        template <typename T> static Match classifyFusion(const T &f, SequinID &id, double fuzzy)
        {
            Match match;
            match.query = f;
            
            if (f.cID_1 != ChrT || f.cID_2 != ChrT)
            {
                if (f.cID_1 != ChrT && f.cID_2 != ChrT)
                {
                    match.state = Negative; // TODO: What to do with endogenous?
                }
                else
                {
                    match.state = EndoChrT;
                }
            }
            else
            {
                match.state = Negative;

                const auto min = std::min(f.l1, f.l2);
                const auto max = std::max(f.l1, f.l2);
                
                const auto &r = Standard::instance().r_fus;

                // Any match by position and fuzzy?
                if ((match.known = r.findFusion(min, max, f.s1, f.s2, fuzzy)))
                {
                    match.state = Positive;
                }
            }

            return match;
        }

        template <typename Options, typename F> static void analyze(const FileName &file, const Options &o, F f)
        {
            o.info("Fuzzy level: " + std::to_string(o.fuzzy));
            o.info("Parsing alignment file");

//            auto positive = [&](const SequinID &id, Reads reads)
//            {
//                assert(!id.empty() && r.match(id));
//
//                //stats.chrT->n_chrT++;
//                //stats.chrT->h.at(id)++;
//
//                f(id, reads);
//                
////                if (shouldMix)
//                {
//                    // Known abundance for the fusion
////                    const auto known = r.match(id)->abund(Mix_1);
//                    
//                    // Measured abundance for the fusion
//  //                  const auto measured = reads;
//                    
//    //                Point p;
//                    
//      //              p.x = known;
//        //            p.y = measured;
//
//                    //(*(stats.chrT))[id] = p;
//                }
//            };

            SequinID id;

            switch (o.caller)
            {
                case Star:
                {
                    ParserStarFusion::parse(Reader(file), [&](const CalledFusion &t, const ParserProgress &)
                    {
                        f(classifyFusion(t, id, o.fuzzy));
                    });

                    break;
                }
                    
                case TopHat:
                {
                    ParserTopFusion::parse(Reader(file), [&](const CalledFusion &t, const ParserProgress &)
                    {
                        f(classifyFusion(t, id, o.fuzzy));
                    });

                    break;
                }
            }
            
            /*
             * Find out all the sequins undetected in the experiment
             */

            //o.info("Detected " + std::to_string(stats.chrT->h.size()) + " sequins in the reference");
            o.info("Checking for missing sequins");
            
//            for (const auto &i : s.seqIDs)
//            {
//                const auto &seqID = i;
//
//                // If the histogram has an entry of zero
//                if (!stats.h.at(seqID))
//                {
//                    if (!s.seqs_1.count(seqID))
//                    {
//                        o.warn(seqID + " defined in the referene but not in the mixture and it is undetected.");
//                        continue;
//                    }
//
//                    o.warn(seqID + " defined in the referene but not detected");
//
//                    const auto seq = s.seqs_1.at(seqID);
//
//                    // Known abundance for the fusion
//                    const auto known = seq.abund() / seq.length;
//
//                    //stats.y.push_back(0); // TODO: We shouldn't even need to add those missing sequins!
//                    //stats.z.push_back(seqID);
//                    //stats.x.push_back(log2f(known));
//
//                    stats.miss.push_back(MissingSequin(seqID, known));
//                }
//            }

            // The references are simply the known fusion points
            //stats.chrT->m.nr() = r.countFusion();

            o.info("Calculating limit of sensitivity");

//            if (shouldMix)
  //          {
    //            stats.chrT->s = r.limit(stats.chrT->h);
      //      }
        }
    };
}

#endif