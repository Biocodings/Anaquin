#ifndef T_REPLICATE_HPP
#define T_REPLICATE_HPP

#include "stats/analyzer.hpp"
#include "TransQuin/t_express.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    struct TReplicate : public Analyzer
    {
        typedef TExpress::TestData TestData;
        
        typedef TExpress::Metrics  Metrics;
        typedef TExpress::Options  Options;
        typedef TExpress::Software Software;
        
        typedef TExpress::Stats Stats;
        
        static std::vector<Stats> analyze(const std::vector<FileName> &, const Options &o);
        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);

        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
