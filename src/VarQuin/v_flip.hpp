#ifndef V_FLIP_HPP
#define V_FLIP_HPP

#include <map>
#include "tools/tools.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct VFlip
    {
        typedef AnalyzerOptions Options;
        
        enum class Status
        {
            ReverseReverse,
            ReverseNotMapped,
            ForwardForward,
            ForwardReverse,
            ForwardNotMapped,
            NotMappedNotMapped,
            RevHang,
            ForHang
        };
        
        struct Stats : public MappingStats
        {
            std::map<VFlip::Status, Counts> counts;
            
            inline Proportion prop(Status status) const
            {
                return ((Proportion) counts.at(status)) / sum(counts);
            }
        };

        struct Impl
        {
            virtual bool isReverse(const ChrID &) = 0;
            virtual void process(const ParserSAM::Data &, const ParserSAM::Data &, Status) = 0;
        };

        static Stats analyze(const FileName &, const Options &, Impl &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
