#ifndef V_FLIP_2_HPP
#define V_FLIP_2_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct VFlip
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public MappingStats
        {
            Counts nHang   = 0;
            Counts nCross  = 0;
            Counts nSingle = 0;
            Counts nPaired = 0;

            Proportion pHang;
            Proportion pCross;
            Proportion pPaired;
            Proportion pSingle;
        };

        struct Impl
        {
            virtual bool isReverse(const ChrID &) = 0;

            // Paired-end reads both mapped and complemented
            virtual void paired(const ParserSAM::Data &, const ParserSAM::Data &) = 0;

            // Paired-end reads on both forward and reverse genome
            virtual void cross(const ParserSAM::Data &, const ParserSAM::Data &) = 0;

            // Single-end reads
            virtual void single(const ParserSAM::Data &) = 0;

            // Paired-end reads that the other mate not found
            virtual void hanging(const ParserSAM::Data &) = 0;
        };

        static Stats analyze(const FileName &, const Options &, Impl &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
