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
            Counts nAmbig   = 0;
            Counts nHang    = 0;
            Counts nCross   = 0;
            Counts nSingle  = 0;
            Counts nPaired  = 0;
            Counts nReverse = 0;

            Proportion pHang;
            Proportion pAmbig;
            Proportion pCross;
            Proportion pPaired;
            Proportion pSingle;
            Proportion pReverse;
        };

        struct Impl
        {
            virtual bool isReverse(const ChrID &) = 0;

            // Both mapped
            virtual void paired(const ParserSAM::Data &, const ParserSAM::Data &) = 0;

            // Mapped to reverse and the other not mapped
            virtual void ambig(const ParserSAM::Data &, const ParserSAM::Data &) = 0;
            
            // Mapped to both forward and reverse
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
