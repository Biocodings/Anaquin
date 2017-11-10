#include <fstream>
#include "VarQuin.hpp"
#include "VarQuin/v_split.hpp"
#include "writers/bam_writer.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

typedef VSplit::Stats Stats;

Stats VSplit::analyze(const FileName &file, const Options &o)
{
    Stats stats;
//    BAMWriter w1, w2, w3, w4, w5, w6, w7;
    BAMWriter w2;

//    w1.open(o.work + "/VarSplit_sample.bam");
    w2.open(o.work + "/VarSplit_germline.bam");
//    w3.open(o.work + "/VarSplit_somatic.bam");
//    w4.open(o.work + "/VarSplit_structural.bam");
//    w5.open(o.work + "/VarSplit_ladder.bam");
//    w6.open(o.work + "/VarSplit_ambiguous.bam");
//    w7.open(o.work + "/VarSplit_unmapped.bam");

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(toString(i.p.i));
        }
        
        if (!x.mapped && !x.isMateAligned)
        {
            stats.nNMap++;
           // w7.write(x);
        }
        else if ((x.mapped && !x.isMateAligned) || (!x.mapped &&  x.isMateAligned))
        {
            stats.nAmbig++;
           // w6.write(x);
        }
        else
        {
            const auto isG1 = isGerm(x.cID);
            const auto isG2 = isGerm(x.rnext);
            const auto isC1 = isCancer(x.cID);
            const auto isC2 = isCancer(x.rnext);
            const auto isS1 = isStruct(x.cID);
            const auto isS2 = isStruct(x.rnext);
            const auto isL1 = isLadQuin(x.cID);
            const auto isL2 = isLadQuin(x.rnext);
            const auto isH1 = !isG1 && !isC1 && !isS1 && !isL1;
            const auto isH2 = !isG2 && !isC2 && !isS2 && !isL2;

            if ((isG1 != isG2) || (isC1 != isC2) || (isS1 != isS2) || (isL1 != isL2) || (isH1 != isH2))
            {
                stats.nAmbig++;
         //       w6.write(x);
            }
            else if (isH1 && isH2)
            {
                stats.nSample++;
        //        w1.write(x);
            }
            else if (isG1 && isG2)
            {
                stats.nGerm++;
                w2.write(x);
                w2.close();
                w2.close();
            }
            else if (isC1 && isC2)
            {
                stats.nSom++;
//                w3.write(x);
            }
            else if (isS1 && isS2)
            {
                stats.nStru++;
//                w4.write(x);
            }
            else if (isL1 && isL2)
            {
                stats.nLad++;
       //         w5.write(x);
            }
            else
            {
                throw std::runtime_error(x.name + " is unknown");
            }
        }
    }, true);
    
//    w1.close();
    w2.close();
//    w3.close();
//    w4.close();
//    w5.close();
//    w6.close();
//    w7.close();
//    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const Stats &stats, const VSplit::Options &o)
{
    const auto summary = "-------VarSplit Summary Statistics\n\n"
                         "-------Input file\n\n"
                         "       Input alignment file: %1%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped:  %2%  (%3$.2f%%)\n"
                         "       Sample:    %4%  (%5$.2f%%)\n"
                         "       Germline:  %6%  (%7$.2f%%)\n"
                         "       Somatic:   %8% (%9$.2f%%)\n"
                         "       Structual: %10% (%11$.2f%%)\n"
                         "       Ladders:   %12% (%13$.2f%%)\n"
                         "       Ambiguous: %14% (%15$.2f%%)\n";

    const auto total = stats.nAmbig + stats.nNMap + stats.nGerm + stats.nSom + stats.nSample + stats.nStru + stats.nLad;

    #define S(x) ((float) x / total)
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file             // 1
                                            % stats.nNMap      // 2
                                            % S(stats.nNMap)   // 3
                                            % stats.nSample    // 4
                                            % S(stats.nSample) // 5
                                            % stats.nGerm      // 6
                                            % S(stats.nGerm)   // 7
                                            % stats.nSom       // 8
                                            % S(stats.nSom)    // 8
                                            % stats.nStru      // 10
                                            % S(stats.nStru)   // 11
                                            % stats.nLad       // 12
                                            % S(stats.nLad)    // 13
                                            % stats.nAmbig     // 14
                                            % S(stats.nAmbig)  // 15
                     ).str());
}

void VSplit::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarSplit_summary.stats
     */

    writeSummary("VarSplit_summary.stats", file, stats, o);
}
