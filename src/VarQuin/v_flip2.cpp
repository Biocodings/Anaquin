#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip2.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName o1 = "VarFlip_sequins_1.fq";
static const FileName o2 = "VarFlip_sequins_2.fq";

VFlip2::Stats VFlip2::analyze(const FileName &align, const Options &o)
{
    Stats stats;

    FileWriter f1(o.work);
    FileWriter f2(o.work);

    f1.open(o1);
    f2.open(o2);

    // Required for pooling paired-end reads
    std::map<ReadName, ParserSAM::Data> seenMates;
    
    ParserSAM::parse(align, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        if (isReverseGenome(x.cID))
        {
            stats.countSyn++;
            
            if (!x.isPassed || x.isSecondary)
            {
                return;
            }
            else if (!seenMates.count(x.name))
            {
                seenMates[x.name] = x;
            }
            else
            {
                const auto &seen  = seenMates[x.name];
                const auto first  = seen.isFirstPair ? &seen : &x;
                const auto second = seen.isFirstPair ? &x : &seen;
                
                f1.write("@" + first->name + "/1");
                f1.write(first->seq);
                f1.write("+");
                f1.write(first->qual);

                f2.write("@" + second->name + "/2");
                f2.write(second->seq);
                f2.write("+");
                f2.write(second->qual);
                
                seenMates.erase(x.name);
            }
        }
        else if (!x.mapped)
        {
            stats.countNA++;
        }
        else
        {
            stats.countGen++;
        }
    }, true);

    f1.close();
    f2.close();
    
    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &align,
                         const FileName &f1,
                         const FileName &f2,
                         const VFlip2::Stats &stats,
                         const VFlip2::Options &o)
{
    const auto summary = "-------VarFlip Output Results\n\n"
                         "-------VarFlip Inputs\n\n"
                         "       Alignment file: %1%\n\n"
                         "-------VarFlip Outputs\n\n"
                         "       Sequin alignments: %2%\n"
                         "                          %3%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %4% (%5%%%)\n"
                         "       Forward:  %6% (%7%%%)\n"
                         "       Reverse:  %8% (%9%%%)\n"
                         "       Dilution: %10$.4f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % f1               // 2
                                            % f2               // 3
                                            % stats.countNA    // 4
                                            % stats.propNA()   // 5
                                            % stats.countGen   // 6
                                            % stats.propGen()  // 7
                                            % stats.countSyn   // 8
                                            % stats.propSyn()  // 9
                                            % stats.dilution() // 10
                     ).str());
}
    
void VFlip2::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats",
                 file,
                 o1,
                 o2,
                 stats,
                 o);
}
