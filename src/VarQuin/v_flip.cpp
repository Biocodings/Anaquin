/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include <fstream>
#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "parsers/parser_fq.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

static const FileName revGenome_1  = "VarFlip_sequins_1.fq";
static const FileName revGenome_2  = "VarFlip_sequins_2.fq";
static const FileName forwGenome_1 = "VarFlip_genome_1.fq";
static const FileName forwGenome_2 = "VarFlip_genome_2.fq";

VFlip::Stats VFlip::analyze(const FileName &seq1,
                            const FileName &seq2,
                            const FileName &align,
                            const Options  &o)
{
    Stats stats;
    
    std::set<ReadName> names;

    ParserSAM::parse(align, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (!x.mapped)
        {
            stats.countNA++;
        }
        else if (Standard::isSynthetic(x.cID))
        {
            stats.countSyn++;
            
            /*
             * Eg: @GV_IDEL_011863_R-1400/1
             * Eg: @GV_IDEL_011863_R-1400/2
             */
            
            names.insert("@" + x.name + "/1");
            names.insert("@" + x.name + "/2");
        }
        else
        {
            stats.countGen++;
        }
    });

    auto f = [&](const FileName &f, const FileName &forw, const FileName &rev)
    {
        std::ofstream fs, fg;
 
        fs.open(rev,  std::ios_base::app);
        fg.open(forw, std::ios_base::app);
        
        ParserFQ::parse(Reader(f), [&](ParserFQ::Data &x, const ParserProgress &)
        {
            if (names.count(x.name))
            {
                // DNA complement of the sequence
                complement(x.seq);
                
                fs << x.name << " " << x.info << std::endl;
                fs << x.seq  << std::endl;
                fs << x.opt  << std::endl;
                fs << x.qual << std::endl;
            }
            else
            {
                fg << x.name << " " << x.info << std::endl;
                fg << x.seq  << std::endl;
                fg << x.opt  << std::endl;
                fg << x.qual << std::endl;
            }
        });
        
        fs.close();
        fg.close();
    };

    o.info("Generating " + forwGenome_1 + " & " + revGenome_1);
    f(seq1, o.work + "/" + forwGenome_1, o.work + "/" + revGenome_1);

    o.info("Generating " + forwGenome_2 + " & " + revGenome_2);
    f(seq2, o.work + "/" + forwGenome_2, o.work + "/" + revGenome_2);

    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &seq1,
                         const FileName &seq2,
                         const FileName &align,
                         const FileName &forw1,
                         const FileName &forw2,
                         const FileName &frev1,
                         const FileName &frev2,
                         const VFlip::Stats &stats,
                         const VFlip::Options &o)
{
    const auto summary = "-------VarFlip Output Results\n\n"
                         "-------VarFlip Inputs\n\n"
                         "       Input sequence files: %1%\n"
                         "                             %2%\n"
                         "       Alignment file:       %3%\n\n"
                         "-------VarFlip Outputs\n\n"
                         "       Forward genome reads: %4%\n"
                         "                             %5%\n\n"
                         "       Flipped reverse genome reads: %6%\n"
                         "                                     %7%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %8% (%9%)\n"
                         "       Forward:  %10% (%11%)\n"
                         "       Reverse:  %12% (%13%)\n"
                         "       Dilution: %14$.2f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % seq1             // 1
                                            % seq2             // 2
                                            % align            // 3
                                            % forw1            // 4
                                            % forw2            // 5
                                            % frev1            // 6
                                            % frev2            // 7
                                            % stats.countNA    // 8
                                            % stats.propNA()   // 9
                                            % stats.countGen   // 10
                                            % stats.propGen()  // 11
                                            % stats.countSyn   // 12
                                            % stats.propSyn()  // 13
                                            % stats.dilution() // 14
                     ).str());
}
    
void VFlip::report(const std::vector<FileName> &files, const Options &o)
{
    const auto seq1  = files[0];
    const auto seq2  = files[1];
    const auto align = files[2];
    
    const auto stats = analyze(seq1, seq2, align, o);
    
    /*
     * Generating VarFlip_summary.stats
     */
    
    writeSummary("VarFlip_summary.stats",
                 seq1,
                 seq2,
                 align,
                 forwGenome_1,
                 forwGenome_2,
                 revGenome_1,
                 revGenome_2,
                 stats,
                 o);
}