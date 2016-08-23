/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include <thread>
#include <fstream>
#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "parsers/parser_fq.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

#ifdef ANAQUIN_DEBUG
#include <mutex>
#endif

using namespace Anaquin;

static const FileName revGenome_1 = "VarFlip_sequins_1.fq";
static const FileName revGenome_2 = "VarFlip_sequins_2.fq";
static const FileName forGenome_1 = "VarFlip_genome_1.fq";
static const FileName forGenome_2 = "VarFlip_genome_2.fq";

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
            if (x.mapq <= 10)
            {
                stats.countLow++;

#ifdef ANAQUIN_DEBUG
                o.logInfo("Mapping: " + x.name + "\t" + x.cID);
#endif
            }
            else
            {
                stats.countSyn++;
                names.insert(x.name);

#ifdef ANAQUIN_DEBUG
                o.logInfo("Added: " + x.name + "\t" + x.cID);
#endif
            }
        }
        else
        {
            stats.countGen++;
        }
    });

#ifdef ANAQUIN_DEBUG
    std::mutex m;
#endif
    
    auto flip = [&](const FileName &f, const FileName &forw, const FileName &rev)
    {
        FileWriter ff(o.work);
        FileWriter rf(o.work);
 
        ff.open(forw);
        rf.open(rev);
        
        std::string n1, n2;
        
        ParserFQ::parse(Reader(f), [&](ParserFQ::Data &x, const ParserProgress &)
        {
            n1 = x.name;

            if (n1.length() >= 2 &&
                ((n1[n1.length()-1] == '1' && n1[n1.length()-2] == '/') ||
                 (n1[n1.length()-1] == '2' && n1[n1.length()-2] == '/')))
            {
                n1 = n1.substr(0, n1.size()-2);
            }
            
            if (n1[0] == '@')
            {
                n2 = n1.substr(1, n1.size()-1);
            }

            if (names.count(n1) || names.count(n2))
            {
#ifdef ANAQUIN_DEBUG
                m.lock();
                o.logInfo("Found: " + x.name + "\t" + std::to_string(names.count(n1)) + "\t" + std::to_string(names.count(n2)) + "\t" + n1 + "\t" + n2);
                m.unlock();
#endif

                // DNA complement of the sequence
                complement(x.seq);
                
                rf.write(x.name + " " + x.info);
                rf.write(x.seq);
                rf.write(x.opt);
                rf.write(x.qual);
            }
            else
            {
                ff.write(x.name + " " + x.info);
                ff.write(x.seq);
                ff.write(x.opt);
                ff.write(x.qual);
            }
        });
        
        ff.close();
        rf.close();
    };

    o.info("Generating " + forGenome_1 + " & " + revGenome_1);
    o.info("Generating " + forGenome_2 + " & " + revGenome_2);

    std::thread t1(flip, seq1, forGenome_1, revGenome_1);
    std::thread t2(flip, seq2, forGenome_2, revGenome_2);
    
    t1.join();
    t2.join();

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
                         "       Forward reads: %4%\n"
                         "                      %5%\n\n"
                         "       Flipped reverse reads: %6%\n"
                         "                              %7%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %8% (%9%%%)\n"
                         "       Forward:  %10% (%11%%%)\n"
                         "       Reverse:  %12% (%13%%%)\n"
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
                 forGenome_1,
                 forGenome_2,
                 revGenome_1,
                 revGenome_2,
                 stats,
                 o);
}