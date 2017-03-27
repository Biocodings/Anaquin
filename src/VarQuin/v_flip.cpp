#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName Paired_1    = "VarFlip_paired_1.fq";
static const FileName Paired_2    = "VarFlip_paired_2.fq";
static const FileName Crossed_1   = "VarFlip_crossed_1.fq";
static const FileName Crossed_2   = "VarFlip_crossed_2.fq";
static const FileName HangingFile = "VarFlip_hanging.fq";
static const FileName SingleFile  = "VarFlip_single.fq";

VFlip::Stats VFlip::analyze(const FileName &file, const Options &o, Impl &impl)
{
    Stats stats;

    // Required for pooling paired-end reads
    std::map<ReadName, ParserSAM::Data> seenMates;
    
    ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (impl.isReverse(x.cID))
        {
            stats.nSyn++;
            
            if (!x.isPassed || x.isSecondary || x.isSupplement)
            {
                return;
            }
            
            if (x.isPaired)
            {
                stats.nPaired++;
                
                if (!seenMates.count(x.name))
                {
                    seenMates[x.name] = x;
                }
                else
                {
                    auto &seen  = seenMates[x.name];
                    auto first  = seen.isFirstPair ? &seen : &x;
                    auto second = seen.isFirstPair ? &x : &seen;

                    if (first->isForward)
                    {
                        complement(first->seq);
                    }
                    else
                    {
                        std::reverse(first->seq.begin(), first->seq.end());
                    }
                    
                    if (second->isForward)
                    {
                        complement(second->seq);
                    }
                    else
                    {
                        std::reverse(second->seq.begin(), second->seq.end());
                    }
                    
                    /*
                     * Cross-alignment occurs when paired-end read mapped to both genomes.
                     */
                    
                    if (!impl.isReverse(first->cID) || !impl.isReverse(second->cID))
                    {
                        stats.nCross++;
                        impl.cross(*first, *second);
                    }
                    else
                    {
                        impl.paired(*first, *second);
                    }

                    seenMates.erase(x.name);
                }
            }
            else
            {
                stats.nSingle++;
                
                // Compute the complement (but not reverse)
                complement(x.seq);

                impl.single(x);
            }
        }
        else if (!x.mapped)
        {
            stats.nNA++;
        }
        else
        {
            stats.nGen++;
        }
    }, true);

    o.info("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);

        impl.hanging(i.second);
    }

    stats.nHang = seenMates.size();
    
    const auto total = stats.nPaired + stats.nSingle + stats.nHang + stats.nCross;
    
    stats.pHang   = static_cast<Proportion>(stats.nHang)   / total;
    stats.pCross  = static_cast<Proportion>(stats.nCross)  / total;
    stats.pPaired = static_cast<Proportion>(stats.nPaired) / total;
    stats.pSingle = static_cast<Proportion>(stats.nSingle) / total;

    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &align,
                         const VFlip::Stats &stats,
                         const VFlip::Options &o)
{
    const auto summary = "-------VarFlip Output Results\n\n"
                         "-------VarFlip Inputs\n\n"
                         "       Alignment file: %1%\n\n"
                         "-------VarFlip Outputs\n\n"
                         "       Paired-end alignments: %2%\n"
                         "                              %3%\n"
                         "       Crossed alignments:    %4%\n"
                         "                              %5%\n"
                         "       Hanging alignments:    %6%\n\n"
                         "       Single alignments:     %7%\n\n"
                         "-------Alignments\n\n"
                         "       Paired:   %8% (%9%%%)\n"
                         "       Crossed:  %10% (%11%%%)\n"
                         "       Hanging:  %12% (%13%%%)\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %14% (%15%%%)\n"
                         "       Forward:  %16% (%17%%%)\n"
                         "       Reverse:  %18% (%19%%%)\n"
                         "       Dilution: %20$.4f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % Paired_1         // 2
                                            % Paired_2         // 3
                                            % Crossed_1        // 4
                                            % Crossed_2        // 5
                                            % HangingFile      // 6
                                            % SingleFile       // 7
                                            % stats.nPaired    // 8
                                            % stats.pPaired    // 9
                                            % stats.nCross     // 10
                                            % stats.pCross     // 11
                                            % stats.nHang      // 12
                                            % stats.pHang      // 13
                                            % stats.nNA        // 14
                                            % stats.propNA()   // 15
                                            % stats.nGen       // 16
                                            % stats.propGen()  // 17
                                            % stats.nSyn       // 18
                                            % stats.propSyn()  // 19
                                            % stats.dilution() // 20
                     ).str());
}
    
void VFlip::report(const FileName &file, const Options &o)
{
    struct Impl : public VFlip::Impl
    {
        Impl(const Options &o)
        {
            p1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            p2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            hg = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            c1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            c2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            sing = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            
            p1->open(Paired_1);
            p2->open(Paired_2);
            c1->open(Crossed_1);
            c2->open(Crossed_2);
            hg->open(HangingFile);
            sing->open(SingleFile);
        }

        ~Impl()
        {
            c1->close();
            c2->close();
            p1->close();
            p2->close();
            hg->close();
            sing->close();
        }
        
        bool isReverse(const ChrID &cID)
        {
            return isReverseGenome(cID);
        }

        void paired(const ParserSAM::Data &x, const ParserSAM::Data &y)
        {
            p1->write("@" + x.name + "/1");
            p1->write(x.seq);
            p1->write("+");
            p1->write(x.qual);
            p2->write("@" + y.name + "/2");
            p2->write(y.seq);
            p2->write("+");
            p2->write(y.qual);
        }

        void cross(const ParserSAM::Data &x, const ParserSAM::Data &y)
        {
            c1->write("@" + x.name + "/1");
            c1->write(x.seq);
            c1->write("+");
            c1->write(x.qual);
            c2->write("@" + y.name + "/2");
            c2->write(y.seq);
            c2->write("+");
            c2->write(y.qual);
        }
        
        void single(const ParserSAM::Data &x)
        {
            if (isReverseGenome(x.cID))
            {
                sing->write("@" + x.name);
                sing->write(x.seq);
                sing->write("+");
                sing->write(x.qual);
            }
        }

        void hanging(const ParserSAM::Data &x)
        {
            if (x.isFirstPair)
            {
                hg->write("@" + x.name + "/1");
            }
            else
            {
                hg->write("@" + x.name + "/2");
            }
            
            hg->write(x.seq);
            hg->write("+");
            hg->write(x.qual);
        }
        
        std::shared_ptr<FileWriter> p1, p2, hg, c1, c2, sing;
    };
    
    Impl impl(o);

    const auto stats = analyze(file, o, impl);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
