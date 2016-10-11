#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName SINGLE   = "VarFlip_single.fq";
static const FileName HANGING  = "VarFlip_hanging.fq";
static const FileName PAIRED_1 = "VarFlip_paired_1.fq";
static const FileName PAIRED_2 = "VarFlip_paired_2.fq";

VFlip::Stats VFlip::analyze(const FileName &align, const Options &o, Impl &impl)
{
    Stats stats;

    // Required for pooling paired-end reads
    std::map<ReadName, ParserSAM::Data> seenMates;
    
    ParserSAM::parse(align, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (impl.isReverse(x.cID))
        {
            stats.countSyn++;
            
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
                    
                    impl.paired(*first, *second);
                    seenMates.erase(x.name);
                }
            }
            else
            {
                stats.nSingle++;
                
                // Compute the complement (but not reverse)
                complement(x.seq);

                impl.nonPaired(x);
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

    o.info("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);

        impl.unknownPaired(i.second);
    }

    stats.nHanging = seenMates.size();
    
    const auto total = stats.nPaired + stats.nSingle + stats.nHanging;
    
    stats.pPaired  = static_cast<Proportion>(stats.nPaired)  / total;
    stats.pSingle  = static_cast<Proportion>(stats.nSingle)  / total;
    stats.pHanging = static_cast<Proportion>(stats.nHanging) / total;

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
                         "       Hanging alignments: %4%\n"
                         "       Single-end alignments: %5%\n\n"
                         "-------Alignments\n\n"
                         "       Paired:   %6% (%7%%%)\n"
                         "       Hanging:  %8% (%9%%%)\n"
                         "       Unpaired: %10% (%11%%%)\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %12% (%13%%%)\n"
                         "       Forward:  %14% (%15%%%)\n"
                         "       Reverse:  %16% (%17%%%)\n"
                         "       Dilution: %18$.4f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % PAIRED_1         // 2
                                            % PAIRED_2         // 3
                                            % HANGING          // 4
                                            % SINGLE           // 5
                                            % stats.nPaired    // 6
                                            % stats.pPaired    // 7
                                            % stats.nHanging   // 8
                                            % stats.pHanging   // 9
                                            % stats.nSingle    // 10
                                            % stats.pSingle    // 11
                                            % stats.countNA    // 12
                                            % stats.propNA()   // 13
                                            % stats.countGen   // 14
                                            % stats.propGen()  // 15
                                            % stats.countSyn   // 16
                                            % stats.propSyn()  // 17
                                            % stats.dilution() // 18
                     ).str());
}
    
void VFlip::report(const FileName &file, const Options &o)
{
    struct Impl : public VFlip::Impl
    {
        Impl(const Options &o)
        {
            up = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            p1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            p2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            hg = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            
            up->open(SINGLE);
            hg->open(HANGING);
            p1->open(PAIRED_1);
            p2->open(PAIRED_2);
        }

        ~Impl()
        {
            up->close();
            p1->close();
            p2->close();
            hg->close();
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

        void nonPaired(const ParserSAM::Data &x)
        {
            up->write("@" + x.name + "/1");
            up->write(x.seq);
            up->write("+");
            up->write(x.qual);
        }

        void unknownPaired(const ParserSAM::Data &x)
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
        
        std::shared_ptr<FileWriter> p1, p2, up, hg;
    };
    
    Impl impl(o);

    const auto stats = analyze(file, o, impl);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
