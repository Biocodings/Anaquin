#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName AMBIG_1     = "VarFlip_ambig_1.fq";
static const FileName AMBIG_2     = "VarFlip_ambig_2.fq";
static const FileName Paired_1    = "VarFlip_paired_1.fq";
static const FileName Paired_2    = "VarFlip_paired_2.fq";
static const FileName Crossed_1   = "VarFlip_crossed_1.fq";
static const FileName Crossed_2   = "VarFlip_crossed_2.fq";
static const FileName HangingFile = "VarFlip_hanging.fq";

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
            stats.nSeqs++;
            
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
                    
                    // Anything mapped to the reverse?
                    const auto anyReverse = impl.isReverse(first->cID) || impl.isReverse(second->cID);
                    
                    // Anything mapped to the forward?
                    const auto anyForward = !impl.isReverse(first->cID) || !impl.isReverse(second->cID);
                    
                    // Anything not mapped?
                    const auto anyNMapped = !first->mapped || !second->mapped;
                    
                    // Crossed alignment?
                    const auto crossed = !anyNMapped && anyReverse && anyForward;
                    
                    // Ambigious alignment?
                    const auto ambig = anyReverse && !anyForward && anyNMapped;
                    
                    if (ambig)
                    {
                        stats.nAmbig++;
                        impl.ambig(*first, *second);
                    }
                    else if (crossed)
                    {
                        stats.nCross++;
                        impl.cross(*first, *second);
                    }
                    else
                    {
                        stats.nReverse++;
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

                // Single-ended alignments
                impl.single(x);
            }
        }
        else if (!x.mapped)
        {
            stats.nNA++;
        }
        else
        {
            stats.nEndo++;
        }
    }, true);

    o.logInfo("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);

        // Alignments without the other mate
        impl.hanging(i.second);
    }

    stats.nHang = seenMates.size();

    stats.pSingle = static_cast<Proportion>(stats.nSingle) / (stats.nSingle + stats.nPaired);
    stats.pPaired = static_cast<Proportion>(stats.nPaired) / (stats.nSingle + stats.nPaired);
    
    const auto total = stats.nAmbig + stats.nReverse + stats.nHang + stats.nCross;
    
    stats.pHang    = static_cast<Proportion>(stats.nHang)    / total;
    stats.pAmbig   = static_cast<Proportion>(stats.nAmbig)   / total;
    stats.pCross   = static_cast<Proportion>(stats.nCross)   / total;
    stats.pReverse = static_cast<Proportion>(stats.nReverse) / total;

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
                         "       Ambiguous alignments:  %6%\n"
                         "                              %7%\n"
                         "       Hanging alignments:    %8%\n\n"
                         "-------Alignments\n\n"
                         "       Paired:   %9% (%10%%%)\n"
                         "       Crossed:  %11% (%12%%%)\n"
                         "       Hanging:  %13% (%14%%%)\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %15% (%16%%%)\n"
                         "       Forward:  %17% (%18%%%)\n"
                         "       Reverse:  %19% (%20%%%)\n"
                         "       Dilution: %21$.4f\n";

    #define P(x) (isnan(x) ? "0" : (std::to_string(x)))
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % Paired_1         // 2
                                            % Paired_2         // 3
                                            % Crossed_1        // 4
                                            % Crossed_2        // 5
                                            % HangingFile      // 6
                                            % AMBIG_1          // 7
                                            % AMBIG_2          // 8
                                            % stats.nPaired    // 9
                                            % P(stats.pPaired) // 10
                                            % stats.nCross     // 11
                                            % P(stats.pCross)  // 12
                                            % stats.nHang      // 13
                                            % P(stats.pHang)   // 14
                                            % stats.nNA        // 15
                                            % stats.propNA()   // 16
                                            % stats.nEndo       // 17
                                            % stats.propGen()  // 18
                                            % stats.nSeqs       // 19
                                            % stats.propSyn()  // 20
                                            % stats.dilution() // 21
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
            a1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            
            a1->open(AMBIG_1);
            a2->open(AMBIG_2);
            p1->open(Paired_1);
            p2->open(Paired_2);
            c1->open(Crossed_1);
            c2->open(Crossed_2);
            hg->open(HangingFile);
        }

        ~Impl()
        {
            a1->close();
            a2->close();
            c1->close();
            c2->close();
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
        
        void ambig(const ParserSAM::Data &x, const ParserSAM::Data &y)
        {
            a1->write("@" + x.name + "/1");
            a1->write(x.seq);
            a1->write("+");
            a1->write(x.qual);
            a2->write("@" + y.name + "/2");
            a2->write(y.seq);
            a2->write("+");
            a2->write(y.qual);
        }

        void single(const ParserSAM::Data &x) {}

        void hanging(const ParserSAM::Data &x)
        {
            if (isReverse(x.cID) && x.mapped)
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
        }
        
        std::shared_ptr<FileWriter> p1, p2, hg, c1, c2, a1, a2;
    };
    
    Impl impl(o);

    const auto stats = analyze(file, o, impl);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
