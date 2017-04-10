#include <algorithm>
#include "data/biology.hpp"
#include "tools/errors.hpp"
#include "VarQuin/v_flip.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName NMNM_1    = "VarFlip_NMapNMap_1.fq";
static const FileName NMNM_2    = "VarFlip_NMapNMap_2.fq";
static const FileName ForNM_1   = "VarFlip_ForNMap_1.fq";
static const FileName ForNM_2   = "VarFlip_ForNMap_2.fq";
static const FileName RevNM_1   = "VarFlip_RevNMap_1.fq";
static const FileName RevNM_2   = "VarFlip_RevNMap_2.fq";
static const FileName RevRev_1  = "VarFlip_RevRev_1.fq";
static const FileName RevRev_2  = "VarFlip_RevRev_2.fq";
static const FileName ForFor_1  = "VarFlip_ForFor_1.fq";
static const FileName ForFor_2  = "VarFlip_ForFor_2.fq";
static const FileName ForRev_1  = "VarFlip_ForRev_1.fq";
static const FileName ForRev_2  = "VarFlip_ForRev_2.fq";
static const FileName ForHang_1 = "VarFlip_ForHang.fq";
static const FileName RevHang_1 = "VarFlip_RevHang.fq";

VFlip::Stats VFlip::analyze(const FileName &file, const Options &o, Impl &impl)
{
    Stats stats;

    stats.counts[Status::RevHang]            = 0;
    stats.counts[Status::ForHang]            = 0;
    stats.counts[Status::ReverseReverse]     = 0;
    stats.counts[Status::ForwardForward]     = 0;
    stats.counts[Status::ForwardReverse]     = 0;
    stats.counts[Status::ReverseNotMapped]   = 0;
    stats.counts[Status::ForwardNotMapped]   = 0;
    stats.counts[Status::NotMappedNotMapped] = 0;
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserSAM::Data> seenMates;
    
    ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (x.mapped)
        {
            stats.nNA++;
        }
        else if (impl.isReverse(x.cID))
        {
            stats.nSeqs++;
        }
        else
        {
            stats.nEndo++;
        }

        if (!x.isPassed || x.isSecondary || x.isSupplement)
        {
            return;
        }
        
        A_CHECK(x.isPaired, x.name + " is not pair-ended. Singled-ended not supported.");
        
        if (!seenMates.count(x.name))
        {
            seenMates[x.name] = x;
        }
        else
        {
            auto &seen  = seenMates[x.name];
            auto first  = seen.isFirstPair ? &seen : &x;
            auto second = seen.isFirstPair ? &x : &seen;
            
            /*
             * Only complement reads aligned to the reverse genome
             */
            
            if (impl.isReverse(first->cID))
            {
                if (first->isForward)
                {
                    complement(first->seq);
                }
                else
                {
                    std::reverse(first->seq.begin(), first->seq.end());
                }
            }
            
            if (impl.isReverse(second->cID))
            {
                if (second->isForward)
                {
                    complement(second->seq);
                }
                else
                {
                    std::reverse(second->seq.begin(), second->seq.end());
                }
            }
            
            const auto bothRev =  impl.isReverse(first->cID) && impl.isReverse(second->cID);
            const auto bothFor = !impl.isReverse(first->cID) && !impl.isReverse(second->cID);
            const auto anyRev  =  impl.isReverse(first->cID) || impl.isReverse(second->cID);
            const auto anyFor  = !impl.isReverse(first->cID) || !impl.isReverse(second->cID);
            const auto anyMap  =  first->mapped ||  second->mapped;
            const auto anyNMap = !first->mapped || !second->mapped;
            
            VFlip::Status status;
            
            if (bothRev && !anyNMap)
            {
                status = Status::ReverseReverse;
            }
            else if (bothFor && !anyNMap)
            {
                status = Status::ForwardForward;
            }
            else if (anyRev && anyFor && !anyNMap)
            {
                status = Status::ForwardReverse;
            }
            else if (anyNMap && anyMap && anyRev)
            {
                status = Status::ReverseNotMapped;
            }
            else if (anyNMap && anyMap && anyFor)
            {
                status = Status::ForwardNotMapped;
            }
            else
            {
                status = Status::NotMappedNotMapped;
            }
            
            stats.counts[status]++;
            impl.process(*first, *second, Status::ReverseReverse);
            seenMates.erase(x.name);
        }
    }, true);

    o.logInfo("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);

        if (impl.isReverse(i.second.cID))
        {
            stats.counts[Status::RevHang]++;
            impl.process(i.second, i.second, Status::RevHang);
        }
        else
        {
            stats.counts[Status::ForHang]++;
            impl.process(i.second, i.second, Status::ForHang);
        }
    }

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
                         "       RevRev:   %9% (%10%%%)\n"
                         "       Crossed:  %11% (%12%%%)\n"
                         "       Hanging:  %13% (%14%%%)\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %15% (%16%%%)\n"
                         "       Forward:  %17% (%18%%%)\n"
                         "       Reverse:  %19% (%20%%%)\n"
                         "       Dilution: %21$.4f\n";

    #define C(x) stats.counts.at(x)
    #define P(x) stats.prop(x)

    typedef VFlip::Status S;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % "????"          // 2
                                            % "????"          // 3
                                            % "????"        // 4
                                            % "????"        // 5
                                            % "????"      // 6
                                            % "????"          // 7
                                            % "????"          // 8
                                            % C(S::ReverseReverse) // 9
                                            % P(S::ReverseReverse) // 10
                                            % C(S::ForwardReverse) // 11
                                            % P(S::ForwardReverse) // 12
                                            % C(S::RevHang)        // 13
                                            % P(S::RevHang)    // 14
                                            % stats.nNA        // 15
                                            % stats.propNA()   // 16
                                            % stats.nEndo      // 17
                                            % stats.propGen()  // 18
                                            % stats.nSeqs      // 19
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
            auto init = [&](Status status, FileName f1, FileName f2)
            {
                w1[status] = std::shared_ptr<FileWriter>(new FileWriter(o.work));
                w1[status]->open(f1);
                
                if (status != Status::ForHang && status != Status::RevHang)
                {
                    w2[status] = std::shared_ptr<FileWriter>(new FileWriter(o.work));
                    w2[status]->open(f2);
                }
            };
            
            init(Status::ForHang, ForHang_1, ForHang_1);
            init(Status::RevHang, RevHang_1, RevHang_1);
            init(Status::ReverseReverse, RevRev_1, RevRev_2);
            init(Status::ForwardForward, ForFor_1, ForFor_2);
            init(Status::ForwardReverse, ForRev_1, ForRev_2);
            init(Status::ReverseNotMapped, RevNM_1, RevNM_2);
            init(Status::ForwardNotMapped, ForNM_1, ForNM_2);
            init(Status::NotMappedNotMapped, NMNM_1, NMNM_2);
        }

        ~Impl()
        {
            for (auto &i : w1) { i.second->close(); }
            for (auto &i : w2) { i.second->close(); }
        }
        
        bool isReverse(const ChrID &cID)
        {
            return isReverseGenome(cID);
        }

        void process(const ParserSAM::Data &x, const ParserSAM::Data &y, VFlip::Status status)
        {
            auto p1 = w1.at(status);
            auto p2 = w2.count(status) ? w2.at(status) : nullptr;

            auto writePaired = [&]()
            {
                p1->write("@" + x.name + "/1");
                p1->write(x.seq);
                p1->write("+");
                p1->write(x.qual);
                p2->write("@" + y.name + "/2");
                p2->write(y.seq);
                p2->write("+");
                p2->write(y.qual);
            };
            
            auto writeSingle = [&]()
            {
                if (x.mapped)
                {
                    if (x.isFirstPair)
                    {
                        p1->write("@" + x.name + "/1");
                    }
                    else
                    {
                        p1->write("@" + x.name + "/2");
                    }
                    
                    p1->write(x.seq);
                    p1->write("+");
                    p1->write(x.qual);
                }
            };
            
            switch (status)
            {
                case Status::ReverseReverse:
                case Status::ForwardReverse:
                case Status::ForwardForward:
                case Status::ReverseNotMapped:
                case Status::ForwardNotMapped:
                case Status::NotMappedNotMapped:
                {
                    writePaired();
                    break;
                }

                case Status::RevHang:
                case Status::ForHang:
                {
                    writeSingle();
                    break;
                }
            }
        }
        
        std::map<Status, std::shared_ptr<FileWriter>> w1, w2;
    };
    
    Impl impl(o);

    const auto stats = analyze(file, o, impl);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
