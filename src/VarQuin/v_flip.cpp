#include <algorithm>
#include "data/biology.hpp"
#include "tools/errors.hpp"
#include "VarQuin/v_flip.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName HANG_1    = "VarFlip_hanging.fq";
static const FileName FLIPPED_1 = "VarFlip_flipped_1.fq";
static const FileName FLIPPED_2 = "VarFlip_flipped_2.fq";
static const FileName AMBIG_1   = "VarFlip_ambiguous_1.fq";
static const FileName AMBIG_2   = "VarFlip_ambiguous_2.fq";

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

        if (!x.mapped)
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
                         "-------Alignments\n\n"
                         "       Unmapped: %2% (%3%%%)\n"
                         "       Forward:  %4% (%5%%%)\n"
                         "       Reverse:  %6% (%7%%%)\n"
                         "       Dilution: %8$.4f\n\n"
                         "-------VarFlip Outputs\n\n"
                         "       Flipped reads:   %9% (%10%%%)\n"
                         "       Ambiguous reads: %11% (%12%%%)\n"
                         "       Hanging reads:   %13% (%14%%%)\n";

    #define C(x) stats.counts.at(x)

    typedef VFlip::Status S;
    
    const auto cf = C(S::ReverseReverse) + C(S::ReverseNotMapped);
    const auto ca = C(S::ForwardForward) + C(S::ForwardReverse) + C(S::ForwardNotMapped) + C(S::NotMappedNotMapped);
    const auto ch = C(S::RevHang) + C(S::ForHang);
    const auto pf = 100.0 * cf / (cf + ca + ch);
    const auto pa = 100.0 * ca / (cf + ca + ch);
    const auto ph = 100.0 * ch / (cf + ca + ch);

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % align            // 1
                                            % stats.nNA        // 2
                                            % stats.propNA()   // 3
                                            % stats.nEndo      // 4
                                            % stats.propGen()  // 5
                                            % stats.nSeqs      // 6
                                            % stats.propSyn()  // 7
                                            % stats.dilution() // 8
                                            % cf               // 9
                                            % pf               // 10
                                            % ca               // 11
                                            % pa               // 12
                                            % ch               // 13
                                            % ph               // 14
                     ).str());
}
    
void VFlip::report(const FileName &file, const Options &o)
{
    struct Impl : public VFlip::Impl
    {
        Impl(const Options &o)
        {
            h1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            f1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            f2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));

            h1->open(HANG_1);
            a1->open(AMBIG_1);
            a2->open(AMBIG_2);
            f1->open(FLIPPED_1);
            f2->open(FLIPPED_2);
        }

        ~Impl()
        {
            h1->close();
            a1->close();
            a2->close();
            f1->close();
            f2->close();
        }
        
        bool isReverse(const ChrID &cID) { return isReverseChr(cID); }

        void process(const ParserSAM::Data &x, const ParserSAM::Data &y, VFlip::Status status)
        {
            auto writePaired = [&](std::shared_ptr<FileWriter> p1, std::shared_ptr<FileWriter> p2)
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
            
            auto writeSingle = [&](std::shared_ptr<FileWriter> p1)
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
                case Status::ReverseNotMapped:
                {
                    writePaired(f1, f2);
                    break;
                }

                case Status::ForwardReverse:
                case Status::ForwardForward:
                case Status::ForwardNotMapped:
                case Status::NotMappedNotMapped:
                {
                    writePaired(a1, a2);
                    break;
                }

                case Status::RevHang:
                case Status::ForHang:
                {
                    writeSingle(h1);
                    break;
                }
            }
        }

        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> f1, f2;
    };
    
    Impl impl(o);

    const auto stats = analyze(file, o, impl);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
