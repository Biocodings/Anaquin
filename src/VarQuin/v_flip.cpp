#include <thread>
#include <algorithm>
#include "tools/errors.hpp"
#include "VarQuin/v_flip.hpp"
#include "VarQuin/v_split.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

static const FileName HANG_1    = "VarFlip_hanging.fq";
static const FileName FLIPPED_1 = "VarFlip_flipped_1.fq";
static const FileName FLIPPED_2 = "VarFlip_flipped_2.fq";
static const FileName AMBIG_1   = "VarFlip_ambiguous_1.fq";
static const FileName AMBIG_2   = "VarFlip_ambiguous_2.fq";

// Multi-threading...
static VFlip::Impl *__impl__;

// Multi-threading...
VFlip::Stats __stats__;

static void VarSplit(const FileName &file, const VFlip::Options &o)
{
    VSplit::Options o2;
    
    o2.work   = o.work;
    o2.logger = o.logger;
    
    // Generate derived alignment files
    VSplit::report(file, o2);
}

static void VarFlip(const FileName &file, const VFlip::Options &o)
{
    typedef VFlip::Stats Stats;
    typedef VFlip::Status Status;
    
    __stats__.counts[Status::RevHang]            = 0;
    __stats__.counts[Status::ForHang]            = 0;
    __stats__.counts[Status::ReverseReverse]     = 0;
    __stats__.counts[Status::ForwardForward]     = 0;
    __stats__.counts[Status::ForwardReverse]     = 0;
    __stats__.counts[Status::ReverseNotMapped]   = 0;
    __stats__.counts[Status::ForwardNotMapped]   = 0;
    __stats__.counts[Status::NotMappedNotMapped] = 0;
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;
    
    o.logInfo("Parsing: " + o.work + "/VarFlip_sequins.bam");
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        if (!x.mapped)
        {
            __stats__.nNA++;
        }
        else if (__impl__->isReverse(x.cID))
        {
            __stats__.nSeqs++; // Reverse genome
        }
        else
        {
            __stats__.nEndo++; // Forward genome
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
            
            if (__impl__->isReverse(first->cID))
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
            
            if (__impl__->isReverse(second->cID))
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
            
            const auto bothRev =  __impl__->isReverse(first->cID) && __impl__->isReverse(second->cID);
            const auto bothFor = !__impl__->isReverse(first->cID) && !__impl__->isReverse(second->cID);
            const auto anyRev  =  __impl__->isReverse(first->cID) || __impl__->isReverse(second->cID);
            const auto anyFor  = !__impl__->isReverse(first->cID) || !__impl__->isReverse(second->cID);
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
            
            __stats__.counts[status]++;
            __impl__->process(*first, *second, status);
            seenMates.erase(x.name);
        }
    }, true);
    
    o.logInfo("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);
        
        if (__impl__->isReverse(i.second.cID))
        {
            __stats__.counts[Status::RevHang]++;
            __impl__->process(i.second, i.second, Status::RevHang);
        }
        else
        {
            __stats__.counts[Status::ForHang]++;
            __impl__->process(i.second, i.second, Status::ForHang);
        }
    }
}

VFlip::Stats VFlip::analyze(const FileName &file, const Options &o)
{
    std::thread t1(VarFlip, file, o);
    std::thread t2(VarSplit, file, o);
    
    t1.join();
    t2.join();

    // Generated by flipping
    return __stats__;
}

static void writeSummary(const FileName &file,
                         const FileName &align,
                         const VFlip::Stats &stats,
                         const VFlip::Options &o)
{
    const auto summary = "-------VarFlip Summary Statistics\n\n"
                         "-------VarFlip Inputs\n\n"
                         "       Alignment file: %1%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %2% (%3$.2f%%)\n"
                         "       Forward:  %4% (%5$.2f%%)\n"
                         "       Reverse:  %6% (%7$.2f%%)\n\n"
                         "-------VarFlip Outputs\n\n"
                         "       Flipped reads:   %8% (%9$.2f%%)\n"
                         "       Ambiguous reads: %10% (%11$.2f%%)\n"
                         "       Hanging reads:   %12% (%13$.2f%%)\n";

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
    o.writer->write((boost::format(summary) % align         // 1
                                            % stats.nNA     // 2
                                            % stats.pNA()   // 3
                                            % stats.nEndo   // 4
                                            % stats.pEndo() // 5
                                            % stats.nSeqs   // 6
                                            % stats.pSyn()  // 7
                                            % cf            // 8
                                            % pf            // 9
                                            % ca            // 10
                                            % pa            // 11
                                            % ch            // 12
                                            % ph            // 13
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
        
        bool isReverse(const ChrID &cID) { return isRevChr(cID); }

        void process(const ParserBAM::Data &x, const ParserBAM::Data &y, VFlip::Status status)
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
                //case Status::ForwardForward:
                //case Status::ForwardNotMapped:
                //case Status::NotMappedNotMapped:
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
                    
                default: { break; }
            }
        }

        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> f1, f2;
    };
    
    Impl impl(o);
    __impl__ = &impl;

    const auto stats = analyze(file, o);
    
    /*
     * Generating VarFlip_summary.stats
     */

    writeSummary("VarFlip_summary.stats", file, stats, o);
}
