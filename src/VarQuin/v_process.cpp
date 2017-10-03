#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_process.hpp"
#include "parsers/parser_bam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

template <typename T, typename F> VProcess::Stats parse(const FileName &file, T o, F f)
{
    typedef VProcess::Stats Stats;
    typedef VProcess::Status Status;

    Stats stats;
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        if (!x.mapped)
        {
            stats.nNA++;
        }
        else if (isRevChr(x.cID))
        {
            stats.nSeqs++; // Reverse genome
        }
        else
        {
            stats.nEndo++; // Forward genome
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
//            auto &seen  = seenMates[x.name];
//            auto first  = seen.isFirstPair ? &seen : &x;
//            auto second = seen.isFirstPair ? &x : &seen;
//
//            /*
//             * Only complement reads aligned to the reverse genome
//             */
//
//            if (__impl__->isReverse(first->cID))
//            {
//                if (first->isForward)
//                {
//                    complement(first->seq);
//                }
//                else
//                {
//                    std::reverse(first->seq.begin(), first->seq.end());
//                }
//            }
//
//            if (__impl__->isReverse(second->cID))
//            {
//                if (second->isForward)
//                {
//                    complement(second->seq);
//                }
//                else
//                {
//                    std::reverse(second->seq.begin(), second->seq.end());
//                }
//            }
//
//            const auto bothRev =  __impl__->isReverse(first->cID) && __impl__->isReverse(second->cID);
//            const auto bothFor = !__impl__->isReverse(first->cID) && !__impl__->isReverse(second->cID);
//            const auto anyRev  =  __impl__->isReverse(first->cID) || __impl__->isReverse(second->cID);
//            const auto anyFor  = !__impl__->isReverse(first->cID) || !__impl__->isReverse(second->cID);
//            const auto anyMap  =  first->mapped ||  second->mapped;
//            const auto anyNMap = !first->mapped || !second->mapped;
//
//            VFlip::Status status;
//
//            if (bothRev && !anyNMap)
//            {
//                status = Status::ReverseReverse;
//            }
//            else if (bothFor && !anyNMap)
//            {
//                status = Status::ForwardForward;
//            }
//            else if (anyRev && anyFor && !anyNMap)
//            {
//                status = Status::ForwardReverse;
//            }
//            else if (anyNMap && anyMap && anyRev)
//            {
//                status = Status::ReverseNotMapped;
//            }
//            else if (anyNMap && anyMap && anyFor)
//            {
//                status = Status::ForwardNotMapped;
//            }
//            else
//            {
//                status = Status::NotMappedNotMapped;
//            }
//
//            __stats__.counts[status]++;
//            __impl__->process(*first, *second, status);
//            seenMates.erase(x.name);
        }
    }, true);

    o.logInfo("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);
        
        if (isRevChr(i.second.cID))
        {
            stats.counts[Status::RevHang]++;
            f(i.second, i.second, Status::RevHang);
        }
        else
        {
            stats.counts[Status::ForHang]++;
            f(i.second, i.second, Status::ForHang);
        }
    }

    return stats;
}

VProcess::Stats VProcess::analyze(const FileName &file, const Options &o)
{
    /*
     * For efficiency, this tool writes output files directly in the analyze() function.
     */
    
    struct Impl
    {
        Impl(const Options &o)
        {
            h1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            f1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            f2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            
            static const FileName HANG_1   = "VarProcess_hanging.fq";
            static const FileName SEQS_1   = "VarProcess_sequins_1.fq";
            static const FileName SEQS_2   = "VarProcess_sequins_2.fq";
            static const FileName GENOME_1 = "VarProcess_genome_1.fq";
            static const FileName GENOME_2 = "VarProcess_genome_2.fq";
            static const FileName AMBIG_1  = "VarProcess_ambiguous_1.fq";
            static const FileName AMBIG_2  = "VarProcess_ambiguous_2.fq";

            h1->open(HANG_1);
            f1->open(SEQS_1);
            f2->open(SEQS_2);
            a1->open(AMBIG_1);
            a2->open(AMBIG_2);
            g1->open(GENOME_1);
            g2->open(GENOME_2);
        }
        
        inline void writeHang(const ParserBAM::Data &x)
        {
            if (x.mapped)
            {
                if (x.isFirstPair)
                {
                    h1->write("@" + x.name + "/1");
                }
                else
                {
                    h1->write("@" + x.name + "/2");
                }
                
                h1->write(x.seq);
                h1->write("+");
                h1->write(x.qual);
            }
        }
        
        inline void writeBefore(const ParserBAM::Data &x1, const ParserBAM::Data &x2)
        {
            s1.push_back(x1);
            s2.push_back(x2);
        }
        
        inline void writePaired(std::shared_ptr<FileWriter> p1, std::shared_ptr<FileWriter> p2, const ParserBAM::Data &x1, const ParserBAM::Data &x2)
        {
            p1->write("@" + x1.name + "/1");
            p1->write(x1.seq);
            p1->write("+");
            p1->write(x1.qual);
            p2->write("@" + x2.name + "/2");
            p2->write(x2.seq);
            p2->write("+");
            p2->write(x2.qual);
        };

        inline void writeAmb(const ParserBAM::Data &x, const ParserBAM::Data &y)
        {
            writePaired(a1, a2, x, y);
        }

        inline void writeGeno(const ParserBAM::Data &x, const ParserBAM::Data &y)
        {
            writePaired(g1, g2, x, y);
        }

        // Alignment records for sequins (we'll need them for sampling)
        std::vector<ParserBAM::Data> s1, s2;
        
        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> f1, f2;
        std::shared_ptr<FileWriter> g1, g2;
    };

    const auto &r = Standard::instance().r_var;

    // Regions without edge effects
    const auto r1 = r.regs1();

    // Region with edge effects
    const auto r2 = r.regs2();
    
    Impl impl(o);
    
    return parse(file, o, [&](const ParserBAM::Data &x1, const ParserBAM::Data &x2, Status status)
    {
        switch (status)
        {
            case Status::ReverseReverse:
            case Status::ReverseNotMapped:
            {
                impl.writeBefore(x1, x2);
                
                auto trim = [&](const ParserBAM::Data &x)
                {
                    std::vector<DInter *> multi;
                    const auto m = x.mapped && r1.count(x.cID) ? r1.at(x.cID).contains(x.l, &multi) : nullptr;
                    
                    if (m)
                    {
                        std::sort(multi.begin(), multi.end(), [&](const DInter * x, const DInter * y)
                        {
                            return x->l().length() < y->l().length();
                        });

                        // The smallest region
                        const auto m = multi.front();
                        
                        const auto lTrim = std::abs(x.l.start - m->l().start) <= o.trim;
                        const auto rTrim = std::abs(x.l.end - m->l().end) <= o.trim;
                        
                        return lTrim || rTrim;
                    }
                    
                    return false;
                };
                
                /*
                 * Perform edge trimming and calculate alignment coverage for sequins
                 */
                
                if (!trim(x1))
                {
                    
                }
                
                if (!trim(x2))
                {
                    
                }
                
                break;
            }

            case Status::ForwardForward:
            {
                impl.writeGeno(x1, x2);
                
                /*
                 * Calculate alignment coverage for genomic
                 */

                break;
            }
                
            case Status::ForwardReverse:
            case Status::ForwardNotMapped:
            case Status::NotMappedNotMapped:
            {
                impl.writeAmb(x1, x2);
                break;
            }

            case Status::RevHang:
            case Status::ForHang:
            {
                impl.writeHang(x1);
                break;
            }
        }
    });
}

void VProcess::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
}
