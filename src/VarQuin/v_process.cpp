#include "tools/random.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_process.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

template <typename T, typename F> VProcess::Stats &parse(const FileName &file, VProcess::Stats &stats, T o, F f)
{
    typedef VProcess::Stats Stats;
    typedef VProcess::Status Status;

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

typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;

template <typename Stats> Coverage stats2cov(const VProcess::Method meth, const Stats &stats)
{
    switch (meth)
    {
        case VProcess::Method::Mean:   { return stats.mean; }
        case VProcess::Method::Median: { return stats.p50;  }
            
        /*
         * Prop and Reads specifies a fixed proportion to subsample. It's not actually a measure
         * to report coverage.
         */

        case VProcess::Method::Prop:
        case VProcess::Method::Reads: { return stats.mean; }
    }
}

static void sample(const Chr2DInters &r2,
                   VProcess::Stats &stats,
                   const Chr2DInters &sampled,
                   const VProcess::Options &o)
{
    typedef VProcess::Method Method;
    
    // For each chromosome...
    for (auto &i : r2)
    {
        const auto cID = i.first;
        
        // For each region...
        for (auto &j : i.second.data())
        {
            const auto &l = j.second.l();
            
            // Genomic statistics for the region
            const auto gs = stats.gInters.at(cID).find(l.key())->stats();
            
            // Sequins statistics for the region
            const auto ss = stats.sInters.at(cID).find(l.key())->stats();
            
            o.info("Calculating coverage for " + j.first);
            
            /*
             * Now we have the data, we'll need to compare coverage and determine the fraction that
             * the synthetic alignments needs to be sampled.
             */
            
            const auto endoC = stats2cov(o.meth, gs);
            const auto seqsC = stats2cov(o.meth, ss);
            
            Proportion norm;
            
            switch (o.meth)
            {
                case Method::Mean:
                case Method::Median: { norm = seqsC ? std::min(endoC / seqsC, 1.0) : NAN; break; }
                case Method::Prop:
                {
                    A_ASSERT(!isnan(o.p));
                    norm = o.p;
                    break;
                }
                    
                case Method::Reads:
                {
                    A_ASSERT(!isnan(o.reads));
                    
                    const auto gAligns = gs.aligns;
                    const auto sAligns = ss.aligns;
                    
                    /*
                     * Nothing to sample if no sequin alignments. Subsample everything if the
                     * genomic region has higher coverage.
                     */
                    
                    norm = sAligns == 0 ? 0 : gAligns >= sAligns ? 1 : ((Proportion) gAligns) / sAligns;
                    
                    break;
                }
            }
            
            stats.allBeforeEndoC.push_back(endoC);
            stats.allBeforeSeqsC.push_back(seqsC);
            
            if (isnan(norm))
            {
                o.logWarn((boost::format("Normalization is NAN for %1%:%2%-%3%") % i.first
                                                                                 % l.start
                                                                                 % norm).str());
                
                // We can't just use NAN...
                norm = 0.0;
            }
            else if (norm == 1.0)
            {
                o.logWarn((boost::format("Normalization is 1 for %1%:%2%-%3% (%4%)") % i.first
                                                                                     % l.start
                                                                                     % l.end
                                                                                     % j.first).str());
            }
            else
            {
                o.logInfo((boost::format("Normalization is %1% for %2%:%3%-%4% (%5%)") % norm
                                                                                       % i.first
                                                                                       % l.start
                                                                                       % l.end
                                                                                       % j.first).str());
            }
            
            stats.c2v[cID][l].nEndo   = gs.aligns;
            stats.c2v[cID][l].nBefore = ss.aligns;

            stats.c2v[cID][l].rID    = j.first;
            stats.c2v[cID][l].endo   = endoC;
            stats.c2v[cID][l].before = seqsC;
            stats.c2v[cID][l].norm   = stats.norms[i.first][l] = norm;

            stats.allNorms.push_back(norm);
        }
    }
    
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;

    /*
     * Initalize independnet random generators for every sampling region
     */

    Selection select;

    for (const auto &i : stats.norms)
    {
        for (const auto &j : i.second)
        {
            A_ASSERT(j.second >= 0 && j.second <= 1.0 && !isnan(j.second));

            // Create independent random generator for each region
            select[i.first][j.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - j.second));
        }
    }

    A_ASSERT(select.size() == stats.norms.size());

    std::shared_ptr<FileWriter> f1, f2;
    
    static const FileName SEQS_1 = "VarProcess_sequins_1.fq";
    static const FileName SEQS_2 = "VarProcess_sequins_2.fq";

    f1->open(SEQS_1);
    f2->open(SEQS_2);

    auto __sample__ = [&](const std::vector<ParserBAM::Data> &v)
    {
        /*
         * We should sample :
         *
         *   - Anything that is not mapped
         *   - Anything outside the sampling regions
         *   - Inside the region with probability
         */

        for (const auto &x : v)
        {
            auto shouldSampled = !x.mapped;
            
            if (!shouldSampled)
            {
                DInter *inter;
                
                /*
                 * Should that be contains or overlap? We prefer overlaps because any read that is overlapped
                 * into the regions still give valuable information and sequencing depth.
                 */
                
                if (sampled.count(x.cID) && (inter = sampled.at(x.cID).overlap(x.l)))
                {
                    // Transform for the trimming region
                    auto l = Locus(inter->l().start + o.edge, inter->l().end - o.edge);
                    
                    if (select.at(x.cID).at(l)->select(x.name))
                    {
                        shouldSampled = true;
                    }
                }
                else
                {
                    // Never throw away reads outside the regions
                    shouldSampled = true;
                }
            }
            
            if (shouldSampled)
            {
                std::cout << x.cID << std::endl;
//                if (trimmed.count(x.cID) && trimmed.at(x.cID).overlap(x.l))
//                {
//                    trimmed.at(x.cID).overlap(x.l)->map(x.l);
//                }
            }
        }
    };
    
    __sample__(stats.s1);
    __sample__(stats.s2);
    
    f1->close();
    f2->close();
}

VProcess::Stats VProcess::analyze(const FileName &file, const Options &o)
{
    Stats stats;

    struct Impl
    {
        Impl(Stats &stats, const Options &o) : stats(stats)
        {
            h1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            g1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            g2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            
            static const FileName HANG_1   = "VarProcess_hanging.fq";
            static const FileName GENOME_1 = "VarProcess_genome_1.fq";
            static const FileName GENOME_2 = "VarProcess_genome_2.fq";
            static const FileName AMBIG_1  = "VarProcess_ambiguous_1.fq";
            static const FileName AMBIG_2  = "VarProcess_ambiguous_2.fq";

            h1->open(HANG_1);
            a1->open(AMBIG_1);
            a2->open(AMBIG_2);
            g1->open(GENOME_1);
            g2->open(GENOME_2);
        }
        
        ~Impl()
        {
            h1->close();
            a1->close();
            a2->close();
            g1->close();
            g2->close();
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
            stats.s1.push_back(x1);
            stats.s2.push_back(x2);
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
        
        Stats &stats;
        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> g1, g2;
    };

    Impl impl(stats, o);
    
    const auto &r = Standard::instance().r_var;

    // Regions without edge effects
    const auto r1 = r.regs1();

    // Region with edge effects
    const auto r2 = r.regs2();

    auto initInts = [&](ID2Intervals &inters)
    {
        // For each chromosome...
        for (const auto &i : r2)
        {
            DIntervals<> x;

            for (const auto &inter : i.second.data())
            {
                const auto &l = inter.second.l();
                x.add(DInter(l.key(), l));
            }

            inters[i.first] = x;
            inters[i.first].build();
        }
    };
    
    initInts(stats.gInters);
    initInts(stats.sInters);

    auto gCov = [&](const ParserBAM::Data &x)
    {
        if (x.mapped && stats.gInters.count(x.cID))
        {
            stats.gInters[x.cID].overlap(x.l)->map(x.l);
        }
    };

    auto sCov = [&](const ParserBAM::Data &x)
    {
        if (x.mapped && stats.sInters.count(x.cID))
        {
            stats.sInters[x.cID].overlap(x.l)->map(x.l);
        }
    };

    parse(file, stats, o, [&](const ParserBAM::Data &x1, const ParserBAM::Data &x2, Status status)
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
                 * Perform edge trimming and calculate alignment coverage for sequin regions
                 */

                if (!trim(x1)) { sCov(x1); }
                if (!trim(x2)) { sCov(x2); }

                break;
            }

            case Status::ForwardForward:
            {
                impl.writeGeno(x1, x2);

                /*
                 * Calculate alignment coverage for genomic regions
                 */

                gCov(x1);
                gCov(x2);

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
    
    /*
     * Subsample sequin reads
     */
    
    //sample(r2, stats, r1, o);

    return stats;
}

void VProcess::report(const FileName &file, const Options &o)
{
    /*
     * For efficiency, this tool writes output files directly in the analyze() function.
     */

    analyze(file, o);
}
