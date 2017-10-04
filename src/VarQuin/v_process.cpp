#include "tools/random.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_process.hpp"
#include "writers/bam_writer.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

typedef VProcess::Stats  Stats;
typedef VProcess::Status Status;
typedef VProcess::Method Method;

typedef std::map<ChrID, Base> Headers;
typedef std::map<ChrID, DIntervals<>> Regions;

/*
 * Implement trimming for sequins (reverse genome and ladders)
 */

template <typename O> bool shouldTrim(const ParserBAM::Data &x, const Headers &heads, const O &o)
{
    assert(heads.count(x.cID));
    
    if (!o.shouldTrim)
    {
        return false;
    }
    
    // Length of the sequin
    const auto len = heads.at(x.cID);
    
    const auto lTrim = std::abs(x.l.start) <= o.trim;
    const auto rTrim = std::abs(x.l.end - len) <= o.trim;

    return lTrim || rTrim;
}

template <typename T, typename F> VProcess::Stats &parse(const FileName &file, VProcess::Stats &stats, T o, F f)
{
    const auto &r = Standard::instance().r_var;
    
    // Regions without edge effects
    const auto r1 = r.r1()->inters();
    
    // Regions with edge effects
    const auto r2 = r.r2()->inters();

    // Sequin names
    const auto seqs = r.r1()->seqs();
    
    // Mapping from sequins to chromosomes
    const auto s2c = r.r1()->s2c();
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;
    
    const auto heads = ParserBAM::header(file);
    
    /*
     * Update alignment coverage for endogenous and sequins
     */
    
    auto coverage = [&](const ParserBAM::Data &x, const ChrID &cID, ID2Intervals &inters)
    {
        if (x.mapped && inters.count(cID))
        {
            const auto matched = inters.at(cID).overlap(x.l);
            
            if (matched)
            {
                matched->map(x.l);
            }
        }
    };

    // Eg: chrev1, LAD_18 ...
    auto isVarQuin = [&](const ChrID &x)
    {
        return isLadQuin(x) || seqs.count(x);
    };
    
    ParserBAM::parse(file, [&](const ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        if (!x.mapped)              { stats.nNA++;  }
        else if (seqs.count(x.cID)) { if (isLadQuin(x.cID)) { stats.nLad++; } stats.nSeqs++; }
        else                        { stats.nEndo++; }
        
        /*
         * Directly write out endogenous alignments (pairing not required)
         */
        
        if (!isLadQuin(x.cID) && (!seqs.count(x.cID) && !seqs.count(x.rnext)))
        {
            // Calculate alignment coverage for endogenous regions
            coverage(x, x.cID, stats.gInters);

            // Write out the read
            f(&x, nullptr, Status::ForwardForward);
            
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

            if (isLadQuin(first->cID) && isLadQuin(second->cID))
            {
                /*
                 * Ladder alignments don't require calibration and flipping
                 */

                if (!shouldTrim(*first, heads, o) && !shouldTrim(*second, heads, o))
                {
                    f(first, second, Status::LadQuin);
                }
            }
            else if (!s2c.count(first->cID) || !s2c.count(second->cID))
            {
                f(first, second, Status::ReverseLadQuin);
            }
            else
            {
                /*
                 * Reverse alignments require trimming, calibration and flipping
                 */

                if (!shouldTrim(*first, heads, o) && !shouldTrim(*second, heads, o))
                {
                    assert(s2c.count(x.cID));
                    
                    // Calculate alignment coverage for sequin regions
                    coverage(x, s2c.at(x.cID), stats.sInters);
                
                    const auto bothRev =  isVarQuin(first->cID) && isVarQuin(second->cID);
                    const auto bothFor = !isVarQuin(first->cID) && !isVarQuin(second->cID);
                    const auto anyRev  =  isVarQuin(first->cID) || isVarQuin(second->cID);
                    const auto anyFor  = !isVarQuin(first->cID) || !isVarQuin(second->cID);
                    const auto anyMap  =  first->mapped ||  second->mapped;
                    const auto anyNMap = !first->mapped || !second->mapped;
                    
                    Status status;
                    
                    if (bothRev && !anyNMap)
                    {
                        status = Status::ReverseReverse;
                    }
                    else if (bothFor && !anyNMap)
                    {
                        throw std::runtime_error("Status::ForwardForward is invalid for sequins");
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
                    
                    // Write out the paired-end reads
                    f(first, second, status);
                }
            }
            
            seenMates.erase(x.name);
        }
    }, true);
    
    for (auto &i: seenMates)
    {
        o.logWarn("Unpaired mate: " + i.first);
        
        // Compute the complement (but not reverse)
        complement(i.second.seq);
        
        if (isVarQuin(i.second.cID))
        {
            stats.counts[Status::RevHang]++;
            f(&i.second, nullptr, Status::RevHang);
        }
        else
        {
            stats.counts[Status::ForHang]++;
            f(&i.second, nullptr, Status::ForHang);
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

/*
 * Calculate calibration factors (but not peforming it)
 */

static void calibrate(VProcess::Stats &stats,
                      const Chr2DInters &r2,
                      const VProcess::Options &o)
{
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
    
    assert(!stats.norms.empty());
}

static void sample(VProcess::Stats &stats,
                   const Chr2DInters &sampled,
                   const VProcess::Options &o)
{
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;
    
    /*
     * Initalize independnet random generators for each sampling region
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
    
    f1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
    f2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));

    f1->open("VarProcess_sequins_1.fq");
    f2->open("VarProcess_sequins_2.fq");
    
    /*
     * We should sample :
     *
     *   - Anything that is not mapped
     *   - Anything outside the sampling regions
     *   - Inside the region with probability
     */

    assert(stats.s1.size() == stats.s2.size());
    
    for (auto i = 0; i < stats.s1.size(); i++)
    {
        auto &x1 = stats.s1.at(i);
        auto &x2 = stats.s2.at(i);

        auto shouldSampled = !x1.mapped || !x2.mapped;
        
        if (!shouldSampled)
        {
            DInter *inter;
            
            /*
             * Should that be contains or overlap? We prefer overlaps because any read that is overlapped
             * into the regions still give valuable information and sequencing depth.
             */
            
            if (sampled.count(x1.cID) && (inter = sampled.at(x1.cID).overlap(x1.l)))
            {
                // Transform for the trimming region
                auto l = Locus(inter->l().start + o.edge, inter->l().end - o.edge);
                
                if (select.at(x1.cID).at(l)->select(x1.name))
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
        
        /*
         * Only complement reads aligned to the reverse genome
         */
        
        auto __reverse__ = [&](ParserBAM::Data &x)
        {
            if (x.isForward)
            {
                complement(x.seq);
            }
            else
            {
                std::reverse(x.seq.begin(), x.seq.end());
            }
        };
        
        if (shouldSampled)
        {
            __reverse__(x1);
            __reverse__(x2);
            
            f1->write("@" + x1.name + "/1");
            f1->write(x1.seq);
            f1->write("+");
            f1->write(x1.qual);
            
            f2->write("@" + x2.name + "/2");
            f2->write(x2.seq);
            f2->write("+");
            f2->write(x2.qual);
        }
    }

    f1->close();
    f2->close();
}

VProcess::Stats VProcess::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    // Regions without edge effects
    const auto r1 = r.r1()->inters();
    
    // Regions with edge effects
    const auto r2 = r.r2()->inters();
    
    Stats stats;
    
    struct Impl
    {
        Impl(Stats &stats, const Options &o) : stats(stats)
        {
            h1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            a2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            l1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
            l2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));

            h1->open("VarProcess_hanging.fq");
            a1->open("VarProcess_ambiguous_1.fq");
            a2->open("VarProcess_ambiguous_2.fq");
            l1->open("VarProcess_ladder_1.fq");
            l2->open("VarProcess_ladder_2.fq");

            // Write to BAM alignment file
            geno.open(o.work + "/VarProcess_genome.bam");
        }
        
        ~Impl()
        {
            h1->close();
            a1->close();
            a2->close();
            l1->close();
            l2->close();
            geno.close();
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
        
        inline void writeLad(const ParserBAM::Data &x, const ParserBAM::Data &y)
        {
            writePaired(l1, l2, x, y);
        }

        inline void writeAmb(const ParserBAM::Data &x, const ParserBAM::Data &y)
        {
            writePaired(a1, a2, x, y);
        }
        
        inline void writeEndo(const ParserBAM::Data &x)
        {
            geno.write(x);
        }
        
        Stats &stats;
        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> l1, l2;

        // Writing BAM for genomical alignments
        BAMWriter geno;
    };
    
    Impl impl(stats, o);
    
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

    parse(file, stats, o, [&](const ParserBAM::Data *x1, const ParserBAM::Data *x2, Status status)
    {
        switch (status)
        {
            case Status::LadQuin:
            {
                impl.writeLad(*x1, *x2);
                break;
            }

            case Status::ReverseReverse:
            case Status::ReverseNotMapped:
            {
                impl.writeBefore(*x1, *x2);
                break;
            }
                
            case Status::ForwardForward:
            {
                // We should directly write out endogenous reads
                assert(x2 == nullptr);
                
                // Write out endogenous alignments
                impl.writeEndo(*x1);
                
                break;
            }
                
            case Status::ForwardReverse:
            case Status::ReverseLadQuin:
            case Status::ForwardNotMapped:
            case Status::NotMappedNotMapped:
            {
                impl.writeAmb(*x1, *x2);
                break;
            }
                
            case Status::RevHang:
            case Status::ForHang:
            {
                assert(x2 == nullptr);                
                impl.writeHang(*x1);
                break;
            }
        }
    });

    /*
     * Checking calibration
     */

    calibrate(stats, r2, o);
    
    /*
     * Calibrating sequin alignments
     */
    
    sample(stats, r1, o);
    
    return stats;
}

void VProcess::report(const FileName &file, const Options &o)
{
    /*
     * For efficiency, this tool writes output files directly in the analyze() function.
     */
    
    analyze(file, o);
}

