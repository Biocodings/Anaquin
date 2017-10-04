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

template <typename O> bool shouldTrim(Stats &stats, const ParserBAM::Data &x, const Headers &heads, const O &o)
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

    if (lTrim) { stats.trim.left++;  }
    if (rTrim) { stats.trim.right++; }

    stats.trim.before++;
    if (!lTrim && !rTrim) { stats.trim.after++; }
    
    return lTrim || rTrim;
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
                      const std::set<SequinID> &seqs,
                      const std::map<SequinID, std::pair<ChrID, std::string>> &s2e,
                      const std::map<SequinID, std::pair<ChrID, std::string>> &s2s,
                      const VProcess::Options &o)
{
    for (const auto &sID : seqs)
    {
        assert(s2s.count(sID) && s2e.count(sID));
        
        const auto &p1 = s2e.at(sID);
        const auto &p2 = s2s.at(sID);
        
        assert(stats.gInters.count(p1.first));
        assert(stats.gInters.at(p1.first).find(p1.second));
        assert(stats.sInters.count(p2.first));
        assert(stats.sInters.at(p2.first).find(p2.second));

        // Locus for the sequin
        const auto &l = stats.gInters.at(p1.first).find(p1.second)->l();
        
        // Chromosome for the sequin
        const auto &cID = p1.first;
        
        // Genomic statistics for the region
        const auto gs = stats.gInters.at(p1.first).find(p1.second)->stats();
        
        // Sequins statistics for the region
        const auto ss = stats.sInters.at(p2.first).find(p2.second)->stats();
        
        o.info("Calculating coverage for " + sID);
        
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
            o.logWarn((boost::format("Normalization is NAN for %1%:%2%-%3%") % cID
                                                                             % l.start
                                                                             % norm).str());
            
            // We can't just use NAN...
            norm = 0.0;
        }
        else if (norm == 1.0)
        {
            o.logWarn((boost::format("Normalization is 1 for %1%:%2%-%3% (%4%)") % cID
                                                                                 % l.start
                                                                                 % l.end
                                                                                 % sID).str());
        }
        else
        {
            o.logInfo((boost::format("Normalization is %1% for %2%:%3%-%4% (%5%)") % norm
                                                                                   % cID
                                                                                   % l.start
                                                                                   % l.end
                                                                                   % sID).str());
        }
        
        stats.c2v[cID][l].nEndo   = gs.aligns;
        stats.c2v[cID][l].nBefore = ss.aligns;
        
        stats.c2v[cID][l].rID    = sID;
        stats.c2v[cID][l].endo   = endoC;
        stats.c2v[cID][l].before = seqsC;
        stats.c2v[cID][l].norm   = stats.norms[cID][l] = norm;
        
        stats.allNorms.push_back(norm);
    }
    
    assert(!stats.norms.empty());
}

static void sample(VProcess::Stats &stats,
                   const Chr2DInters &sampled,
                   const Chr2DInters &r2,
                   const std::map<SequinID, std::pair<ChrID, std::string>> &s2e,
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
            
            const auto cID1 = s2e.at(x1.cID).first;
            const auto cID2 = s2e.at(x2.cID).first;

            if (r2.count(cID1) && r2.at(cID1).overlap(x1.l))
            {
                r2.at(cID1).overlap(x1.l)->map(x1.l);
            }

            if (r2.count(cID2) && r2.at(cID2).overlap(x2.l))
            {
                r2.at(cID2).overlap(x2.l)->map(x2.l);
            }
        }
    }

    f1->close();
    f2->close();
}

struct SampledInfo
{
    SequinID rID;
    
    // Alignment coverage for the endogenous sample
    Coverage endo;
    
    // Alignment coverage before subsampling
    Coverage before;
    
    // Alignment coverage after subsampling
    Coverage after;
    
    // Number of alignments before and after
    Counts nEndo, nBefore, nAfter;
    
    // Normalization factor
    Proportion norm;
};

template <typename O> double afterSeqsC(const Chr2DInters &tRegs, Stats &stats, const O &o)
{
    std::vector<double> allAfterSeqsC;
    
    // For each chromosome...
    for (auto &i : tRegs)
    {
        // For each region...
        for (auto &j : i.second.data())
        {
            // Coverage after subsampling
            const auto cov = stats2cov(o.meth, j.second.stats());
            
            // Alignments after subsampling
            const auto aligns = j.second.stats().aligns;
            
            stats.c2v[i.first].at(j.second.l()).after  = cov;
            stats.c2v[i.first].at(j.second.l()).nAfter = aligns;
            
            // Required for generating summary statistics
            allAfterSeqsC.push_back(cov);
        }
    }
    
    return SS::mean(allAfterSeqsC);
}

struct GenomeSequins
{
    Counts nEndo = 0;
    Counts nSeqs = 0;
};

GenomeSequins tBefore(const Stats &x)
{
    GenomeSequins r;
    
    r.nEndo = x.nEndo;
    r.nSeqs = x.nSeqs;
    
    return r;
}

GenomeSequins sAfter(const Stats &x)
{
    GenomeSequins r;
    
    r.nEndo = 9999; //x.es.nMap;
    r.nSeqs = 9999; //y.nMap;
    
    return r;
}

GenomeSequins sBefore(const Stats &x)
{
    GenomeSequins r;
    
    r.nEndo = 9999; //x.es.nMap;
    r.nSeqs = 9999; //x.ss.nMap;
    
    return r;
}

GenomeSequins tAfter(const Stats &x)
{
    GenomeSequins r;
    
    r.nEndo = x.nEndo;
    r.nSeqs = 9999; //y.nNA + y.nMap;
    
    return r;
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
    
    stats.nRegs = countMap(r1, [&](ChrID, const DIntervals<> &x) { return x.size(); });
    
    /*
     * Initalize genomic and sequin regions
     */
    
    // Mapping between sequins to their endogenous regions
    std::map<SequinID, std::pair<ChrID, std::string>> s2e;
    
    // Mapping between sequins to their sequin regions
    std::map<SequinID, std::pair<ChrID, std::string>> s2s;
    
    auto initRegs = [&]()
    {
        // For each chromosome...
        for (const auto &i : r2)
        {
            DIntervals<> x1;
            
            for (const auto &inter : i.second.data())
            {
                assert(!inter.second.name().empty());
                
                // Sequin for the region
                const auto &sID = inter.second.name();
                
                const auto &l1 = inter.second.l();
                x1.add(DInter(l1.key(), l1));
                
                // Update mapping for endogenous
                s2e[sID] = std::pair<ChrID, std::string>(i.first, l1.key());
                
                // Locus for the sequin region
                const auto l2 = Locus(0, inter.second.l().end - inter.second.l().start);
                
                DIntervals<> x2;
                x2.add(DInter(l2.key(), l2));
                
                // Only single chromosome for each sequin
                stats.sInters[sID] = x2;
                stats.sInters[sID].build();
                
                // Update mapping for endogenous
                s2s[sID] = std::pair<ChrID, std::string>(sID, l2.key());
            }
            
            stats.gInters[i.first] = x1;
            stats.gInters[i.first].build();
        }
    };
    
    initRegs();
    
    assert(!s2e.empty() && !s2s.empty());
    assert(s2e.size() == s2s.size());

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
        else if (seqs.count(x.cID)) { if (isLadQuin(x.cID)) { stats.lad.nLad++; } stats.nSeqs++; }
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
                
                if (!shouldTrim(stats, *first, heads, o) && !shouldTrim(stats, *second, heads, o))
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
                
                if (!shouldTrim(stats, *first, heads, o) && !shouldTrim(stats, *second, heads, o))
                {
                    assert(s2c.count(x.cID));
                    assert(stats.sInters.count(x.cID));

                    // Calculate alignment coverage for sequin regions
                    coverage(x, x.cID, stats.sInters);
                    
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
    
    /*
     * Checking calibration before sampling
     */
    
    calibrate(stats, r2, seqs, s2e, s2s, o);
    
    /*
     * Calibrating sequin alignments
     */
    
    sample(stats, r1, r2, s2e, o);
    
    /*
     * Checking calibration after sampling
     */

    stats.afterSeqs = afterSeqsC(r2, stats, o);
    
//    stats.tBefore = tBefore(stats.cStats, after);
//    stats.tAfter  = tAfter (stats.cStats, after);
//    stats.sBefore = sBefore(stats.cStats, after);
//    stats.sAfter  = sAfter (stats.cStats, after);

    return stats;
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

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VProcess::Stats &stats, const VProcess::Options &o)
{
    extern FileName BedRef();

    const auto summary = "-------VarProcess Summary Statistics\n\n"
                         "-------VarFlip Inputs\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Input alignment file: %2%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %3% regions\n"
                         "       Method:  %4%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %5% (64$.2f%%)\n"
                         "       Genome:   %7% (%8$.2f%%)\n"
                         "       Sequins:  %9% (%10$.2f%%)\n\n"
                         "-------Trimming\n\n"
                         "       Left:  %11% reads\n"
                         "       Right: %12% reads\n\n"
                         "-------Before trimming\n\n"
                         "       Number of alignments: %13%\n\n"
                         "-------After trimming\n\n"
                         "       Number of alignments: %14%\n\n"
                         "-------Sequin Outputs\n\n"
                         "       Flipped reads:   %15% (%16$.2f%%)\n"
                         "       Ambiguous reads: %17% (%18$.2f%%)\n"
                         "       Hanging reads:   %19% (%20$.2f%%)\n\n"
                         "-------Before calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %21$.2f\n"
                         "       Sequin coverage (average): %22$.2f\n\n"
                         "-------After calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %23$.2f\n"
                         "       Sequin coverage (average): %24$.2f\n\n"
                         "       Scaling Factor: %25% \u00B1 %26%\n\n"
                         "-------Total alignments (before subsampling)\n\n"
                         "       Sample: %27%\n"
                         "       Sequin: %28%\n\n"
                         "-------Total alignments (after subsampling)\n\n"
                         "       Sample: %29%\n"
                         "       Sequin: %30%\n\n"
                         "-------Alignments within specified regions (before subsampling)\n\n"
                         "       Sample: %31%\n"
                         "       Sequin: %32%\n\n"
                         "-------Alignments within specified regions (after subsampling)\n\n"
                         "       Sample: %33%\n"
                         "       Sequin: %34%\n\n";

    #define C(x) stats.counts.at(x)
    
    const auto cf = C(Status::ReverseReverse) + C(Status::ReverseNotMapped);
    const auto ca = C(Status::ForwardForward) + C(Status::ForwardReverse) + C(Status::ForwardNotMapped) + C(Status::NotMappedNotMapped);
    const auto ch = C(Status::RevHang) + C(Status::ForHang);
    const auto pf = 100.0 * cf / (cf + ca + ch);
    const auto pa = 100.0 * ca / (cf + ca + ch);
    const auto ph = 100.0 * ch / (cf + ca + ch);
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()          // 1
                                            % src               // 2
                                            % stats.nRegs       // 3
                                            % "LeftRight"       // 4
                                            % stats.nNA         // 5
                                            % stats.pNA()       // 6
                                            % stats.nEndo       // 7
                                            % stats.pEndo()     // 8
                                            % stats.nSeqs       // 9
                                            % stats.pSyn()      // 10
                                            % stats.trim.left   // 11
                                            % stats.trim.right  // 12
                                            % stats.trim.before // 13
                                            % stats.trim.after  // 14
                                            % cf                // 15
                                            % pf                // 16
                                            % ca                // 17
                                            % pa                // 18
                                            % ch                // 19
                                            % ph                // 20
                                            % "????"                // 21
                                            % "????"                // 22
                                            % "????"                // 23
                                            % "????"                // 24
                                            % "????"                // 25
                                            % "????"                // 26
                                            % "????"                // 27
                                            % "????"                // 28
                                            % "????"                // 29
                                            % "????"                // 30
                                            % "????"                // 31
                                            % "????"                // 32
                                            % "????"                // 33
                                            % "????"                // 34
                     ).str());
}

static void writeSequins(const FileName &file, const Stats &stats, const VProcess::Options &o)
{
    o.generate(file);
    
    const auto format = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8$.2f\t%9$.2f\t%10$.2f\t%11$.2f");
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Start"
                                           % "End"
                                           % "SampleReads"
                                           % "BeforeReads"
                                           % "AfterReads"
                                           % "GenomeCov"
                                           % "BeforeCov"
                                           % "AfterCov"
                                           % "Norm").str());
    
    // For each chromosome...
    for (const auto &i : stats.c2v)
    {
        // For each variant...
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % j.second.rID
                                                   % i.first
                                                   % j.first.start
                                                   % j.first.end
                                                   % j.second.nEndo
                                                   % j.second.nBefore
                                                   % j.second.nAfter
                                                   % j.second.endo
                                                   % j.second.before
                                                   % j.second.after
                                                   % j.second.norm).str());
        }
    }
    
    o.writer->close();
}

void VProcess::report(const FileName &file, const Options &o)
{
    // For efficiency, this tool writes some of the output files directly in the analyze() function.
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarProcess_summary.stats
     */
    
    writeSummary("VarProcess_summary.stats", file, stats, o);

    /*
     * Generating VarProcess_sequins.tsv
     */
    
    writeSequins("VarProcess_sequins.tsv", stats, o);
}

