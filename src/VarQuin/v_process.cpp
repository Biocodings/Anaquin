#include "tools/random.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_process.hpp"
#include "writers/bam_writer.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

typedef VProcess::Stats   Stats;
typedef VProcess::Status  Status;
typedef VProcess::Method  Method;
typedef VProcess::Options Options;

typedef std::map<ChrID, Base> Headers;

/*
 * Implement trimming for sequins (reverse genome and ladders)
 */

static bool shouldTrim(const ParserBAM::Data &x, const Headers &heads, const Options &o, bool &lTrim, bool &rTrim)
{
    assert(heads.count(x.cID));

    // Unable to trim if it's not mapped
    if (!o.shouldTrim || !x.mapped)
    {
        return false;
    }
    
    // Length of the sequin from BAM header
    const auto len = heads.at(x.cID);
    
    lTrim = std::abs(x.l.start) <= o.trim;
    rTrim = std::abs(x.l.end - len) <= o.trim;
 
    return lTrim || rTrim;
}

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
                      const VProcess::Options &o)
{
    auto &mStats = stats.mStats;
    
    for (const auto &sID : seqs)
    {
        assert(stats.mStats.s2s.count(sID) && stats.mStats.s2e.count(sID));
        
        const auto &p1 = mStats.s2e.at(sID);
        const auto &p2 = mStats.s2s.at(sID);
        
        assert(mStats.eInters.count(p1.first));
        assert(mStats.eInters.at(p1.first).find(p1.second));
        assert(mStats.bInters.count(p2.first));
        assert(mStats.bInters.at(p2.first).find(p2.second));

        // Chromosome for the sequin
        const auto &cID = p1.first;
        
        // Endogenous statistics for the region
        const auto es = mStats.eInters.at(p1.first).find(p1.second)->stats();
        
        // Sequins statistics for the region
        const auto ss = mStats.bInters.at(p2.first).find(p2.second)->stats();
        
        /*
         * Now we have the data, we'll need to compare coverage and determine the fraction that
         * the synthetic alignments needs to be sampled.
         */
        
        const auto endoC = stats2cov(o.meth, es);
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
                
                const auto gAligns = es.aligns;
                const auto sAligns = ss.aligns;
                
                /*
                 * Nothing to sample if no sequin alignments. Subsample everything if the
                 * genomic region has higher coverage.
                 */
                
                norm = sAligns == 0 ? 0 : gAligns >= sAligns ? 1 : ((Proportion) gAligns) / sAligns;
                
                break;
            }
        }
        
        stats.cStats.allBeforeEndoC.push_back(endoC);
        stats.cStats.allBeforeSeqsC.push_back(seqsC);
        
        // Forward locus for the region
        const auto &l = stats.mStats.eInters.at(p1.first).find(p1.second)->l();

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
        
        stats.c2v[sID].nEndo   = es.aligns;
        stats.c2v[sID].nBefore = ss.aligns;
        
        stats.c2v[sID].rID    = sID;
        stats.c2v[sID].endo   = endoC;
        stats.c2v[sID].before = seqsC;
        stats.c2v[sID].norm   = norm;
        
        stats.cStats.covs[sID]  = endoC;
        stats.cStats.norms[sID] = norm;
        
        stats.cStats.allNorms.push_back(norm);
    }
    
    A_ASSERT(!stats.cStats.norms.empty());
}

static Counts sample(Stats &stats, const Chr2DInters &r1, std::set<ReadName> &sampled, const Options &o)
{
    A_ASSERT(!stats.cStats.covs.empty());
    A_ASSERT(!stats.cStats.norms.empty());

    typedef std::map<SequinID, std::shared_ptr<RandomSelection>> Selection;
    
    /*
     * Initalize independnet random generators for each sampling region
     */
    
    Selection select;
    
    for (const auto &i : stats.cStats.norms)
    {
        A_ASSERT(i.second >= 0 && i.second <= 1.0 && !isnan(i.second));
        
        // Create independent random generator for each region
        select[i.first]= std::shared_ptr<RandomSelection>(new RandomSelection((1.0 - i.second), 100));
    }
    
    A_ASSERT(select.size() == stats.cStats.norms.size());
    
    auto f1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
    auto f2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
    auto s1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
    auto s2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));

    f1->open("VarProcess_flipped_1.fq");
    f2->open("VarProcess_flipped_2.fq");

    /*
     * We should sample :
     *
     *   - Anything that is not mapped
     *   - Anything outside the sampling regions
     *   - Inside the region with probability
     */

    assert(stats.s1.size() == stats.s2.size());
    
    Counts nSeqs = 0;
    
    for (auto i = 0; i < stats.s1.size(); i++)
    {
        auto &x1 = stats.s1.at(i);
        auto &x2 = stats.s2.at(i);

        auto &after = stats.mStats.aInters;
        
        auto shouldKeep = [&](const ParserBAM::Data &x)
        {
            if (!x.mapped) { return true; }
            else
            {
                // Targeted coverage
                const auto exp = stats.cStats.covs.at(x.cID);
                
                // Latest coverage
                const auto obs = after.at(x.cID).stats().mean;
                
                // Stop if the target coverage reached
                if (obs >= exp)
                {
                    o.info("Stopped for " + x.cID);
                    return false;
                }

                return select.at(x.cID)->select(x.name);
            }
        };

        const auto k1 = shouldKeep(x1);
        const auto k2 = shouldKeep(x2);

        if (k1 || k2)
        {
            auto addCov = [&](const ParserBAM::Data &x)
            {
                if (after.at(x.cID).overlap(x.l))
                {
                    after.at(x.cID).overlap(x.l)->map(x.l);
                }
            };
            
            addCov(x1);
            addCov(x2);

            sampled.insert(x1.name);
            sampled.insert(x2.name);

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

            __reverse__(x1);
            __reverse__(x2);
            
            nSeqs++;
            nSeqs++;

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

    return nSeqs;
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

template <typename O> double checkAfter(Stats &stats, const O &o)
{
    std::vector<double> all;
    
    for (const auto &seq : stats.mStats.seqs)
    {
        const auto x = stats.mStats.aInters.at(seq).stats();

        // Coverage after sampling
        const auto cov = stats2cov(o.meth, x);
        
        // Alignments after sampling
        const auto aligns = x.aligns;
        
        stats.c2v[seq].after  = cov;
        stats.c2v[seq].nAfter = aligns;
        
        // Required for generating summary statistics
        all.push_back(cov);
    }

    return SS::mean(all);
}

template <typename T, typename F> VProcess::Stats &parse(const FileName &file, VProcess::Stats &stats, T o, F f)
{
    o.info("Edge: " + std::to_string(o.edge));

    const auto &r = Standard::instance().r_var;
    
    // Regions without edge effects
    const auto r1 = r.r1()->inters();
    
    // Regions with edge effects
    const auto r2 = r.r2()->inters();
    
    // Sequin names
    stats.mStats.seqs = r.r1()->seqs();
    
    // Mapping from sequins to their chromosomes
    const auto s2c = r.r1()->s2c();
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;
    
    // Cached for proprocessing
    std::set<ReadName> kept;

    stats.gStats.nRegs = countMap(r1, [&](ChrID, const DIntervals<> &x) { return x.size(); });
    
    /*
     * Initalize genomic and sequin regions
     */
    
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
                stats.mStats.s2e[sID] = std::pair<ChrID, std::string>(i.first, l1.key());
                
                /*
                 * Construct interval for sequins and consider edge effects
                 */
                
                // Locus for the sequin region
                auto l2 = Locus(0, inter.second.l().end - inter.second.l().start);

                l2.end   -= o.edge;
                l2.start += o.edge;
                
                A_ASSERT(l2.end > l2.start);
                
                DIntervals<> x2;
                x2.add(DInter(l2.key(), l2));
                
                stats.mStats.bInters[sID] = x2;
                stats.mStats.bInters[sID].build();
                stats.mStats.aInters[sID] = x2;
                stats.mStats.aInters[sID].build();
                
                stats.mStats.s2s[sID] = std::pair<ChrID, std::string>(sID, l2.key());
            }
            
            stats.mStats.eInters[i.first] = x1;
            stats.mStats.eInters[i.first].build();
        }
    };
    
    initRegs();

    assert(!stats.mStats.s2e.empty() && !stats.mStats.s2s.empty());
    assert(stats.mStats.s2e.size() == stats.mStats.s2s.size());

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
        return isLadQuin(x) || stats.mStats.seqs.count(x);
    };
    
    bool t1l, t2l; // Left trimming for paired-end?
    bool t1r, t2r; // Right trimming for paired-end?

    ParserBAM::parse(file, [&](const ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        A_CHECK(x.isPaired, x.name + " is not pair-ended. Singled-ended not supported.");

        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        if (!x.mapped)             { stats.nNA++;  }
        else if (isLadQuin(x.cID)) { stats.lad.nLad++; stats.nSeqs++; }
        else if (isVarQuin(x.cID)) { stats.nSeqs++; }
        else                       { stats.nEndo++; }
        
        /*
         * Write out endogenous alignments
         */
        
        if (!isVarQuin(x.cID) && !isVarQuin(x.rnext))
        {
            // Calculate alignment coverage for endogenous regions
            coverage(x, x.cID, stats.mStats.eInters);
            
            f(&x, nullptr, Status::ForwardForward);
            
            if (stats.mStats.eInters.count(x.cID) && stats.mStats.eInters.at(x.cID).overlap(x.l))
            {
                stats.gStats.bREndo++;
            }
            
            return;
        }
        
        /*
         * Write out ambiguous alignments (VarQuin and endeogenous)
         */
        
        else if ((!isVarQuin(x.cID) && isVarQuin(x.rnext)) || (isVarQuin(x.cID) && !isVarQuin(x.rnext)))
        {
            f(&x, nullptr, Status::ForwardVarQuin);
            return;
        }
        else if (!x.isPassed || x.isSecondary || x.isSupplement)
        {
            return;
        }

        // Relevant sequin alignments before trimming
        f(&x, nullptr, Status::Original);

        if (!seenMates.count(x.name))
        {
            seenMates[x.name] = x;
        }
        else
        {
            auto &seen  = seenMates[x.name];
            auto first  = seen.isFirstPair ? &seen : &x;
            auto second = seen.isFirstPair ? &x : &seen;

            if (first->mapped)  { stats.trim.before++; }
            if (second->mapped) { stats.trim.before++; }

            if (shouldTrim(*first, heads, o, t1l, t1r) || shouldTrim(*second, heads, o, t2l, t2r))
            {
                if (t1l || t2l) { stats.trim.left  += 2; }
                if (t1r || t2r) { stats.trim.right += 2; }
                seenMates.erase(x.name);
                return;
            }

            kept.insert(first->name);
            kept.insert(second->name);

            if (first->mapped)  { stats.trim.after++; }
            if (second->mapped) { stats.trim.after++; }

            if (isLadQuin(first->cID) && isLadQuin(second->cID))
            {
                f(first, second, Status::LadQuinLadQuin);
            }
            else if (!s2c.count(first->cID) || !s2c.count(second->cID))
            {
                stats.counts[Status::ReverseLadQuin]++;
                f(first, second, Status::ReverseLadQuin);
            }
            else
            {
                const auto bothVar =  isVarQuin(first->cID) && isVarQuin(second->cID);
                const auto bothFor = !isVarQuin(first->cID) && !isVarQuin(second->cID);
                const auto anyVar  =  isVarQuin(first->cID) || isVarQuin(second->cID);
                const auto anyFor  = !isVarQuin(first->cID) || !isVarQuin(second->cID);
                const auto anyMap  =  first->mapped ||  second->mapped;
                const auto anyNMap = !first->mapped || !second->mapped;
                
                assert(s2c.count(x.cID));
                assert(stats.mStats.bInters.count(x.cID));
                
                // Calculate alignment coverage for sequin regions before calibration
                coverage(*first,  first->cID,  stats.mStats.bInters);
                coverage(*second, second->cID, stats.mStats.bInters);

                Status status;
                
                if (bothVar && !anyNMap)
                {
                    status = Status::ReverseReverse;
                }
                else if (bothFor && !anyNMap)
                {
                    throw std::runtime_error("Status::ForwardForward is invalid for sequins");
                }
                else if (anyVar && anyFor && !anyNMap)
                {
                    throw std::runtime_error("Status::ForReverse is invalid for sequins");
                }
                else if (anyNMap && anyMap && anyVar)
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
                stats.counts[status]++;
                
                f(first, second, status);
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
            stats.counts[Status::ReverseHang]++;
            f(&i.second, nullptr, Status::ReverseHang);
        }
        else
        {
            stats.counts[Status::ForwardHang]++;
            f(&i.second, nullptr, Status::ForwardHang);
        }
    }
    
    /*
     * Checking calibration before sampling
     */ 
    
    calibrate(stats, r2, stats.mStats.seqs, o);
    
    /*
     * Calibrating sequin alignments
     */

    std::set<ReadName> sampled;
    stats.gStats.aTSeqs = sample(stats, r1, sampled, o);
    
    /*
     * Checking calibration after sampling
     */

    stats.afterSeqs = checkAfter(stats, o);
    
    stats.gStats.bTEndo = stats.nEndo;
    stats.gStats.bTSeqs = stats.nSeqs;
    stats.gStats.aTEndo = stats.nEndo;
    stats.gStats.aREndo = stats.gStats.bREndo;
    
    /*
     * Write out alignments after trimming
     */
    
    ParserBAM::parse(file, [&](const ParserBAM::Data &x, const ParserBAM::Info &)
    {
        if (kept.count(x.name))    { f(&x, nullptr, Status::Passed);  }
        if (sampled.count(x.name)) { f(&x, nullptr, Status::Sampled); }
    }, false);

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

            orig.open(o.work + "/VarProcess_original.bam");
            endo.open(o.work + "/VarProcess_sample.bam");
            pass.open(o.work + "/VarProcess_passed.bam");
            sample.open(o.work + "/VarProcess_sampled.bam");
        }
        
        ~Impl()
        {
            h1->close();
            a1->close();
            a2->close();
            l1->close();
            l2->close();
            endo.close();
            pass.close();
            orig.close();
            sample.close();
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
        
        inline void writeBefore1(const ParserBAM::Data &x)
        {
            stats.s1.push_back(x);
        }

        inline void writeBefore2(const ParserBAM::Data &x)
        {
            stats.s2.push_back(x);
        }

        inline void writeFASTQ(std::shared_ptr<FileWriter> p, const ParserBAM::Data &x)
        {
            p->write("@" + x.name + "/1");
            p->write(x.seq);
            p->write("+");
            p->write(x.qual);
        };

        inline void writeLad1(const ParserBAM::Data &x)
        {
            writeFASTQ(l1, x);
        }

        inline void writeLad2(const ParserBAM::Data &x)
        {
            writeFASTQ(l2, x);
        }

        inline void writeAmb1(const ParserBAM::Data &x)
        {
            writeFASTQ(a1, x);
        }

        inline void writeAmb2(const ParserBAM::Data &x)
        {
            writeFASTQ(a2, x);
        }

        inline void writeEndo(const ParserBAM::Data &x)
        {
            endo.write(x);
        }

        inline void writeSample(const ParserBAM::Data &x)
        {
            sample.write(x);
        }
        
        inline void writePass(const ParserBAM::Data &x)
        {
            pass.write(x);
        }

        inline void writeOrig(const ParserBAM::Data &x)
        {
            orig.write(x);
        }

        Stats &stats;
        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> l1, l2;

        BAMWriter endo, sample, pass, orig;
    };
    
    Impl impl(stats, o);
    
    parse(file, stats, o, [&](const ParserBAM::Data *x1, const ParserBAM::Data *x2, Status status)
    {
        switch (status)
        {
            case Status::Original:
            {
                impl.writeOrig(*x1);
                break;
            }
                
            case Status::Sampled:
            {
                assert(x2 == nullptr);
                impl.writeSample(*x1);
                break;
            }
                
            case Status::Passed:
            {
                assert(x2 == nullptr);
                impl.writePass(*x1);
                break;
            }
                
            case Status::LadQuinLadQuin:
            {
                impl.writeLad1(*x1);
                impl.writeLad2(*x2);
                break;
            }

            case Status::ReverseReverse:
            case Status::ReverseNotMapped:
            {
                impl.writeBefore1(*x1);
                impl.writeBefore2(*x2);
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
                
            case Status::ReverseLadQuin:
            case Status::ForwardVarQuin:
            case Status::ForwardNotMapped:
            case Status::NotMappedNotMapped:
            {
                impl.writeAmb1(*x1);
                if (x2) { impl.writeAmb2(*x2); }
                break;
            }
                
            case Status::ReverseHang:
            case Status::ForwardHang:
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
                         "-------Input files\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Input alignment file:      %2%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %3% regions\n"
                         "       Edge:    %4% regions\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %5% (%6$.2f%%)\n"
                         "       Genome:   %7% (%8$.2f%%)\n"
                         "       Sequins:  %9% (%10$.2f%%)\n"
                         "       Ladders:  %11%\n\n"
                         "-------Trimming\n\n"
                         "       Left:  %12% alignments\n"
                         "       Right: %13% alignments\n\n"
                         "-------Before trimming\n\n"
                         "       Number of alignments: %14% (only primary alignments)\n\n"
                         "-------After trimming\n\n"
                         "       Number of alignments: %15%\n\n"
                         "-------Sequin Outputs\n\n"
                         "       Flipped reads:   %16% (%17$.2f%%)\n"
                         "       Ambiguous reads: %18% (%19$.2f%%)\n"
                         "       Hanging reads:   %20% (%21$.2f%%)\n\n"
                         "-------Before calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %22$.2f\n"
                         "       Sequin coverage (average): %23$.2f\n\n"
                         "-------After calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %24$.2f\n"
                         "       Sequin coverage (average): %25$.2f\n\n"
                         "       Scaling Factor: %26% \u00B1 %27%\n\n"
                         "-------Alignments within reference regions (before subsampling)\n\n"
                         "       Sample: %28%\n"
                         "       Sequin: %29%\n\n"
                         "-------Alignments within reference regions (after subsampling)\n\n"
                         "       Sample: %30%\n"
                         "       Sequin: %31%\n";

    #define C(x) stats.counts.at(x)
    
    const auto cf = C(Status::ReverseReverse) + C(Status::ReverseNotMapped) + C(Status::LadQuinLadQuin);
    const auto ca = C(Status::ForwardVarQuin) + C(Status::ForwardNotMapped) + C(Status::NotMappedNotMapped) + C(Status::ReverseLadQuin);
    const auto ch = C(Status::ReverseHang) + C(Status::ForwardHang);
    const auto pf = 100.0 * cf / (cf + ca + ch);
    const auto pa = 100.0 * ca / (cf + ca + ch);
    const auto ph = 100.0 * ch / (cf + ca + ch);
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()                 // 1
                                            % src                      // 2
                                            % stats.gStats.nRegs       // 3
                                            % o.edge                   // 4
                                            % stats.nNA                // 5
                                            % stats.pNA()              // 6
                                            % stats.nEndo              // 7
                                            % stats.pEndo()            // 8
                                            % stats.nSeqs              // 9
                                            % stats.pSyn()             // 10
                                            % stats.lad.nLad           // 11
                                            % stats.trim.left          // 12
                                            % stats.trim.right         // 13
                                            % stats.trim.before        // 14
                                            % stats.trim.after         // 15
                                            % cf                       // 16
                                            % pf                       // 17
                                            % ca                       // 18
                                            % pa                       // 19
                                            % ch                       // 20
                                            % ph                       // 21
                                            % stats.cStats.meanBEndo() // 22
                                            % stats.cStats.meanBSeqs() // 23
                                            % stats.cStats.meanBEndo() // 24
                                            % stats.afterSeqs          // 25
                                            % stats.cStats.normMean()  // 26
                                            % stats.cStats.normSD()    // 27
                                            % stats.gStats.bREndo      // 28
                                            % stats.gStats.bTSeqs      // 29
                                            % stats.gStats.aREndo      // 30
                                            % stats.gStats.aTSeqs      // 31
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
    
    for (const auto &i : stats.c2v)
    {
        const auto &j = stats.mStats.s2e.at(i.first);

        const auto x = stats.mStats.eInters.at(j.first).find(j.second);
        assert(x != nullptr);
        
        o.writer->write((boost::format(format) % i.first
                                               % j.first
                                               % x->l().start
                                               % x->l().end
                                               % i.second.nEndo
                                               % i.second.nBefore
                                               % i.second.nAfter
                                               % i.second.endo
                                               % i.second.before
                                               % i.second.after
                                               % i.second.norm).str());
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

