#include "tools/random.hpp"
#include "VarQuin/VarQuin.hpp"
#include "writers/bam_writer.hpp"
#include "VarQuin/v_partition.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

typedef VPartition::Stats     Stats;
typedef VPartition::Paired    Paired;
typedef VPartition::Method    Method;
typedef VPartition::Options   Options;
typedef std::map<ChrID, Base> Headers;

static ParserBAM::Data flip(const ParserBAM::Data &x)
{
    auto r = x;
    
    if (r.isForward)
    {
        complement(r.seq);
    }
    else
    {
        std::reverse(r.seq.begin(), r.seq.end());
    }
    
    return r;
};

static bool shouldTrim(const ParserBAM::Data &x,
                       const Headers &heads,
                       const std::map<ChrID, DIntervals<>> &r3,
                       const Options &o,
                       bool &lTrim,
                       bool &rTrim)
{
    A_ASSERT(heads.count(x.cID));

    // Unable to trim if it's not mapped
    if (!o.shouldTrim || !x.mapped)
    {
        return false;
    }
    
    auto lrTrim = [&](Base start, Base end)
    {
        // Trimming by left?
        lTrim = x.l.start < start ? false : (x.l.start - start) <= o.trim;
        
        // Trimming by right?
        rTrim = x.l.end > end ? false : (end - x.l.end) <= o.trim;
        
        return lTrim || rTrim;
    };
    
    // Trimming by looking at BAM headers?
    auto shouldTrim = lrTrim(0, heads.at(x.cID));

    // Is this a structural variant region?
    if (!shouldTrim && r3.count(x.cID))
    {
        A_ASSERT(r3.at(x.cID).data().size() == 1);
        
        // Trimming interval
        const auto &l = r3.at(x.cID).data().begin()->second.l();
        
        // Should we trim by structural variant?
        shouldTrim = lrTrim(l.start, l.end);
    }

    return shouldTrim;
}

template <typename Stats> Coverage stats2cov(const Method meth, const Stats &stats)
{
    switch (meth)
    {
        case Method::Mean:   { return stats.mean; }
        case Method::Median: { return stats.p50;  }
            
        /*
         * Prop and Reads specifies a fixed proportion to subsample. It's not actually a measure
         * to report coverage.
         */
            
        case Method::Prop:
        case Method::Reads: { return stats.mean; }
    }
}

static SequinID trimSID(const SequinID &sID)
{
    if (isSubstr(sID, "_R") || isSubstr(sID, "_V"))
    {
        // Eg: DEL_02_R to DEL_02
        return noLast(sID, "_");
    }
    
    return sID;
}

/*
 * Calculate calibration factors (but not peforming it)
 */

static void calibrate(Stats &stats,
                      const Chr2DInters &r2,
                      const std::set<SequinID> &seqs,
                      const Options &o)
{
    auto &mStats = stats.mStats;
    
    for (const auto &sID : seqs)
    {
        auto sID_ = trimSID(sID);
        A_ASSERT(stats.mStats.s2s.count(sID_) && stats.mStats.s2e.count(sID_));
        
        const auto &p1 = mStats.s2e.at(sID_);
        const auto &p2 = mStats.s2s.at(sID_);
        
        A_ASSERT(mStats.eInters.count(p1.first));
        A_ASSERT(mStats.eInters.at(p1.first).find(p1.second));
        A_ASSERT(mStats.bInters.count(p2.first));
        A_ASSERT(mStats.bInters.at(p2.first).find(p2.second));

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
                                                                                 % sID_).str());
        }
        else
        {
            o.logInfo((boost::format("Normalization is %1% for %2%:%3%-%4% (%5%)") % norm
                                                                                   % cID
                                                                                   % l.start
                                                                                   % l.end
                                                                                   % sID_).str());
        }
        
        stats.c2v[sID_].nEndo   = es.aligns;
        stats.c2v[sID_].nBefore = ss.aligns;
        
        stats.c2v[sID_].rID    = sID_;
        stats.c2v[sID_].endo   = endoC;
        stats.c2v[sID_].before = seqsC;
        stats.c2v[sID_].norm   = norm;
        
        stats.cStats.covs[sID_]  = endoC;
        stats.cStats.norms[sID_] = norm;
        
        stats.cStats.allNorms.push_back(norm);
    }
    
    A_ASSERT(!stats.cStats.norms.empty());
}

static Counts sample(Stats &stats, std::set<ReadName> &sampled, const Options &o)
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

    f1->open("VarPartition_flipped_1.fq");
    f2->open("VarPartition_flipped_2.fq");

    /*
     * We should sample :
     *
     *   - Anything that is not mapped
     *   - Anything outside the sampling regions
     *   - Inside the region with probability
     */

    A_ASSERT(stats.s1.size() == stats.s2.size());
    
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
                auto sID = trimSID(x.cID);
                
                // Targeted coverage
                const auto exp = stats.cStats.covs.at(sID);
                
                // Latest coverage
                const auto obs = after.at(sID).stats().mean;
                
                // Stop if the target coverage reached
                if (obs >= exp)
                {
                    o.logInfo("Stopped for " + x.cID);
                    return false;
                }

                return select.at(sID)->select(x.name);
            }
        };

        const auto k1 = shouldKeep(x1);
        const auto k2 = shouldKeep(x2);

        if (k1 || k2)
        {
            auto addCov = [&](const ParserBAM::Data &x)
            {
                if (after.at(trimSID(x.cID)).overlap(x.l))
                {
                    after.at(trimSID(x.cID)).overlap(x.l)->map(x.l);
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
            
            stats.nFlip++;

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

    return 2 * stats.nFlip;
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
    
    for (const auto &sID : stats.mStats.seqs)
    {
        const auto x = stats.mStats.aInters.at(trimSID(sID)).stats();

        // Coverage after sampling
        const auto cov = stats2cov(o.meth, x);
        
        // Alignments after sampling
        const auto aligns = x.aligns;
        
        stats.c2v[trimSID(sID)].after  = cov;
        stats.c2v[trimSID(sID)].nAfter = aligns;
        
        // Required for generating summary statistics
        all.push_back(cov);
    }

    return SS::mean(all);
}

template <typename T, typename F> Stats &parse(const FileName &file, Stats &stats, T o, F f)
{
    o.info("Edge: " + std::to_string(o.edge));
    o.info("Trimming: " + std::to_string(o.trim));

    const auto &r = Standard::instance().r_var;
    
    // Regions without edge effects
    const auto r1 = r.r1()->inters();
    
    // Regions with edge effects
    const auto r2 = r.r2()->inters();

    // Structuring variant trimming region
    const auto r3 = r.r3()->inters();

    // Sequin names
    stats.mStats.seqs = r.r1()->seqs();
    
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;
    
    // Cached for proprocessing
    std::set<ReadName> kept;

    stats.gStats.nRegs = countMap(r1, [&](ChrID, const DIntervals<> &x) { return x.size(); });
    
    /*
     * Initalize genomic and sequin regions
     */
    
    auto initR = [&]()
    {
        std::map<ChrID, std::map<SequinID, Locus>> x;

        for (const auto &i : r1)
        {
            for (const auto &j : i.second.data())
            {
                const auto sID = trimSID(j.second.name());
                const auto &l  = j.second.l();

                if (x.count(i.first) && x[i.first].count(sID))
                {
                    x[i.first][sID].start = std::min(x[i.first][sID].start, l.start);
                    x[i.first][sID].end   = std::max(x[i.first][sID].end, l.end);
                }
                else
                {
                    x[i.first][sID] = l;
                }
            }
        }

        // For each chromosome...
        for (const auto &i : x)
        {
            DIntervals<> x1;
            
            for (const auto &j : i.second)
            {
                A_ASSERT(!j.first.empty());
                
                // Sequin for the region
                const auto &sID = j.first;
                
                const auto &l1 = j.second;
                x1.add(DInter(l1.key(), l1));
                
                // Update mapping for endogenous
                stats.mStats.s2e[sID] = std::pair<ChrID, std::string>(i.first, l1.key());
                
                /*
                 * Construct interval for sequins and consider edge effects
                 */
                
                // Locus for the sequin region
                auto l2 = Locus(0, j.second.end - j.second.start);

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
    
    initR();

    A_ASSERT(!stats.mStats.s2e.empty() && !stats.mStats.s2s.empty());
    A_ASSERT(stats.mStats.s2e.size() == stats.mStats.s2s.size());

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
    
    auto isReverse = [&](const ChrID &x) { return stats.mStats.seqs.count(x);   };
    auto isVarQuin = [&](const ChrID &x) { return isLadQuin(x) || isReverse(x); };
    
    bool t1l, t2l; // Left trimming for paired-end?
    bool t1r, t2r; // Right trimming for paired-end?

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        A_CHECK(x.isPaired, x.name + " is not pair-ended. Singled-ended not supported.");

        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        /*
         * 1. Forward    -> Forward   (ForwardForward)
         * 2. Forward    -> NotMapped (ForwardNotMapped)
         * 3. NotMappted -> Forward   (ForwardNotMapped)
         * 4. NotMapped  -> NotMapped (NotMappedNotMapped)
         * 5. Forward    -> VarQuin   (ForwardVarQuin)
         * 6. VarQuin    -> Forward   (ForwardVarQuin)
         */

        const auto isMap1 = x.mapped;
        const auto isMap2 = !x.rnext.empty() && x.rnext != "*";
        const auto isVar1 = isMap1 && isVarQuin(x.cID);
        const auto isVar2 = isMap2 && isVarQuin(x.rnext);
        const auto isLad1 = isMap1 && isLadQuin(x.cID);
        //const auto isLad2 = isMap2 && isLadQuin(x.rnext);
        const auto isRev1 = isMap1 && isReverse(x.cID);
        //const auto isRev2 = isMap2 && isReverse(x.rnext);

        if      (!isMap1) { stats.nNA++;   }
        else if (isLad1)  { stats.nLad++;  }
        else if (isRev1)  { stats.nRev++;  }
        else              { stats.nEndo++; }
        
        if (isMap1 && !isVar1 && isMap2 && !isVar2)
        {
            // Calculate alignment coverage for endogenous regions
            coverage(x, x.cID, stats.mStats.eInters);
            
            f(&x, nullptr, Paired::ForwardForward);
            
            // Within sequin regions?
            if (stats.mStats.eInters.count(x.cID) && stats.mStats.eInters.at(x.cID).overlap(x.l))
            {
                stats.gStats.bREndo++;
            }
            
            return;
        }
        else if ((isMap1 && isVar1 && !isMap2))
        {
            // Calculate alignment coverage for endogenous regions
            coverage(x, x.cID, stats.mStats.eInters);

            f(&x, nullptr, Paired::ForwardNotMapped);
            
            // Within sequin regions?
            if (stats.mStats.eInters.count(x.cID) && stats.mStats.eInters.at(x.cID).overlap(x.l))
            {
                stats.gStats.bREndo++;
            }

            return;
        }
        else if ((!isMap1 && isMap2 && !isVar2))
        {
            f(&x, nullptr, Paired::ForwardNotMapped);
            return;
        }
        else if (!isMap1 && !isMap2)
        {
            f(&x, nullptr, Paired::NotMappedNotMapped);
            return;
        }
        
        // Don't trim and calibrate unless it's primary
        else if (!x.isPassed || x.isSecondary || x.isSupplement)
        {
            return;
        }
        
        x.lSeq();
        x.lQual();
        x.lName();

        if (!seenMates.count(x.name))
        {
            seenMates[x.name] = x;
        }
        else
        {
            auto &seen  = seenMates[x.name];
            auto first  = seen.isFirstPair ? &seen : &x;
            auto second = seen.isFirstPair ? &x : &seen;

            if ((!isVarQuin(first->cID) && isVarQuin(second->cID)) || (isVarQuin(first->cID) && !isVarQuin(second->cID)))
            {
                f(&x, nullptr, Paired::ForwardVarQuin);
                stats.pairs[Paired::ForwardVarQuin]++;
                seenMates.erase(x.name);
                return;
            }
            
            if (first->mapped)  { stats.trim.before++; }
            if (second->mapped) { stats.trim.before++; }

            if (shouldTrim(*first, heads, r3, o, t1l, t1r) || shouldTrim(*second, heads, r3, o, t2l, t2r))
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
                stats.pairs[Paired::LadQuinLadQuin]++;
                f(first, second, Paired::LadQuinLadQuin);
            }
            else if (isLadQuin(first->cID) || isLadQuin(second->cID))
            {
                stats.pairs[Paired::ReverseLadQuin]++;
                f(first, second, Paired::ReverseLadQuin);
            }
            else
            {
                const auto bothVar =  isVarQuin(first->cID) && isVarQuin(second->cID);
                const auto bothFor = !isVarQuin(first->cID) && !isVarQuin(second->cID);
                const auto anyVar  =  isVarQuin(first->cID) || isVarQuin(second->cID);
                const auto anyFor  = !isVarQuin(first->cID) || !isVarQuin(second->cID);
                const auto anyMap  =  first->mapped ||  second->mapped;
                const auto anyNMap = !first->mapped || !second->mapped;
                
                if (first->isFirstPair == second->isFirstPair)
                {
                    throw std::runtime_error(x.name + " is invalid. No paied-end mate detected for this read.");
                }
                
                A_ASSERT(stats.mStats.bInters.count(trimSID(x.cID)));
                A_ASSERT(stats.mStats.bInters.count(trimSID(x.cID)));
                
                /*
                 * Calculate alignment coverage for sequins before calibration
                 */
                
                coverage(*first,  trimSID(first->cID),  stats.mStats.bInters);
                coverage(*second, trimSID(second->cID), stats.mStats.bInters);

                Paired status;
                
                if (bothVar && !anyNMap)
                {
                    status = Paired::ReverseReverse;
                }
                else if (bothFor && !anyNMap)
                {
                    throw std::runtime_error("Invalid ForwardForward");
                }
                else if (anyVar && anyFor && !anyNMap)
                {
                    throw std::runtime_error("Invalid ForReverse");
                }
                else if (anyNMap && anyMap && anyVar)
                {
                    status = Paired::ReverseNotMapped;
                }
                else if (anyNMap && anyMap && anyFor)
                {
                    throw std::runtime_error("Invalid ForwardNotMapped");
                }
                else
                {
                    throw std::runtime_error("Invalid NotMappedNotMapped");
                }

                stats.pairs[status]++;
                f(first, second, status);
                
                if (o.notCalib && status == Paired::ReverseReverse)
                {
                    f(first, second, Paired::TrimmedNotCalibrated);
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
            stats.pairs[Paired::ReverseHang]++;
            f(&i.second, nullptr, Paired::ReverseHang);
        }
        else
        {
            stats.pairs[Paired::ForwardHang]++;
            f(&i.second, nullptr, Paired::ForwardHang);
        }
    }
    
    /*
     * Checking calibration before sampling
     */ 
    
    o.info("Checking before calibration");
    calibrate(stats, r2, stats.mStats.seqs, o);
    
    /*
     * Calibrating sequin alignments
     */

    std::set<ReadName> sampled;

    o.info("Performing calibration");
    stats.gStats.aTSeqs = sample(stats, sampled, o);
    
    /*
     * Checking calibration after sampling
     */

    o.info("Checking after calibration");
    stats.afterSeqs = checkAfter(stats, o);
    
    stats.gStats.bTEndo = stats.nEndo;
    stats.gStats.bTSeqs = stats.trim.after;
    stats.gStats.aTEndo = stats.nEndo;
    stats.gStats.aREndo = stats.gStats.bREndo;
    
    if (o.notCalib)
    {
        ParserBAM::parse(file, [&](const ParserBAM::Data &x, const ParserBAM::Info &)
         {
#ifdef DEBUG
             if (sampled.count(x.name)) { f(&x, nullptr, Paired::TrimmedCalibrated); }
#endif
         }, false);
    }

    return stats;
}

Stats VPartition::analyze(const FileName &file, const Options &o)
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

            h1->open("VarPartition_hanging.fq");
            a1->open("VarPartition_ambiguous_1.fq");
            a2->open("VarPartition_ambiguous_2.fq");
            l1->open("VarPartition_ladder_1.fq");
            l2->open("VarPartition_ladder_2.fq");

            endo.open(o.work + "/VarPartition_sample.bam");
            
            if (o.notCalib)
            {
                n1 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
                n2 = std::shared_ptr<FileWriter>(new FileWriter(o.work));
                n1->open("VarPartition_notCalibrated_1.fq");
                n2->open("VarPartition_notCalibrated_2.fq");
            }

#ifdef DEBUG
            calib.open(o.work + "/VarPartition_calibrated.bam");
#endif
        }
        
        ~Impl()
        {
            h1->close();
            a1->close();
            a2->close();
            l1->close();
            l2->close();
            endo.close();
            
            if (n1) { n1->close(); }
            if (n2) { n2->close(); }

#ifdef DEBUG
            calib.close();
#endif
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

        inline void writeFASTQ(std::shared_ptr<FileWriter> p, const ParserBAM::Data &x, bool isFirst)
        {
            p->write("@" + x.name + (isFirst ? "/1" : "/2"));
            p->write(x.seq);
            p->write("+");
            p->write(x.qual);
        };

        inline void writeLad1(const ParserBAM::Data &x)
        {
            writeFASTQ(l1, x, true);
        }

        inline void writeLad2(const ParserBAM::Data &x)
        {
            writeFASTQ(l2, x, false);
        }

        inline void writeAmb1(const ParserBAM::Data &x)
        {
            writeFASTQ(a1, x, true);
        }

        inline void writeAmb2(const ParserBAM::Data &x)
        {
            writeFASTQ(a2, x, false);
        }

        inline void writeEndo(const ParserBAM::Data &x)
        {
            endo.write(x);
        }

#ifdef DEBUG
        inline void writeCalib(const ParserBAM::Data &x)
        {
            calib.write(x);
        }
#endif

        inline void writeTrimmedNotCalibrated1(const ParserBAM::Data &x)
        {
            writeFASTQ(n1, flip(x), true);
        }

        inline void writeTrimmedNotCalibrated2(const ParserBAM::Data &x)
        {
            writeFASTQ(n2, flip(x), false);
        }

        Stats &stats;
        std::shared_ptr<FileWriter> h1;
        std::shared_ptr<FileWriter> a1, a2;
        std::shared_ptr<FileWriter> l1, l2;
        std::shared_ptr<FileWriter> n1, n2;

        BAMWriter endo;

#ifdef DEBUG
        BAMWriter calib;
#endif
    };
    
    Impl impl(stats, o);
    
    parse(file, stats, o, [&](const ParserBAM::Data *x1, const ParserBAM::Data *x2, Paired status)
    {
        switch (status)
        {
            case Paired::TrimmedCalibrated:
            {
#ifdef DEBUG
                A_ASSERT(x2 == nullptr);
                impl.writeCalib(*x1);
#endif
                break;
            }
                
            case Paired::TrimmedNotCalibrated:
            {
                impl.writeTrimmedNotCalibrated1(*x1);
                impl.writeTrimmedNotCalibrated2(*x1);
                break;
            }
                
            case Paired::LadQuinLadQuin:
            {
                impl.writeLad1(*x1);
                impl.writeLad2(*x2);
                break;
            }

            case Paired::ReverseReverse:
            case Paired::ReverseNotMapped:
            {
                impl.writeBefore1(*x1);
                impl.writeBefore2(*x2);
                break;
            }
                
            case Paired::ForwardForward:
            case Paired::ForwardNotMapped:
            case Paired::NotMappedNotMapped:
            {
                A_ASSERT(x2 == nullptr);
                impl.writeEndo(*x1);
                break;
            }

            case Paired::ReverseLadQuin:
            case Paired::ForwardVarQuin:
            {
                impl.writeAmb1(*x1);
                if (x2) { impl.writeAmb2(*x2); }
                break;
            }
                
            case Paired::ReverseHang:
            case Paired::ForwardHang:
            {
                A_ASSERT(x2 == nullptr);
                impl.writeHang(*x1);
                break;
            }
        }
    });

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const Stats &stats, const Options &o)
{
    extern FileName Bed1Ref();

    const auto summary = "-------VarPartition Summary Statistics\n\n"
                         "-------Input files\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Input alignment file:      %2%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %3% regions\n"
                         "       Edge:    %4%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped: %5% (%6$.2f%%)\n"
                         "       Sample:   %7% (%8$.2f%%)\n"
                         "       Reverse:  %9% (%10$.2f%%)\n"
                         "       Ladders:  %11% (%12$.2f%%)\n\n"
                         "-------Trimming\n\n"
                         "       Left:  %13% alignments\n"
                         "       Right: %14% alignments\n\n"
                         "-------Before trimming\n\n"
                         "       Number of alignments: %15% (only primary alignments)\n\n"
                         "-------After trimming\n\n"
                         "       Number of alignments: %16%\n\n"
                         "-------Sequin Outputs\n\n"
                         "       Flipped pairs:   %17% (%18$.2f%%)\n"
                         "       Ladder pairs:    %19% (%20$.2f%%)\n"
                         "       Ambiguous pairs: %21% (%22$.2f%%)\n"
                         "       Hanging reads:   %23% (%24$.2f%%)\n\n"
                         "-------Before calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %25$.2f\n"
                         "       Sequin coverage (average): %26$.2f\n\n"
                         "-------After calibration (within sampling regions)\n\n"
                         "       Sample coverage (average): %27$.2f\n"
                         "       Sequin coverage (average): %28$.2f\n\n"
                         "       Scaling Factor: %29% \u00B1 %30%\n\n"
                         "-------Alignments within reference regions (before calibration but after trimming)\n\n"
                         "       Sample: %31%\n"
                         "       Sequin: %32%\n\n"
                         "-------Alignments within reference regions (after calibration)\n\n"
                         "       Sample: %33%\n"
                         "       Sequin: %34%\n";

    #define C(x) stats.pairs.at(x)
    
    const auto cf = stats.nFlip;
    const auto ca = C(Paired::ForwardVarQuin) + C(Paired::ForwardNotMapped) + C(Paired::NotMappedNotMapped) + C(Paired::ReverseLadQuin);
    const auto ch = C(Paired::ReverseHang) + C(Paired::ForwardHang);
    const auto cl = C(Paired::LadQuinLadQuin);

    const auto tot1 = cf + ca + ch + cl;
    const auto pf   = 100.0 * cf / tot1;
    const auto pa   = 100.0 * ca / tot1;
    const auto ph   = 100.0 * ch / tot1;
    const auto pl   = 100.0 * cl / tot1;

    const auto pNA   = 100.0 * stats.nNA   / (stats.nNA + stats.nEndo + stats.nLad + stats.nRev);
    const auto pEndo = 100.0 * stats.nEndo / (stats.nNA + stats.nEndo + stats.nLad + stats.nRev);
    const auto pLad  = 100.0 * stats.nLad  / (stats.nNA + stats.nEndo + stats.nLad + stats.nRev);
    const auto pRev  = 100.0 * stats.nRev  / (stats.nNA + stats.nEndo + stats.nLad + stats.nRev);

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % Bed1Ref()                 // 1
                                            % src                      // 2
                                            % stats.gStats.nRegs       // 3
                                            % o.edge                   // 4
                                            % stats.nNA                // 5
                                            % pNA                      // 6
                                            % stats.nEndo              // 7
                                            % pEndo                    // 8
                                            % stats.nRev               // 9
                                            % pRev                     // 10
                                            % stats.nLad               // 11
                                            % pLad                     // 12
                                            % stats.trim.left          // 13
                                            % stats.trim.right         // 14
                                            % stats.trim.before        // 15
                                            % stats.trim.after         // 16
                                            % cf                       // 17
                                            % pf                       // 18
                                            % cl                       // 19
                                            % pl                       // 20
                                            % ca                       // 21
                                            % pa                       // 22
                                            % ch                       // 23
                                            % ph                       // 24
                                            % stats.cStats.meanBEndo() // 25
                                            % stats.cStats.meanBSeqs() // 26
                                            % stats.cStats.meanBEndo() // 27
                                            % stats.afterSeqs          // 28
                                            % stats.cStats.normMean()  // 29
                                            % stats.cStats.normSD()    // 30
                                            % stats.gStats.bREndo      // 31
                                            % stats.gStats.bTSeqs      // 32
                                            % stats.gStats.aREndo      // 33
                                            % stats.gStats.aTSeqs      // 34
                     ).str());
}

static void writeSequins(const FileName &file, const Stats &stats, const Options &o)
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

static void writeRegions(const FileName &file, const Stats &, const Options &o)
{
    extern FileName Bed2Ref();
    System::copy(Bed2Ref(), o.work + "/" + file);
}

void VPartition::report(const FileName &file, const Options &o)
{
    // For efficiency, this tool writes some of the output files directly in the analyze() function.
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarPartition_summary.stats
     */
    
    writeSummary("VarPartition_summary.stats", file, stats, o);

    /*
     * Generating VarPartition_sequins.tsv
     */
    
    writeSequins("VarPartition_sequins.tsv", stats, o);
    
    /*
     * Generating VarPartition_sequins.bed
     */
    
    writeRegions("VarPartition_regions.bed", stats, o);
}

