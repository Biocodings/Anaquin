#include <ss/stats.hpp>
#include "tools/random.hpp"
#include "VarQuin/v_sample.hpp"
#include "writers/writer_sam.hpp"

extern Anaquin::FileName BedRef();

using namespace Anaquin;

ParserBAMBED::Stats VSample::sample(const FileName    &file,
                                    const NormFactors &norms,
                                    const Chr2DInters &sampled,
                                    const Chr2DInters &trimmed,
                                    const VSample::Options &o)
{
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;
    
    /*
     * Initalize independnet random generators for every sampling region
     */
    
    Selection select;
    
    for (const auto &i : norms)
    {
        for (const auto &j : i.second)
        {
            A_ASSERT(j.second >= 0 && j.second <= 1.0 && !isnan(j.second));
            
            // Create independent random generator for each region
            select[i.first][j.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - j.second));
        }
    }

    A_ASSERT(select.size() == norms.size());

    o.info("Sampling: " + file);
    
    WriterSAM writer;
    writer.open("");

    return ParserBAMBED::parse(file, sampled, [&](const ParserBAM::Data &x,
                                                  const ParserBAM::Info &info,
                                                  const DInter *inter)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logInfo(std::to_string(info.p.i));
        }
        
        /*
         * We should sample :
         *
         *   - Anything that is not mapped
         *   - Anything outside the sampling regions
         *   - Inside the region with probability
         */

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
            // Write SAM read to console
            writer.write(x);
            
            if (trimmed.count(x.cID) && trimmed.at(x.cID).overlap(x.l))
            {
                trimmed.at(x.cID).overlap(x.l)->map(x.l);
            }
            
            return ParserBAMBED::Response::OK;
        }

        return ParserBAMBED::Response::SKIP_EVERYTHING;
    });
}

template <typename Stats> Coverage stats2cov(const VSample::Method meth, const Stats &stats)
{
    switch (meth)
    {
        case VSample::Method::Mean:   { return stats.mean; }
        case VSample::Method::Median: { return stats.p50;  }

        /*
         * Prop and Reads specifies a fixed proportion to subsample. It's not actually a measure
         * to report coverage.
         */

        case VSample::Method::Prop:
        case VSample::Method::Reads: { return stats.mean; }
    }
}

VSample::CalibrateStats VSample::check(const FileName &endo,
                                       const FileName &seqs,
                                       const Chr2DInters &tRegs,
                                       const Chr2DInters &regs,
                                       const VSample::Options &o)
{
    A_ASSERT(!tRegs.empty());
    A_ASSERT(tRegs.size() == regs.size());
 
    VSample::CalibrateStats stats;
    
    // Checking endogenous alignments before sampling
    stats.es = ParserBAMBED::parse(endo, tRegs, [&](const ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logWait(std::to_string(info.p.i));
        }
        
        if (x.mapped)
        {
            stats.nEndo++;
        }
        
        return ParserBAMBED::Response::OK;
    });
    
    // Checking sequin alignments before sampling
    stats.ss = ParserBAMBED::parse(seqs, tRegs, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *inter)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logWait(std::to_string(info.p.i));
        }
        
        if (x.mapped)
        {
            stats.nSeqs++;
        }
        
        return ParserBAMBED::Response::OK;
    });
    
    // For each chromosome...
    for (auto &i : tRegs)
    {
        const auto cID = i.first;
        
        // For each region...
        for (auto &j : i.second.data())
        {
            const auto &l = j.second.l();
            
            // Endogenous statistics for the region
            const auto gs = stats.es.inters.at(cID).find(l.key())->stats();
            
            // Sequins statistics for the region
            const auto ss = stats.ss.inters.at(cID).find(l.key())->stats();
            
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
                case VSample::Method::Mean:
                case VSample::Method::Median: { norm = seqsC ? std::min(endoC / seqsC, 1.0) : NAN; break; }
                case VSample::Method::Prop:
                {
                    A_ASSERT(!isnan(o.p));
                    norm = o.p;
                    break;
                }
                    
                case VSample::Method::Reads:
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
    
    return stats;
}

VSample::Stats VSample::analyze(const FileName &endo, const FileName &seqs, const Options &o)
{
    o.analyze(endo);
    o.analyze(seqs);

    o.logInfo("Edge: " + std::to_string(o.edge));
    
    const auto &r = Standard::instance().r_var;
    
    VSample::Stats stats;
    
    // Regions to subsample after trimming
    const auto tRegs = r.regs2();
    
    // Regions without trimming
    const auto regs = r.regs1();
    
    // Check calibration statistics
    stats.cStats = check(endo, seqs, tRegs, regs, o);
    
    const auto norms = stats.cStats.norms;
    stats.c2v = stats.cStats.c2v;
    
    const auto allBeforeEndoC = stats.cStats.allBeforeEndoC;
    const auto allBeforeSeqsC = stats.cStats.allBeforeSeqsC;
    
    // We have the normalization factors so we can proceed with subsampling.
    const auto after = VSample::sample(seqs, norms, regs, tRegs, o);
    
    stats.afterSeqs = VSample::afterSeqsC(tRegs, stats.c2v, o);
    
    stats.tBefore = tBefore(stats.cStats, after);
    stats.tAfter  = tAfter (stats.cStats, after);
    stats.sBefore = sBefore(stats.cStats, after);
    stats.sAfter  = sAfter (stats.cStats, after);

    return stats;
}

double VSample::afterSeqsC(const Chr2DInters &tRegs, std::map<ChrID, std::map<Locus, SampledInfo>> &c2v, VSample::Options o)
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
            
            c2v[i.first].at(j.second.l()).after  = cov;
            c2v[i.first].at(j.second.l()).nAfter = aligns;

            // Required for generating summary statistics
            allAfterSeqsC.push_back(cov);
        }
    }
    
    return SS::mean(allAfterSeqsC);
}

VSample::GenomeSequins VSample::tBefore(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    VSample::GenomeSequins r;
    
    r.nEndo = x.nEndo;
    r.nSeqs = x.nSeqs;
    
    return r;
}

VSample::GenomeSequins VSample::sAfter(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    VSample::GenomeSequins r;
    
    r.nEndo = x.es.nMap;
    r.nSeqs = y.nMap;
    
    return r;
}

VSample::GenomeSequins VSample::sBefore(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    VSample::GenomeSequins r;
    
    r.nEndo = x.es.nMap;
    r.nSeqs = x.ss.nMap;

    return r;
}

VSample::GenomeSequins VSample::tAfter(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    VSample::GenomeSequins r;
    
    r.nEndo = x.nEndo;
    r.nSeqs = y.nNA + y.nMap;
    
    return r;
}

static void generateCSV(const FileName &file, const VSample::Stats &stats, const VSample::Options &o)
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

static void generateSummary(const FileName &file,
                            const FileName &endo,
                            const FileName &seqs,
                            const VSample::Stats &stats,
                            const VSample::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);
    
    const auto summary = "-------VarSubsample Summary Statistics\n\n"
                         "-------VarSubsample Inputs\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Alignment file (sample):   %2%\n"
                         "       Alignment file (sequin):   %3%\n\n"
                         "       Method: %5%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %4% regions\n\n"
                         "-------Before subsampling (within sampling regions)\n\n"
                         "       Sample coverage (average): %16$.2f\n"
                         "       Sequin coverage (average): %17$.2f\n\n"
                         "-------After subsampling (within sampling regions)\n\n"
                         "       Sample coverage (average): %18$.2f\n"
                         "       Sequin coverage (average): %19$.2f\n\n"
                         "       Scaling Factor: %14% \u00B1 %15%\n\n"
                         "-------Total alignments (before subsampling)\n\n"
                         "       Sample: %6%\n"
                         "       Sequin: %7%\n\n"
                         "-------Total alignments (after subsampling)\n\n"
                         "       Sample: %8%\n"
                         "       Sequin: %9%\n\n"
                         "-------Alignments within specified regions (before subsampling)\n\n"
                         "       Sample: %10%\n"
                         "       Sequin: %11%\n\n"
                         "-------Alignments within specified regions (after subsampling)\n\n"
                         "       Sample: %12%\n"
                         "       Sequin: %13%\n\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()                 // 1
                                            % endo                     // 2
                                            % seqs                     // 3
                                            % r.nRegs()                // 4
                                            % meth2Str(o.meth)         // 5
                                            % stats.tBefore.nEndo      // 6
                                            % stats.tBefore.nSeqs      // 7
                                            % stats.tAfter.nEndo       // 8
                                            % stats.tAfter.nSeqs       // 9
                                            % stats.sBefore.nEndo      // 10
                                            % stats.sBefore.nSeqs      // 11
                                            % stats.sAfter.nEndo       // 12
                                            % stats.sAfter.nSeqs       // 13
                                            % stats.cStats.normMean()  // 14
                                            % stats.cStats.normSD()    // 15
                                            % stats.cStats.meanBEndo() // 16
                                            % stats.cStats.meanBSeqs() // 17
                                            % stats.cStats.meanBEndo() // 18
                                            % stats.afterSeqs          // 19
                     ).str());
    o.writer->close();
}

void VSample::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto stats = analyze(endo, seqs, o);
    
    /*
     * Generating VarSubsample_summary.stats
     */
    
    generateSummary("VarSubsample_summary.stats", endo, seqs, stats, o);

    /*
     * Generating VarSubsample_sequins.csv
     */

    generateCSV("VarSubsample_sequins.csv", stats, o);
}
