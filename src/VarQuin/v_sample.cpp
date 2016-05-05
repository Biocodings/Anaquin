#include "tools/sampling.hpp"
#include "VarQuin/v_sample.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/writer_sam.hpp"

using namespace Anaquin;

/*
 * Histogram manipulations and operations
 */

template <typename T> static Counts sums(const std::map<T, Counts> &m)
{
    Counts c = 0;
    
    for (const auto &i : m)
    {
        if (i.second == 0)
        {
            c++;
        }
        else
        {
            c += i.second;
        }
    }
    
    return c;
}

static bool checkAlign(const ChrID &genoID, const ChrID &id, const Locus &l)
{
    const auto &r = Standard::instance().r_var;

    if (id == ChrT)
    {
        return r.match(l, MatchRule::Contains);
    }
    else if (id == genoID)
    {
        return r.findGeno(genoID, l);
    }

    return false;
}

VSample::Stats VSample::stats(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    assert(!r.genoID().empty());
    
    o.info("Query: " + r.genoID());
    o.analyze(file);

    Stats stats;

    stats.cov = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && !(p.i % 1000000))
        {
            //o.wait(std::to_string(p.i));
        }
        
        if (align.cID == ChrT)
        {
            std::cout << checkAlign(r.genoID(), align.cID, align.l) << std::endl;
        }

        return checkAlign(r.genoID(), align.cID, align.l);
    });

    if (!stats.cov.hist.count(ChrT))
    {
        throw std::runtime_error("Failed to find any alignment for " + ChrT);
    }
    else if (!stats.cov.hist.count(r.genoID()))
    {
        throw std::runtime_error("Failed to find any alignment for " + r.genoID());
    }

    o.info(std::to_string(sums(stats.cov.hist)) + " alignments in total");
    o.info(std::to_string(stats.cov.hist.at(ChrT)) + " alignments to chrT");
    o.info(std::to_string(stats.cov.hist.at(r.genoID())) + " alignments to " + r.genoID());
    o.info(std::to_string(stats.cov.inters.size()) + " intervals generated");

    o.info("Generating statistics for " + std::string(ChrT));
    stats.chrT = stats.cov.inters.find(ChrT)->stats([&](const ChrID &id, Base i, Base j, Coverage cov)
    {
        return static_cast<bool>(r.match(Locus(i, j), MatchRule::Contains));
    });

    o.info("Generating statistics for " + r.genoID());
    stats.endo = stats.cov.inters.find(r.genoID())->stats([&](const ChrID &id, Base i, Base j, Coverage cov)
    {
        return static_cast<bool>(r.findGeno(r.genoID(), Locus(i, j)));
    });

    assert(stats.chrT.mean && stats.endo.mean);

    o.info("Calculating coverage for " + ChrT + " and " + r.genoID());

    /*
     * Now we have the data, we'll need to compare the coverage and determine what fraction that
     * the synthetic chromosome needs to be subsampled.
     */
    
    switch (o.method)
    {
        case ArithAverage:
        {
            stats.chrTC = stats.chrT.mean;
            stats.endoC = stats.endo.mean;
            break;
        }

        case Maximum:
        {
            stats.chrTC = stats.chrT.max;
            stats.endoC = stats.endo.max;
            break;
        }

        case Median:
        {
            stats.chrTC = stats.chrT.p50;
            stats.endoC = stats.endo.p50;
            break;
        }

        case Percentile75:
        {
            stats.chrTC = stats.chrT.p75;
            stats.endoC = stats.endo.p75;
            break;
        }
    }

    assert(stats.chrTC && stats.endoC);

    if (stats.endoC >= stats.chrTC)
    {
        // Not an error because it could happen in a simulation
        o.warn("The coverage for the genome is higher than the synthetic chromosome. This is unexpected because the genome is much wider.");
    }

    return stats;
}

void VSample::sample(const FileName &src, const FileName &dst, const Stats &stats, const Options &o)
{
    assert(!src.empty() && !dst.empty());
    
    /*
     * Generating a subsampled alignment. It's expected that coverage would roughly match between
     * the genome and synthetic chromosome.
     */
    
    o.info("Sampling the alignment");
    
    WriterSAM writer;
    writer.open(dst);

    if (stats.sample() == 0.0)
    {
        o.warn("Sampling fraction is zero. This could be an error in the inputs.");
    }
    else if (stats.sample() == 1.0)
    {
        o.warn("Sampling fraction is one. This could be an error in the inputs.");
    }
    
    SamplingTool sampler(1 - stats.sample());
    
    ParserSAM::parse(src, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
                         
        const bam1_t *b    = reinterpret_cast<bam1_t *>(info.data);
        const bam_hdr_t *h = reinterpret_cast<bam_hdr_t *>(info.header);
        
        if (!align.i)
        {
            /*
             * This is the key, randomly write the read with certain probability
             */

            if (align.cID != ChrT || sampler.select(bam_get_qname(b)))
            {
                writer.write(h, b);
            }
        }
    });
    
    writer.close();
}

// Generate and report statistics for subsampling
void VSample::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    auto meth2Str = [&]()
    {
        switch (o.method)
        {
            case Percentile75: { return "75th";    }
            case ArithAverage: { return "Mean";    }
            case Maximum:      { return "Maximum"; }
            case Median:       { return "Median";  }
        }
    };
    
    o.info(meth2Str());
    
    // Statistics before alignment
    const auto before = VSample::stats(file, o);

    // Subsample the alignment
    VSample::sample(file, o.work + "/VarSubsample_sampled.sam", before);

    // Statistics after alignment
    const auto after = VSample::stats(o.work + "/VarSubsample_sampled.sam", o);

    /*
     * Generating bedgraph for the pre-statistics
     */
    
    auto pre = CoverageTool::CoverageBedGraphOptions();
    
    pre.writer = o.writer;
    pre.file   = "VarSubsample_before.bedgraph";

    CoverageTool::bedGraph(before.cov, pre, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        return checkAlign(r.genoID(), id, Locus(i, j));
    });

    /*
     * Generating statistics for the post-statistics
     */
    
    auto post = CoverageTool::CoverageBedGraphOptions();

    post.writer = o.writer;
    post.file   = "VarSubsample_after.bedgraph";

    CoverageTool::bedGraph(after.cov, post, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        return checkAlign(r.genoID(), id, Locus(i, j));
    });
    
    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Proportion of alignments mapped to the synthetic and genome\n"
                         "   ***\n\n"
                         "   Unmapped:  %2% aligns\n"
                         "   Synthetic: %3% aligns\n"
                         "   Genome:    %4% aligns\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %5%\n\n"
                         "   Reference: %6% sequins\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Genome)\n"
                         "   ***\n\n"
                         "   File: %7%\n\n"
                         "   Reference:  %8% intervals\n\n"
                         "   ********************************************************************\n"
                         "   ***                                                              ***\n"
                         "   ***      How the coverage is calculated? Possibilities:          ***\n"
                         "   ***                                                              ***\n"
                         "   ***           - Mean                                             ***\n"
                         "   ***           - Quartile                                         ***\n"
                         "   ***           - Median                                           ***\n"
                         "   ***           - Maximum                                          ***\n"
                         "   ***                                                              ***\n"
                         "   ***      Please refer to the online documentation for details    ***\n"
                         "   ***                                                              ***\n"
                         "   ********************************************************************\n\n"
                         "   Method: %9%\n\n"
                         "   *******************************************\n"
                         "   ***                                     ***\n"
                         "   ***    Statistics before subsampling    ***\n"
                         "   ***                                     ***\n"
                         "   *******************************************\n\n"
                         "   Aligns: %10%\n"
                         "   Depth (Synthetic): %11%\n"
                         "   Depth (Genome):    %12%\n\n"
                         "   *******************************************\n"
                         "   ***                                     ***\n"
                         "   ***    Statistics after subsampling     ***\n"
                         "   ***                                     ***\n"
                         "   *******************************************\n\n"
                         "   Aligns: %13%\n"
                         "   Depth (Synthetic): %14%\n"
                         "   Depth (Genome):    %15%\n";

    o.writer->open("VarSubsample_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % before.cov.unmapped
                                            % before.cov.n_chrT
                                            % before.cov.n_geno
                                            % o.rChrT
                                            % r.countSeqs()
                                            % o.rGeno
                                            % r.countInters()
                                            % meth2Str()
                                            % sums(before.cov.hist)
                                            % before.chrTC
                                            % before.endoC
                                            % sums(after.cov.hist)
                                            % after.chrTC
                                            % after.endoC).str());
    o.writer->close();
}