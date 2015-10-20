#include "tools/sampling.hpp"
#include "variant/v_sample.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/writer_sam.hpp"

using namespace Anaquin;

VSample::Stats VSample::stats(const FileName &file, const Options &o)
{
    assert(!o.queryID.empty());

    Stats stats;

    const auto &r = Standard::instance();
    
    o.info("Query: " + o.queryID);
    o.analyze(file);

    /*
     * Generating coverage on both chromosomes
     */
     
    stats.cov = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && !(p.i % 1000000))
        {
            o.wait(std::to_string(p.i));
        }
        
        if (align.id == r.id)
        {
            return static_cast<bool>(r.r_var.findGeno(align.l));
        }
        else if (align.id == o.queryID)
        {
            return static_cast<bool>(r.r_var.findQuery(o.queryID, align.l));
        }

        return false;
    });

    if (!stats.cov.hist.count(r.id))
    {
        throw std::runtime_error("Failed to find any alignment for " + r.id);
    }
    else if (!stats.cov.hist.count((o.queryID)))
    {
        throw std::runtime_error("Failed to find any alignment for " + o.queryID);
    }
    
    o.info(std::to_string(sums(stats.cov.hist)) + " alignments in total");
    o.info(std::to_string(stats.cov.hist.at(r.id)) + " alignments to chrT");
    o.info(std::to_string(stats.cov.hist.at(o.queryID)) + " alignments to " + o.queryID);
    o.info(std::to_string(stats.cov.inters.size()) + " intervals generated");

    o.info("Generating statistics for " + r.id);
    stats.chrT = stats.cov.inters.find(r.id)->stats([&](const ChromoID &id, Base i, Base j, Coverage cov)
    {
        return static_cast<bool>(r.r_var.findGeno(Locus(i, j)));
    });

    o.info("Generating statistics for " + o.queryID);
    stats.query = stats.cov.inters.find(o.queryID)->stats([&](const ChromoID &id, Base i, Base j, Coverage cov)
    {
        return static_cast<bool>(r.r_var.findQuery(o.queryID, Locus(i, j)));
    });

    assert(stats.chrT.mean && stats.query.mean);

    o.info("Calculating coverage for " + r.id + " and " + o.queryID);

    /*
     * Now we have the data, we'll need to compare the coverage and determine what fraction that
     * the synthetic chromosome needs to be subsampled.
     */
    
    switch (o.method)
    {
        case ArithAverage:
        {
            stats.chrTC  = stats.chrT.mean;
            stats.queryC = stats.query.mean;
            break;
        }

        case Maximum:
        {
            stats.chrTC  = stats.chrT.max;
            stats.queryC = stats.query.max;
            break;
        }

        case Median:
        {
            stats.chrTC  = stats.chrT.p50;
            stats.queryC = stats.query.p50;
            break;
        }

        case Percentile75:
        {
            stats.chrTC  = stats.chrT.p75;
            stats.queryC = stats.query.p75;
            break;
        }
    }

    assert(stats.chrTC && stats.queryC);

    if (stats.queryC >= stats.chrTC)
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
            
            if (align.id != Standard::instance().id || sampler.select(bam_get_qname(b)))
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
    VSample::sample(file, o.working + "/VarSample_sampled.bam", before);

    // Statistics after alignment
    const auto after = VSample::stats(o.working + "/VarSample_sampled.bam", o);
    
    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    
    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:    %2% alignments\n"
                         "   Genome:      %3% alignments\n"
                         "   Synthetic:   %4% alignments\n\n"
                         "   Reference:   %5% sequins\n\n"
                         "   Method: %6%\n\n"
                         "   ***********************\n"
                         "   Before subsampling:\n"
                         "   ***********************\n\n"
                         "   Alignments: %7%\n"
                         "   Query Coverage: %8%\n"
                         "   Synthetic Coverage: %9%\n\n"
                         "   ***********************\n"
                         "   After subsampling:\n"
                         "   ***********************\n\n"
                         "   Alignments: %10%\n"
                         "   Query Coverage: %11%\n"
                         "   Synthetic Coverage: %12%\n";
    
    o.writer->open("VarSubsample_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % before.cov.unmapped
                                            % before.cov.n_expT
                                            % before.cov.n_chrT
                                            % Standard::instance().r_var.countRefGenes()
                                            % meth2Str()
                                            % sums(before.cov.hist)
                                            % before.queryC
                                            % before.chrTC
                                            % sums(after.cov.hist)
                                            % after.queryC
                                            % after.chrTC).str());
    o.writer->close();
}