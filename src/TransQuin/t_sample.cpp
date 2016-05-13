#include "tools/sample.hpp"
#include <ss/maths/stats.hpp>
#include "TransQuin/t_sample.hpp"
#include "writers/writer_sam.hpp"
#include "parsers/parser_cufflink.hpp"

using namespace Anaquin;

TSample::Stats TSample::stats(const FileName &file, const Options &o)
{
    TSample::Stats stats;
    
    o.info(file);
    
    switch (o.soft)
    {
        case Software::None:
        {
            break;
        }
            
        case Software::Cufflinks:
        {
            ParserCufflink::parse(file, [&](const ParserCufflink::Data &x, const ParserProgress &)
            {
                if (x.cID == ChrT)
                {
                    stats.n_chrT++;
                    stats.chrT.push_back(x.abund);
                }
                else
                {
                    stats.n_geno++;
                    stats.geno.push_back(x.abund);
                }
            });
            
            break;
        }
    }

    /*
     * Calculate the normalization factor
     */
     
    switch (o.meth)
    {
        case Method::OneP:
        {
            stats.prop = 0.01;
            break;
        }
            
        case Method::FiveP:
        {
            stats.prop = 0.05;
            break;
        }

        case Method::TenP:
        {
            stats.prop = 0.10;
            break;
        }

        case Method::Mean:
        {
            stats.genoBefore = SS::mean(stats.geno);
            stats.chrTBefore = SS::mean(stats.chrT);
            
            if (stats.genoBefore >= stats.chrTAfter)
            {
                throw std::runtime_error("Sequencing depth for genomic transcripts is at least as high as sequins");
            }
            
            stats.prop = stats.genoBefore / stats.chrTBefore;
            
            stats.chrTAfter  = stats.genoBefore;
            stats.chrTBefore = stats.genoBefore;

            break;
        }
    }

    o.info("Proportion: " + toString(stats.prop));
    
    assert(stats.prop > 0 && stats.prop < 1.0);
    return stats;
}

static void generateSummary(const FileName &file, const TSample::Stats &stats, const TSample::Options &o)
{
    o.generate(file);
    
    auto meth2Str = [&]()
    {
        switch (o.meth)
        {
            case TSample::Method::OneP:  { return "1%";   }
            case TSample::Method::FiveP: { return "5%";   }
            case TSample::Method::TenP:  { return "10%";  }
            case TSample::Method::Mean:  { return "Mean"; }
        }
    };
    
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
                         "   Method: %7%\n\n"
                         "   *******************************************\n"
                         "   ***                                     ***\n"
                         "   ***    Statistics before subsampling    ***\n"
                         "   ***                                     ***\n"
                         "   *******************************************\n\n"
                         "   Aligns: %8%\n"
                         "   Depth (Synthetic): %9%\n"
                         "   Depth (Genome):    %10%\n\n"
                         "   *******************************************\n"
                         "   ***                                     ***\n"
                         "   ***    Statistics after subsampling     ***\n"
                         "   ***                                     ***\n"
                         "   *******************************************\n\n"
                         "   Aligns: %11%\n"
                         "   Depth (Synthetic): %12%\n"
                         "   Depth (Genome):    %13%\n";
    
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % "-"
                                            % stats.n_chrT
                                            % stats.n_geno
                                            % o.rChrT
                                            % "??"
                                            % meth2Str()
                                            % "??"
                                            % stats.chrTBefore
                                            % stats.genoBefore
                                            % "??"
                                            % stats.chrTAfter
                                            % stats.genoAfter).str());
    o.writer->close();
}

void TSample::report(const FileName &file, const Options &o)
{
    const auto stats = TSample::stats(file, o);
    
    /*
     * Generating TransSubsample_summary.stats
     */
    
    generateSummary("TransSubsample_summary.stats", stats, o);
    
    /*
     * Generating TransSubsample_sampled.sam
     */
    
    o.generate("TransSubsample_sampled.sam");
    
    WriterSAM writer;
    writer.open(o.work + "/TransSubsample_sampled.sam");

    SamplingTool sampler(1 - stats.prop);

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        const auto *b = reinterpret_cast<bam1_t *>(info.data);
        const auto *h = reinterpret_cast<bam_hdr_t *>(info.header);
        
        if (!align.i)
        {
            /*
             * This is the key, randomly write the reads with certain probability
             */
            
            if (align.cID != ChrT || sampler.select(bam_get_qname(b)))
            {
                writer.write(h, b);
            }
        }
    });
    
    writer.close();
}