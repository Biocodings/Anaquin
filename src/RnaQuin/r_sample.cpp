#include "tools/sample.hpp"
#include <ss/maths/stats.hpp>
#include "RnaQuin/r_sample.hpp"
#include "writers/writer_sam.hpp"
#include "parsers/parser_cufflink.hpp"

using namespace Anaquin;

RSample::Stats RSample::stats(const FileName &file, const Options &o)
{
    RSample::Stats stats;
    
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
                    stats.n_syn++;
                    stats.chrT.push_back(x.abund);
                }
                else
                {
                    stats.n_gen++;
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
        case Method::_1:  { stats.prop = 0.01; break; }
        case Method::_5:  { stats.prop = 0.05; break; }
        case Method::_10: { stats.prop = 0.10; break; }
        case Method::_15: { stats.prop = 0.15; break; }
        case Method::_20: { stats.prop = 0.20; break; }
        case Method::_50: { stats.prop = 0.50; break; }
            
            
/*
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
 */
    }

    assert(stats.prop > 0 && stats.prop < 1.0);
    
    o.info("Proportion: " + toString(stats.prop));
    
    return stats;
}

static void generateSummary(const FileName &file, const RSample::Stats &stats, const RSample::Options &o)
{
    o.generate(file);
    
    auto meth2Str = [&]()
    {
        switch (o.meth)
        {
            case RSample::Method::_1:  { return "1%";  }
            case RSample::Method::_5:  { return "5%";  }
            case RSample::Method::_10: { return "10%"; }
            case RSample::Method::_15: { return "15%"; }
            case RSample::Method::_20: { return "20%"; }
            case RSample::Method::_50: { return "50%"; }
        }
    };
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Proportion of alignments mapped to the synthetic and genome\n"
                         "   ***\n\n"
                         "   Unmapped:  %2% reads\n"
                         "   Synthetic: %3% reads\n"
                         "   Genome:    %4% reads\n\n"
                         "   ***\n"
                         "   *** Reference annotation\n"
                         "   ***\n\n"
                         "   File: %5%\n\n"
                         "   Reference: %6% sequins\n\n"
                         "   ***                                               \n"
                         "   *** How the coverage is calculated?               \n"
                         "   ***                                               \n"
                         "   *** Please refer to the documentation for details \n"
                         "   ***                                               \n\n"
                         "   Method: %7%\n\n"
                         "   ***                                     \n"
                         "   ***    Statistics before subsampling    \n"
                         "   ***                                     \n\n"
                         "   Aligns: %8%\n"
                         "   Coverage (Synthetic): %9%\n"
                         "   Coverage (Genome):    %10%\n\n"
                         "   ***                                     \n"
                         "   ***    Statistics after subsampling     \n"
                         "   ***                                     \n\n"
                         "   Aligns: %11%\n"
                         "   Coverage (Synthetic): %12%\n"
                         "   Coverage (Genome):    %13%\n";
    
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % "-"
                                            % stats.n_syn
                                            % stats.n_gen
                                            % o.rAnnot
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

void RSample::report(const FileName &file, const Options &o)
{
    const auto stats = RSample::stats(file, o);
    
    /*
     * Generating RnaSubsample_summary.stats
     */
    
    generateSummary("RnaSubsample_summary.stats", stats, o);
    
    /*
     * Generating RnaSubsample_sampled.sam
     */
    
    o.generate("RnaSubsample_sampled.sam");
    
    WriterSAM writer;
    writer.open(o.work + "/RnaSubsample_sampled.sam");

    SamplingTool sampler(1 - stats.prop);

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &info)
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