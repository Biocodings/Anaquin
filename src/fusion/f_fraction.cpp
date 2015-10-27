#include "fusion/f_fraction.hpp"
#include "fusion/f_discover.hpp"
#include "fusion/f_classify.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

FFraction::Stats FFraction::stats(const FileName &chim, const FileName &splice, const Options &o)
{
    FFraction::Stats stats;
    const auto &r = Standard::instance().r_fus;

    // Measured abundance for the normal genes
    std::map<SequinID, Counts> normals;
    
    // Measured abundance for the fusion genes
    std::map<SequinID, Counts> fusions;
    
    /*
     * Parse the normal junctions
     */
    
    ParserSTab::parse(Reader(splice), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
    {
        const SequinData *match;

        if ((match = r.findSplice(c.l)))
        {
            normals[match->id] = c.unique;
        }
    });

    /*
     * Parse the chimeric junctions
     */
    
    ParserStarFusion::parse(Reader(chim), [&](const ParserStarFusion::Fusion &f, const ParserProgress &)
    {
        const auto r = FClassify::classifyFusion(f, o);

        if (r.code == FClassify::Code::Positive)
        {
            fusions[r.match->id] = f.reads;
        }
    });

    o.info("Found " + std::to_string(normals.size()) + " introns");
    o.info("Found " + std::to_string(fusions.size()) + " fusions");
    
    /*
     * Compare those related genes that are detected in both conditions
     */
    
    for (const auto &i : normals)
    {
        for (const auto &j : fusions)
        {
            if (r.normalToFusion(i.first) == j.first)
            {
                // Either the normal ID or fusion ID can be used
                const auto expected = r.findSpliceChim(i.first);
                
                // Measured fold change between normal and fusion gene
                const auto measured = static_cast<double>(i.second) / j.second;
                
                stats.add(i.first + " - " + j.first, expected->fold(), measured);
            }
        }
    }

    return stats;
}

void FFraction::report(const FileName &splice, const FileName &chim, const Options &o)
{
    const auto stats = FFraction::stats(splice, chim, o);

    /*
     * Generating summary statistics
     */

    o.info("Generating summary statistics");
    AnalyzeReporter::linear("FusionFraction_summary.stats", splice + " & " + chim, stats, "fusions", o.writer);

    /*
     * Generating Bioconductor
     */
    
    o.info("Generating Bioconductor");
    AnalyzeReporter::scatter(stats, "", "FusionFraction", "Expected Fold", "Measured Fold", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
}