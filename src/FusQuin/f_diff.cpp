#include "FusQuin/f_diff.hpp"
#include "FusQuin/FUSQuin.hpp"
#include "FusQuin/f_discover.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

FDiff::Stats FDiff::analyze(const FileName &normal, const FileName &fusion, const Options &o)
{
    FDiff::Stats stats;
    
    const auto &r = Standard::instance().r_fus;

    // Measured abundance for the normal genes
    std::map<SequinID, Counts> normals;
    
    // Measured abundance for the fusion genes
    std::map<SequinID, Counts> fusions;
    
    /*
     * Parse normal junctions
     */
    
    ParserSTab::parse(Reader(normal), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
    {
        //if (c.id == ChrT) { stats.chrT->n_chrT++; }
        //else              { stats.chrT->n_geno++; }

        const SequinData *match;

        if (c.id == ChrT && (match = r.findSplice(c.l)))
        {
            normals[match->id] = c.unique;
            //stats.chrT->h.at(match->id)++;
        }
    });

    /*
     * Parse chimeric junctions
     */
    
    ParserStarFusion::parse(Reader(fusion), [&](const CalledFusion &f, const ParserProgress &)
    {
/*
        const auto r = FClassify::classifyFusion(f, o);

        switch (r.code)
        {
            case FUSQuin::Label::Genome:
            case FUSQuin::Label::GenomeChrT: { stats.chrT->n_geno++; }
            case FUSQuin::Label::Positive:
            case FUSQuin::Label::Negative:   { stats.chrT->n_chrT++; }
        }
        
        if (r.code == FUSQuin::Label::Positive)
        {
            fusions[r.match->id] = f.reads;
            //stats.chrT->h.at(r.match->id)++;
        }
*/
    });

    o.info("Found " + std::to_string(normals.size()) + " introns");
    o.info("Found " + std::to_string(fusions.size()) + " fusions");
    
    /*
     * Compare the genes that are detected in both conditions.
     */
    
    for (const auto &i : normals)
    {
        for (const auto &j : fusions)
        {
            if (r.normalToFusion(i.first) == j.first)
            {
//                // Either the normal ID or fusion ID can be used
//                const auto expected = r.findSpliceChim(i.first);
//                
//                // Measured fold change between normal and fusion gene
//                const auto measured = static_cast<double>(i.second) / j.second;
//                
//                stats.chrT->add(i.first + " - " + j.first, expected->fold(), measured);
            }
        }
    }

    //stats.chrT->ss = Standard::instance().r_fus.limit(stats.chrT->h);

    return stats;
}

void FDiff::report(const FileName &splice, const FileName &chim, const Options &o)
{
    const auto stats = FDiff::analyze(splice, chim, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    //AnalyzeReporter::linear("FusDiff_summary.stats", splice + " & " + chim, stats, "fusions", o.writer);

    /*
     * Generating scatter plot
     */
    
    //AnalyzeReporter::scatter(stats, "", "FusDiff", "Expected Fold", "Measured Fold", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
}