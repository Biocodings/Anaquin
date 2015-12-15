#include "trans/t_diffs.hpp"
#include "parsers/parser_cdiffs.hpp"

using namespace SS;
using namespace Anaquin;

template <typename T> void update(TDiffs::Stats &stats, const T &t, const GenericID &id, bool isoform, const TDiffs::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    
    if (t.cID != Standard::chrT)
    {
        stats.chrT->n_expT++;
        return;
    }
    
    stats.chrT->n_chrT++;
    
    // The known and observed fold-change
    Fold known = NAN;
    
    // It's NAN if the sequin defined in reference but not in mixture
    Fold measured = NAN;
    
    /*
     * Differential expression at the gene level
     */
    
    auto g = [&](const GeneID &id, double fpkm_1, double fpkm_2)
    {
        // The known and observed fold-change
        Fold known = NAN;
        
        // It's NAN if the sequin defined in reference but not in mixture
        Fold measured = NAN;
        
        const auto *g = r.findGene(id);
        
        if (g)
        {
            // Calculate the known fold-change between B and A
            known = (g->abund(Mix_2) / g->abund(Mix_1));
        }
        
        if (g && !isnan(fpkm_1) && !isnan(fpkm_2) && fpkm_1 && fpkm_2)
        {
            stats.chrT->h.at(id)++;
            
            // Measured fold-change between the two mixtures
            measured = fpkm_2 / fpkm_1;
        }
        
        stats.chrT->add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
    };

    switch (o.level)
    {
        case TDiffs::Gene:
        {
            if ((t.status != NoTest) && stats.chrT->h.count(t.id))
            {
                g(t.id, t.fpkm_1, t.fpkm_2);
            }
            
            break;
        }
            
        case TDiffs::Isoform:
        {
            if ((t.status == NoTest) || !stats.chrT->h.count(id))
            {
                return;
            }
            
            const auto *seq = r.match(id);
            
            if (seq)
            {
                // Known fold-change between the two mixtures
                known = seq->abund(Mix_2) / seq->abund(Mix_1);
            }
            
            if ((t.status != NoTest) && t.fpkm_1 && t.fpkm_2)
            {
                stats.chrT->h.at(id)++;
                
                // Measured fold-change between the two mixtures
                measured = t.fpkm_2 / t.fpkm_1;
            }
            
            //stats.chrT->add(t.testID, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            stats.chrT->add(id, !isnan(known) ? known : NAN, !isnan(measured) ? measured : NAN);
            
            break;
        }
    }
}

template <typename Functor> TDiffs::Stats calculate(const TDiffs::Options &o, Functor f)
{
    TDiffs::Stats stats;
    const auto &r = Standard::instance().r_trans;
    
    stats.chrT = std::shared_ptr<TDiffs::Stats::ChrT>(new TDiffs::Stats::ChrT());
    
    const auto isoform = o.level == TDiffs::Isoform;
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    // Construct for a histogram at the appropriate level
    stats.chrT->h = isoform ? r.hist() : r.geneHist();
    
    o.info("Parsing tracking file");

    f(stats);
    
    o.info("Calculating limit of sensitivity");
    
    stats.chrT->ss = isoform ? r.limit(stats.chrT->h) : r.limitGene(stats.chrT->h);
    
    return stats;
}

TDiffs::Stats TDiffs::analyze(const std::vector<DiffTest> &tests, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        for (auto &test : tests)
        {
            update(stats, test, test.id, o.level == TDiffs::Isoform, o);
        }
    });
}

TDiffs::Stats TDiffs::analyze(const FileName &file, const Options &o)
{
    return calculate(o, [&](TDiffs::Stats &stats)
    {
        switch (o.soft)
        {
            case Cuffdiffs:
            {
                const auto isIsoform = o.level == TDiffs::Isoform;
                
                ParserCDiffs::parse(file, [&](const TrackingDiffs &t, const ParserProgress &)
                {
                    update(stats, t, isIsoform ? t.testID : t.id, isIsoform, o);
                });

                break;
            }

            case EdgeR:
            case DESeq2:
            {
                throw "Not implemented";
            }
        }
    });
}

std::vector<TDiffs::Stats> TDiffs::analyze(const std::vector<FileName> &files, const Options &o)
{
    std::vector<TDiffs::Stats> stats;
    
    for (const auto &file : files)
    {
        stats.push_back(analyze(file, o));
    }

    return stats;
}

void TDiffs::report(const FileName &file, const Options &o)
{
    const auto stats = TDiffs::analyze(file, o);
    const auto units = (o.level == Isoform) ? "isoforms" : "genes";
    
    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    AnalyzeReporter::linear("TransDiff_summary.stats", file, stats, units, o.writer);
    
    /*
     * Generating scatter plot
     */
    
    o.info("Generating scatter plot");
    AnalyzeReporter::scatter(stats, "", "TransDiff", "Expected fold change of mixture A and B", "Measured fold change of mixture A and B", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
}

void TDiffs::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = analyze(files, o);
    const auto units = (o.level == Isoform) ? "isoforms" : "genes";

    /*
     * Generating summary statistics for each replicate
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        const auto file = (boost::format("TransDiff_%1%_summary.stats") % files[i]).str();
        AnalyzeReporter::linear(file, files[i], stats[i], units, o.writer);
    }
    
    /*
     * Generating scatter plots for each replicate
     */

    for (auto i = 0; i < files.size(); i++)
    {
        const auto file = (boost::format("TransDiff_%1%_summary.stats") % files[i]).str();
        AnalyzeReporter::scatter(stats[i], "", file, "Expected fold change of mixture A and B", "Measured fold change of mixture A and B", "Expected log2 fold change of mixture A and B", "Expected log2 fold change of mixture A and B", o.writer);
    }
    
    /*
     * Generating summary statistics for all replicates
     */
    
}