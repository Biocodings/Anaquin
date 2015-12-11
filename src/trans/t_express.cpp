#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tracking.hpp"
#include "parsers/parser_stringtie.hpp"

using namespace SS;
using namespace Anaquin;

template <typename T> void update(TExpress::Stats &stats, const T &t, const GenericID &id, bool isoform, const TExpress::Options &o)
{
    if (t.cID != Standard::chrT)
    {
        stats.n_expT++;
        return;
    }
    
    stats.n_chrT++;
    
    /*
     * There're two possibilities here, comparing at the isoform or gene level. While the workflow is similar, the underlying
     * data-structure is different.
     */

    const auto &r = Standard::instance().r_trans;

    if (isoform)
    {
        const TransData *m = nullptr;
        
        // Try to match by name if possible
        m = r.match(id);
        
        if (!m)
        {
            // Try to match by locus (de-novo assembly)
            m = r.match(t.l, Overlap);
        }
        
        if (!m)
        {
            o.logWarn((boost::format("%1% not found. Unknown isoform.") % id).str());
        }
        else
        {
            stats.h.at(m->id)++;
            
            if (t.fpkm)
            {
                stats.add(id, m->abund(Mix_1), t.fpkm);
            }
        }
    }
    else
    {
        const TransRef::GeneData *m = nullptr;
        
        // Try to match by name if possible
        m = r.findGene(id);
        
        if (!m)
        {
            // Try to match by locus (de-novo assembly)
            m = r.findGene(t.l, Contains);
        }
        
        if (!m)
        {
            o.logWarn((boost::format("%1% not found. Unknown gene.") % id).str());
        }
        else
        {
            stats.h.at(m->id)++;
            
            if (t.fpkm)
            {
                stats.add(id, m->abund(Mix_1), t.fpkm);
            }
        }
    }
}

template <typename Functor> TExpress::Stats calculate(const TExpress::Options &o, Functor f)
{
    TExpress::Stats stats;
    
    const bool isoform = o.level == TExpress::Isoform;
    o.logInfo(isoform ? "Isoform tracking" : "Gene tracking");
    
    f(stats);
    
    return stats;
}

TExpress::Stats TExpress::analyze(const std::vector<Expression> &exps, const Options &o)
{
    return calculate(o, [&](TExpress::Stats &stats)
    {
        const auto &r = Standard::instance().r_trans;
        
        // Construct for a histogram at the appropriate level
        stats.h = o.level == TExpress::Isoform ? r.hist() : r.geneHist();
        
        for (const auto &i : exps)
        {
            update(stats, i, i.id, o.level == TExpress::Isoform, o);
        }
        
        stats.ss = (o.level == TExpress::Isoform) ? r.limit(stats.h) : r.limitGene(stats.h);
    });
}

TExpress::Stats TExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);

    return calculate(o, [&](TExpress::Stats &stats)
    {
        const auto &r = Standard::instance().r_trans;
        
        // Construct for a histogram at the appropriate level
        stats.h = o.level == TExpress::Isoform ? r.hist() : r.geneHist();

        const auto isIsoform = o.level == TExpress::Isoform;
        
        switch (o.tool)
        {
            case Cufflinks:
            {
                ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
                {
                    update(stats, t, isIsoform ? t.trackID : t.id, isIsoform, o);
                });
                
                break;
            }
                
            case StringTie:
            {
                if (isIsoform)
                {
                    ParserStringTie::parseIsoforms(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                    {
                        update(stats, t, t.id, isIsoform, o);
                    });
                }
                else
                {
                    ParserStringTie::parseGenes(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                    {
                        update(stats, t, t.id, isIsoform, o);
                    });
                }

                break;
            }
        }
        
        stats.ss = isIsoform ? r.limit(stats.h) : r.limitGene(stats.h);
    });
}

void TExpress::report(const FileName &file, const Options &o)
{
    const auto stats = TExpress::analyze(file, o);
    const auto units = (o.level == Isoform) ? "isoforms" : "genes";
    
    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    AnalyzeReporter::linear("TransExpress_summary.stats", file, stats, units, o.writer);
    
    /*
     * Generating scatter plot
     */
    
    o.info("Generating scatter plot");
    AnalyzeReporter::scatter(stats, "", "TransExpress", "Expected concentration (attomol/ul)", "Measured coverage (FPKM)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 FPKM)", o.writer);
}