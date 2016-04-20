#include "MetaQuin/m_diff.hpp"
#include "MetaQuin/m_assembly.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMFold();

MDiff::Stats MDiff::analyze(const FileName &file_1, const FileName &file_2, const Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MDiff::Stats stats;
    
    // Initalize the sequins
    stats.hist = r.hist();
    
    assert(!o.pA.empty() && !o.pB.empty());
    
    o.info((boost::format("Analyzing: %1%") % o.pA).str());
    stats.align_1 = MBlat::analyze(o.pA);
    
    o.info((boost::format("Analyzing: %1%") % o.pB).str());
    stats.align_2 = MBlat::analyze(o.pB);
    
    /*
     * The implementation is very similar to a single sample. The only difference is that
     * we're interested in the log-fold change.
     */
    
    o.info((boost::format("Analyzing: %1%") % file_1).str());
    const auto dStats_1 = Velvet::analyze<MAssembly::Stats, Contig>(file_1, &stats.align_1);
    
    o.info((boost::format("Analyzing: %1%") % file_2).str());
    const auto dStats_2 = Velvet::analyze<MAssembly::Stats, Contig>(file_2, &stats.align_2);
    
    /*
     * Plot the coverage relative to the known concentration (in attamoles/ul) of each assembled contig.
     */
    
    // Marginal for mixture A
    std::map<SequinID, Coverage> y1;
    
    // Marginal for mixture B
    std::map<SequinID, Coverage> y2;
    
    /*
     * Analyzing for the first sample
     */
    
    o.info("Analyzing for the first sample");
    
    for (const auto &meta : stats.align_1.metas)
    {
        const auto &align = meta.second;
        
        const auto p = MExpress::calculate(stats, stats.align_1, dStats_1, align->id(), *meta.second, o, o.coverage);
        y1[align->id()] = p.y;
    }
    
    o.info((boost::format("Detected %1% sequins in the first sample") % y1.size()).str());
    
    /*
     * Analyzing for the second sample
     */
    
    o.info("Analyzing for the second sample");
    
    for (const auto &meta : stats.align_2.metas)
    {
        const auto &align = meta.second;
        
        const auto p = MExpress::calculate(stats, stats.align_2, dStats_2, align->id(), *meta.second, o, o.coverage);
        y2[align->id()] = p.y;
    }
    
    o.info((boost::format("Detected %1% sequins in the second sample") % y1.size()).str());
    
    /*
     * Merging data, note that we can only do a differential comparison if the sequin appears
     * detected in both samples.
     */
    
    for (const auto &meta : stats.align_1.metas)
    {
        const auto &align = meta.second;
        
        if (!align->contigs.empty())
        {
            // Only when the alignment is detected in both samples
            if (y2.at(align->id()) && y1.at(align->id()))
            {
                // Known concentration
                const auto known = align->seq->abund(Mix_2, false) / align->seq->abund(Mix_1, false);
                
                // Ratio of the marginal concentration
                const auto measured = y2.at(align->id()) / y1.at(align->id());
                
                stats.add(align->id(), known, measured);
                
                SequinDiff d;
                
                d.id    = align->id();
                d.e1    = align->seq->abund(Mix_1, false);
                d.e2    = align->seq->abund(Mix_2, false);
                d.m1    = y1.at(align->id());
                d.m2    = y2.at(align->id());
                d.eFold = known;
                d.mFold = measured;
                
                stats.diffs.insert(d);
            }
        }
    }
    
    stats.n_chrT = dStats_1.contigs.size() + dStats_2.contigs.size();
    stats.n_geno = (dStats_1.n + dStats_2.n) - stats.n_chrT;
    
    // Calculating the absolute detection limit
    stats.limit = r.absolute(stats.hist);
    
    o.info((boost::format("Detected %1% sequin pairs in estimating differential") % stats.size()).str());

    return stats;
}

static void writeCSV(const FileName &file, const MDiff::Stats &stats, const MDiff::Options &o)
{
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    /*
     * Generating detailed statistics for each sequin
     */
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Expected 1 (attomol/ul)"
                                           % "Expected 2 (attomol/ul)"
                                           % "Measured 1 (k-mer average)"
                                           % "Measured 2 (k-mer average)"
                                           % "Expected Fold"
                                           % "Measured Fold"
                                           % "Expected Log-Fold"
                                           % "Measured Log-Fold").str());
    
    for (const auto &diff : stats.diffs)
    {
        o.writer->write((boost::format(format) % diff.id
                                               % diff.e1
                                               % diff.m1
                                               % diff.e2
                                               % diff.m2
                                               % diff.eFold
                                               % diff.mFold
                                               % log2(diff.eFold)
                                               % log2(diff.mFold)).str());
    }
    
    o.writer->close();
}

void MDiff::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto stats = MDiff::analyze(file1, file2, o);

    /*
     * Generating differential comparisons for both samples
     */
    
    o.info("Generating MetaDiff_summary.stats");
    o.writer->open("MetaDiff_summary.stats");
    //o.writer->write(StatsWriter::linearSummary(file1 + " & " + file2, o.rChrT, stats, stats.hist));
    o.writer->close();

    /*
     * 2. Generating differential results
     */
    
    writeCSV("MetaDiff_quin.stats", stats, o);

    /*
     * 3. Generating log-fold plot
     */
    
    o.info("Generating MetaDiff_fold.R");
    o.writer->open("MetaDiff_fold.R");
    o.writer->write(RWriter::createScript("MetaDiff_quins.csv", PlotMFold()));
    o.writer->close();
}
