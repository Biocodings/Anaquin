#include "LadQuin/l_diff.hpp"

using namespace Anaquin;

LDiff::Stats LDiff::analyze(const FileName &fileA, const FileName &fileB, const Options &o)
{    
    LDiff::Stats stats;
    
    auto opt = LNorm::Options();
    
    opt.writer = o.writer;
    opt.logger = o.logger;
    opt.output = o.output;
    
    opt.mix = Mix_1;
    stats.a = LNorm::analyze(fileA, opt);
    
    opt.mix = Mix_2;
    stats.b = LNorm::analyze(fileB, opt);

    /*
     * Print a warning message for each sequin detected in B but not in A
     */
    
    o.info("Checking for sequins in mix B but not in mix A");
    
    for (const auto &i : stats.b.data.normalized)
    {
        const auto &id = i.first;
        
        if (!stats.a.data.normalized.at(id))
        {
            o.warn((boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str());
        }
    }
    
    o.info((boost::format("%1% sequins in mix A") % stats.a.data.normalized.size()).str());
    o.info((boost::format("%1% sequins in mix B") % stats.b.data.normalized.size()).str());

    return stats;
}

void LDiff::report(const FileName &fileA, const FileName &fileB, const Options &o)
{
    const auto stats = LDiff::analyze(fileA, fileB, o);
    
    /*
     * Try for each detected sequin. But only if it's detected in both mixtures.
     */
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";

    o.writer->open("LadDiff_quins.stats");
    o.writer->write((boost::format(format) % "ID"
                                           % "Expected A"
                                           % "Expected B"
                                           % "Expected B/A"
                                           % "Measured A"
                                           % "Measured B"
                                           % "Measured B/A"
                                           % "Normalized A"
                                           % "Normalized B"
                                           % "Normalized B/A"
                                           % "Adjusted A"
                                           % "Adjusted B"
                                           % "Adjusted B/A").str());

    for (const auto &i : stats.a.data.normalized)
    {
        // Eg: C_02_C
        const auto &seqID = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!stats.b.data.normalized.at(seqID))
        {
            o.writer->write((boost::format(format) % seqID
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA").str());
            continue;
        }

        // Known fold change between mixture A and B
        const auto known = (stats.b.data.expect.at(seqID) / stats.a.data.expect.at(seqID));

        // Measured fold change between mixture A and B
        const auto measured = (stats.b.data.measured.at(seqID) / stats.a.data.measured.at(seqID));

        // Normalized fold change between mixture A and B
        const auto normalized = (stats.b.data.normalized.at(seqID) / stats.a.data.normalized.at(seqID));

        // Measured adjusted fold change between mixture A and B
        const auto adjusted = (stats.b.data.adjusted.at(seqID) / stats.a.data.adjusted.at(seqID));

        assert(stats.a.data.measured.count(seqID) && stats.b.data.measured.count(seqID));

        //stats.data[ChrT].add(seqID, known, normalized);
        //stats.add(seqID, known, adjusted);

        o.writer->write((boost::format(format) % seqID
                                               % stats.a.data.expect.at(seqID)
                                               % stats.b.data.expect.at(seqID)
                                               % known
                                               % stats.a.data.measured.at(seqID)
                                               % stats.b.data.measured.at(seqID)
                                               % measured
                                               % stats.a.data.normalized.at(seqID)
                                               % stats.b.data.normalized.at(seqID)
                                               % normalized
                                               % stats.a.data.adjusted.at(seqID)
                                               % stats.b.data.adjusted.at(seqID)
                                               % adjusted).str());
    }
    
    o.writer->close();

    /*
     * Generating summary statistics
     */

    o.info("Generating summary statistics");
    //AnalyzeReporter::linear("LadderDifferent_summary.stats", fileA, fileB, a, b, stats, "sequins", o.writer);

    /*
     * Generating LadDiff_diff.R
     */
    
    o.info("Generating LadDiff_diff.R");
    ////AnalyzeReporter::scatter(stats, "", "LadderDifferent", "Expected concentration (attomol/ul)", "Measured coverage (q)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer, true, false);
}