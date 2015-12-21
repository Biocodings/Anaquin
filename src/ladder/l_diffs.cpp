#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"

using namespace Anaquin;

LDiffs::Stats LDiffs::report(const FileName &fileA, const FileName &fileB, const Options &o)
{
    LDiffs::Stats stats;

    // Copy the pointers across
    auto opt = LAbund::Options();
    
    opt.writer = o.writer;
    opt.logger = o.logger;
    opt.output = o.output;

    opt.mix = Mix_1;
    const auto a = LAbund::analyze(fileA, opt);

    opt.mix = Mix_2;
    const auto b = LAbund::analyze(fileB, opt);

    /*
     * Print a warning message for each sequin detected in B but not in A
     */

    o.info("Checking for sequins in mix B but not in mix A");

    for (const auto &i : b.data.at(ChrT).normalized)
    {
        const auto &id = i.first;
        
        if (!a.data.at(ChrT).normalized.at(id))
        {
            o.warn((boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str());
        }
    }
    
    o.info((boost::format("%1% sequins in mix A") % a.data.at(ChrT).normalized.size()).str());
    o.info((boost::format("%1% sequins in mix B") % b.data.at(ChrT).normalized.size()).str());

    /*
     * Try for each detected sequin. But only if it's detected in both mixtures.
     */
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";

    o.writer->open("LadderDifferent_quin.csv");
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

    for (const auto &i : a.data.at(ChrT).normalized)
    {
        // Eg: C_02_C
        const auto &seqID = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!b.data.at(ChrT).normalized.at(seqID))
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
        const auto known = (b.data.at(ChrT).expect.at(seqID) / a.data.at(ChrT).expect.at(seqID));

        // Measured fold change between mixture A and B
        const auto measured = (b.data.at(ChrT).measured.at(seqID) / a.data.at(ChrT).measured.at(seqID));

        // Normalized fold change between mixture A and B
        const auto normalized = (b.data.at(ChrT).normalized.at(seqID) / a.data.at(ChrT).normalized.at(seqID));

        // Measured adjusted fold change between mixture A and B
        const auto adjusted = (b.data.at(ChrT).adjusted.at(seqID) / a.data.at(ChrT).adjusted.at(seqID));

        assert(a.data.at(ChrT).measured.count(seqID) && b.data.at(ChrT).measured.count(seqID));

        stats.data[ChrT].add(seqID, known, normalized);
        //stats.add(seqID, known, adjusted);

        o.writer->write((boost::format(format) % seqID
                                               % a.data.at(ChrT).expect.at(seqID)
                                               % b.data.at(ChrT).expect.at(seqID)
                                               % known
                                               % a.data.at(ChrT).measured.at(seqID)
                                               % b.data.at(ChrT).measured.at(seqID)
                                               % measured
                                               % a.data.at(ChrT).normalized.at(seqID)
                                               % b.data.at(ChrT).normalized.at(seqID)
                                               % normalized
                                               % a.data.at(ChrT).adjusted.at(seqID)
                                               % b.data.at(ChrT).adjusted.at(seqID)
                                               % adjusted).str());
    }
    
    o.writer->close();

    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");
    //AnalyzeReporter::linear("LadderDifferent_summary.stats", fileA, fileB, a, b, stats, "sequins", o.writer);

    /*
     * Generating Bioconductor
     */
    
    o.info("Generating Bioconductor");
    ////AnalyzeReporter::scatter(stats, "", "LadderDifferent", "Expected concentration (attomol/ul)", "Measured coverage (q)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer, true, false);

	return stats;
}