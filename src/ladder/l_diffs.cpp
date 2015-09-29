#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"
#include <ss/regression/lm.hpp>

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
    o.analyze(fileA);
    const auto a = LAbund::report(fileA, opt);

    opt.mix = Mix_2;
    o.analyze(fileB);
    const auto b = LAbund::report(fileB, opt);

    /*
     * Print a warning message for each sequin detected in B but not in A
     */

    o.info("Checking for sequins in mix B but not in mix A");

    for (const auto &i : b.normalized)
    {
        const auto &id = i.first;
        
        if (!a.normalized.at(id))
        {
            o.warn((boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str());
        }
    }
    
    o.info("Merging mixtures");
    o.info((boost::format("%1% sequins in mix A") % a.normalized.size()).str());
    o.info((boost::format("%1% sequins in mix B") % b.normalized.size()).str());

    /*
     * Try for each detected sequin. But only if it's detected in both mixtures.
     */
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";

    o.writer->open("LadderDifferent_quin.csv");
    o.writer->write((boost::format(format) % "ID"
                                           % "Expect A"
                                           % "Expect B"
                                           % "Expect D"
                                           % "Measure A"
                                           % "Measure B"
                                           % "Norm A"
                                           % "Norm B"
                                           % "Adjust A"
                                           % "Adjust B"
                                           % "Adjust_D").str());

    for (const auto &i : a.normalized)
    {
        // Eg: C_02_C
        const auto &seqID = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!b.normalized.at(seqID))
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
                                                   % "NA").str());
            continue;
        }

        // Known fold change between mixture A and B
        const auto known = (b.expect.at(seqID) / a.expect.at(seqID));

        // Measured fold change between mixture A and B
        const auto adjusted = (b.adjusted.at(seqID) / a.adjusted.at(seqID));

        assert(a.measured.count(seqID) && b.measured.count(seqID));

        stats.add(seqID, log2(known), log2(adjusted));
        
        o.writer->write((boost::format(format) % seqID
                                               % a.expect.at(seqID)
                                               % b.expect.at(seqID)
                                               % known
                                               % a.measured.at(seqID)
                                               % b.measured.at(seqID)
                                               % a.normalized.at(seqID)
                                               % b.normalized.at(seqID)
                                               % a.adjusted.at(seqID)
                                               % b.adjusted.at(seqID)
                                               % adjusted).str());
    }
    
    o.writer->close();

	return stats;
}