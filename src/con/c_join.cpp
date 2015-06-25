#include "con/c_join.hpp"
#include "stats/expression.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_sam.hpp"

using namespace Spike;

std::vector<double> create(Counts rA, Counts rB, Counts rC, Counts rD, double fold, Counts size)
{
    const auto nA = rA * (fold / size);
    const auto nB = rB * (fold / size);
    const auto nC = rC * (fold / size);
    const auto nD = rD * (fold / size);

    return std::vector<double> { nA, nB, nC, nD };
}

CJoin::Stats CJoin::analyze(const std::string &file, const Options &options)
{
    CJoin::Stats stats;

    // We'll need it to construct expected library size
    std::set<BaseID> baseIDs;
    
    // Construct a histogram of the aligned sequins
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        // Don't repeat the same read if it's spliced
        if (align.i == 0)
        {
            stats.actTotal++;
            stats.raw[align.id]++;
            baseIDs.insert(align.id.substr(0, align.id.size() - 2));
        }
    });

    assert(!baseIDs.empty());
    assert(stats.actTotal);

    const auto &s = Standard::instance();

    /*
     * Calculate for the expected library size. The size depends on the detected sequins.
     */
    
    for (const auto &id : baseIDs)
    {
        const auto fold = s.c_seqs_A.at(id).abund();

        // Expected size for this particular sequin
        const auto size = 1.0 * fold + 2.0 * fold + 4.0 * fold + 8.0 * fold;
     
        stats.expTotal += size;
    }
    
    assert(stats.expTotal);
    
    for (const auto &i : s.c_seqs_A)
    {
        const std::string base = i.first;

        const auto baseA = base + "_A";
        const auto baseB = base + "_B";
        const auto baseC = base + "_C";
        const auto baseD = base + "_D";
    
        if (!stats.raw.count(baseA) || !stats.raw.count(baseB) || !stats.raw.count(baseC) || !stats.raw.count(baseD))
        {
            continue;
        }
        
        #define COUNT(x) stats.raw.count(x) ? stats.raw.at(x) : 0

        // Create a vector for normalized measured coverage
        const auto normal = create(COUNT(baseA), COUNT(baseB), COUNT(baseC), COUNT(baseD), 1.0, stats.actTotal);

        // Create a vector for normalized expected coverage
        const auto expect = create(1.0, 2.0, 4.0, 8.0, s.c_seqs_A.at(base).abund(), stats.expTotal);

        // Fit a linear regression model
        const auto lm = SS::lm("y ~ x", SS::data.frame(SS::c(normal), SS::c(expect)));

        // Regression slope that we'll correct to 1
        const auto slope = lm.coeffs[1].v;
        
        std::vector<double> correct;
        correct.resize(normal.size());

        std::transform(normal.begin(), normal.end(), correct.begin(), [&](double c)
        {
            return c / slope;
        });

        assert(expect[0] && expect[1] && expect[2] && expect[3]);

        stats.expect[baseA]  = expect[0];
        stats.expect[baseB]  = expect[1];
        stats.expect[baseC]  = expect[2];
        stats.expect[baseD]  = expect[3];

        stats.normal[baseA]  = normal[0];
        stats.normal[baseB]  = normal[1];
        stats.normal[baseC]  = normal[2];
        stats.normal[baseD]  = normal[3];
        
        stats.correct[baseA] = correct[0];
        stats.correct[baseB] = correct[1];
        stats.correct[baseC] = correct[2];
        stats.correct[baseD] = correct[3];
    }
    
    /*
     * Write out histogram
     */

    auto writeHist = [&](const std::string &file,
                         const std::map<SequinID, Coverage> &expect,
                         const std::map<SequinID, Counts> &raw,
                         const std::map<SequinID, Coverage> &correct)
    {
        const std::string format = "%1%\t%2%\t%3%\t%4%";
        
        options.writer->open(file);
        options.writer->write((boost::format(format) % "ID" % "expect" % "raw" % "correct").str());

        assert(raw.size() == correct.size());
        
        for (const auto &i : stats.raw)
        {
            assert(correct.count(i.first));
            options.writer->write((boost::format(format) % i.first
                                                         % expect.at(i.first)
                                                         % raw.at(i.first)
                                                         % correct.at(i.first)).str());
        }

        options.writer->close();
    };

    writeHist("conjoint_histogram.stats", stats.expect, stats.raw, stats.correct);

	return stats;
}