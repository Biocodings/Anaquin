#include "con/c_single.hpp"
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

CSingle::Stats CSingle::analyze(const std::string &file, const Options &options)
{
    CSingle::Stats stats;

    // We'll need it to construct expected library size
    std::set<BaseID> baseIDs;
    
    // Construct a histogram of the aligned sequins
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        // Don't repeat the same read if it's spliced
        if (align.i == 0)
        {
            stats.actTotal++;
            stats.abund[align.id]++;
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
    
        if (!stats.abund.count(baseA) ||
            !stats.abund.count(baseB) ||
            !stats.abund.count(baseC) ||
            !stats.abund.count(baseD))
        {
            continue;
        }
        
        #define COUNT(x) stats.abund.count(x) ? stats.abund.at(x) : 0

        // Create a vector for normalized measured coverage
        const auto actual = create(COUNT(baseA), COUNT(baseB), COUNT(baseC), COUNT(baseD), 1.0, stats.actTotal);

        // Create a vector for normalized expected coverage
        const auto expect = create(1.0, 2.0, 4.0, 8.0, s.c_seqs_A.at(base).abund(), stats.expTotal);

        // Fit a linear regression model
        const auto lm = SS::lm("y ~ x", SS::data.frame(SS::c(actual), SS::c(expect)));

        // Regression slope that we'll correct to 1
        const auto slope = lm.coeffs[1].v;
        
        std::vector<double> correct;
        correct.resize(actual.size());

        std::transform(actual.begin(), actual.end(), correct.begin(), [&](double c)
        {
            return c / slope;
        });

        assert(expect[0] && expect[1] && expect[2] && expect[3]);

        stats.expect[baseA]  = expect[0];
        stats.expect[baseB]  = expect[1];
        stats.expect[baseC]  = expect[2];
        stats.expect[baseD]  = expect[3];

        stats.actual[baseA]  = actual[0];
        stats.actual[baseB]  = actual[1];
        stats.actual[baseC]  = actual[2];
        stats.actual[baseD]  = actual[3];
        
        stats.correct[baseA] = correct[0];
        stats.correct[baseB] = correct[1];
        stats.correct[baseC] = correct[2];
        stats.correct[baseD] = correct[3];

        stats.s_correct[base] = correct[0] + correct[1] + correct[2] + correct[3];
    }
    
    /*
     * Write out histogram
     */

    auto writeHist = [&](const std::string &file,
                         const std::map<SequinID, Counts>   &abund,
                         const std::map<SequinID, Coverage> &expect,
                         const std::map<SequinID, Coverage> &actual,
                         const std::map<SequinID, Coverage> &correct)
    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        options.writer->open(file);
        options.writer->write((boost::format(format) % "ID" % "abund" % "expect" % "actual" % "correct").str());

        std::cout << abund.size() << " " << correct.size() << " " << expect.size() << " " << actual.size() << std::endl;

        assert(expect.size() == correct.size());
        assert(expect.size() == actual.size());

        for (const auto &i : stats.abund)
        {
            if (!correct.count(i.first))
            {
                continue;
            }
            
            assert(correct.count(i.first));
            options.writer->write((boost::format(format) % i.first
                                                         % abund.at(i.first)
                                                         % expect.at(i.first)
                                                         % actual.at(i.first)
                                                         % correct.at(i.first)).str());
        }

        options.writer->close();
    };

    writeHist("conjoint.stats", stats.abund, stats.expect, stats.actual, stats.correct);

	return stats;
}