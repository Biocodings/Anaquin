#include "ladder/l_abund.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

static std::vector<double> create(Counts rA, Counts rB, Counts rC, Counts rD, double fold, Counts size)
{
    const auto nA = rA * (fold / size);
    const auto nB = rB * (fold / size);
    const auto nC = rC * (fold / size);
    const auto nD = rD * (fold / size);

    return std::vector<double> { nA, nB, nC, nD };
}

LAbund::Stats LAbund::report(const std::string &file, const Options &o)
{
    LAbund::Stats stats;
    const auto &r = Standard::instance().r_lad;

    // The sequins detected
    std::set<SequinID> seqIDs;

    o.info("Parsing alignment file");

    /*
     * Construct a histogram or distribution
     */
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && (p.i % 1000000) == 0)
        {
            o.wait(std::to_string(p.i));
        }
        
        // Don't repeat the same read if it's spliced
        if (align.i == 0)
        {
            stats.obsTotal++;
            stats.measured[align.id]++;

            // Eg: C_16_A to C_16
            const auto baseID = align.id.substr(0, align.id.find_last_of("_"));

            seqIDs.insert(align.id);
            o.logger->write((boost::format("%1%: %2%") % p.i % baseID).str());
        }
    });

    if (seqIDs.empty())
    {
        o.warn("Empty histogram. Please check the inputs and try again.");
        return stats;
    }

    /*
     * Histogram has been created for detected sequins, such as C_14_C.
     */
    
    assert(stats.obsTotal);
    o.info("Histogram created. Calculating the expected library size.");

    /*
     * Calculate for the expected library size. The size depends on the sequins detected.
     */

    for (const auto &seqID : seqIDs)
    {
        o.info((boost::format("Calculating for sequin: %1%") % seqID).str());

        if (!r.match(seqID))
        {
            o.warn(seqID + " is in alignment but not found in the mixture file");
        }
        else
        {
            const auto typeID = seqID.substr(seqID.size() - 1);

            const static std::map<std::string, double> fold =
            {
                { "A", 1.0 },
                { "B", 2.0 },
                { "C", 4.0 },
                { "D", 8.0 },
            };

            stats.expTotal += fold.at(typeID) * r.match(seqID)->abund(o.mix, false);
        }
    }

    if (!stats.expTotal)
    {
        o.error("stats.expTotal == 0");

        // Report a common and useful error message
        throw std::runtime_error("Unable to find anything in the alignment that matches with the mixture. Usually this is caused by an incorrect mixture file. Please check your mixture file.");
    }

    /*
     * Adjusting the observed abundance
     */

    o.info("Adjusting the measured abundance");

    for (const auto &jID : r.joinIDs())
    {
        const auto A = jID + "_A";
        const auto B = jID + "_B";
        const auto C = jID + "_C";
        const auto D = jID + "_D";

        // Skip over any ladder defined in the mixture but totally undetected
        if (!seqIDs.count(A) && !seqIDs.count(B) && !seqIDs.count(C) && !seqIDs.count(D))
        {
            continue;
        }

        //stats.h[jID]++;

        #define COUNT(x) stats.measured.count(x) ? stats.measured.at(x) : 0
        
        // Create a vector for normalized measured coverage
        const auto normalize = create(COUNT(A), COUNT(B), COUNT(C), COUNT(D), 1.0, stats.obsTotal);
        
        // Create a vector for normalized expected coverage
        const auto expect = create(r.match(A)->mixes.at(Mix_1),
                                   r.match(B)->mixes.at(Mix_1),
                                   r.match(C)->mixes.at(Mix_1),
                                   r.match(D)->mixes.at(Mix_1), 1.0, stats.expTotal);

        // There must be something detected, otherwise we wouldn't have been here
        assert(SS::sum(normalize));

        // Fit a linear regression model
        const auto lm = SS::lm("y ~ x", SS::R::data.frame(SS::R::c(normalize), SS::R::c(expect)));

        // Regression slope that we'll correct to 1
        const auto slope = lm.coeffs[1].value;

        std::vector<double> adjusted;
        adjusted.resize(normalize.size());

        /*
         * This is the key in ladder analysis. Adjust the normalized abundance.
         */
        
        std::transform(normalize.begin(), normalize.end(), adjusted.begin(), [&](double c)
        {
            return c / slope;
        });

        assert(expect[0] && expect[1] && expect[2] && expect[3]);
        
        stats.expect[A]  = expect[0];
        stats.expect[B]  = expect[1];
        stats.expect[C]  = expect[2];
        stats.expect[D]  = expect[3];

        stats.normalized[A] = normalize[0];
        stats.normalized[B] = normalize[1];
        stats.normalized[C] = normalize[2];
        stats.normalized[D] = normalize[3];
        
        stats.adjusted[A] = adjusted[0];
        stats.adjusted[B] = adjusted[1];
        stats.adjusted[C] = adjusted[2];
        stats.adjusted[D] = adjusted[3];

        stats.sequinAdjusted[jID] = adjusted[0] + adjusted[1] + adjusted[2] + adjusted[3];
    }

    o.info("Comparing expected with measured");

    // Try for each detected sequin to form an abundance plot
    for (const auto &i : stats.normalized)
    {
        const auto &seqID = i.first;
        //const auto baseID = s.seq2base.at(seqID);
        const auto known  = stats.expect.at(seqID);
        const auto actual = stats.normalized.at(seqID);

        assert(!isnan(known)  && !isinf(known));
        assert(!isnan(actual) && !isinf(actual));

        o.logInfo((boost::format("0x1234 - %1% %2% %3%") % seqID % known % actual).str());

        stats.add(seqID, log2(known), actual ? log2(actual) : 0);
    }

    o.info("Calculating sensitivity");
    //stats.s = Expression::analyze(stats.h, mix);

    o.info("Generating statistics");
    
    auto writeHist = [&](const std::string &file,
                         const std::map<SequinID, Counts>   &abund,
                         const std::map<SequinID, Coverage> &expect,
                         const std::map<SequinID, Coverage> &actual,
                         const std::map<SequinID, Coverage> &adjust)
    {
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->open(file);
        o.writer->write((boost::format(format) % "ID" % "abund" % "expect" % "observed" % "adjusted").str());

        /*
         * The argument abund is a histogram of abundance before normalization. It's directly taken off from
         * the alignment file. Not all sequins would be detected, in fact anything could be in the histogram.
         */
        
        assert(expect.size() == actual.size() && expect.size() == adjust.size());
        
        for (const auto &i : adjust)
        {
            // Eg: GA322_B
            const auto id = i.first;
            
            if (abund.count(id))
            {
                o.writer->write((boost::format(format) % id
                                                       % abund.at(id)
                                                       % expect.at(id)
                                                       % actual.at(id)
                                                       % adjust.at(id)).str());
            }
            else
            {
                o.writer->write((boost::format(format) % id
                                                       % "NA"
                                                       % "NA"
                                                       % "NA"
                                                       % "NA").str());
            }
        }

        o.writer->close();
    };

    //AnalyzeReporter::linear(stats, "LadderAbundance", "FPKM", o.writer);
    writeHist("LadderAbundance_hist.csv", stats.measured, stats.expect, stats.normalized, stats.adjusted);

	return stats;
}