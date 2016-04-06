#include "LadQuin/l_abund.hpp"
#include "parsers/parser_sam.hpp"
#include <ss/regression/linear.hpp>

extern std::string PlotLadderAbund();

using namespace Anaquin;

static std::vector<double> create(Coverage rA, Coverage rB, Coverage rC, Coverage rD, double fold, Counts size)
{
    const auto nA = rA * (fold / size);
    const auto nB = rB * (fold / size);
    const auto nC = rC * (fold / size);
    const auto nD = rD * (fold / size);

    return std::vector<double> { nA, nB, nC, nD };
}

LAbund::Stats LAbund::analyze(const FileName &file, const Options &o)
{
    LAbund::Stats stats;
    const auto &r = Standard::instance().r_lad;

    // Sequins detected in the experiment
    std::set<SequinID> seqIDs;

    o.analyze(file);

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        else if (!align.mapped || !r.match(align.cID))
        {
            return;
        }
        else if (!align.i)
        {
            stats.data.obsTotal++;
            stats.data.measured[align.cID]++;
            seqIDs.insert(align.cID);
        }

        stats.n_chrT++;
    });

    if (seqIDs.empty())
    {
        o.warn("No sequin detected. Please check and try again.");
        return stats;
    }

    assert(stats.data.obsTotal);

    o.info((boost::format("Detected %1% sequins in reference") % r.data().size()).str());
    o.info((boost::format("Detected %1% sequins in query")     % seqIDs.size()).str());

    o.info("Estimating the expected library size.");

    /*
     * Estimating expected library size. The size depends on the detected sequins.
     */

    o.info("Estimating the expected library size.");

    for (const auto &seqID : seqIDs)
    {
        stats.data.expTotal += r.match(seqID)->abund(o.mix);
    }

    if (!stats.data.expTotal)
    {
        o.error("Expected library size == 0");
        throw std::runtime_error("Error in mixture input. Please check and try again.");
    }

    /*
     * Adjusting the observed abundance
     */

    o.info("Adjusting observations");

    for (const auto &jID : r.joinIDs())
    {
        const auto A = jID + "_A";
        const auto B = jID + "_B";
        const auto C = jID + "_C";
        const auto D = jID + "_D";

        // Skip over any ladder defined but not undetected at all
        if (!seqIDs.count(A) && !seqIDs.count(B) && !seqIDs.count(C) && !seqIDs.count(D))
        {
            continue;
        }

        stats.h_joined[jID]++;
        
        if (seqIDs.count(A)) { stats.hist[A]++; }
        if (seqIDs.count(B)) { stats.hist[B]++; }
        if (seqIDs.count(C)) { stats.hist[C]++; }
        if (seqIDs.count(D)) { stats.hist[D]++; }

        #define COUNT(x) stats.data.measured.count(x) ? stats.data.measured.at(x) : 0
        
        // Create a vector for normalized measured coverage
        const auto normalize = create(COUNT(A), COUNT(B), COUNT(C), COUNT(D), 1.0, stats.data.obsTotal);

        // Create a vector for normalized expected coverage
        const auto expect = create(r.match(A)->abund(o.mix),
                                   r.match(B)->abund(o.mix),
                                   r.match(C)->abund(o.mix),
                                   r.match(D)->abund(o.mix), 1.0, stats.data.expTotal);

        assert(SS::sum(expect));
        
        // There must be something detected, otherwise we wouldn't have been here
        assert(SS::sum(normalize));

        // Fit a linear regression model
        const auto lm = SS::linearModel(normalize, expect);

        // Regression slope that we'll correct to 1
        const auto slope = lm.coeffs[1].est;

        std::vector<double> adjusted;
        adjusted.resize(normalize.size());

        /*
         * This is the key in ladder analysis, adjusting the normalized abundance
         */
        
        std::transform(normalize.begin(), normalize.end(), adjusted.begin(), [&](double c)
        {
            return c / slope;
        });

        assert(expect[0] && expect[1] && expect[2] && expect[3]);

        stats.data.expect[A]  = expect[0];
        stats.data.expect[B]  = expect[1];
        stats.data.expect[C]  = expect[2];
        stats.data.expect[D]  = expect[3];

        stats.data.normalized[A] = normalize[0];
        stats.data.normalized[B] = normalize[1];
        stats.data.normalized[C] = normalize[2];
        stats.data.normalized[D] = normalize[3];
        
        stats.data.adjusted[A] = adjusted[0];
        stats.data.adjusted[B] = adjusted[1];
        stats.data.adjusted[C] = adjusted[2];
        stats.data.adjusted[D] = adjusted[3];

        stats.data.joinAdjusted[jID] = adjusted[0] + adjusted[1] + adjusted[2] + adjusted[3];
    }

    o.info("Comparing expected with measured");

    // Try for each detected sequin to form an abundance plot
    for (const auto &i : stats.data.normalized)
    {
        const auto &seqID = i.first;

        const auto known  = stats.data.expect.at(seqID);
        const auto actual = stats.data.normalized.at(seqID);

        if (!known || !actual)
        {
            continue;
        }
        
        assert(!isnan(known)  && !isinf(known));
        assert(!isnan(actual) && !isinf(actual));

        stats.data.add(seqID, known, actual);
    }

    o.info("Calculating detection limit (joined level)");
    stats.s_joined = r.limitJoin(stats.h_joined);

    o.info("Calculating detection limit (unjoined level)");
    stats.absolute = r.absolute(stats.hist);

  	return stats;
}

static void writeCSV(const FileName &file, const LAbund::Stats &stats, const LAbund::Options &o)
{
    const auto &abund  = stats.data.measured;
    const auto &expect = stats.data.expect;
    const auto &actual = stats.data.normalized;
    const auto &adjust = stats.data.adjusted;

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Abund"    // "Abund (counts)"
                                           % "Expected" // "Expected (attomol/ul)"
                                           % "Observed" // "Observed (normalized counts)"
                                           % "Adjusted" // "Adjusted (normalized counts)"
                                           % "Ratio").str());

    /*
     * The argument abund is a histogram of abundance before normalization. It's directly taken off from
     * the alignment file. Not all sequins would be detected, in fact anything could be in the histogram.
     */
    
    assert(expect.size() == actual.size() && expect.size() == adjust.size());
    
    for (const auto &i : adjust)
    {
        // Eg: GA322_B
        const auto &id = i.first;
        
        if (abund.count(id))
        {
            o.writer->write((boost::format(format) % id
                                                   % abund.at(id)
                                                   % expect.at(id)
                                                   % actual.at(id)
                                                   % adjust.at(id)
                                                   % (adjust.at(id) / actual.at(id))).str());
        }
        else
        {
            o.writer->write((boost::format(format) % id
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA").str());
        }
    }
    
    o.writer->close();
};


void LAbund::report(const FileName &file, const Options &o)
{
    const auto stats = LAbund::analyze(file, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.writer->open("LadderAbundance_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats.data,
                                                "reads"));
    o.writer->close();

    /*
     * Generating CSV for the abundance
     */
    
    writeCSV("LadderAbundance_quins.csv", stats, o);

    /*
     * Generating scatter plot
     */

    o.writer->open("LadderAbundance_scatter.R");
    o.writer->write(RWriter::createScript("LadderAbundance_quins.csv", PlotLadderAbund()));
    o.writer->close();
}