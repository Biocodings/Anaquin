#include "VARQuin/v_align.hpp"
#include "VARQuin/v_sample.hpp"
#include "parsers/parser_sam.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

//template <typename T> static void sums(const std::map<T, Counts> &m, Counts &c)
//{
//    for (const auto &i : m)
//    {
//        if (i.second == 0)
//        {
//            c++;
//        }
//        else
//        {
//            c += i.second;
//        }
//    }
//    
//    assert(c);
//}

static VAlign::Stats init()
{
    VAlign::Stats stats;
    
    stats.data[ChrT];
    stats.data["chr21"];

    // Distribution for the sequins
    stats.h = Standard::instance().r_var.hist();
    
    return stats;
}

static void classifyChrT(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    const SequinData * match;

    if (classify(stats.data.at(ChrT).m, align, [&](const Alignment &)
    {
        return (match = Standard::instance().r_var.match(align.l, MatchRule::Contains)) ? Positive : Negative;
    }))
    {
        /*
         * It's important to map the position relative to the beginning of the sequin, as interval
         * has been designed for chromosome and thus starts from position 0.
         */
        
        const auto l = inters.find(match->id);

        assert(l);
        assert(l->l() == match->l);

        const auto t = Locus(align.l.start - l->l().start, align.l.end - l->l().start);
        assert(t.length() <= l->l().length());

        inters.find(match->id)->add(t);
        stats.h.at(match->id)++;
    }
}

static void classifyEndo(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    // Empty Implementation
}

VAlign::Stats VAlign::analyze(const FileName &file, const Options &o)
{
    auto stats = init();
    
    o.analyze(file);

    const auto &r = Standard::instance().r_var;

    /*
     * We'll need the intervals for measuring at the base level.
     */
    
    Intervals<> inters;
    
    /*
     * Construcing an interval each sequin. This enables us to calculate sensitivity for each sequin.
     */

    for (auto &i : r.data())
    {
        if (boost::algorithm::ends_with(i.first, "_V"))
        {
            continue;
        }
        
        inters.add(Interval(i.first, i.second.l));
    }

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (align.spliced)
        {
            o.warn("Splice read: " + align.name + " detected");
            return;
        }
        
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        stats.update(align);
        
        if (!align.mapped)
        {
            return;
        }
        
        if (align.cID == ChrT)
        {
            classifyChrT(align, stats, inters);
        }
        else if (align.cID == "chr21")
        {
            classifyEndo(align, stats, inters);
        }
    });

    /*
     * Calculating limit of sensitivity
     */
    
    // Calculate for the sensitivity
    stats.limit = r.limit(stats.h);

    /*
     * Calculating base statistics for both synthetic and endogenous
     */

    o.info("Calculating base statistics");

    Base totCov = 0;
    Base totLen = 0;
    
    for (const auto &i: inters.data())
    {
        Base covered = 0;
        
        i.second.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                covered += j - i;
            }
        });

        stats.data[ChrT].covered[i.first] = covered;
        stats.data[ChrT].length [i.first] = i.second.l().length();

        totCov += covered;
        totLen += i.second.l().length();

        assert(totCov <= totLen);
    }
    
    assert(totLen == inters.length());
    return stats;
}

static void writeSummary(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;

    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Fraction of reads mapped to the synthetic and experimental chromosomes\n"
                         "   ***\n\n"
                         "   Unmapped:   %2% reads\n"
                         "   Experiment: %3% (%4%%%) reads\n"
                         "   Synthetic:  %5% (%6%%%) reads\n\n"
                         "   Dilution:   %7%\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %8%\n\n"
                         "   Synthetic: %9% sequins\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Experiment)\n"
                         "   ***\n\n"
                         "   File: %10%\n\n"
                         "   Experiment: %11% genes\n\n"
                         "   *************************************************\n"
                         "   ***                                           ***\n"
                         "   ***    Comparison with synthetic annotation   ***\n"
                         "   ***                                           ***\n"
                         "   *************************************************\n\n"
                         "   Sensitivity:  %12%\n"
                         "   Specificity:  %13%\n\n"
                         "   Base Covered: %14%\n\n"
                         "   Detection limit: %15% (%16%)\n\n"
                         "   *************************************************\n"
                         "   ***                                           ***\n"
                         "   ***    Comparison with endogenous annotation   ***\n"
                         "   ***                                           ***\n"
                         "   *************************************************\n\n"
                         "   Sensitivity:  %17%\n"
                         "   Specificity:  %18%\n\n"
                         "   Base Covered: %19%\n\n";
    
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.unmapped
                                            % stats.n_chrT
                                            % stats.chrTProp()
                                            % stats.n_endo
                                            % stats.endoProp()
                                            % stats.dilution()
                                            % o.rChrT()
                                            % r.countSeqs()
                                            % o.rEndo()         // 10
                                            % "NA"              // 11
                                            % "NA"              // 12
                                            % "NA"              // 13
                                            % stats.bSN(ChrT)   // 14
                                            % stats.limit.abund // 15
                                            % stats.limit.id    // 16
                                            % "NA"              // 17
                                            % "NA"
                                            % stats.bSN("chr21")).str());
    o.writer->close();
}

void VAlign::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Write out summary statistics
     */
    
    writeSummary("VarAlign_summary.stats", stats, o);
    

    /*
     * Generating alignment statistics
     */
    
    
    /*
     * Generating sequin statistics
     */
    
    o.writer->open("VarAlign_quins.stats");
    o.writer->write((boost::format("Summary for input: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "ID" % "Counts (reads)").str());
    
//    for (const auto &i : stats.h)
//    {
//        o.writer->write((boost::format(format) % i.first % stats.h.at(i.first)).str());
//    }
    
    o.writer->close();
}
