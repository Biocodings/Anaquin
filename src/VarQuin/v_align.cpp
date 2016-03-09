#include "VARQuin/v_align.hpp"
#include "VARQuin/v_sample.hpp"
#include "parsers/parser_sam.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

template <typename T> struct AlignmentMatch
{
    inline bool contains() const { return !cMatch; }
    inline bool overlaps() const { return !oMatch; }

    T * oMatch = nullptr;
    T * cMatch = nullptr;

    // Number of bases to the left of the reference (overlaps == true)
    Base lGaps = 0;

    // Number of bases to the right of the reference (overlaps == true)
    Base rGaps = 0;
};

static AlignmentMatch<Interval> __match__;

template <typename T, typename F> const AlignmentMatch<T> * matchT(const Alignment &align, AlignmentMatch<T> &match, F f)
{
    match = AlignmentMatch<T>();
    
    if ((match.cMatch = f(align.l, MatchRule::Contains)))
    {
        match.oMatch = match.cMatch;
    }
    else if ((match.oMatch = f(align.l, MatchRule::Overlap)))
    {
        // Empty Implementation
    }

    auto x = __match__.cMatch ? __match__.cMatch : __match__.oMatch;
    
    if (x)
    {
        x->map(align.l, &match.lGaps, &match.rGaps);
    }

    return match.cMatch ? &match : nullptr;
}

static VAlign::Stats init()
{
    VAlign::Stats stats;

    const auto &r = Standard::instance().r_var;

    stats.data[ChrT].hist = r.hist();
    stats.data[Endo].hist = r.endoHist();
    
    return stats;
}

static void classifyChrT(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    const auto &r = Standard::instance().r_var;
    const SequinData * match;

    if ((match = r.match(align.l, MatchRule::Contains)))
    {
        stats.data[ChrT].tp++;
        
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
        stats.data[ChrT].hist.at(match->id)++;
    }
    else
    {
        stats.data[ChrT].fp++;
    }
}

static void classifyEndo(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    const auto &r = Standard::instance().r_var;

    if (matchT(align, __match__, [&](const Locus &l, MatchRule rule)
    {
        return r.findEndo(align.cID, align.l);
    }))
    {
        stats.data[Endo].tp++;
        stats.data[Endo].hist.at(__match__.cMatch->id())++;
    }
    else
    {
        stats.data[Endo].fp++;
    }
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
        else if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        stats.update(align);
        
        if (!align.mapped)
        {
            return;
        }        
        else if (align.cID == ChrT)
        {
            classifyChrT(align, stats, inters);
        }
        else
        {
            classifyEndo(align, stats, inters);
        }
    });

    stats.limit = r.limit(stats.data[ChrT].hist);

    /*
     * Calculating interval statistics
     */

    o.info("Calculating interval statistics");

    auto f = [&](Stats::Data &data, const Intervals<> &inters)
    {
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
            
            data.covered[i.first] = covered;
            data.length [i.first] = i.second.l().length();

            totCov += covered;
            totLen += i.second.l().length();
            
            assert(totCov <= totLen);
        }
        
        assert(totLen == inters.length());
    };
    
    f(stats.data[ChrT], inters);
    f(stats.data[Endo], r.endoInters());

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VAlign::Stats &stats, const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Fraction of reads mapped to the synthetic and experimental chromosomes\n"
                         "   ***\n\n"
                         "   Unmapped:   %2% reads\n"
                         "   Synthetic:  %3% (%4%%%) reads\n"
                         "   Experiment: %5% (%6%%%) reads\n\n"
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
                         "   Detection limit: %14% (%15%)\n\n"
                         "   *************************************************\n"
                         "   ***                                           ***\n"
                         "   ***    Comparison with endogenous annotation   ***\n"
                         "   ***                                           ***\n"
                         "   *************************************************\n\n"
                         "   Sensitivity:  %16%\n"
                         "   Specificity:  %17%\n\n";
    
    o.writer->open(file);
    o.writer->write((boost::format(summary) % src
                                            % stats.unmapped
                                            % stats.n_chrT
                                            % stats.chrTProp()
                                            % stats.n_endo
                                            % stats.endoProp()
                                            % stats.dilution()
                                            % o.rChrT
                                            % r.countSeqs()
                                            % o.rEndo            // 10
                                            % r.countInters()    // 11
                                            % stats.sn(ChrT)     // 12
                                            % stats.pc(ChrT)     // 13
                                            % stats.limit.abund  // 14
                                            % stats.limit.id     // 15
                                            % stats.sn(ChrT)     // 16
                                            % stats.pc(ChrT)     // 17
                     ).str());
    o.writer->close();
}

void VAlign::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Write out summary statistics
     */
    
    writeSummary("VarAlign_summary.stats", file, stats, o);

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
