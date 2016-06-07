#include "VarQuin/v_align.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"

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
    const auto &r = Standard::instance().r_var;

    VAlign::Stats stats;
    stats.data[ChrT].hist = r.baseHist();
    
    if (!r.genoID().empty())
    {
        stats.data[r.genoID()].hist = r.genomeHist();
    }

    return stats;
}

static void classifySynth(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    assert(Standard::isSynthetic(align.cID));

    const auto &r = Standard::instance().r_var;
    const SequinData * match;

    if ((match = r.match(align.l, MatchRule::Contains)))
    {
        stats.data[ChrT].tp++;

        const auto bID = baseID(match->id);
        
        /*
         * It's important to map the position relative to the beginning of the sequin, as interval
         * has been designed for chromosome and thus starts from position 0.
         */
        
        const auto l = inters.find(bID);

        assert(l);
        assert(l->l() == match->l);

        const auto t = Locus(align.l.start - l->l().start, align.l.end - l->l().start);
        assert(t.length() <= l->l().length());

        inters.find(bID)->add(t);
        stats.data[ChrT].hist.at(bID)++;
    }
    else
    {
        stats.data[ChrT].fp++;
        stats.data[ChrT].afp.push_back(align.name);
    }
}

static void classifyGenome(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    const auto &r = Standard::instance().r_var;

    if (matchT(align, __match__, [&](const Locus &l, MatchRule rule)
    {
        return r.findGeno(align.cID, align.l);
    }))
    {
        stats.data[align.cID].tp++;
        stats.data[align.cID].hist.at(__match__.cMatch->id())++;
    }
    else
    {
        stats.data[align.cID].fp++;
    }
}

VAlign::Stats VAlign::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    auto stats = init();
    
    o.analyze(file);

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
        
        inters.add(Interval(baseID(i.first), i.second.l));
    }

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &info)
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
        else if (Standard::isSynthetic(align.cID))
        {
            classifySynth(align, stats, inters);
        }
        else
        {
            classifyGenome(align, stats, inters);
        }
    });

    stats.limit = r.absoluteBase(stats.data[ChrT].hist);

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
            
            i.second.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
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
    
    for (auto &i : stats.data)
    {
        if (i.first == ChrT)
        {
            f(stats.data[ChrT], inters);
        }
        else
        {
            f(stats.data[i.first], r.genoInters());
        }
    }
    
    /*
     * Mapping from sequins to reads
     */
    
    for (const auto &i : stats.data.at(ChrT).hist)
    {
        stats.s2r[i.first] = i.second;
    }
    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VAlign::Stats &stats, const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Number of alignments aligned to the synthetic and genome\n"
                         "   ***\n\n"
                         "   Unmapped:  %2%\n"
                         "   Synthetic: %3% (%4%%%)\n"
                         "   Genome:    %5% (%6%%%)\n\n"
                         "   Dilution:  %7%\n\n"
                         "   ***\n"
                         "   *** Reference annotation\n"
                         "   ***\n\n"
                         "   File: %8%\n\n"
                         "   Synthetic: %9% sequins\n"
                         "   Detected:  %10% sequins\n\n"
                         "   Genome:    %11% genes\n"
                         "   Detected:  - genes\n\n"
                         "   ***                                      \n"
                         "   *** Comparison with synthetic annotation \n"
                         "   ***                                      \n"
                         "   Sensitivity: %12%\n"
                         "   Precision:   %13%\n\n"
                         "   Detection limit: %14% (%15%)\n\n"
                         "   ***                                    \n"
                         "   *** Comparison with genomic annotation \n"
                         "   ***                                    \n"
                         "   Sensitivity: %16%\n"
                         "   Precision:   %17%\n\n";
    
    const auto hasGeno = !r.genoID().empty();

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % src
                                            % stats.n_unmap
                                            % stats.n_chrT
                                            % (100 * stats.chrTProp())
                                            % stats.n_geno
                                            % (100 * stats.endoProp())
                                            % stats.dilution()
                                            % o.rChrT
                                            % stats.data.at(ChrT).hist.size()
                                            % count(stats.data.at(ChrT).hist)
                                            % (!hasGeno ? "-" : toString(r.countInters()))      // 11
                                            % stats.sn(ChrT)                                    // 12
                                            % stats.pc(ChrT)                                    // 13
                                            % stats.limit.abund                                 // 14
                                            % stats.limit.id                                    // 15
                                            % (!hasGeno ? "-" : toString(stats.sn(r.genoID()))) // 16
                                            % (!hasGeno ? "-" : toString(stats.pc(r.genoID()))) // 17
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4%";
    o.writer->write((boost::format(format) % "seq"
                                           % "input"
                                           % "reads"
                                           % "sn"
                                           % "pc").str());

    for (const auto &i : stats.data.at(ChrT).hist)
    {
        o.writer->write((boost::format(format) % i.first
                                               % r.findGene(i.first)->concent()
                                               % stats.s2r.at(i.first)
                                               % stats.sn(ChrT, i.first)
                                               % "-").str());
    }
    
    o.writer->close();
}

static void writeQueries(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "read" % "label").str());

    for (const auto &i : stats.data)
    {
        for (const auto &j : i.second.afp)
        {
            o.writer->write((boost::format(format) % j % "FP").str());
        }
    }

    o.writer->close();
}

void VAlign::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating VarAlign_summary.stats
     */
    
    writeSummary("VarAlign_summary.stats", file, stats, o);

    /*
     * Generating VarAlign_quins.stats
     */
    
    writeQuins("VarAlign_quins.stats", stats, o);

    /*
     * Generating VarAlign_queries.stats
     */
    
    writeQueries("VarAlign_queries.stats", stats, o);
    
    /*
     * Generating VarAlign_report.pdf
     */

    o.report->open("VarAlign_report.pdf");
    o.report->addTitle("VarAlign");
    o.report->addFile("VarAlign_summary.stats");
    o.report->addFile("VarAlign_quins.stats");
    o.report->addFile("VarAlign_queries.stats");
}