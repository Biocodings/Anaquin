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
        stats.data[ChrT].gtp[bID]++;
    }
    else
    {
        stats.data[ChrT].fp++;
        stats.data[ChrT].afp.push_back(align.name);

        if ((match = r.match(align.l, MatchRule::Overlap)))
        {
            stats.data[ChrT].gfp[baseID(match->id)]++;
        }
    }
}

static void classifyGenome(const Alignment &align, VAlign::Stats &stats, Intervals<> &inters)
{
    const auto &r = Standard::instance().r_var;

    if (Standard::isGenomic(align.cID))
    {
        if (matchT(align, __match__, [&](const Locus &l, MatchRule rule)
        {
            return r.findGeno(align.cID, align.l, rule);
        }))
        {
            const auto gID = __match__.cMatch->id();
            stats.data[align.cID].tp++;
            stats.data[align.cID].hist.at(gID)++;
            stats.data[align.cID].gtp[gID]++;
        }
        else
        {
            /*
             * Here, we know we have the reference chromosome. But what happens if it's outside
             * the regions? It can be TP because it could align to a gene outside the regions.
             * It could also be FP because it might fail to align correctly. We simply don't have
             * the information. Thus, it's only FP if it aligns at least partially with our regions.
             */
            
            if (__match__.oMatch)
            {
                stats.data[align.cID].fp++;
                stats.data[align.cID].gfp[__match__.oMatch->id()]++;
            }
        }
    }
}

VAlign::Stats VAlign::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    auto stats = init();
    o.analyze(file);

    /*
     * We'll need the intervals for measuring the bases. Construcing an interval for each sequin.
     * This enables us to calculate sensitivity.
     */
    
    Intervals<> inters;
    
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
        if (!align.i && info.p.i && !(info.p.i % 1000000))
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
        if (Standard::isSynthetic(i.first))
        {
            f(stats.data[i.first], inters);
        }
        else
        {
            f(stats.data[i.first], r.genoInters());
        }
    }
    
    /*
     * -------------------- Calculating genomic statistics --------------------
     */
    
    /*
     * 1. Covered and length for each sequin
     */

    for (const auto &i : stats.data)
    {
        if (!Standard::isSynthetic(i.first))
        {
            stats.s2c[i.first] = 0;
            stats.s2l[i.first] = 0;
            
            for (const auto &j : i.second.covered)
            {
                stats.s2c[i.first] += j.second;
            }
            
            for (const auto &j : i.second.length)
            {
                stats.s2l[i.first] += j.second;
            }

            assert(stats.s2l[i.first] >= stats.s2c[i.first]);
            
            stats.gtp += i.second.tp;
            stats.gfp += i.second.fp;            
        }
    }
    
    stats.gc = sum(stats.s2c);
    stats.gl = sum(stats.s2l);
    assert(stats.gl >= stats.gl);
    
    /*
     * -------------------- Calculating sequin statistics --------------------
     */
    
    /*
     * 1. Reads aligned to each sequin
     */
    
    for (const auto &i : stats.data.at(ChrT).hist)
    {
        stats.s2r[i.first] = i.second;
    }
    
    /*
     * 2. Covered and length for each sequin
     */

    stats.s2l = stats.data.at(ChrT).length;
    stats.s2c = stats.data.at(ChrT).covered;
    assert(stats.s2l.size() == stats.s2c.size());
    
    /*
     * 3. Sensitivity for each sequin
     */
    
    for (const auto &i : stats.s2l)
    {
        stats.s2s[i.first] = static_cast<Proportion>(stats.s2c.at(i.first)) / stats.s2l.at(i.first);
    }

    assert(stats.s2s.size() == stats.s2l.size());

    /*
     * 4: Precision for each sequin
     */
    
    for (const auto &i : stats.data.at(ChrT).hist)
    {
        #define FROM_MAP(x) x.count(i.first) ? x.at(i.first) : 0
        
        const auto tp = FROM_MAP(stats.data.at(ChrT).gtp);
        const auto fp = FROM_MAP(stats.data.at(ChrT).gfp);
        const auto pc = (tp + fp) ? static_cast<Proportion>(tp) / (tp + fp) : NAN;
        
        assert(pc == NAN || (stats.s2p[i.first] >= 0 && stats.s2p[i.first] <= 1.0));
        stats.s2p[i.first] = pc;
    }

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VAlign::Stats &stats, const VAlign::Options &o)
{
    const auto sums2c = sum(stats.s2c);
    const auto sums2l = sum(stats.s2l);
    
    const auto summary = "-------VarAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       User alignment file: %2%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped:  %3%\n"
                         "       Synthetic: %4% (%5%)\n"
                         "       Genome:    %6% (%7%)\n"
                         "       Dilution:  %8%\n\n"
                         "-------Comparison of alignments to annotation (Synthetic)\n\n"
                         "       *Region level\n"
                         "       Covered:     %9%\n"
                         "       Uncovered:   %10%\n"
                         "       Total:       %11%\n"
                         "       Sensitivity: %12%\n"
                         "       Precision:   %13%\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %14%\n"
                         "       Uncovered:   %15%\n"
                         "       Total:       %16%\n"
                         "       Sensitivity: %17%\n"
                         "       Precision:   %18%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % o.rAnnot
                                            % src
                                            % stats.n_unmap
                                            % stats.n_syn
                                            % (100 * stats.synProp())
                                            % stats.n_gen
                                            % (100 * stats.genProp())
                                            % stats.dilution()
                                            % stats.gc
                                            % (stats.gl - stats.gc)
                                            % stats.gl
                                            % stats.gsn()
                                            % stats.gpc()
                                            % sums2c
                                            % (sums2l - sums2c)
                                            % sums2l
                                            % stats.sn(ChrT)
                                            % stats.pc(ChrT)).str());
    o.writer->close();
}

static void writeQuins(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    o.writer->write((boost::format(format) % "Seq"
                                           % "Length"
                                           % "Reads"
                                           % "Sn"
                                           % "Pc").str());

    for (const auto &i : stats.data.at(ChrT).hist)
    {
        assert(stats.s2s.at(i.first) == stats.sn(ChrT, i.first));
        o.writer->write((boost::format(format) % i.first
                                               % stats.s2l.at(i.first)
                                               % stats.s2r.at(i.first)
                                               % stats.sn(ChrT, i.first)
                                               % stats.s2p.at(i.first)).str());
    }

    o.writer->close();
}

static void writeQueries(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "reads" % "label").str());

    for (const auto &j : stats.data.at(ChrT).afp)
    {
        o.writer->write((boost::format(format) % j % "FP").str());
    }

    o.writer->close();
}

void VAlign::report(const FileName &file, const Options &o)
{
    o.info("Genome: [" + Standard::instance().r_var.genoID() + "]");
    
    const auto stats = analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating VarAlign_summary.stats
     */
    
    writeSummary("VarAlign_summary.stats", file, stats, o);

    /*
     * Generating VarAlign_quins.stats
     */
    
    writeQuins("VarAlign_quins.csv", stats, o);

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
    o.report->addFile("VarAlign_quins.csv");
    o.report->addFile("VarAlign_queries.stats");
}