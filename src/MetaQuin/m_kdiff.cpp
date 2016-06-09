#include "MetaQuin/m_kdiff.hpp"
#include "MetaQuin/m_abund.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "parsers/parser_stamp.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMFold();

// Defined in resources.cpp // TODO:...
extern Scripts PlotVProb();

MKDiff::Stats MKDiff::analyze(const std::vector<FileName> &files, const Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MKDiff::Stats stats;
    
    // Initalize the sequins
    stats.hist = r.hist();

//    switch (o.soft)
//    {
//        case Software::BWA:
//        case Software::Bowtie:
//        {
//            const auto statsA = MAbund::analyze(files[0]);
//            const auto statsB = MAbund::analyze(files[1]);
//            
//            const auto dataA = statsA.data(false);
//            const auto dataB = statsB.data(false);
//
//            for (const auto &i : r.data())
//            {
//                if (dataA.id2x.count(i.first) && dataB.id2x.count(i.first))
//                {
//                    const auto expected = i.second.concent(Mix_2) / i.second.concent(Mix_1);
//                    const auto measured = dataB.id2y.at(i.first) / dataA.id2y.at(i.first);
//                    
//                    stats.add(i.first, expected, measured);
//                }
//                else
//                {
//                    std::cout << i.first << std::endl;
//                }
//            }
//
//            break;
//        }
//            
//        case Software::STAMP:
//        {
//            ParserStamp::parse(Reader(files[0]), [&](const ParserStamp::Data &x, const ParserProgress &)
//            {
//                const auto m = r.match(x.id);
//                
//                if (m)
//                {
//                    stats.s2p[m->id] = x.p;
//                    stats.s2q[m->id] = x.q;
//                    
//                    const auto expected = m->concent(Mix_2) / m->concent(Mix_1);
//                    const auto measured = exp(x.effect);
//                    
//                    stats.add(m->id, expected, measured);
//                }
//            });
//            
//            break;
//        }
//    }
//    
//    for (auto &seq : stats.hist)
//    {
//        if (!stats.s2p.count(seq.first))
//        {
//            stats.miss.insert(seq.first);
//        }
//    }
    
    return stats;
}

static void writeCSV(const FileName &file, const MKDiff::Stats &stats, const MKDiff::Options &o)
{
//    const auto format = "%1%\t%2%\t%3%";
//    
//    /*
//     * Generating detailed statistics for each sequin
//     */
//    
//    o.writer->open(file);
//    o.writer->write((boost::format(format) % "seq"
//                                           % "input"
//                                           % "measured").str());
//    
//    for (const auto &diff : stats.diffs)
//    {
//        o.writer->write((boost::format(format) % diff.id
//                                               % log2(diff.eFold)
//                                               % log2(diff.mFold)).str());
//    }
//    
//    o.writer->close();
}


static Scripts writeQuery(const MKDiff::Stats &stats, const MKDiff::Options &o)
{
    const auto &r = Standard::instance().r_meta;

    std::stringstream ss;
    
    const auto format = "%1%\t%2%\t%3%";

    for (auto &seq : stats.hist)
    {
        if (stats.s2p.count(seq.first))
        {
            const auto m = r.match(seq.first);
            const auto a = m->concent(Mix_2) / m->concent(Mix_1);
            
            
            
            ss << (boost::format(format) % seq.first % stats.s2p.at(seq.first)
                   
                   % a
                   ).str();
            ss << std::endl;
        }
        else
        {
        }
    }
    
    return ss.str();
}



void MKDiff::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = MKDiff::analyze(files, o);

    /*
     * Generating MetaKDiff_summary.stats
     */
    
    o.generate("MetaKDiff_summary.stats");
    o.writer->open("MetaKDiff_summary.stats");
    o.writer->write(StatsWriter::linearSummary("????", o.rAnnot, stats, stats, stats.hist, "sequins"));
    o.writer->close();

    /*
     * Generating MetaKDiff_sequins.stats
     */
    
    o.generate("MetaKDiff_sequins.stats");
    //writeCSV("MetaKDiff_sequins.stats", stats, o);
    o.writer->open("MetaKDiff_sequins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaKDiff_queries.stats
     */
    
    o.generate("MetaKDiff_queries.stats");
    o.writer->open("MetaKDiff_queries.stats");
    o.writer->write(writeQuery(stats, o));
    o.writer->close();
    
    /*
     * Generating MetaDiff_fold.R
     */
    
    o.generate("MetaDiff_fold.R");
    o.writer->open("MetaDiff_fold.R");
    o.writer->write(RWriter::createScript("MetaDiff_sequins.stats", PlotMFold()));
    o.writer->close();
}
