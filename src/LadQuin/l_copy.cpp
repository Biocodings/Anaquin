#include "LadQuin/l_copy.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_cnvnator.hpp"

using namespace Anaquin;

LCopy::Stats LCopy::analyze(const FileName &file, const Options &o)
{
    LCopy::Stats stats;

  	return stats;
}

static void writeCSV(const FileName &file, const LCopy::Stats &stats, const LCopy::Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Expected"
                                           % "Observed").str());

    for (auto i = 0; i < stats.ids.size(); i++)
    {
        o.writer->write((boost::format(format) % stats.ids[i]
                                               % stats.expected[i]
                                               % stats.measured[i]).str());
    }

    o.writer->close();
};

static void writeBed(const FileName &file, const LCopy::Stats &stats, const LCopy::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%2%-%3%";
    
    o.writer->open(file);
    
    for (const auto &i : stats.detected)
    {
        o.writer->write((boost::format(format) % "chrR"
                                               % i.start
                                               % i.end).str());
    }
    
    o.writer->close();
}

void LCopy::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_lad;
    
    std::vector<ParserBed::Data> x;

    LCopy::Stats stats;
    
    ParserBed::parse(Reader("/Users/tedwong/Desktop/CNVSource/chrF.bed"), [&](const ParserBed::Data &data, const ParserProgress &)
    {
        x.push_back(data);
    });

    ParserCNV::parse(Reader(file), [&](const ParserCNV::Data &data, const ParserProgress &)
    {
        stats.detected.insert(data.l);
        
        int m = std::numeric_limits<int>::max();;
        
        for (auto &i : x)
        {
            const auto lRes = std::abs(i.l.start - data.l.start);
            const auto rRes = std::abs(i.l.end - data.l.end);

            m = std::min((int)m, (int)(lRes + rRes));

            if ((lRes + rRes) < 1000)
            {
                stats.ids.push_back(i.id);
                stats.expected.push_back(r.match(i.id)->abund(o.mix));
                stats.measured.push_back(data.fold);
            }
        }
    });
    
    std::cout << stats.detected.size() << std::endl;
    
    writeCSV("LadCopy_quins.csv", stats, o);
    
    /*
     * Generating a BED file for all duplicates
     */

    writeBed("LadCopy_detected.bed", stats, o);
}



