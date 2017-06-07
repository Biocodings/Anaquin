#include "tools/tools.hpp"
#include "VarQuin/v_trim.hpp"
#include "writers/writer_sam.hpp"
#include "parsers/parser_bambed.hpp"

extern Anaquin::FileName BedRef();

using namespace Anaquin;

VTrim::Stats VTrim::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    o.analyze(file);
    VTrim::Stats stats;
    
    // Regions without edge effects
    const auto regs = r.regs1();

    const auto shouldL = o.meth == Method::Left  || o.meth == Method::LeftRight;
    const auto shouldR = o.meth == Method::Right || o.meth == Method::LeftRight;
    
    std::vector<DInter *> multi;
    
    // Check trimming reads...
    ParserBAMBED::parse(file, regs, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *inter)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logWait(std::to_string(info.p.i));
        }

        multi.clear();
        const auto m = x.mapped && regs.count(x.cID) ? regs.at(x.cID).contains(x.l, &multi) : nullptr;
        
        if (m)
        {
            std::sort(multi.begin(), multi.end(), [&](const DInter * x, const DInter * y)
            {
                return x->l().length() < y->l().length();
            });

            // The smallest region
            const auto m = multi.front();
            
            const auto lTrim = std::abs(x.l.start - m->l().start) <= o.trim;
            const auto rTrim = std::abs(x.l.end - m->l().end) <= o.trim;
            
            if (shouldL && lTrim) { stats.lTrim.insert(x.name); }
            if (shouldR && rTrim) { stats.rTrim.insert(x.name); }
        }

        stats.before++;
        return ParserBAMBED::Response::OK;
    });
    
    stats.left  = stats.lTrim.size();
    stats.right = stats.rTrim.size();
    
    WriterSAM writer;
    writer.openTerm();
    
    // Triming away the paired reads ...
    ParserBAMBED::parse(file, regs, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *inter)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logWait(std::to_string(info.p.i));
        }

        if (!stats.lTrim.count(x.name) && !stats.rTrim.count(x.name))
        {
            stats.after++;
            writer.write(x);
        }
        
        return ParserBAMBED::Response::OK;
    });
    
    stats.nRegs = countMap(regs, [&](ChrID, const DIntervals<> &x)
    {
        return x.size();
    });

    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &src,
                         const VTrim::Stats &stats,
                         const VTrim::Options &o)
{
    o.generate(file);
    
    auto meth2Str = [&]()
    {
        switch (o.meth)
        {
            case VTrim::Method::Left:      { return "Left";      }
            case VTrim::Method::Right:     { return "Right";     }
            case VTrim::Method::LeftRight: { return "LeftRight"; }
        }
    };

    const auto summary = "-------VarTrim Summary Statistics\n\n"
                         "-------VarTrim Inputs\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Input alignment file: %2%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %3% regions\n"
                         "       Method:  %4%\n\n"
                         "-------Trimming\n\n"
                         "       Left:  %5% reads\n"
                         "       Right: %6% reads\n\n"
                         "-------Before trimming\n\n"
                         "       Number of alignments: %7%\n\n"
                         "-------After trimming\n\n"
                         "       Number of alignments: %8%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()     // 1
                                            % src         // 2
                                            % stats.nRegs  // 3
                                            % meth2Str()   // 4
                                            % stats.left   // 5
                                            % stats.right  // 6
                                            % stats.before // 7
                                            % stats.after  // 8
                     ).str());
    o.writer->close();
}

void VTrim::report(const FileName &file, const Options &o)
{
    /*
     * Generating VarTrim_summary.stats
     */
    
    writeSummary("VarTrim_summary.stats", file, analyze(file, o), o);
}
