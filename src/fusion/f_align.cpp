#include "fusion/f_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

FAlign::Stats FAlign::report(const FileName &file, const Options &o)
{
    FAlign::Stats stats;
//    const auto &r = Standard::instance().r_fus;
//
//    o.analyze(file);
//
//    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
//    {
//        if (!align.i && !(info.p.i % 1000000))
//        {
//            o.wait(std::to_string(info.p.i));
//        }
//        
//        stats.update(align);
//        
//        const SequinData *data = nullptr;
//        
//        if ((data = r.findNormal(align.l)) || (data = r.findFusion(align.l)))
//        {
//            if ((data = r.findFusion(align.l)))
//            {
//                std::cout << "Fusion" << std::endl;
//            }
//            else
//            {
//                std::cout << "Normal" << std::endl;
//            }
//        }
//    });
//
    return stats;
}