#include "fusion/f_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

FAlign::Stats FAlign::report(const FileName &file, const Options &options)
{
    FAlign::Stats stats;
    //const auto &s = Standard::instance();

    std::vector<Alignment> exons, introns;
    
    options.info("Parsing alignment file");

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        
    });

    return stats;
}