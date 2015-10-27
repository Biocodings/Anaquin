#include "fusion/f_fraction.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

FFraction::Stats FFraction::stats(const FileName &splice, const FileName &chim, const Options &o)
{
    FFraction::Stats stats;
    const auto &r = Standard::instance().r_fus;

    
    
    
    
    /*
     * Read the normal junctions
     */
    
    ParserSTab::parse(Reader(chim), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
    {
            const auto rr = r.match(c.l, MatchRule::Contains);
            std::cout << c.l.start << " " << c.l.end << " " << c.unique << std::endl;
    });
    
    
    /*
     * Read the splicing junctions (eg: star-fusion.fusion_candidates.txt)
     */
    
    ParserStarFusion::parse(Reader(splice), [&](const ParserStarFusion::Fusion &f, const ParserProgress &)
    {
//        std::cout << f. f.reads << std::endl;
    });

    return stats;
}

void FFraction::report(const FileName &splice, const FileName &chim, const Options &o)
{
    const auto stats = FFraction::stats(splice, chim, o);

    /*
     * Generating summary statistics
     */
    
    
    
    
}