#include "fusion/f_fraction.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

FFraction::Stats FFraction::stats(const FileName &chim, const FileName &splice, const Options &o)
{
    FFraction::Stats stats;
    const auto &r = Standard::instance().r_fus;

    // Measured abundance for the normal genes
    std::map<SequinID, Counts> normals;
    
    // Measured abundance for the fusion genes
    std::map<SequinID, Counts> fusions;
    
    /*
     * Read the normal junctions
     */
    
    const SequinData *match;

    ParserSTab::parse(Reader(splice), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
    {
        if (c.unique == 878034)
        {
            match = match;
        }
        
        if ((match = r.findSplice(c.l)))
        {
            //std::cout << match->id << std::endl;
        }
    });

    /*
     * Read the chimeric junctions
     */
    
    ParserStarFusion::parse(Reader(chim), [&](const ParserStarFusion::Fusion &f, const ParserProgress &)
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