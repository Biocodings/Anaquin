#include "fusion/f_fusion.hpp"
#include "parsers/parser_fusion.hpp"

using namespace Spike;

FFusion::Stats FFusion::analyze(const std::string &file, const Options &options)
{
    FFusion::Stats stats;
    const auto &s = Standard::instance();
    
    ParserFusion::parse(Reader(file), [&](const ParserFusion::Fusion &f, const ParserProgress &)
    {
        if (f.chr_1 == s.id && f.chr_2 == s.id)
        {
            
//            if (classify(stats.p.m, f, [&](const ParserFusion::Fusion &)
  //          {
                
                
                
                
    //            return true;
      //      }));
        }
    });
    
    /*
     if (align.id == s.id)
     {
     stats.n_chrT++;
     
     if (classify(stats.p.m, align, [&](const Alignment &)
     {
     matched = findMap(s.d_vars, align, MatchRule::Contains);
     
     if (options.filters.count(matched->id))
     {
     return Ignore;
     }
     
     return matched ? Positive : Negative;
     }))
     {
     stats.c.at(matched->id)++;
     }
     }
     else
     {
     stats.n_samps++;
     }
*/
    
    
    
    
    return stats;
}