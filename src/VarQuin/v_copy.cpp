#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

VCopy::Stats VCopy::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
 
    VCopy::Stats stats;
    
    /*
     * Filter all sequins with CNV equal to the option
     */
    
    for (auto &i : r.data())
    {
        auto m = r.match(i.first);
        
        if (m->concent() == o.copy)
        {
            
            
        }
    }
    
    /*
     * Use the sequins from last section for calculating normalization factor
     */
 
    
    /*
     * Adjust all other sequins relative to the calculated sequins
     */
    
    
    /*
     * Perform normalization for all sequins together
     */
    
    
    return stats;
}
