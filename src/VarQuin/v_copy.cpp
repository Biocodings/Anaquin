#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

VCopy::Stats VCopy::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
 
    VCopy::Stats stats;
    
    /*
     * Filter all sequins with CNV equal to the input option
     */
    
    for (auto &i : r.data())
    {
        auto m = r.match(i.first);
        
        if (m->concent() == o.copy)
        {
            
            
        }
    }
    
    /*
     * Check calibration statistics for all sequins
     */
 
    
    /*
     * Adjust non-reference sequins relative to the references
     */
    
    
    /*
     * Perform normalization for all sequins
     */
    
    
    return stats;
}
