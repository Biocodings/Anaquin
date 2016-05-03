#ifndef HIST_HPP
#define HIST_HPP

#include <map>
#include "data/types.hpp"

namespace Anaquin
{
    typedef std::map<long, Counts> HashHist;
    
    typedef std::map<SequinID, Counts> GenomeHist;
    typedef std::map<SequinID, Counts> SequinHist;
}

#endif