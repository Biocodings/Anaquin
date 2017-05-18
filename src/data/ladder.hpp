#ifndef LADDER_HPP
#define LADDER_HPP

#include <map>
#include "data/data.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    struct Ladder
    {
        void add(const SequinID &seq, Mixture m, Concent c)
        {
            seqs.insert(seq);
            
            switch (m)
            {
                case Mix_1: { m1[seq] = c; break; }
                case Mix_2: { m2[seq] = c; break; }
            }
        }
        
        inline Counts count() const { return seqs.size(); }
        
        std::set<SequinID> seqs;
        
        // Ladder for mixture 1
        std::map<SequinID, Concent> m1;
        
        // Ladder for mixture 2
        std::map<SequinID, Concent> m2;
    };
}

#endif
