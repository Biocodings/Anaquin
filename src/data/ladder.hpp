#ifndef LADDER_HPP
#define LADDER_HPP

#include <map>
#include "data/data.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    class Translate
    {
        public:
        
            inline void add(const Name &from, const Name &to)
            {
                x[from] = to;
            }
        
            inline Name translate(const Name &from)
            {
                return x.at(from);
            }

            inline std::size_t size() const { return x.size(); }

        private:
            std::map<Name, Name> x;
    };
    
    struct Ladder
    {
        inline void add(const SequinID &id, Mixture m, Concent c)
        {
            seqs.insert(id);
            
            switch (m)
            {
                case Mix_1: { m1[id] = c; break; }
                case Mix_2: { m2[id] = c; break; }
            }
        }
        
        inline Counts count() const { return seqs.size(); }

        inline Concent input(const SequinID &id, Mixture m)
        {
            return m == Mix_1 ? m1.at(id) : m2.at(id);
        }

        inline void remove(const SequinID &id)
        {
            seqs.erase(id);
            m1.erase(id);
            m2.erase(id);
        }

        std::set<SequinID> seqs;
        
        // Ladder for mixture 1
        std::map<SequinID, Concent> m1;
        
        // Ladder for mixture 2
        std::map<SequinID, Concent> m2;
    };
}

#endif
