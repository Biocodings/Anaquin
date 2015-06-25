#ifndef GI_SEQUIN_HPP
#define GI_SEQUIN_HPP

#include <numeric>
#include "data/types.hpp"
#include "data/locus.hpp"

namespace Spike
{
    class Sequin
    {
        public:
            enum Group { A, B, C, D, E, F, G, H, I, J, K };

            operator Locus()    const { return l;  }
            operator SequinID() const { return id; }

            // Eg: R2_7_1, R2_7_2 etc
            SequinID id;

            // Eg: R2_7, R2_72 etc
            BaseID baseID;
        
            // Eg: 1, 2, 3 etc
            TypeID typeID;

            Locus l;

            BasePair length;

            // Amount of abundance, a non-const method
            Concentration &abund() { return _abund; }
        
            // Amount of abundance, a const method
            const Concentration &abund() const { return _abund; }

        private:
            Concentration _abund;
    };

    struct Base
    {
        Concentration abund() const
        {
            return std::accumulate(sequins.begin(), sequins.end(), 0,
                [&](int sum, const std::pair<TypeID, Sequin> &p)
                {
                    return sum + p.second.abund();
                });
        }

        std::map<TypeID, Sequin> sequins;
    };
    
    struct VariantBase : public Base
    {
        inline double alleleFreq() const
        {
            assert(sequins.size() == 2);
            
            const auto ref = sequins.begin()->first;
            const auto var = sequins.rbegin()->first;

            // Abundance for the reference
            const auto r = sequins.at(ref).abund();
            
            // Abundance for the variant
            const auto v = sequins.at(var).abund();
            
            // Abundance ratio of reference to variant DNA standard
            return v / (r + v);
        }
    };
}

#endif