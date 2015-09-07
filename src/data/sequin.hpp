#ifndef GI_SEQUIN_HPP
#define GI_SEQUIN_HPP

#include <numeric>
#include "data/locus.hpp"

namespace Anaquin
{
    class Sequin
    {
        public:
            operator Locus()    const { return l;  }
            operator SequinID() const { return id; }

            // Eg: R2_7_1, R2_7_2 etc
            SequinID id;

            // Eg: R2_7, R2_72 etc
            BaseID baseID;
        
            // Eg: 1, 2, 3 etc
            TypeID typeID;

            Locus l;

            Base length;

            // Abundance spiked, a non-const method
            Concentration &abund() { return _abund; }
        
            // Abundance spiked, a const method
            const Concentration &abund() const { return _abund; }

        private:
            Concentration _abund;
    };
}

#endif