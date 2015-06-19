#ifndef GI_SEQUIN_HPP
#define GI_SEQUIN_HPP

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

            // The group of log-fold changes between samples
            Group grp;

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
        inline Concentration abund() const
        {
            return r->abund() + (v ? v->abund() : 0);
        }

        // Reference mixture, always valid
        const Sequin *r = NULL;

        // Variant mixture, could be NULL
        const Sequin *v = NULL;

        std::map<TypeID, Sequin> sequins;
    };
}

#endif