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

            SequinID id;
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

    struct GeneSequin
    {
        inline Concentration abund() const
        {
            return r.abund() + v.abund();
        }

        // Reference and variant mixtures
        Sequin r, v;
    };
}

#endif