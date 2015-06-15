#ifndef GI_SEQUIN_HPP
#define GI_SEQUIN_HPP

#include "data/types.hpp"
#include "data/locus.hpp"

namespace Spike
{
    class Sequin
    {
        public:
            operator Locus()     const { return l;  }
            operator IsoformID() const { return id; }

            SequinID id;
            Locus l;

            BasePair length;

            // Amount of abundance, a non-const method
            Concentration &abund() { return _abund; }
        
            // Amount of abundance, a const method
            const Concentration &abund() const { return _abund; }

        private:
            Concentration _abund;
    };

    struct Sequins
    {
        inline Concentration abund() const
        {
            return r.abund() + v.abund();
        }
        
        // Each mixture represents a transcript for a gene
        GeneID geneID;
        
        // Reference and variant mixtures
        Sequin r, v;
    };
}

#endif