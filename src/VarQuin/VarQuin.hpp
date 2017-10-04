#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/data.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    inline std::string gt2str(Genotype x)
    {
        switch (x)
        {
            case Genotype::Somatic:     { return "Somatic";      }
            case Genotype::Homozygous:  { return "Homozygous";   }
            case Genotype::Heterzygous: { return "Heterozygous"; }
        }
    }

    inline std::string var2str(Variation x)
    {
        switch (x)
        {
            case Variation::SNP:         { return "SNP";         }
            case Variation::Deletion:    { return "Deletion";    }
            case Variation::Insertion:   { return "Insertion";   }
            case Variation::Duplication: { return "Duplication"; }
            case Variation::Inversion:   { return "Inversion";   }
        }
    }
    
    // Eg: LAD_18
    inline bool isLadQuin(const ChrID &x)
    {
        A_ASSERT(!x.empty());
        return x.find("LAD_") != std::string::npos;
    }
    
    // Eg: chrev1, chrev10 etc...
    inline bool isRevChr(const ChrID &x)
    {
        A_ASSERT(!x.empty());
        return x.find("rev") != std::string::npos;
    }    
}

#endif
