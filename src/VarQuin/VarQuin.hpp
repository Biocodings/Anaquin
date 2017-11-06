#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/data.hpp"
#include "tools/tools.hpp"
#include "data/biology.hpp"

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

    template <typename T> T noRV(const T &x)
    {
        return (isSubstr(x, "_R") || isSubstr(x, "_V")) ? noLast(x, "_") : x;
    }
    
    template <typename T> bool isCancer(const T &x)
    {
        return x.find("CS_") != std::string::npos;
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
