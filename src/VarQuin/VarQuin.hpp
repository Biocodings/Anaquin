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
    
    template <typename T> bool isSomatic(const T &x)
    {
        A_ASSERT(!x.empty());
        return isSubstr(x, "CS_") || isSubstr(x, "CI_");
    }

    template <typename T> bool isStruct(const T &x)
    {
        A_ASSERT(!x.empty());
        return isSubstr(x, "DEL_") ||
               isSubstr(x, "DUP_") ||
               isSubstr(x, "INS_") ||
               isSubstr(x, "INV_") ||
               isSubstr(x, "MEI_");
    }

    template <typename T> bool isGerm(const T &x)
    {
        A_ASSERT(!x.empty());
        return isSubstr(x, "GS_")    ||
               isSubstr(x, "GI_")    ||
               isSubstr(x, "LoGC_")  ||
               isSubstr(x, "VLGC_")  ||
               isSubstr(x, "HiGC_")  ||
               isSubstr(x, "VHGC_")  ||
               isSubstr(x, "SH_")    ||
               isSubstr(x, "LH_")    ||
               isSubstr(x, "SD_")    ||
               isSubstr(x, "LD_")    ||
               isSubstr(x, "ST_")    ||
               isSubstr(x, "LT_")    ||
               isSubstr(x, "SQ_")    ||
               isSubstr(x, "MS_")    ||
               isSubstr(x, "LQ_");
    }

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

    enum class VCFFilter
    {
        NotFiltered,
        Passed,
    };
    
    // Hask key for mapping a variant
    typedef long VarHashKey;
    
    typedef std::map<VarHashKey, Counts> VarHashTable;
}

#endif
