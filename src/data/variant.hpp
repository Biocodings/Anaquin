#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <map>
#include <cmath>
#include "data/locus.hpp"
#include "data/biology.hpp"
#include <boost/format.hpp>

namespace Anaquin
{
    inline long var2hash(const SequinID &id, Variation type, const Locus &l)
    {
        const auto str = (boost::format("%1%_%2%_%3%_%4%") % id
                                                           % type
                                                           % l.start
                                                           % l.end).str();
        return std::hash<std::string>{}(str);
    }

    struct Variant
    {
        operator const Locus &() const { return l; }
        
        inline bool operator<(const Locus &x) const { return l < x; }

        inline bool isSV() const
        {
            return alt == "<DUP>" || alt == "<INV>" || alt == "<INS>" || alt == "<DEL>";
        }
        
        inline Variation type() const
        {
            if (alt == "<DUP>" || alt[0] == '=')
            {
                return Variation::Duplication;
            }
            else if (alt == "<INV>" || alt[0] == '^')
            {
                return Variation::Inversion;
            }
            else if (alt == "<INS>" || alt[0] == '+')
            {
                return Variation::Insertion;
            }
            else if (alt == "<DEL>" || alt[0] == '-')
            {
                return Variation::Deletion;
            }
            else if (ref.size() == alt.size())
            {
                return Variation::SNP;
            }
            else if (ref.size() > alt.size())
            {
                return Variation::Deletion;
            }
            else
            {
                return Variation::Insertion;
            }
        }

        inline long key() const
        {
            return var2hash(name, type(), l);
        }
        
        ChrID cID;

        // Eg: GI_005
        std::string name;

        // The reference position, with the 1st base having position 1
        Locus l;
        
        // Reference and alternative allele
        Sequence ref, alt;
        
        // Allelle frequency
        Proportion allF = NAN;
        
        // Quality score (eg: QUAL in VCF (62.74))
        double qual = NAN;
        
        // P-value (not always provided)
        Probability p = NAN;
        
        // Eg: AD for VCF and REF for VarScan
        Counts readR = NAN;
        
        // Eg: AD for VCF and ALT for VarScan
        Counts readV = NAN;
        
        // Depth coverage (eg: DP for VCF and REF+ALT for VarScan)
        Counts depth = NAN;
        
        // Optional data
        std::map<std::string, std::string> opts;
    };
}

#endif
