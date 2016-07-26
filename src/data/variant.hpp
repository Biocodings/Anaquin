#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <cmath>
#include "data/data.hpp"
#include "data/locus.hpp"
#include <boost/format.hpp>

namespace Anaquin
{
    inline long var2hash(const SequinID &id, Mutation type, const Locus &l)
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

        inline Mutation type() const
        {
            if (alt[0] == '-')
            {
                return Deletion;
            }
            else if (alt[0] == '+')
            {
                return Insertion;
            }
            else if (ref.size() == alt.size())
            {
                return SNP;
            }
            else if (ref.size() > alt.size())
            {
                return Deletion;
            }
            else
            {
                return Insertion;
            }
        }

        inline Proportion alleleFreq() const
        {
            return static_cast<Proportion>(readV) / (readR + readV);
        }
        
        inline long key() const
        {
            return var2hash(id, type(), l);
        }
        
        // Eg: chrIS
        ChrID cID;

        // Eg: D_1_10
        SequinID id;

        // The reference position, with the 1st base having position 1
        Locus l;
        
        Sequence ref, alt;
        
        // Allelle frequency
        Proportion allF = NAN;
        
        // Quality score (eg: QUAL in VCF (62.74))
        double qual = NAN;
        
        // Quality of the reference (eg: provided by VarScan)
        double qualR = NAN;
        
        // Quality of the variant (eg: provided by VarScan)
        double qualV = NAN;

        // P-value (not always provided)
        Probability p = NAN;
        
        // Eg: AD for VCF and REF for VarScan
        Counts readR = NAN;
        
        // Eg: AD for VCF and ALT for VarScan
        Counts readV = NAN;
        
        // Depth coverage (eg: DP for VCF and REF+ALT for VarScan)
        Counts depth = NAN;
    };
    
    typedef Variant CalledVariant;
}

#endif