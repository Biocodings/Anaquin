#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <cmath>
#include <data/locus.hpp>
#include <data/biology.hpp>
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
        
        // Eg: chrT
        ChrID cID;

        // Eg: D_1_10
        SequinID id;

        // The reference position, with the 1st base having position 1
        Locus l;
        
        Sequence ref, alt;
        
        // Allelle frequency
        Proportion allF = NAN;
        
        // Base quality of the reference (not always provided)
        int qualR = NAN;
        
        // Base quality of the variant (not always provided)
        int qualV = NAN;

        // P-value (not always provided)
        Probability p = NAN;
        
        // Eg: AD for VCF and REF for VarScan
        Counts readR;
        
        // Eg: AD for VCF and ALT for VarScan
        Counts readV;
        
        // Total coverage (eg: DP for VCF and REF+ALT for VarScan)
        Counts cov = NAN;
    };
    
    typedef Variant CalledVariant;
}

#endif