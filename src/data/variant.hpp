#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <map>
#include <cmath>
#include "data/locus.hpp"
#include "data/biology.hpp"
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
                return Mutation::Deletion;
            }
            else if (alt[0] == '+')
            {
                return Mutation::Insertion;
            }
            else if (ref.size() == alt.size())
            {
                return Mutation::SNP;
            }
            else if (ref.size() > alt.size())
            {
                return Mutation::Deletion;
            }
            else
            {
                return Mutation::Insertion;
            }
        }

        inline Proportion alleleFreq() const
        {
            return static_cast<Proportion>(readV) / (readR + readV);
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
        
        // Optional data
        std::map<std::string, std::string> opts;
    };
}

#endif
