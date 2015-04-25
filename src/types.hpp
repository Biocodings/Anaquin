#ifndef GI_TYPES_HPP
#define GI_TYPES_HPP

#include <string>

namespace Spike
{
    typedef unsigned long Counts;
    typedef double Percentage;
    
    typedef std::string Sequence;
    
    // Defined as long long because there could be many reads
    typedef long long Reads;
    
    typedef float FPKM;
    typedef float Concentration;
    
    // Number of lines in a file (most likely a large file)
    typedef long long Lines;
    
    typedef double Fold;
    
    typedef std::string GeneID;
    typedef std::string ChromoID;
    typedef std::string SequinID;
    typedef std::string IsoformID;
    typedef std::string VariantID;
    typedef std::string FeatureName;
    typedef std::string TranscriptID;

    // Eg: 388488 from the first matching base
    typedef long long BasePair;    
}

#endif