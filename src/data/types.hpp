#ifndef GI_TYPES_HPP
#define GI_TYPES_HPP

#include <string>

namespace Anaquin
{
    typedef std::string ContigID;
    typedef std::string MetaQuinID;

    typedef unsigned long Counts;
    typedef double Percentage;
    
    typedef std::string Sequence;
    
    // Defined as long long because there could be many reads
    typedef long long Reads;
    
    typedef double FPKM;
    typedef double Coverage;
    typedef double Concentration;

    // Number of lines in a file (most likely a large file)
    typedef long long Lines;
    
    typedef double Fold;

    typedef std::string TypeID;
    typedef std::string BaseID;
    typedef std::string GeneID;
    typedef std::string ChromoID;
    typedef std::string SequinID;
    typedef std::string VariantID;
    typedef std::string TranscriptID;

    typedef std::string FileName;
    typedef std::string FeatureName;

    // Eg: 388488 from the first matching base
    typedef long long BasePair;    
}

#endif