#ifndef GI_TYPES_HPP
#define GI_TYPES_HPP

#include <string>

namespace Anaquin
{
    typedef unsigned long Counts;

    typedef std::string Sequence;
    
    // Defined as long long because there could be many reads
    typedef long long Reads;
    
    typedef double FPKM;
    typedef double Coverage;
    typedef double Concentration;

    // Number of lines in a file (most likely a large file)
    typedef long long Lines;
    
    typedef double Fold;
    typedef double Percentage;

    typedef std::string BaseID;
    typedef std::string GeneID;
    typedef std::string ChromoID;
    typedef std::string SequinID;
    typedef std::string ContigID;
    typedef std::string IsoformID;
    typedef std::string TranscriptID;

    typedef std::string Units;
    typedef std::string FileName;
    typedef std::string FilePath;
    typedef std::string AxisLabel;
    typedef std::string FeatureName;

    typedef long long Base;
}

#endif