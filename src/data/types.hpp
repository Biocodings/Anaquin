#ifndef TYPES_HPP
#define TYPES_HPP

#include <string>

namespace Anaquin
{
    const std::string ChrT = "chrT";
    const std::string ExpT = "expT";
    const std::string Endo = "endo";

    typedef std::string Sequence;
    
    typedef long long KMers;
    typedef long long Reads;
    
    typedef double FPKM;
    typedef double Concent;
    typedef double Coverage;
    typedef double Concentration;

    // Number of lines in a file (most likely a large file)
    typedef long long Lines;
    
    typedef double Fold;
    typedef double LogFold;
    typedef double Percentage;

    typedef std::string BinID;
    typedef std::string GeneID;
    typedef std::string ExonID;
    typedef std::string GenoID;
    typedef std::string GenomeID;
    typedef std::string ChromoID;
    typedef std::string SequinID;
    typedef std::string ContigID;
    typedef std::string IntronID;
    typedef std::string IsoformID;
    typedef std::string GenericID;
    typedef std::string FeatureID;
    typedef std::string SampleName;
    typedef std::string TranscriptID;

    typedef std::string Path;
    typedef std::string Label;
    typedef std::string Units;
    typedef std::string Scripts;
    typedef std::string FileName;
    typedef std::string FilePath;
    typedef std::string AxisLabel;
    typedef std::string FeatureName;

    typedef long long Base;
    typedef long long Counts;
}

#endif