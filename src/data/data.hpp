#ifndef DATA_HPP
#define DATA_HPP

#include <map>
#include <cmath>
#include <string>

namespace Anaquin
{
    /*
     * Represents a mathced element that can be identified
     */
    
    struct Matched
    {
        virtual std::string name() const = 0;
    };
    
    typedef std::string Sequence;
    
    typedef long long KMers;
    typedef long long Reads;
    
    typedef double FPKM;
    typedef double Concent;
    typedef double Coverage;
    typedef double Measured;
    
    typedef double Fold;
    typedef double LogFold;
    typedef double Proportion;
    typedef long double Probability;
    
    typedef std::string ChrID;
    typedef std::string GeneID;
    typedef std::string ExonID;
    typedef std::string GenoID;
    typedef std::string TransID;
    typedef std::string GenomeID;
    typedef std::string SequinID;
    typedef std::string IntronID;
    typedef std::string RegionID;
    typedef std::string IsoformID;
    typedef std::string GenericID;

    typedef std::string Path;
    typedef std::string Label;
    typedef std::string Token;
    typedef std::string Units;
    typedef std::string Command;
    typedef std::string Scripts;
    typedef std::string ContigID;
    typedef std::string FileName;
    typedef std::string ReadName;
    
    typedef long long Base;
    typedef long long Depth;
    typedef long long Counts;
    typedef long long Quality;
    
    typedef std::map<SequinID, Counts> Hist;

    enum Mixture
    {
        Mix_1,
        Mix_2,
    };

    struct Limit
    {
        SequinID id;
        
        // Expected concentration for the sequin
        Concent abund = NAN;
    };
}

#endif
