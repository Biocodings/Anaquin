#ifndef AS_STANDARD_HPP
#define AS_STANDARD_HPP

#include <set>
#include <map>
#include <list>
#include <vector>
#include "Locus.hpp"
#include "Feature.hpp"
#include "Sequence.hpp"

typedef std::string GeneID;
typedef std::string TranscriptID;

struct Gene
{
    GeneID id;
    
    // Location of the gene relative to the chromosome
    Locus l;
    
    // List of exons in this gene sorted by positions
    std::vector<Feature> exons;

    // List of known spliced junctions
    std::vector<Feature> js;
};

struct Standard
{
    typedef std::map<std::string, Sequence> SequinMap;
    
	//inline bool matchChromo(const Feature &q) const
	//{
    //    return loc.contains(q.loc);
	//}

	inline bool matchFeature(const Feature &q) const
	{
		for (auto r : fs)
		{
			if (r.type == q.type && r.l.contains(q.l))
			{
				return true;
			}
		}

		return false;
	}

    std::string id;

    // The location of the standard
    Locus l;

    // List of genes mixed in the standard sorted by positions
    std::vector<Gene> genes;
    
    // List of sequins added to the experiment
    SequinMap sequins;
    
    // List of known features
    std::list<Feature> fs;
};

#endif