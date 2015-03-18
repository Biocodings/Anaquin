#ifndef AS_STANDARD_HPP
#define AS_STANDARD_HPP

#include <set>
#include <map>
#include <list>
#include <vector>
#include "Locus.hpp"
#include "Feature.hpp"
#include "Sequence.hpp"
#include "ConfusionMatrix.hpp"

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

struct Concentration
{
    GeneID id;
    
    // Amount of concentration be added
    Amount amounts;
};

struct MixturePair
{
    Concentration x, y;
};

struct Standard
{
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

    ChromoID id;

    // The location of the chromosome
    Locus l;

    // Genes mixed in the standard sorted by positions
    std::vector<Gene> genes;

    // Sequins added for mixture A
    std::list<MixturePair> mixA;

    // Known features (exons, introns etc)
    std::vector<Feature> fs;

    // Known introns (also known as junctions)
    std::vector<Feature> introns;
};

#endif