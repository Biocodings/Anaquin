#ifndef AS_STANDARD_HPP
#define AS_STANDARD_HPP

#include <map>
#include <list>
#include "Types.hpp"
#include "Feature.hpp"
#include "Sequence.hpp"

typedef std::map<std::string, Sequence> SequinMap;

struct Standard
{
	inline bool matchChromo(const Feature &q) const
	{
		return (q.start >= start && q.end <= end);
	}

	inline bool matchFeature(const Feature &q) const
	{
		for (auto r : fs)
		{
			if (r.type == q.type && q.start >= r.start && q.end <= r.end)
			{
				return true;
			}
		}

		return false;
	}

    std::string id;
    
	Locus end;
	Locus start;

    // List of sequins added to the experiment
    SequinMap sequins;
    
    // List of known spliced junctions
    std::list<Feature> js;

    // List of known features
    std::list<Feature> fs;
};

#endif