#include <fstream>
#include <assert.h>
#include <algorithm>
#include "ParserFA.hpp"
#include "ParserGTF.hpp"
#include "StandardFactory.hpp"
#include <iostream>
using namespace std;

Standard StandardFactory::reference()
{
	std::ifstream in("/Users/user1/Sources/QA/Data/Standards/ChrT.5.10.fa");
	//std::ifstream in("C://Sources//QA//Data//Standards//ChrT.5.10.fa");
	std::string line;

	// Assume that the first line contains only the name of the chromosome
	std::getline(in, line);

    Standard r;
    
	// Remove the '<' prefix
	r.id = line.substr(1, line.size());
    
    /*
     * The region occupied by the chromosome is the smallest area contains all the features.
     */
    
	r.loc.end = std::numeric_limits<BasePair>::min();
	r.loc.start = std::numeric_limits<BasePair>::max();

    /*
     * Data-structures required to build up the chromosome. Orders of the features are not guarenteed.
     */
    
    std::set<GeneID> g_ids;
    std::map<TranscriptID, GeneID> t_ids;

    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f)
	//ParserGTF::parse("C://Sources//QA//Data//Standards//RNAstandards.gtf", [&](const Feature &f)
	{
		r.loc.end = std::max(r.loc.end, f.loc.end);
		r.loc.start = std::min(r.loc.start, f.loc.start);
        r.fs.push_back(f);
        
        assert(f.options.size() == 2);
        assert(f.options.count("gene_id"));
        assert(f.options.count("transcript_id"));

        // TODO: Linking error in Xcode
        Options t = f.options;
        
        const auto gid = t["gene_id"];
        
        g_ids.insert(gid);
        t_ids[gid] = gid;
    });

    /*
     * Construct the data-structure for each gene.
     */

    for (auto g_id : g_ids)
    {
        Gene g;
        
        g.id = g_id;
        g.loc.end = std::numeric_limits<BasePair>::min();
        g.loc.start = std::numeric_limits<BasePair>::max();

        /*
         * Add all features for this gene.
         */
        
        for (auto f : r.fs)
        {
            if (f.options["gene_id"] == g.id)
            {
                g.loc.end = std::max(g.loc.end, f.loc.end);
                g.loc.start = std::min(g.loc.start, f.loc.start);
                
                if (f.type == Exon)
                {
                    g.exons.push_back(f);
                }
            }
        }
        
        assert(g.loc.end > g.loc.start);
        
        // Sort the exons by starting positions
        std::sort(g.exons.begin(), g.exons.end(), [](const Feature& x, const Feature& y)
        {
            return (x.loc.start < y.loc.start);
        });
        
        assert(!g.exons.empty());
        r.genes.push_back(g);
    }
    
    assert(!r.genes.empty());
    
    // Sort the genes by starting positions
    std::sort(r.genes.begin(), r.genes.end(), [](const Gene& x, const Gene& y)
    {
        return (x.loc.start < y.loc.start);
    });
    
    /*
     * Create data-structure for each junction between exons. Only possible once the exons
     * have created and sorted.
     */
/*
    std::for_each(r.genes.begin(), r.genes.end(), [&](Gene &g)
    {
        for (auto i = 0; i < g.exons.size(); i++)
        {
            if (i)
            {
                const auto p = g.exons[i - 1];
                const auto e = g.exons[i];
                
                // Assume the exons have been sorted
                assert(p.loc.end < e.loc.start);
                
                Feature j;
                
                j.update(p.end, p.)
                
                j.start = p.end;
                j.end = e.start;
                j.type = Junction;
                j.length = j.end - j.start;
                
                assert(j.end > j.start);
                g.js.push_back(j);
            }
        }

        assert(g.js.size() == g.exons.size() - 1);
    });
*/
	ParserFA::parse("/Users/user1/Sources/QA/Data/Standards/RNAsequins.fa", [&](const Sequence &s)
	//ParserFA::parse("C://Sources//QA//Data//Standards//RNAsequins.fa", [&](const Sequence &s)
	{
		r.sequins[s.id] = s;
	});

    return r;
}