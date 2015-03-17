#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "ParserFA.hpp"
#include "ParserBED.hpp"
#include "ParserGTF.hpp"
#include "StandardFactory.hpp"

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
    
	r.l.end = std::numeric_limits<BasePair>::min();
	r.l.start = std::numeric_limits<BasePair>::max();

    /*
     * Data-structures required to build up the chromosome. Orders of the features are not guarenteed.
     */
    
    std::set<GeneID> gids;

    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f, ParserProgress &p)
	//ParserGTF::parse("C://Sources//QA//Data//Standards//RNAstandards.gtf", [&](const Feature &f)
	{
		r.l.end = std::max(r.l.end, f.l.end);
		r.l.start = std::min(r.l.start, f.l.start);
        r.fs.push_back(f);
        
        assert(f.options.size() == 2);
        assert(f.options.count("gene_id"));
        assert(f.options.count("transcript_id"));

        // TODO: Linking error in Xcode
        Options t = f.options;
        
        const auto gid = t["gene_id"];
        gids.insert(gid);
    });

    /*
     * Construct the data-structure for each gene.
     */

    for (auto gid : gids)
    {
        Gene g;
        
        g.id = gid;
        g.l.end = std::numeric_limits<BasePair>::min();
        g.l.start = std::numeric_limits<BasePair>::max();

        /*
         * Add all features for this gene.
         */
        
        for (auto f : r.fs)
        {
            if (f.options["gene_id"] == g.id)
            {
                g.l.end = std::max(g.l.end, f.l.end);
                g.l.start = std::min(g.l.start, f.l.start);
                
                if (f.type == Exon)
                {
                    g.exons.push_back(f);
                }
            }
        }
        
        assert(g.l.end > g.l.start);
        
        // Sort the exons by starting positions
        std::sort(g.exons.begin(), g.exons.end(), [](const Feature& x, const Feature& y)
        {
            return (x.l.start < y.l.start);
        });
        
        assert(!g.exons.empty());
        r.genes.push_back(g);
    }
    
    assert(!r.genes.empty());
    
    // Sort the genes by starting positions
    std::sort(r.genes.begin(), r.genes.end(), [](const Gene& x, const Gene& y)
    {
        return (x.l.start < y.l.start);
    });

    /*
     * Create data-structure for the known junctions between exons.
     */

    ParserBED::parse("/Users/user1/Sources/QA/Data/Standards/RNAstandards.bed", [&](const BedFeature &f)
    {
        /*
         * In this context, a block is simply an exon. The name of a BED line would be the name of the gene.
         */
        
        const auto iter = std::find_if(r.genes.begin(), r.genes.end(), [&](const Gene &g)
        {
            return (g.id == f.name);
        });

        assert(iter != r.genes.end());
        
        for (std::size_t i = 0; i < f.blocks.size(); i++)
        {
            if (i)
            {
                Feature j;

                j.chromo = r.id;
                j.type = Junction;

                // Junction (intron) is a region between exons that have been spliced
                j.l = Locus(f.blocks[i - 1].end, f.blocks[i].start);

                // TODO: Fix this
                r.introns.push_back(j);
                
                iter->js.push_back(j);
            }
        }

        assert(iter->exons.size() == iter->js.size() + 1);
    });

    /*
     * Create data-structure for the known sequins.
     */
    
	ParserFA::parse("/Users/user1/Sources/QA/Data/Standards/RNAsequins.fa", [&](const Sequence &s)
	//ParserFA::parse("C://Sources//QA//Data//Standards//RNAsequins.fa", [&](const Sequence &s)
	{
		r.sequins[s.id] = s;
	});

    return r;
}