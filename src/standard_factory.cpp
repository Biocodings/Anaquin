#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "parser_fa.hpp"
#include "parser_bed.hpp"
#include "parser_gtf.hpp"
#include "parser_csv.hpp"
#include "parser_vcf.hpp"
#include "standard_factory.hpp"

const std::string VARIANT_FILE = "/Users/tedwong/Sources/QA/Data/RNA/variant.ChrT51.vcf";

Standard StandardFactory::reference()
{
	std::ifstream in("/Users/tedwong/Sources/QA/Data/RNA/ChrT.5.10.fa");
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
    std::set<IsoformID> iids;

    ParserGTF::parse("/Users/tedwong/Sources/QA/Data/RNA/RNAstandards.gtf", [&](const Feature &f, ParserProgress &p)
	{
		r.l.end = std::max(r.l.end, f.l.end);
		r.l.start = std::min(r.l.start, f.l.start);
        r.fs.push_back(f);
        
        assert(!f.iID.empty());
        assert(!f.geneID.empty());

        iids.insert(f.iID);
        gids.insert(f.geneID);
        
        // Construct a mapping between isoformID to geneID
        r.iso2Gene[f.iID] = f.geneID;
    });

    /*
     * Create data-structure for mutation and variations
     */
    
    ParserVCF::parse(VARIANT_FILE, [&] (const VCFVariant &v)
    {
                         
    });
    
    /*
     * Construct the data-structure for genes
     */

    for (auto gid : gids)
    {
        Gene g;
        
        g.id = gid;
        g.l.end = std::numeric_limits<BasePair>::min();
        g.l.start = std::numeric_limits<BasePair>::max();

        /*
         * Add all features for this gene
         */
        
        for (auto f : r.fs)
        {
            if (f.geneID == g.id)
            {
                g.l.end = std::max(g.l.end, f.l.end);
                g.l.start = std::min(g.l.start, f.l.start);
                
                if (f.type == Exon)
                {
                    g.exons.push_back(f);
                    r.exons.push_back(f);
                }
            }
        }

        assert(g.l.end > g.l.start);

        assert(!r.exons.empty());
        r.genes.push_back(g);
    }
    
    assert(!r.genes.empty());

    // Sort the exons by starting position
    std::sort(r.exons.begin(), r.exons.end(), [](const Feature& x, const Feature& y)
    {
        return (x.l.start < y.l.start);
    });

    // Sort the genes by starting position
    std::sort(r.genes.begin(), r.genes.end(), [](const Gene& x, const Gene& y)
    {
        return (x.l.start < y.l.start);
    });

    /*
     * Create data-structure for known introns
     */

    ParserBED::parse("/Users/tedwong/Sources/QA/Data/RNA/RNAstandards.bed", [&](const BedFeature &t)
    {
        /*
         * In this context, we're given a transcript. Each block in the transcript is an exon.
         */

        Feature j;

        for (auto i = 0; i < t.blocks.size(); i++)
        {
            if (i)
            {
                j.id = r.id;
                j.type = Junction;

                // Intron is a region between exons that have been spliced
                j.l = Locus(t.blocks[i - 1].end, t.blocks[i].start);

                r.introns.push_back(j);
            }
        }

        const auto iter = std::find_if(r.genes.begin(), r.genes.end(), [&](const Gene &g)
        {
            return (g.id == r.iso2Gene[t.name]);
        });
        
        assert(iter != r.genes.end());
    });

    std::map<std::string, Group> gs =
    {
        { "A", A }, { "B", B }, { "C", C }, { "D", D }
    };

    /*
     * Create data-structure for the sequins in each of the mix. Refer to user-manual for more details.
     */
    
    GMixture g;
    
    ParserCSV::parse("/Users/tedwong/Sources/QA/Data/RNA/Standard_A.csv", [&](const Fields &fields)
    {
        /*
         * Create data-structure for isoforms
         */
        
        IMixture i;
        
        i.id   = fields[3];
        i.fold = stoi(fields[4]);
        i.exp  = stod(fields[5]);

        r.isoA[i.id] = i;
        
        /*
         * Create data-structure for genes
         */

        if (fields[0] == "1")
        {
            g.gr     = gs[fields[1]];
            g.id     = fields[2];
            g.r.id   = i.id;
            g.r.exp  = i.exp;
            g.r.fold = i.fold;
            
            assert(g.id != g.r.id);
        }
        else
        {
            assert(g.id == fields[2]);
            assert(g.gr == gs[fields[1]]);
            
            g.v.id   = i.id;
            g.v.exp  = i.exp;
            g.v.fold = i.fold;

            const auto fold  = g.r.fold + g.v.fold;
            const auto total = g.r.exp  + g.v.exp;
            
            assert((g.r.fold / fold * total) == g.r.exp);
            assert((g.v.fold / fold * total) == g.v.exp);
            
            r.mixA[g.id] = g;
        }
    });
    
    return r;
}