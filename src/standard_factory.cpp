#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "file.hpp"
#include "tokens.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

Standard StandardFactory::reference()
{
	std::ifstream in("data/silico.fa");
	std::string line;

    if (!in.good())
    {
        std::cerr << "Error: Failed to load the reference chromosome" << std::endl;
        throw std::runtime_error("Error: Failed to load the reference chromosome");
    }

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

    ParserGTF::parse("data/RNA/standards.gtf", [&](const Feature &f, ParserProgress &p)
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

    assert(!r.iso2Gene.empty());
    std::vector<std::string> ts;

    /*
     * Create data-structure for DNA mutations & variations
     */
    
    ParserBED::parse("data/DNA/ChrT.5.8.Variation.bed", [&](const BedFeature &f)
    {
        ts.clear();

        Variation v;
        v.l = f.l;
        
        /*
         * Example: D_3_3_R_C/A
         */
        
        Tokens::split(f.name, "_/", ts);

        v.r = ts[ts.size() - 2];
        v.m = ts[ts.size() - 1];
        
        r.vars.push_back(v);
    });

    assert(!r.vars.empty());
    
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

    ParserBED::parse("data/RNA/standards.bed", [&](const BedFeature &t)
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
     * Create data-structure for the sequins. Refer to the documentation for more details.
     */
    
    auto read_mixture = [&](const std::string &file, std::map<GeneID, GMixture> &mg, std::map<IsoformID, IMixture> &ig)
    {
        GMixture gm;
        IMixture im;

        int n = 0;
        int n2 = 0;
        
        ParserCSV::parse(file, [&](const Fields &fields)
        {
            if (!n++)
            {
                return;
            }

            /*
             * REF, VAR, Grp, REF length, VAR length, AIM, RATIO, RATIO, CON, CON, CON_READ, CON_READ, PER KB, ROUND, ROUND
             */

            gm.grp = gs[fields[2]];
            gm.geneID = fields[0].substr(0, fields[0].length() - 2);
            im.tranID = fields[0];

            /*
             * Create data-structure for the reference mixture
             */
            
            im.fold  = stoi(fields[6]);
            im.reads = stof(fields[10]);
            
            assert(ig.count(im.tranID) == 0);
            ig[im.tranID] = (gm.r = im);
            
            /*
             * Create data-structure for the variant mixture
             */

            if (!fields[11].empty())
            {
                im.fold   = stoi(fields[7]);
                im.reads  = stof(fields[11]);
                im.tranID = fields[1];
                
                assert(ig.count(im.tranID) == 0);
                ig[im.tranID] = (gm.v = im);
                assert(gm.r.tranID != gm.v.tranID);
            }
            
            mg[gm.geneID]   = gm;
        });
    };
    
    read_mixture("data/RNA/mixture_A.csv", r.mix_gA, r.mix_iA);
    read_mixture("data/RNA/mixture_B.csv", r.mix_gB, r.mix_iB);
    
    return r;
}