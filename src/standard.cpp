#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "file.hpp"
#include "tokens.hpp"
#include "standard.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

Standard::Standard()
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

	// Remove the '<' prefix
	id = line.substr(1, line.size());
    
    /*
     * The region occupied by the chromosome is the smallest area contains all the features.
     */
    
	l.end = std::numeric_limits<BasePair>::min();
	l.start = std::numeric_limits<BasePair>::max();

    /*
     * Data-structures required to build up the chromosome. Orders of the features are not guarenteed.
     */

    std::set<GeneID> gids;
    std::set<TranscriptID> iids;

    ParserGTF::parse("data/RNA/standards.gtf", [&](const Feature &f, ParserProgress &p)
	{
		l.end = std::max(l.end, f.l.end);
		l.start = std::min(l.start, f.l.start);
        
        fs.push_back(f);
        
        assert(!f.iID.empty());
        assert(!f.geneID.empty());

        iids.insert(f.iID);
        gids.insert(f.geneID);
        
        // Construct a mapping between isoformID to geneID
        iso2Gene[f.iID] = f.geneID;
    });

    assert(!iso2Gene.empty());
    std::vector<std::string> ts;

//    /*
//     * Create data-structure for DNA mutations & variations
//     */
//
//    ParserBED::parse("data/DNA/ChrT.5.8.Variation.bed", [&](const BedFeature &f)
//    {
//        ts.clear();
//
//        Variation v;
//        v.l = f.l;
//        
//        /*
//         * Example: D_3_3_R_C/A
//         */
//        
//        Tokens::split(f.name, "_/", ts);
//
//        v.r = ts[ts.size() - 2];
//        v.m = ts[ts.size() - 1];
//        
//        r.vars.push_back(v);
//    });
//
//    assert(!r.vars.empty());
    
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
        
        for (auto f : fs)
        {
            if (f.geneID == g.id)
            {
                g.l.end = std::max(g.l.end, f.l.end);
                g.l.start = std::min(g.l.start, f.l.start);
                
                if (f.type == Exon)
                {
                    g.exons.push_back(f);
                    exons.push_back(f);
                }
            }
        }

        assert(g.l.end > g.l.start);

        assert(!exons.empty());
        genes.push_back(g);
    }
    
    assert(!genes.empty());

    // Sort the exons by starting position
    std::sort(exons.begin(), exons.end(), [](const Feature& x, const Feature& y)
    {
        return (x.l.start < y.l.start);
    });

    // Sort the genes by starting position
    std::sort(genes.begin(), genes.end(), [](const Gene& x, const Gene& y)
    {
        return (x.l.start < y.l.start);
    });

    /*
     * Create data-structure for RNA standards. An identical GTF file has already been
     * parsed. For simplistic, we prefer to directly extract the locus from a BED file.
     */

    // Required while reading mixtures
    std::map<IsoformID, Locus> temp;

    ParserBED::parse("data/RNA/standards.bed", [&](const BedFeature &t)
    {
        // How easy this is, we'd have to perform unions from a GTF file
        temp[t.name] = t.l;
        
        /*
         * In this context, we're given a transcript. Each block in the transcript is an exon.
         */

        Feature j;

        for (auto i = 0; i < t.blocks.size(); i++)
        {
            if (i)
            {
                j.id = id;
                j.type = Junction;

                // Intron is a region between exons that have been spliced
                j.l = Locus(t.blocks[i - 1].end - 1, t.blocks[i].start); // ????

                introns.push_back(j);
            }
        }

        const auto iter = std::find_if(genes.begin(), genes.end(), [&](const Gene &g)
        {
            return (g.id == iso2Gene[t.name]);
        });

        assert(iter != genes.end());
    });

    assert(!introns.empty());

//    for (auto i = 0; i < r.introns.size(); i++)
  //  {
    //    std::cout << r.introns[i].l.start << " " << r.introns[i].l.end << std::endl;
    //}
    
    static std::map<std::string, Group> gs =
    {
        { "A", A }, { "B", B }, { "C", C }, { "D", D }
    };
    
    /*
     * Create data-structure for sequins. Refer to the documentation for more details.
     */

    auto create_sequin = [&](const IsoformID &id, Group grp, Fold fold, Concentration raw, Concentration fpkm, std::map<IsoformID, Sequin> &ig)
    {
        Sequin seq;
        
        seq.id = id;
        
        // Make sure the sequin is known
        assert(temp.count(seq.id));
        
        // The BED file has given out the position for the sequin
        seq.l = temp[seq.id];
        
        seq.raw  = raw;
        seq.fpkm = fpkm;
        seq.fold = fold;
        
        assert(ig.count(seq.id) == 0);
        assert(seq.l.start != 0 || seq.l.end != 0);
        
        return (ig[seq.id] = seq);
    };
    
    auto read_mixture = [&](const std::string &file, std::map<GeneID, Sequins> &mg, std::map<TranscriptID, Sequin> &ig)
    {
        Sequins seqs;

        unsigned n = 0;

        ParserCSV::parse(file, [&](const Fields &fields)
        {
            if (!n++)
            {
                return;
            }

            /*
             * REF, VAR, Grp, REF_lEN, VAR_LEN, AIM, RATIO, RATIO, CON, CON, CON_READ, CON_READ, PER_KB
             */

            seqs.geneID = fields[0].substr(0, fields[0].length() - 2);
            seqs.r = create_sequin(fields[0], gs[fields[2]], stoi(fields[6]), stof(fields[10]), stof(fields[12]), ig);
            
            /*
             * Create data-structure for the variant mixture
             */

            if (!fields[11].empty())
            {
                seqs.v = create_sequin(fields[1], gs[fields[2]], stoi(fields[7]), stof(fields[11]), stof(fields[13]), ig);
            }

            mg[seqs.geneID] = seqs;
        });
    };

    read_mixture("data/RNA/mixture_A.csv", seqs_gA, seqs_iA);
    read_mixture("data/RNA/mixture_B.csv", seqs_gB, seqs_iB);
}