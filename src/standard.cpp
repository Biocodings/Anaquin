#include <set>
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

    rna();
    dna();
    meta();
}

void Standard::meta()
{
    
}

void Standard::dna()
{
    Variation v;
    std::vector<std::string> toks;
    
    ParserBED::parse("data/dna/ChrT.5.8.Variation.bed", [&](const BedFeature &b)
    {
        v.id  = b.id;
        v.pos = b.l.start;

        // Eg: D_3_3_R_C/A
        Tokens::split(b.name, "_/", toks);
        
        v.r = toks[toks.size() - 2];
        v.m = toks[toks.size() - 1];

        d_vars[v.id] = v;
    });

    assert(!d_vars.empty());
}

void Standard::rna(const std::string &mix)
{
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

    std::map<TranscriptID, GeneID> r_iso2Gene;
    
    ParserGTF::parse("data/RNA/standards.gtf", [&](const Feature &f)
	{
		l.end = std::max(l.end, f.l.end);
		l.start = std::min(l.start, f.l.start);
        
        r_fs.push_back(f);
        
        assert(!f.iID.empty());
        assert(!f.geneID.empty());

        iids.insert(f.iID);
        gids.insert(f.geneID);
        
        // Construct a mapping between isoformID to geneID
        r_iso2Gene[f.iID] = f.geneID;
    });

    assert(!r_iso2Gene.empty());
    std::vector<std::string> toks;

    /*
     * Construct the data-structure for genes
     */

    for (auto gid : gids)
    {
        TGene g;
        
        g.id = gid;
        g.l.end = std::numeric_limits<BasePair>::min();
        g.l.start = std::numeric_limits<BasePair>::max();

        /*
         * Add all features for this gene
         */
        
        for (auto f : r_fs)
        {
            if (f.geneID == g.id)
            {
                g.l.end = std::max(g.l.end, f.l.end);
                g.l.start = std::min(g.l.start, f.l.start);
                
                if (f.type == Exon)
                {
                    g.exons.push_back(f);
                    r_exons.push_back(f);
                }
            }
        }

        assert(g.l.end > g.l.start);

        assert(!r_exons.empty());
        r_genes.push_back(g);
    }
    
    assert(!r_genes.empty());

    // Sort the exons by starting position
    std::sort(r_exons.begin(), r_exons.end(), [](const Feature& x, const Feature& y)
    {
        return (x.l.start < y.l.start);
    });

    // Sort the genes by starting position
    std::sort(r_genes.begin(), r_genes.end(), [](const TGene& x, const TGene& y)
    {
        return (x.l.start < y.l.start);
    });

    /*
     * Create data-structure for RNA standards. An identical GTF file has already been
     * parsed. For simplistic, we prefer to directly extract the locus from a BED file.
     */

    // Required while reading mixtures
    std::map<IsoformID, Locus> temp;

    ParserBED::parse("data/rna/standards.bed", [&](const BedFeature &t)
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

                r_introns.push_back(j);
            }
        }

        const auto iter = std::find_if(r_genes.begin(), r_genes.end(), [&](const TGene &g)
        {
            return (g.id == r_iso2Gene[t.name]);
        });

        assert(iter != r_genes.end());
    });

    assert(!r_introns.empty());

    static std::map<std::string, Group> gs =
    {
        { "A", A }, { "B", B }, { "C", C }, { "D", D }
    };

    ParserCSV::parse(mix, [&](const Fields &fields, unsigned i)
    {
        enum RNAField
        {
            Ref,
            Var,
            Group,
            RLen,
            VLen,
            AimA,
            AimB,
            RatioR,
            RatioV,
        };

        if (i == 0)
        {
            return;
        }

        Sequins seqs;
        
        seqs.r.id   = fields[Ref];
        seqs.v.id   = fields[Var];
        seqs.grp    = gs[fields[Group]];
        seqs.geneID = fields[Ref].substr(0, fields[0].length() - 2);

        seqs.r.l = temp.at(seqs.r.id);
        seqs.v.l = !seqs.v.id.empty() ? temp.at(seqs.v.id) : Locus();

        const auto ratio_r = stoi(fields[RatioR]);
        const auto ratio_v = stoi(fields[RatioV]);
        const auto total_r = ratio_r + ratio_v;

        seqs.fold = static_cast<Fold>(ratio_r) / ratio_v;
        
        const auto r_len = stoi(fields[RLen]);
        const auto v_len = !seqs.v.id.empty() ? stoi(fields[VLen]) : 0;

        auto f = [&](RNAField f, std::map<SequinID, Sequins> &g, std::map<TranscriptID, Sequin> &i)
        {
            assert(f == AimA || f == AimB);

            // Ratio of the reference, reverse between mixtures
            const auto ratio = f == AimA ? ratio_r : ratio_v;

            const auto aim = stof(fields[f]);

            seqs.r.raw  = aim * (static_cast<Fold>(ratio) / total_r);
            seqs.v.raw  = aim - seqs.r.raw;
            seqs.r.fpkm = seqs.r.raw / (1000.0 / r_len);

            // There is always an entry for the reference
            i[seqs.r.id] = seqs.r;
            
            if (!seqs.v.id.empty())
            {
                seqs.v.fpkm  = seqs.v.raw / (1000.0 / v_len);
                i[seqs.v.id] = seqs.v;
            }
            
            g[seqs.geneID] = seqs;
        };

        // Calculate the concentration required by mixture A
        f(AimA, r_seqs_gA, r_seqs_iA);

        // Calculate the concentration required by mixture B
        f(AimB, r_seqs_gB, r_seqs_iB);

        assert(!seqs.geneID.empty() && !seqs.r.id.empty());
    });

    assert(!r_exons.empty() && !r_introns.empty());
}