#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "reader.hpp"
#include "tokens.hpp"
#include "standard.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_gtf.hpp"

extern std::string d_tab_fa();
extern std::string silico_fa();
extern std::string d_ref_bed();
extern std::string standards_bed();
extern std::string standards_gtf();
extern std::string d_variant_vcf();
extern std::string r_standards_txt();

using namespace Spike;

template <typename Iter> BasePair countLocus(const Iter &iter)
{
    BasePair n = 0;
    
    for (const auto &i : iter)
    {
        n += static_cast<Locus>(i).length();
    }

    return n;
}

Standard::Standard()
{
    std::stringstream in(silico_fa());
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

void Standard::meta(const std::string &mix)
{
    // Empty Implementation
}

void Standard::dna(const std::string &mix)
{
    ParserVCF::parse(d_variant_vcf(), [&](const VCFVariant &v, const ParserProgress &)
    {
        d_vars[v.l] = v;
    }, ParserMode::String);

    ParserBED::parse(d_ref_bed(), [&](const BedFeature &f, const ParserProgress &)
    {
        d_exons.push_back(f);
    }, ParserMode::String);

    //ParserCSV::parse(mix, [&](const Fields &fields, const ParserProgress &)
    //{
        // Empty Implementation
    //}, ParserMode::String);

    ParserFA::parse(d_tab_fa(), [&](const FALine &l, const ParserProgress &)
    {
        d_seqs.insert(l.id);
    }, ParserMode::String);

    assert(!d_seqs.empty() && !d_vars.empty() && !d_exons.empty());
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

    ParserGTF::parse(standards_gtf(), [&](const Feature &f, const ParserProgress &)
	{
		l.end = std::max(l.end, f.l.end);
		l.start = std::min(l.start, f.l.start);
        
        r_fs.push_back(f);
        
        assert(!f.tID.empty());
        assert(!f.geneID.empty());

        iids.insert(f.tID);
        gids.insert(f.geneID);
        
        // Construct a mapping between isoformID to geneID
        r_iso2Gene[f.tID] = f.geneID;
    }, ParserGTF::String);

    assert(!r_iso2Gene.empty());
    std::vector<std::string> toks;

    /*
     * Construct the data-structure for genes
     */

    for (auto gid : gids)
    {
        Feature g;

        g.id = g.geneID = gid;
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
                    r_exons.push_back(f);
                    r_l_exons.push_back(RNALocus(f.geneID, f.l));
                }
            }
        }

        assert(g.l.end > g.l.start);
        r_genes.push_back(g);
    }
    
    CHECK_AND_SORT(r_exons);
    CHECK_AND_SORT(r_genes);

    /*
     * An identical GTF file has already been parsed. For simplistic, we prefer to directly
     * extract the locus from a BED file.
     */

    // Required while reading mixtures
    std::map<IsoformID, Locus> temp;

    ParserBED::parse(standards_bed(), [&](const BedFeature &t, const ParserProgress &p)
    {
        assert(!t.name.empty());
        
        // How easy this is, we'd have to perform unions from a GTF file
        temp[t.name] = t.l;

        /*
         * In this context, we're given a transcript. Each block in the transcript is an exon.
         */

        Feature j;

        j.id   = id;
        j.tID  = t.name;
        j.type = Intron;

        for (auto i = 0; i < t.blocks.size(); i++)
        {
            if (i)
            {
                j.l = Locus(t.blocks[i - 1].end + 1, t.blocks[i].start - 1);
                r_introns.push_back(j);
            }
        }

        const auto iter = std::find_if(r_genes.begin(), r_genes.end(), [&](const Feature &g)
        {
            return (g.id == r_iso2Gene[t.name]);
        });

        assert(iter != r_genes.end());
    }, ParserMode::String);

    CHECK_AND_SORT(r_introns);

    static std::map<std::string, Group> gs =
    {
        { "A", A }, { "B", B }, { "C", C }, { "D", D }
    };

    // This definition is the easiest because variant sequin might be unavailable
    std::map<GeneID, std::vector<Fields>> seqs;

    enum RNAField
    {
        ID,
        RNA_ID,
        Group,
        Con_A,
        Con_B,
        Fold,
        LogFold,
    };
    
    ParserCSV::parse(r_standards_txt(), [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }
        
        seqs[fields[RNA_ID].substr(0, fields[RNA_ID].length() - 2)].push_back(fields);
    }, ParserMode::String);

    assert(!seqs.empty());
    
    for (const auto &p : seqs)
    {
        auto parse_sequin = [&](const Fields &f, SequinMap &m, Mixture mix)
        {
            Sequin s;

            s.id  = f[RNA_ID];
            s.l   = temp.at(s.id);
            s.raw = stof(f[mix == MixA ? Con_A : Con_B]);

            return (m[s.id] = s);
        };

        auto build_sequins = [&](const Sequin &r, const Sequin &v, SequinsMap &m)
        {
            Sequins seqs;

            seqs.r   = r;
            seqs.v   = v;
            seqs.grp = gs.at(p.second[0][Group]);
            seqs.geneID = r.id.substr(0, r.id.length() - 2);
            seqs.fold = 0.0; // TODO: This might not even needed

            m[seqs.geneID] = seqs;
        };

        const auto &r_fields = p.second[0];
        
        const auto seq_r_a = parse_sequin(r_fields, r_seqs_iA, MixA);
        const auto seq_r_b = parse_sequin(r_fields, r_seqs_iB, MixB);
        assert(r_fields.size() == 7);
        
        if (p.second.size() > 1)
        {
            const auto &v_fields = p.second[1];
            assert(v_fields.size() == 7);

            const auto seq_v_a = parse_sequin(v_fields, r_seqs_iA, MixA);
            const auto seq_v_b = parse_sequin(v_fields, r_seqs_iB, MixB);

            build_sequins(seq_r_a, seq_r_a, r_seqs_gA);
            build_sequins(seq_r_b, seq_v_b, r_seqs_gB);
        }
        else
        {
            build_sequins(seq_r_a, seq_r_b, r_seqs_gA);
            build_sequins(seq_r_b, seq_r_b, r_seqs_gB);
        }
    }
    
    assert(r_l_exons.size() == r_exons.size());
    assert(!r_exons.empty() && !r_introns.empty());

    /*
     * Merging overlapping regions for the exons
     */

    r_c_exons = countLocus(r_l_exons = Locus::merge<RNALocus, RNALocus>(r_l_exons));
    
    assert(r_c_exons);
    assert(!r_l_exons.empty());
    assert(!Locus::overlap(r_l_exons));
    
    for (const auto &i: r_seqs_iA)
    {
        r_sequins.push_back(i.second);
    }
}