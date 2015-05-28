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

enum CSVField
{
    ID,
    RDM_ID,
    Group,
    Con_A,
    Con_B,
    Fold,
    LogFold,
};

extern std::string silico_f();

extern std::string r_bed_f();
extern std::string r_gtf_f();
extern std::string r_mix_f();

extern std::string d_mix_f();
extern std::string d_vcf_f();
extern std::string d_bed_f();
extern std::string d_tab_f();

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

static Sequins createSequins(const Sequin &r, const Sequin &v)
{
    Sequins seqs;

    seqs.r = r;
    seqs.v = v;
    seqs.grp = Group::A; // TODO: Fix this....
    seqs.geneID = r.id.substr(0, r.id.length() - 2);

    return seqs;
}

static Sequin createSequin(const Fields &f, Mixture mix)
{
    Sequin s;

    s.id  = f[RDM_ID];
    //s.l   = temp.at(s.id);
    s.raw = stof(f[mix == MixA ? Con_A : Con_B]);

    return s;
};

Standard::Standard()
{
    std::stringstream in(silico_f());
    std::string line;
    
    if (!in.good())
    {
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
    // Empty Implementation
}

void Standard::dna()
{
    ParserVCF::parse(d_vcf_f(), [&](const VCFVariant &v, const ParserProgress &)
    {
        d_vars[v.l] = v;
    }, ParserMode::String);

    ParserBED::parse(d_bed_f(), [&](const BedFeature &f, const ParserProgress &)
    {
        d_annot.push_back(f);
    }, ParserMode::String);

    std::map<SequinID, std::vector<Fields>> seqs;

    ParserCSV::parse(d_mix_f(), [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }

        seqs[fields[RDM_ID].substr(0, fields[RDM_ID].length() - 2)].push_back(fields);
    }, ParserMode::String);

    /*
     * Build sequins from the CSV lines
     */

    for (const auto &p : seqs)
    {
        const auto seq_ra = createSequin(p.second[0], MixA);
        const auto seq_va = createSequin(p.second[1], MixA);
        const auto seq_rb = createSequin(p.second[0], MixB);
        const auto seq_vb = createSequin(p.second[1], MixB);

        d_seq_A[seq_ra.id]  = seq_ra;
        d_seq_A[seq_va.id]  = seq_va;
        d_seq_B[seq_rb.id]  = seq_rb;
        d_seq_B[seq_vb.id]  = seq_vb;
        d_pair_A[seq_ra.id] = d_pair_A[seq_va.id] = createSequins(seq_ra, seq_va);
        d_pair_B[seq_rb.id] = d_pair_B[seq_vb.id] = createSequins(seq_rb, seq_vb);
    }

    assert(!d_annot.empty()  && !d_vars.empty());
    assert(!d_seq_A.empty()  && !d_seq_B.empty());
    assert(!d_pair_A.empty() && !d_pair_B.empty());
}

void Standard::rna()
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

    ParserGTF::parse(r_gtf_f(), [&](const Feature &f, const ParserProgress &)
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
    }, String);

    assert(!r_iso2Gene.empty());
    std::vector<std::string> toks;

    /*
     * Construct the data-structure for genes
     */

    for (const auto &gid : gids)
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

    ParserBED::parse(r_bed_f(), [&](const BedFeature &t, const ParserProgress &p)
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

    ParserCSV::parse(r_mix_f(), [&](const Fields &fields, const ParserProgress &p)
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

        auto build_sequins = [&](const Sequin &r, const Sequin &v, PairMap &m)
        {
            Sequins seqs;

            seqs.r      = r;
            seqs.v      = v;
            seqs.grp    = gs.at(p.second[0][Group]);
            seqs.geneID = r.id.substr(0, r.id.length() - 2);

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

            build_sequins(seq_r_a, seq_v_a, r_seqs_gA);
            build_sequins(seq_r_b, seq_v_b, r_seqs_gB);
        }
        else
        {
            build_sequins(seq_r_a, seq_r_b, r_seqs_gA); // TODO: How to fix this??
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