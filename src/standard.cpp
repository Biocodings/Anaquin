#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "standard.hpp"
#include "data/reader.hpp"
#include "data/tokens.hpp"
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

extern std::string m_bed_f();
extern std::string m_mix_f();
extern std::string m_tab_f();

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
    //seqs.grp = Group::A; // TODO: Fix this....
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

static void parseMix__(const std::string &file, Standard::SequinMap &a, Standard::SequinMap &b)
{
    std::map<SequinID, Fields> seqs;
    
    ParserCSV::parse(file, [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0 || fields.size() != 3)
        {
            return;
        }

        assert(seqs.count(fields[0]) == 0);
        seqs[fields[0]] = fields;
    }, DataMode::String, ",");

    for (const auto &seq : seqs)
    {
        Sequin t;

        t.id  = seq.first;
        t.raw = stof(seq.second[1]);

        a[t.id] = t;

        
        
        t.id  = seq.first;
        t.raw = stof(seq.second[2]);
        
        b[t.id] = t;
    }
}

static void parseMix(const std::string &file, Standard::SequinMap &seq_A, Standard::SequinMap &seq_B,
                                              Standard::PairMap  &pair_A, Standard::PairMap  &pair_B)
{
    std::map<SequinID, std::vector<Fields>> seqs;
    
    ParserCSV::parse(file, [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }

        seqs[fields[RDM_ID].substr(0, fields[RDM_ID].length() - 2)].push_back(fields);
    }, DataMode::String);

    /*
     * Build sequins from the CSV lines
     */

    for (const auto &p : seqs)
    {
        const auto seq_ra = createSequin(p.second[0], MixA);
        const auto seq_va = createSequin(p.second[1], MixA);
        const auto seq_rb = createSequin(p.second[0], MixB);
        const auto seq_vb = createSequin(p.second[1], MixB);
        
        seq_A[seq_ra.id]  = seq_ra;
        seq_A[seq_va.id]  = seq_va;
        seq_B[seq_rb.id]  = seq_rb;
        seq_B[seq_vb.id]  = seq_vb;
        pair_A[seq_ra.id] = pair_A[seq_va.id] = createSequins(seq_ra, seq_va);
        pair_B[seq_rb.id] = pair_B[seq_vb.id] = createSequins(seq_rb, seq_vb);
    }

    assert(!seq_A.empty()  && !seq_B.empty());
    assert(!pair_A.empty() && !pair_B.empty());
}

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
    // Parse mixtures
    parseMix__(m_mix_f(), m_seq_A, m_seq_B);

    // Parse annotations
    ParserBED::parse(m_bed_f(), [&](const BedFeature &f, const ParserProgress &)
    {
        m_model.push_back(f);
        
        if (m_seq_A.count(f.name))
        {
            m_seq_A.at(f.name).l = f.l;
            m_seq_B.at(f.name).l = f.l;
        }
        else
        {
            //std::cout << f.name << std::endl;
        }
    }, DataMode::String);

    assert(!m_model.empty());
    assert(!m_seq_A.empty()  && !m_seq_A.empty());
}

void Standard::dna()
{
    // Parse variations for reference
    ParserVCF::parse(d_vcf_f(), [&](const VCFVariant &v, const ParserProgress &)
    {
        d_vars[v.l] = v;
    }, DataMode::String);

    // Parse mixtures
    parseMix(d_mix_f(), d_seq_A, d_seq_B, d_pair_A, d_pair_B);

    // Parse annotation
    ParserBED::parse(d_bed_f(), [&](const BedFeature &f, const ParserProgress &)
    {
        d_annot.push_back(f);
        d_seq_A.at(f.name).l    = d_seq_B.at(f.name).l    = f.l;
        d_pair_A.at(f.name).r.l = d_pair_A.at(f.name).v.l = f.l;
        d_pair_B.at(f.name).r.l = d_pair_B.at(f.name).v.l = f.l;
    }, DataMode::String);

    assert(!d_annot.empty()  && !d_vars.empty());
    assert(!d_seq_A.empty()  && !d_seq_A.empty());
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
    }, DataMode::String);

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
    }, DataMode::String);

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
