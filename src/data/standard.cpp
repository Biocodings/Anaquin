#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "data/standard.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_gtf.hpp"

enum CSVField
{
    CSV_ID,
    CSV_RDM_ID,
    CSV_Group,
    CSV_Con_A,
    CSV_Con_B,
    CSV_Fold,
    CSV_LogFold,
};

extern std::string RNADataFA();
extern std::string RNADataTab();
extern std::string RNADataBed();
extern std::string RNADataGTF();
extern std::string RNADataMix();

extern std::string MetaDataFA();
extern std::string MetaDataBed();
extern std::string MetaDataMix();
extern std::string MetaDataTab();

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

//static Sequins createSequins(const Sequin &r, const Sequin &v)
//{
//    Sequins seqs;
//
//    seqs.r = r;
//    seqs.v = v;
//    //seqs.grp = Group::A; // TODO: Fix this....
//    seqs.geneID = r.id.substr(0, r.id.length() - 2);
//
//    return seqs;
//}

static void parseMix(const Reader &r, Standard::SequinMap &a, Standard::SequinMap &b)
{
    // Define here to detect duplicates
    std::map<SequinID, Fields> seqs;

    ParserCSV::parse(r, [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0 || fields.size() != 5)
        {
            return;
        }

        // Make sure there's no duplicate in the mixture file
        assert(seqs.count(fields[0]) == 0);

        seqs[fields[0]] = fields;
    }, ",");

    for (const auto &seq : seqs)
    {
        Sequin s;

        s.id  = seq.first;
        s.grp = static_cast<Sequin::Group>(seq.second[4][0] - 'A');

        // Length of the sequin
        s.length = stoi(seq.second[1]);

        // Concentration for mixture A
        s.abund() = stof(seq.second[2]);
        
        // Create an entry for mixture A
        a[s.id] = s;

        // Concentration for mixture B
        s.abund() = stof(seq.second[3]);
        
        // Create an entry for mixture B
        b[s.id] = s;
    }
    
    assert(!a.empty() && !b.empty());
}

Standard::Standard()
{
    std::stringstream in(RNADataFA());
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

void Standard::meta_mod(const Reader &r)
{
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        m_model.push_back(f);
                         
        if (m_seq_A.count(f.name))
        {
            m_seq_A.at(f.name).l = f.l;
            m_seq_B.at(f.name).l = f.l;
        }
    });

    assert(!m_model.empty());
}

void Standard::meta_mix(const Reader &r)
{
    m_seq_A.clear();
    m_seq_B.clear();
    
    // Parse a mixture file 
    parseMix(r, m_seq_A, m_seq_B);

    assert(!m_seq_A.empty() && !m_seq_B.empty());
}

void Standard::meta()
{
    // Apply the default mixture file
    meta_mix(Reader(MetaDataMix(), DataMode::String));

    // Apply the default annotation file
    meta_mod(Reader(MetaDataBed(), DataMode::String));
}

void Standard::dna()
{
//    // Parse variations for reference
//    ParserVCF::parse(Reader(d_vcf_f(), DataMode::String), [&](const VCFVariant &v, const ParserProgress &)
//    {
//        d_vars[v.l] = v;
//    });
//
//    // Parse mixtures
//    parseMix(d_mix_f(), d_seq_A, d_seq_B, d_pair_A, d_pair_B);
//    
//    // Parse annotation
//    ParserBED::parse(Reader(d_bed_f(), DataMode::String), [&](const BedFeature &f, const ParserProgress &)
//    {
//        d_annot.push_back(f);
//        d_seq_A.at(f.name).l    = d_seq_B.at(f.name).l    = f.l;
//        d_pair_A.at(f.name).r.l = d_pair_A.at(f.name).v.l = f.l;
//        d_pair_B.at(f.name).r.l = d_pair_B.at(f.name).v.l = f.l;
//    });
//
//    assert(!d_annot.empty()  && !d_vars.empty());
//    assert(!d_seq_A.empty()  && !d_seq_A.empty());
//    assert(!d_pair_A.empty() && !d_pair_B.empty());
}

void Standard::rna_mod(const Reader &r)
{
    /*
     * The region occupied by the chromosome is the smallest area contains all features.
     */
    
    l.end   = std::numeric_limits<BasePair>::min();
    l.start = std::numeric_limits<BasePair>::max();

    std::vector<Feature> fs;
    std::set<GeneID> geneIDs;

    /*
     * The orders in a GTF file is not guaranteed. For simplicity, we'll defer most of the workloads
     * after the parsing.
     */

    ParserGTF::parse(Reader(RNADataGTF(), DataMode::String), [&](const Feature &f, const ParserProgress &)
    {
        assert(!f.tID.empty() && !f.geneID.empty());
        
        l.end   = std::max(l.end, f.l.end);
        l.start = std::min(l.start, f.l.start);
        
        // Automatically filter out the duplicates
        geneIDs.insert(f.geneID);
        
        fs.push_back(f);
        
        // Construct a mapping between isoformID to geneID
        r_isoformToGene[f.tID] = f.geneID;
        
        if (f.type == Exon)
        {
            r_exons.push_back(f);
        }
    });
    
    assert(!r_exons.empty());
    assert(!r_isoformToGene.empty());
    assert(l.end   != std::numeric_limits<BasePair>::min());
    assert(l.start != std::numeric_limits<BasePair>::min());
    
    /*
     * Construct data-structure for each gene
     */
    
    for (const auto &geneID : geneIDs)
    {
        Feature g;

        // The name of the gene is also the it's ID
        g.id = g.geneID = geneID;

        g.l.end   = std::numeric_limits<BasePair>::min();
        g.l.start = std::numeric_limits<BasePair>::max();
        
        /*
         * Add all exons for this gene
         */

        for (const auto &f : fs)
        {
            if (f.geneID == g.id)
            {
                g.l.end   = std::max(g.l.end, f.l.end);
                g.l.start = std::min(g.l.start, f.l.start);
                
                if (f.type == Exon)
                {
                    // TODO: ????
                    r_l_exons.push_back(RNALocus(f.geneID, f.l));
                }
            }
        }

        assert(g.l.end   != std::numeric_limits<BasePair>::min());
        assert(g.l.start != std::numeric_limits<BasePair>::min());
        assert(g.l.end > g.l.start);

        r_genes.push_back(g);
    }

    assert(!r_l_exons.empty() && !r_genes.empty());

    CHECK_AND_SORT(r_exons);
    CHECK_AND_SORT(r_genes);

    /*
     * Extract introns between each successive pair of exons for each gene.
     * Note that it's only possible once the list has been sorted.
     */
    
    for (std::size_t i = 0; i < r_exons.size(); i++)
    {
        if (i && r_exons[i].geneID == r_exons[i-1].geneID)
        {
            Feature f;
            
            f.id   = id;
            f.tID  = r_exons[i].tID;
            f.type = Intron;

            // The locus is simply whatever between the two successive exons
            f.l = Locus(r_exons[i-1].l.end + 1, r_exons[i].l.start - 1);
            
            r_introns.push_back(f);
        }
    }

    assert(!r_introns.empty());
    CHECK_AND_SORT(r_introns);

    assert(!r_geneToIsoform_r.empty() && !r_geneToIsoform_v.empty());
}

void Standard::rna_mix(const Reader &r)
{
    // Parse the mixture file
    parseMix(r, r_seqs_A, r_seqs_B);

    assert(r_l_exons.size() == r_exons.size());
    
    /*
     * Merging overlapping regions for the exons
     */
    
    r_c_exons = countLocus(r_l_exons = Locus::merge<RNALocus, RNALocus>(r_l_exons));
    
    assert(r_c_exons);
    assert(!r_l_exons.empty());
    assert(!Locus::overlap(r_l_exons));

    for (const auto &i: r_seqs_A)
    {
        r_sequins.push_back(i.second);
    }
    
    assert(!r_seqs_gA.empty() && !r_seqs_gB.empty());
}

void Standard::rna()
{
    rna_mod(Reader(RNADataGTF(), DataMode::String));
    rna_mix(Reader(RNADataMix(), DataMode::String));
}