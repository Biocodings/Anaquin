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

extern std::string DNADataBed();
extern std::string DNADataMix();
extern std::string DNADataVCF();

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

typedef std::set<std::string> MixtureFilter;

template <typename SequinMap, typename BaseMap> void parseMix(const Reader &r,
                                                              SequinMap &a,
                                                              SequinMap &b,
                                                              BaseMap   &ba,
                                                              BaseMap   &bb,
                                                              const MixtureFilter &filters)
{
    a.clear();
    b.clear();
    ba.clear();
    bb.clear();
    
    // Used to detect duplicates
    std::set<SequinID> sequinIDs;
    
    // Used to link sequins for each base
    std::map<BaseID, std::set<TypeID>> baseIDs;

    ParserCSV::parse(r, [&](const Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }
        
        Sequin s;
        
        // Make sure there's no duplicate in the mixture file
        assert(sequinIDs.count(fields[0]) == 0);
        
        sequinIDs.insert(s.id = fields[0]);

        // Base ID is simply the ID without the last part
        s.baseID = s.id.substr(0, s.id.find_last_of("_"));

        // Skip over "_"
        s.typeID = s.id.substr(s.id.find_last_of("_") + 1);
        
        assert(s.id != s.baseID && !s.typeID.empty());

        // Length of the sequin
        s.length = stoi(fields[1]);
        
        // Concentration for mixture A
        s.abund() = stof(fields[2]);
        
        // Create an entry for mixture A
        a[s.id] = s;
        
        // Concentration for mixture B
        s.abund() = stof(fields[3]);
        
        // Create an entry for mixture B
        b[s.id] = s;

        baseIDs[s.baseID].insert(s.typeID);
    });

    assert(!a.empty() && !b.empty());

    for (const auto &baseID : baseIDs)
    {
        const auto &typeIDs = baseID.second;
        assert(typeIDs.size() >= 1);
        
        auto f = [&](const SequinMap &m, BaseMap &bm)
        {
            typename BaseMap::mapped_type base;

            for (auto iter = typeIDs.begin(); iter != typeIDs.end(); iter++)
            {
                // Reconstruct the sequinID
                const auto id = baseID.first + "_" + *iter;
                
                base.sequins.insert(std::pair<TypeID, Sequin>(*iter, m.at(id)));
            }

            bm[baseID.first] = base;
        };

        f(a, ba);
        f(b, bb);
    }

    assert(!a.empty() && !b.empty() && !ba.empty() && !bb.empty());
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
    
    // Throw an exception if any inconsistency
    check();
}

void Standard::check() const
{
    for (const auto &i : d_seqs_bA)
    {
        if (i.second.alleleFreq() == 0)
        {
            throw std::runtime_error("Zero allele frequency: " + i.first);
        }
    }
}

void Standard::dna_mod(const Reader &r)
{
    std::vector<std::string> tokens;
    
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        // Eg: D_1_10_R.G/A
        Tokens::split(f.name, ".", tokens);

        // Eg: D_1_10_R and G/A
        assert(tokens.size() == 2);

        // Eg: G/GACTCTCATTC
        const auto var = tokens[1];

        // Eg: D_1_10
        const auto id = tokens[0].substr(0, tokens[0].find_last_of("_"));
        
        Variation v;

        // Eg: D_1_10
        d_seqIDs.insert(v.id = id);
        
        // Eg: G/GACTCTCATTC
        Tokens::split(var, "/", tokens);

        // Eg: G and GACTCTCATTC
        assert(tokens.size() == 2);
        
        v.type = ParserVCF::strToSNP(tokens[0], tokens[1]);
        v.ref  = tokens[0];
        v.alt  = tokens[1];

        d_vars[v.l = f.l] = v;
    });

    assert(!d_vars.empty());
    assert(d_vars.size() >= d_seqIDs.size());
}

void Standard::dna_mix(const Reader &r)
{
    // Parse the mixture file
    parseMix(r, d_seqs_A, d_seqs_B, d_seqs_bA, d_seqs_bB, MixtureFilter());
}

void Standard::meta_mod(const Reader &r)
{
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        m_model.push_back(f);
    });

    assert(!m_model.empty());
}

void Standard::meta_mix(const Reader &r)
{
    m_seqs_A.clear();
    m_seqs_B.clear();
    
    // Parse a mixture file 
    parseMix(r, m_seqs_A, m_seqs_B, m_seqs_bA, m_seqs_bB, MixtureFilter());

    assert(!m_seqs_A.empty()  && !m_seqs_B.empty());
    assert(!m_seqs_bA.empty() && !m_seqs_bB.empty());
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
    dna_mod(Reader(DNADataBed(), DataMode::String));
    dna_mix(Reader(DNADataMix(), DataMode::String));
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
    std::set<GeneID> sequinIDs;

    /*
     * The orders in a GTF file is not guaranteed. For simplicity, we'll defer most of the workloads
     * after parsing.
     */

    ParserGTF::parse(Reader(RNADataGTF(), DataMode::String), [&](const Feature &f, const ParserProgress &)
    {
        // TODO: Please fix me!
        if (f.tID == "R1_140_1" || f.tID == "R1_143_1" || f.tID == "R1_53_2")
        {
            return; // TODO!
        }

        assert(!f.tID.empty() && !f.geneID.empty());
        
        l.end   = std::max(l.end, f.l.end);
        l.start = std::min(l.start, f.l.start);
        
        // Automatically filter out the duplicates
        geneIDs.insert(f.geneID);
        sequinIDs.insert(f.tID);
        
        fs.push_back(f);
        
        // Construct a mapping between sequinID to geneID
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
     * Construct a data-structure that maps from sequinID to it's positon
     */

    for (const auto &sequinID : sequinIDs)
    {
        r_sequins[sequinID] = Locus::expand(r_exons, [&](const Feature &f)
        {
            return f.tID == sequinID;
        });

        // Make sure the position is valid
        assert(r_sequins.at(sequinID) != Locus());
    }

    assert(!r_sequins.empty());

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
}

void Standard::rna_mix(const Reader &r)
{
    // Parse the mixture file
    parseMix(r, r_seqs_A, r_seqs_B, r_seqs_gA, r_seqs_gB, MixtureFilter());

    /*
     * Merging overlapping regions for the exons
     */
    
    r_c_exons = countLocus(r_l_exons = Locus::merge<RNALocus, RNALocus>(r_l_exons));
    
    assert(r_c_exons);
    assert(!r_l_exons.empty());
    assert(!Locus::overlap(r_l_exons));

    /*
     * The mixture file isn't able to tell us the position of a sequin. Here, we'll
     * merge the information with the model. Note that we might have a sequin defined
     * in a mixture file that doesn't not appear in the model.
     */

    for (auto &i: r_seqs_A)
    {
        if (!r_sequins.count(i.first))
        {
            std::cout << "Warning: " << i.first << " defined in mixture but not in the model" << std::endl;
        }
        else
        {
            // Note that we're assigning by reference here
            i.second.l = r_sequins.at(i.first);

            // Make sure it's not an empty range
            assert(i.second.l != Locus());
            
            r_seqIDs.insert(i.first);
        }
    }

    assert(!r_seqIDs.empty());
}

void Standard::rna()
{
    rna_mod(Reader(RNADataGTF(), DataMode::String));
    rna_mix(Reader(RNADataMix(), DataMode::String));
}