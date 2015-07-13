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

extern std::string FusDataMix();
extern std::string FusDataRef();

extern std::string LadderDataMix();

extern std::string RNADataFA();
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

using namespace Anaquin;

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

struct ParseSequinInfo
{
    // Used to detect duplicates
    std::set<SequinID> sequinIDs;
    
    // Used to link sequins for each base
    std::map<BaseID, std::set<TypeID>> baseIDs;
};

template <typename SequinMap, typename BaseMap> void mergeMix(const Reader &r,
                                                              const ParseSequinInfo &info,
                                                              const SequinMap &a,
                                                              const SequinMap &b,
                                                              BaseMap &ba,
                                                              BaseMap &bb)
{
    for (const auto &baseID : info.baseIDs)
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

template <typename SequinMap> ParseSequinInfo parseMix(const Reader &r, SequinMap &m, unsigned column)
{
    m.clear();
    ParseSequinInfo info;
    
    ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }
        
        Sequin s;
        
        // Make sure there's no duplicate in the mixture file
        assert(info.sequinIDs.count(fields[0]) == 0);
        
        info.sequinIDs.insert(s.id = fields[0]);
        
        // Base ID is simply the ID without the last part
        s.baseID = s.id.substr(0, s.id.find_last_of("_"));
        
        // Skip over "_"
        s.typeID = s.id.substr(s.id.find_last_of("_") + 1);
        
        // Length of the sequin
        s.length = stoi(fields[1]);
        
        assert(s.length);
        
        // Concentration for the mixture
        s.abund() = stof(fields[column]);
        
        // Create an entry for the mixture
        m[s.id] = s;
        
        info.baseIDs[s.baseID].insert(s.typeID);
    }, ",");
    
    return info;
}

template <typename SequinMap> ParseSequinInfo parseMix(const Reader &r, SequinMap &a, SequinMap &b)
{
    a.clear();
    b.clear();
    
    ParseSequinInfo info;
    
    ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
    {
        if (p.i == 0)
        {
            return;
        }
        
        Sequin s;
        
        // Make sure there's no duplicate in the mixture file
        assert(info.sequinIDs.count(fields[0]) == 0);
        
        info.sequinIDs.insert(s.id = fields[0]);

        // Base ID is simply the ID without the last part
        s.baseID = s.id.substr(0, s.id.find_last_of("_"));

        // Skip over "_"
        s.typeID = s.id.substr(s.id.find_last_of("_") + 1);
        
        // Length of the sequin
        s.length = stoi(fields[1]);
        
        // Concentration for mixture A
        s.abund() = stof(fields[2]);
        
        assert(s.length && s.abund());
        
        // Create an entry for mixture A
        a[s.id] = s;
        
        // Concentration for mixture B
        s.abund() = stof(fields[3]);
        
        // Create an entry for mixture B
        b[s.id] = s;

        info.baseIDs[s.baseID].insert(s.typeID);

        // Don't assert s.abund() because there might not be mixture B
        assert(s.length);
    });

    assert(!a.empty() && !b.empty());
    
    return info;
}

Standard::Standard()
{
    id = "chrT";
    
    rna();
    meta();
    cancer();
    fusion();
    ladder();
    variant();
    clinical();
}

void Standard::v_ref(const Reader &r)
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

void Standard::v_mix(const Reader &r)
{
    mergeMix(r, parseMix(r, d_seqs_A, d_seqs_B), d_seqs_A, d_seqs_B, d_seqs_bA, d_seqs_bB);
}

void Standard::m_ref(const Reader &r)
{
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        m_model.push_back(f);
    });

    assert(!m_model.empty());
}

void Standard::m_mix(const Reader &r)
{
    mergeMix(r, parseMix(r, m_seqs_A, m_seqs_B), m_seqs_A, m_seqs_B, m_seqs_bA, m_seqs_bB);
}

void Standard::l_mix(const Reader &r)
{
    parseMix(r, l_seqs_A, l_seqs_B);

    for (const auto &i : l_seqs_A)
    {
        const auto &base = i.first;

        l_map[base + "_" + "A"] = base;
        l_map[base + "_" + "B"] = base;
        l_map[base + "_" + "C"] = base;
        l_map[base + "_" + "D"] = base;
    }
}

void Standard::f_mix(const Reader &r)
{
    parseMix(r, f_seqs_A, 2);
}

void Standard::f_ref(const Reader &r)
{
    f_f_fusions.clear();
    f_r_fusions.clear();

    ParserCSV::parse(r, [&](const ParserCSV::Fields &f, const ParserProgress &)
    {
        assert(f[0] == "chrT-chrT");
        
        const auto l  = Locus(stod(f[1]) + 1, stod(f[2]) + 1);
        const auto id = f[4];

        if (f[3] == "fr")
        {
            f_f_fusions[id] = l;
        }
        else
        {
            f_r_fusions[id] = l;
        }
    });

    assert(!f_f_fusions.empty() && !f_r_fusions.empty());
}

void Standard::ladder()
{
    l_mix(Reader(LadderDataMix(), DataMode::String));
}

void Standard::fusion()
{
    f_mix(Reader(FusDataMix(), DataMode::String));
    f_ref(Reader(FusDataRef(), DataMode::String));
}

void Standard::meta()
{
    // Apply the default mixture file
    m_mix(Reader(MetaDataMix(), DataMode::String));

    // Apply the default annotation file
    m_ref(Reader(MetaDataBed(), DataMode::String));
}

void Standard::variant()
{
    v_ref(Reader(DNADataBed(), DataMode::String));
    v_mix(Reader(DNADataMix(), DataMode::String));
}

void Standard::r_ref(const Reader &r)
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
            return; // TODO! To save assembly from crashing!! Defined in the model but not in sequin!
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

void Standard::r_mix(const Reader &r)
{
    mergeMix(r, parseMix(r, r_seqs_A, r_seqs_B), r_seqs_A, r_seqs_B, r_seqs_gA, r_seqs_gB);
    
    /*
     * Merging overlapping regions for the exons
     */
    
    r_c_exons = countLocus(r_l_exons = Locus::merge<RNALocus, RNALocus>(r_l_exons));
    
    assert(r_c_exons);
    assert(!r_l_exons.empty());
    assert(!Locus::overlap(r_l_exons));

    /*
     * The mixture file is unable to give us the position of a sequin. Here, we'll
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

void Standard::cancer()
{
    // Empty Implementation
}

void Standard::clinical()
{
    // Empty Implementation
}

void Standard::rna()
{
    r_ref(Reader(RNADataGTF(), DataMode::String));
    r_mix(Reader(RNADataMix(), DataMode::String));
}