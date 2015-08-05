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

extern std::string FusionDataMix();
extern std::string FusionDataRef();
extern std::string FusionNormalRef();
extern std::string FusionMutatedRef();

extern std::string LadderDataMix();

extern std::string TransDataGTF();
extern std::string TransDataMix();

extern std::string MetaDataFA();
extern std::string MetaDataBed();
extern std::string MetaDataMix();
extern std::string MetaDataTab();

extern std::string VarDataBed();
extern std::string VarDataMix();

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
    std::set<SequinID> seqIDs;
    
    // Used to link sequins for each base
    std::map<BaseID, std::set<TypeID>> baseIDs;
};

/*
 * This function is intended to merge related sequins. For example, merging transcripts for a gene.
 * It should be called after parseMix(). For example, merge(parseMix(..), ... , ...).
 */

template <typename SequinMap, typename BaseMap> void merge(const ParseSequinInfo &info, const SequinMap &m, BaseMap &b)
{
    // We can't merge something that is empty
    assert(!m.empty());
    
    b.clear();

    for (const auto &i : info.baseIDs)
    {
        // Eg: R1_1
        const auto &baseID = i.first;
        
        // Eg: R and V (R1_1_R and R1_1_V)
        const auto &typeIDs = i.second;

        assert(typeIDs.size() >= 1);
        
        typename BaseMap::mapped_type base;

        for (auto iter = typeIDs.begin(); iter != typeIDs.end(); iter++)
        {
            // Reconstruct the sequinID
            const auto seqID = baseID + "_" + *iter;

            base.sequins.insert(std::pair<TypeID, Sequin>(*iter, m.at(seqID)));
        }

        b[baseID] = base;
    }

    assert(!b.empty());
}

template <typename SequinMap, typename BaseMap> void mergeMix__(const Reader &r,
                                                              const ParseSequinInfo &info,
                                                              const SequinMap &a,
                                                              const SequinMap &b,
                                                              BaseMap &ba,
                                                              BaseMap &bb)
{
    ba.clear();
    bb.clear();

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

template <typename SequinMap> ParseSequinInfo parseMix(const Reader &r, SequinMap &m, unsigned column=2)
{
    m.clear();
    ParseSequinInfo info;

    try
    {
        ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
        {
            // Don't bother if this is the first line or an invalid line
            if (p.i == 0 || fields.size() <= 1)
            {
                return; 
            }

            Sequin s;
            
            // Make sure there's no duplicate in the mixture file
            assert(info.seqIDs.count(fields[0]) == 0);
            
            info.seqIDs.insert(s.id = fields[0]);
            
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
    } catch (...)
    {
        std::cerr << "[Warn]: Error in the mixture file" << std::endl;
    }

    return info;
}

template <typename SequinMap> ParseSequinInfo parseMix__(const Reader &r, SequinMap &a, SequinMap &b)
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
        assert(info.seqIDs.count(fields[0]) == 0);
        
        info.seqIDs.insert(s.id = fields[0]);

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
    v_vars.clear();

    std::vector<std::string> tokens;
    
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        // Eg: D_1_10_R_G/A
        Tokens::split(f.name, "_", tokens);

        // Eg: D_1_10_R and G/A
        assert(tokens.size() == 5);

        // Eg: G/GACTCTCATTC
        const auto var = tokens[4];

        // Eg: D_1_10
        const auto id = tokens[0] + "_" + tokens[1] + "_" + tokens[2];

        // Eg: G/GACTCTCATTC
        Tokens::split(var, "/", tokens);

        // Eg: G and GACTCTCATTC
        assert(tokens.size() == 2);
        
        Variation v;

        v.id   = id;
        v.alt  = tokens[1];
        v.ref  = tokens[0];
        v.type = ParserVCF::strToSNP(tokens[0], tokens[1]);

        v_vars[v.l = f.l] = v;
    });

    assert(!v_vars.empty());
}

void Standard::v_mix(const Reader &r)
{
    mergeMix__(r, parseMix__(r, v_seqs_A, v_seqs_B), v_seqs_A, v_seqs_B, v_seqs_bA, v_seqs_bB);
}

void Standard::m_ref(const Reader &r)
{
    m_model.clear();
    
    ParserBED::parse(r, [&](const BedFeature &f, const ParserProgress &)
    {
        m_model.push_back(f);
    });

    assert(!m_model.empty());
}

void Standard::m_mix(const Reader &r)
{
    merge(parseMix(r, m_seqs_A, 2), m_seqs_A, m_seqs_bA);
    merge(parseMix(Reader(r), m_seqs_B, 3), m_seqs_B, m_seqs_bB);
}

void Standard::l_mix(const Reader &r)
{
    seq2base.clear();
    parseMix__(r, l_seqs_A, l_seqs_B);

    for (const auto &i : l_seqs_A)
    {
        const auto base = i.first.substr(0, i.first.size() - 2);
        seq2base[i.first] = base;
        baseIDs.insert(base);
    }
}

void Standard::f_mix(const Reader &r)
{
    parseMix(r, f_seqs_A, 2);
}

void Standard::f_ref(const Reader &r)
{
    ParserCSV::parse(r, [&](const ParserCSV::Fields &f, const ParserProgress &)
    {
        if (f[0] != "chrT-chrT")
        {
            throw std::runtime_error("Invalid reference file. chrT-chrT is expected.");
        }

        FusionBreak b;
        
        b.id = f[4];;
        b.l1 = stod(f[1]) + 1;
        b.l2 = stod(f[2]) + 1;

        if      (f[3] == "ff") { b.s1 = Strand::Forward;  b.s2 = Strand::Forward;  }
        else if (f[3] == "fr") { b.s1 = Strand::Forward;  b.s2 = Strand::Backward; }
        else if (f[3] == "rf") { b.s1 = Strand::Backward; b.s2 = Strand::Forward;  }
        else if (f[3] == "rr") { b.s1 = Strand::Backward; b.s2 = Strand::Backward; }

        // Add a new entry for the known fusion point
        f_breaks.insert(b);

        seqIDs.insert(b.id);
    });
    
    f_n_seq2locus.clear();
    f_f_seq2locus.clear();
    
    /*
     * Read reference file for normal RNA genes
     */

    ParserBED::parse(Reader(FusionNormalRef(), DataMode::String), [&](const BedFeature &f, const ParserProgress &)
    {
        f_n_seq2locus[f.name] = f.l;
    });

    /*
     * Read reference file for fusion genes
     */
    
    ParserBED::parse(Reader(FusionNormalRef(), DataMode::String), [&](const BedFeature &f, const ParserProgress &)
    {
        f_f_seq2locus[f.name] = f.l;
    });

    assert(!seqIDs.empty() && !f_breaks.empty());
    assert(!f_n_seq2locus.empty() && !f_f_seq2locus.empty());
}

void Standard::ladder()
{
    l_mix(Reader(LadderDataMix(), DataMode::String));
}

void Standard::fusion()
{
    f_ref(Reader(FusionDataRef(), DataMode::String));
    f_mix(Reader(FusionDataMix(), DataMode::String));
}

void Standard::meta()
{
    // Apply the default annotation file
    m_ref(Reader(MetaDataBed(), DataMode::String));
    
    // Apply the default mixture file
    m_mix(Reader(MetaDataMix(), DataMode::String));
}

void Standard::variant()
{
    v_ref(Reader(VarDataBed(), DataMode::String));
    v_mix(Reader(VarDataMix(), DataMode::String));
}

void Standard::r_ref(const Reader &r)
{
    r_exons.clear();
    r_genes.clear();
    r_sequins.clear();
    r_introns.clear();
    r_l_exons.clear();
    r_isoformToGene.clear();
    
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

    ParserGTF::parse(r, [&](const Feature &f, const ParserProgress &)
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
    
    r_l_exons = Locus::merge<RNALocus, RNALocus>(r_l_exons);
    
    assert(!r_l_exons.empty());
}

void Standard::r_mix(const Reader &r)
{
    r_seqIDs.clear();
    
    mergeMix__(r, parseMix__(r, r_seqs_A, r_seqs_B), r_seqs_A, r_seqs_B, r_seqs_gA, r_seqs_gB);
    
    /*
     * Merging overlapping regions for the exons
     */
    
    r_c_exons = countLocus(r_l_exons);

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
            // TODO: std::cout << "Warning: " << i.first << " defined in mixture but not in the model" << std::endl;
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
    r_ref(Reader(TransDataGTF(), DataMode::String));
    r_mix(Reader(TransDataMix(), DataMode::String));
}