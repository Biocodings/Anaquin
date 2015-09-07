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
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_feature.hpp"

extern std::string TransStandGTF();

using namespace Anaquin;

struct ParseSequinInfo
{
    // Used to detect duplicates
    std::set<SequinID> seqIDs;
    
    // Used to link sequins for each base
    std::map<BaseID, std::set<TypeID>> baseIDs;
};

template <typename Reference> void readMixture(const Reader &r, Reference &ref, Mixture m, unsigned column=2)
{
    try
    {
        ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
        {
            // Don't bother if this is the first line or an invalid line
            if (p.i == 0 || fields.size() <= 1)
            {
                return;
            }

            ref.add(fields[0], stoi(fields[1]), stof(fields[column]), m);
        });
    }
    catch (...)
    {
        std::cerr << "[Warn]: Error in the mixture file" << std::endl;
    }

    if (!ref.countMixes())
    {
        throw std::runtime_error("Failed to read any sequin in the mixture file. A CSV file format is expected. Please check and try again.");
    }
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
            
            // TODO: Should be this be here?
            //Standard::instance():seqIDs.insert(s.id);
            
            info.baseIDs[s.baseID].insert(s.typeID);
        });
    }
    catch (...)
    {
        std::cerr << "[Warn]: Error in the mixture file" << std::endl;
    }
    
    if (m.empty())
    {
        throw std::runtime_error("Failed to read any sequin in the mixture file. A CSV file format is expected. Please check and try again.");
    }

    return info;
}

void Standard::v_std(const Reader &r)
{
//    // TODO: Fix this
//    seqIDs.clear();
//    
//    ParserFeature::parse(r, [&](const Feature &f, const ParserProgress &)
//    {
//        if (f.type == Exon)
//        {
//            fs_1.push_back(f);
//            
//            /*
//             * TODO: Fix this!!!! Getting sequinIDs should be done somewhere else
//             */
//
//            const auto seqID = f.tID;
//            seqIDs.insert(seqID);
//        }
//    });
//
//    assert(!fs_1.empty());
}

void Standard::v_var(const Reader &r)
{
    std::vector<std::string> toks;

    ParserBED::parse(r, [&](const ParserBED::Annotation &f, const ParserProgress &)
    {
        // Eg: D_1_10_R_G/A
        Tokens::split(f.name, "_", toks);

        // Eg: D_1_10_R and G/A
        assert(toks.size() == 5);

        // Eg: G/GACTCTCATTC
        const auto var = toks[4];

        Variation v;
        
        // Eg: D_1_10
        v.bID = toks[0] + "_" + toks[1] + "_" + toks[2];

        // Eg: D_1_10_R
        v.id  = toks[0] + "_" + toks[1] + "_" + toks[2] + "_" + toks[3];
        
        // Eg: G/GACTCTCATTC
        Tokens::split(var, "/", toks);

        // Eg: G and GACTCTCATTC
        assert(toks.size() == 2);
        
        v.l    = f.l;
        v.alt  = toks[1];
        v.ref  = toks[0];
        v.type = ParserVCF::strToSNP(toks[0], toks[1]);

        r_var.addVar(v);
    });
}

void Standard::v_mix(const Reader &r)
{
    readMixture(r, r_var, Mix_1, 2);
    readMixture(r, r_var, Mix_2, 3);
}

void Standard::m_mix_1(const Reader &r)
{
    readMixture(r, r_meta, Mix_1, 2);
}

void Standard::m_mix_2(const Reader &r)
{
    readMixture(r, r_meta, Mix_2, 3);
}

void Standard::l_mix(const Reader &r)
{
 //TODO   merge(parseMix(r, seqs_1, 2), seqs_1, bases_1);
//TODO    merge(parseMix(Reader(r), seqs_2, 3), seqs_2, bases_2);
}

void Standard::f_mix(const Reader &r)
{
   //TODO parseMix(r, seqs_1, 2);
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

      //TODO  seqIDs.insert(b.id);
    }, "\t");
    
   //TODO assert(!seqIDs.empty() && !f_breaks.empty());
}

void Standard::r_ref(const Reader &r)
{
    ParserGTF::parse(r, [&](const Feature &f, const ParserProgress &)
    {
       if (f.id == id && f.type == Exon)
        {
            r_trans.addRef(f.tID, f.geneID, f.l);
        }
    });
}

void Standard::r_mix(const Reader &r)
{
    readMixture(r, r_trans, Mix_1, 2);
    readMixture(r, r_trans, Mix_2, 3);
}