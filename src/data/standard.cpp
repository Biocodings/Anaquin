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

using namespace Anaquin;

enum MixtureFormat
{
    ID_Length_Mix, // Eg: MG_33  10  60.23529412
    ID_Mix,        // Eg: MG_33  60.23529412
};

static unsigned countColumns(const Reader &r)
{
    std::size_t n = 0;

    ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
    {
        n = std::max(n, fields.size());
    }, ",");

    ParserCSV::parse(Reader(r), [&](const ParserCSV::Fields &fields, const ParserProgress &p)
    {
        n = std::max(n, fields.size());
    }, "\t");
    
    return static_cast<unsigned>(n);
}

template <typename Reference> void readMixture
                (const Reader &r, Reference &ref, Mixture m, MixtureFormat format, unsigned column=2)
{
    auto f = [&](const std::string &delim)
    {
        const auto t = Reader(r);

        try
        {
            bool succceed = false;
            
            ParserCSV::parse(t, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
            {
                // Don't bother if this is the first line or an invalid line
                if (p.i == 0 || fields.size() <= 1)
                {
                    return;
                }
                
                switch (format)
                {
                    case ID_Length_Mix:
                    {
                        succceed = true;
                        ref.add(fields[0], stoi(fields[1]), stof(fields[column]), m); break;
                    }

                    case ID_Mix:
                    {
                        succceed = true;
                        ref.add(fields[0], 0.0, stof(fields[column]), m); break;
                    }
                }
            }, delim);

            return succceed;
        }
        catch (const std::exception &ex)
        {
            std::cout << ex.what() << std::endl;
            std::cerr << "[Error]: Error in the mixture file" << std::endl;
            throw;
        }
    };

    if (!f("\t"))
    {
        if (!f(","))
        {
            throw std::runtime_error("Failed to read any sequin in the mixture file. A CSV file format is expected. Please check and try again.");
        }
    }
}

void Standard::v_std(const Reader &r)
{
    ParserBED::parse(r, [&](const ParserBED::Annotation &f, const ParserProgress &)
    {
        r_var.addStand(f.id, f.l);
    });
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
    readMixture(r, r_var, Mix_1, ID_Length_Mix, 2);
}

void Standard::m_ref(const Reader &r)
{
    /*
     * MG_15   1       4720    MG_15   1       +
     * MG_19   1       4216    MG_19   1       +
     */
    
    ParserBED::parse(r, [&](const ParserBED::Annotation &f, const ParserProgress &)
    {
        r_meta.addStand(f.id, f.l.length());
    });
}

void Standard::m_mix(const Reader &r)
{
    const auto n = countColumns(r);
    readMixture(r, r_meta, Mix_1, ID_Length_Mix);

    if (n >= 4)
    {
        readMixture(Reader(r), r_meta, Mix_2, ID_Length_Mix, 3);
    }
}

void Standard::l_mix(const Reader &r)
{
    readMixture(r, r_lad, Mix_1, ID_Mix, 1);
    readMixture(Reader(r), r_lad, Mix_2, ID_Mix, 2);
}

void Standard::f_mix(const Reader &r)
{
    readMixture(r, r_fus, Mix_1, ID_Length_Mix, 2);
}

void Standard::f_ref(const Reader &r)
{
    ParserCSV::parse(r, [&](const ParserCSV::Fields &f, const ParserProgress &)
    {
        if (f[0] != "chrT-chrT")
        {
            throw std::runtime_error("Invalid reference file. chrT-chrT is expected.");
        }

        FusionRef::FusionPoint b;

        b.id = f[4];;
        b.l1 = stod(f[1]) + 1;
        b.l2 = stod(f[2]) + 1;

        if      (f[3] == "ff") { b.s1 = Strand::Forward;  b.s2 = Strand::Forward;  }
        else if (f[3] == "fr") { b.s1 = Strand::Forward;  b.s2 = Strand::Backward; }
        else if (f[3] == "rf") { b.s1 = Strand::Backward; b.s2 = Strand::Forward;  }
        else if (f[3] == "rr") { b.s1 = Strand::Backward; b.s2 = Strand::Backward; }

        r_fus.addBreak(b);
    }, "\t");
}

void Standard::r_ref(const Reader &r)
{
    ParserGTF::parse(r, [&](const Feature &f, const std::string &, const ParserProgress &)
    {
       if (f.id == id && f.type == Exon)
        {
            r_trans.addRef(f.tID, f.geneID, f.l);
        }
    });
}

void Standard::r_mix(const Reader &r)
{
    const auto n = countColumns(r);
    readMixture(Reader(r), r_trans, Mix_1, ID_Length_Mix, 2);
    
    if (n >= 4)
    {
        readMixture(Reader(r), r_trans, Mix_2, ID_Length_Mix, 3);
    }
}