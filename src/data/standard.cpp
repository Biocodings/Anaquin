#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "data/errors.hpp"
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

template <typename Reference> void readMixture(const Reader &r, Reference &ref, Mixture m, MixtureFormat format, unsigned column=2)
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
        catch (...)
        {
            throw std::runtime_error("[Error]: Error in the mixture file. Please check and try again.");
        }
    };

    if (!f("\t") && !f(","))
    {
        throw std::runtime_error("[Error]: No sequin is found in the mixture file. Please check and try again.");
    }
}

void Standard::addInters(const Reader &r)
{
    ParserBed::parse(r, [&](const ParserBed::Data &f, const ParserProgress &)
    {
        r_var.addRInterval(f.cID, Interval(f.id, f.l));
    });
}

void Standard::addStd(const Reader &r)
{
    ParserBed::parse(r, [&](const ParserBed::Data &f, const ParserProgress &)
    {
        r_var.addStand(f.id, f.l);
    });
}

void Standard::addVar(const Reader &r)
{
    std::vector<std::string> toks;

    ParserVCF::parse(r, [&](const ParserVCF::Data &x, const ParserProgress &)
    {
        Variant v;
        
        // Eg: D_1_3_R
        v.id  = x.id;

        v.l   = x.l;
        v.alt = x.alt;
        v.ref = x.ref;

        r_var.addVar(v);
    });
}

void Standard::addVMix(const Reader &r)
{
    readMixture(r, r_var, Mix_1, ID_Length_Mix, 2);
}

void Standard::addMRef(const Reader &r)
{
    ParserBed::parse(r, [&](const ParserBed::Data &f, const ParserProgress &)
    {
        r_meta.addStand(f.cID, f.l);
    });
}

void Standard::addMMix(const Reader &r)
{
    const auto n = countColumns(r);
    readMixture(r, r_meta, Mix_1, ID_Length_Mix);

    if (n >= 4)
    {
        readMixture(Reader(r), r_meta, Mix_2, ID_Length_Mix, 3);
    }
}

void Standard::addLMix(const Reader &r)
{
    const auto n = countColumns(r);
    readMixture(Reader(r), r_lad, Mix_1, ID_Length_Mix, 1);

    if (n >= 4)
    {
        readMixture(Reader(r), r_lad, Mix_2, ID_Length_Mix, 2);
    }
}

void Standard::addFMix(const Reader &r)
{
    readMixture(r, r_fus, Mix_1, ID_Length_Mix, 2);
}

void Standard::addFStd(const Reader &r)
{
    ParserBed::parse(r, [&](const ParserBed::Data &f, const ParserProgress &)
    {
        r_fus.addStand(f.id, f.l);
    });
}

void Standard::addFSplice(const Reader &r)
{
    ParserBed::parse(r, [&](const ParserBed::Data &f, const ParserProgress &)
    {
        r_fus.addSplice(f.id, f.l);
    });
}

void Standard::addFRef(const Reader &r)
{
    ParserCSV::parse(r, [&](const ParserCSV::Fields &f, const ParserProgress &)
    {
        if (f[0] != "chrT-chrT")
        {
            throw std::runtime_error("Invalid reference file. chrT-chrT is expected.");
        }

        FusionRef::KnownFusion known;

        known.id = f[4];;
        known.l1 = stod(f[1]) + 1;
        known.l2 = stod(f[2]) + 1;

        if      (f[3] == "ff") { known.s1 = Strand::Forward;  known.s2 = Strand::Forward;  }
        else if (f[3] == "fr") { known.s1 = Strand::Forward;  known.s2 = Strand::Backward; }
        else if (f[3] == "rf") { known.s1 = Strand::Backward; known.s2 = Strand::Forward;  }
        else if (f[3] == "rr") { known.s1 = Strand::Backward; known.s2 = Strand::Backward; }

        r_fus.addFusion(known);
    }, "\t");
}

void Standard::addTRef(const Reader &r)
{
    ParserGTF::parse(r, [&](const Feature &f, const std::string &, const ParserProgress &)
    {
        switch (f.type)
        {
            case Gene: { r_trans.addGene(f.cID, f.gID, f.l);        break; }
            case Exon: { r_trans.addExon(f.cID, f.gID, f.tID, f.l); break; }
            default:   { break; }
        }
    });
}

void Standard::addTMix(const Reader &r)
{
    ASSERT(countColumns(r) == 3, "Invalid mixture file. Expected three columns for a single mixture.");

    readMixture(Reader(r), r_trans, Mix_1, ID_Length_Mix, 2);
}

void Standard::addTDMix(const Reader &r)
{
    ASSERT(countColumns(r) == 4, "Invalid mixture file. Expected four columns for a double mixture.");
    
    readMixture(Reader(r), r_trans, Mix_1, ID_Length_Mix, 2);
    readMixture(Reader(r), r_trans, Mix_2, ID_Length_Mix, 3);
}