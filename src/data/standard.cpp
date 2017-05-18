#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "data/bData.hpp"
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "tools/errors.hpp"
#include "data/standard.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_gtf.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

enum MixtureFormat
{
    Name_Len_Mix,
    Name_Len_X_Mix,
    Name_Len_Mix_Mix
};

static unsigned countColumns(const Reader &r)
{
    std::size_t n = 0;

    ParserCSV::parse(r, [&](const ParserCSV::Data &d, const ParserProgress &p)
    {
        n = std::max(n, d.size());
    }, ",");

    ParserCSV::parse(Reader(r), [&](const ParserCSV::Data &d, const ParserProgress &p)
    {
        n = std::max(n, d.size());
    }, "\t");
    
    return static_cast<unsigned>(n);
}

BedData Standard::readBED(const Reader &r, Base trim)
{
    RegionOptions o;
    o.trim = trim;
    return readRegions(Reader(r), [&](const ParserBed::Data &, const ParserProgress &) {}, o);
}

template <typename Reference> Ladder readLadder(const Reader &r, Reference &ref, Mixture m, MixtureFormat format)
{
    Ladder x;
    
    auto parse = [&](const std::string &delim)
    {
        const auto t = Reader(r);
        
        ParserCSV::parse(t, [&](const ParserCSV::Data &d, const ParserProgress &p)
        {
            // Don't bother if this is the first line or an invalid line
            if (p.i == 0 || d.size() <= 1)
            {
                return;
            }
            
            switch (format)
            {
                case Name_Len_Mix:     { x.add(d[0], m, stof(d[2]));              break; }
                case Name_Len_X_Mix:   { x.add(d[0], m, stof(d[3]));              break; }
                case Name_Len_Mix_Mix: { x.add(d[0], m, stoi(d[2]) * stoi(d[3])); break; }
            }
        }, delim);
        
        return x.count();
    };
    
    if (!parse("\t") && !parse(","))
    {
        throw std::runtime_error("No sequin found in the ladder file. Please check and try again.");
    }

    return x;
}

Ladder Standard::addCNV(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file for CNV ladder.");
    return readLadder(Reader(r), r_var, Mix_1, Name_Len_Mix);
}

Ladder Standard::addCon(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readLadder(Reader(r), r_var, Mix_1, Name_Len_Mix_Mix);
}

Ladder Standard::addAF(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file for allele frequnecy ladder.");
    return readLadder(Reader(r), r_var, Mix_1, Name_Len_Mix);
}

void Standard::addMDMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected four columns for a double mixture.");
    
    readLadder(Reader(r), r_meta, Mix_1, Name_Len_Mix);
    readLadder(Reader(r), r_meta, Mix_2, Name_Len_Mix);
}

void Standard::addMMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file. Expected three columns for a single mixture.");
    readLadder(Reader(r), r_meta, Mix_1, Name_Len_Mix);
}

void Standard::addRMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file. Expected three columns for a single mixture.");
    readLadder(Reader(r), r_rna, Mix_1, Name_Len_Mix);
}

void Standard::addRDMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected four columns for a double mixture.");
    
    readLadder(Reader(r), r_rna, Mix_1, Name_Len_Mix);
    readLadder(Reader(r), r_rna, Mix_2, Name_Len_Mix);
}
