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
#include "parsers/parser_bed.hpp"
#include "parsers/parser_gtf.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

enum MixtureFormat
{
    Name_Mix,
    Name_X_Mix,
    Name_Len_X_Mix,
    NUMU,
    UNUM,
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

GTFData Standard::readGTF(const Reader &r)
{
    return gtfData(r);
}

template <typename Reference> Ladder readLadder(const Reader &r, Reference &ref, Mixture m, MixtureFormat format, Ladder x = Ladder())
{
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
                case Name_Mix:       { x.add(d[0], m, stof(d[1])); break; }
                case NUMU:           { x.add(d[0], m, stof(d[2])); break; }
                case UNUM:           { x.add(d[1], m, stof(d[3])); break; }
                case Name_X_Mix:     { x.add(d[0], m, stof(d[2])); break; }
                case Name_Len_X_Mix: { x.add(d[0], m, stof(d[3])); break; }
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
    A_CHECK(countColumns(r) == 2, "Invalid mixture file for CNV ladder.");
    return readLadder(Reader(r), r_var, Mix_1, Name_Mix);
}

Ladder Standard::addCon1(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readLadder(Reader(r), r_var, Mix_1, NUMU);
}

Ladder Standard::addCon2(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readLadder(Reader(r), r_var, Mix_1, UNUM);
}

Ladder Standard::addAF(const Reader &r)
{
    A_CHECK(countColumns(r) == 2, "Invalid mixture file for allele frequnecy ladder.");
    return readLadder(Reader(r), r_var, Mix_1, Name_Mix);
}

Ladder Standard::addMMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected three or more columns.");
    auto l = readLadder(Reader(r), r_meta, Mix_1, Name_Mix);
    return readLadder(Reader(r), r_meta, Mix_2, Name_X_Mix, l);
}

Ladder Standard::addIsoform(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file. Expected three columns.");
    auto l = readLadder(Reader(r), r_rna, Mix_1, Name_Mix);
    return readLadder(Reader(r), r_rna, Mix_2, Name_X_Mix, l);
}

Ladder Standard::addGene(const Reader &r)
{
    auto l = addIsoform(r);
    
    // Ladder for genes
    Ladder genes;
    
    auto aggregate = [&](Mixture m)
    {
        for (const auto &i : l.seqs)
        {
            const auto gene = isoform2Gene(i);
            
            if (genes.seqs.count(gene))
            {
                genes.add(gene, m, l.input(i, m) + genes.input(gene, m));
            }
            else
            {
                genes.add(gene, m, l.input(i, m));
            }
        }
    };
    
    aggregate(Mix_1);
    genes.seqs.clear();
    aggregate(Mix_2);

    return genes;
}


