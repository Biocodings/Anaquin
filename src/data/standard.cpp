#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
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

// Static definition
std::set<ChrID> Standard::genoIDs;

enum MixtureFormat
{
    ID_Length_Mix, // Eg: MG_33  10  60.23529412
    ID_Mix,        // Eg: MG_33  60.23529412
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

template <typename Reference> void readMixture(const Reader &r, Reference &ref, Mixture m, MixtureFormat format, unsigned column=2)
{
    auto f = [&](const std::string &delim)
    {
        const auto t = Reader(r);

        bool succceed = false;
        
        ParserCSV::parse(t, [&](const ParserCSV::Data &d, const ParserProgress &p)
        {
            // Don't bother if this is the first line or an invalid line
            if (p.i == 0 || d.size() <= 1)
            {
                return;
            }
            
            const auto seq = d[0];
            const auto con = stof(d[column]);
            
            switch (format)
            {
                case ID_Length_Mix:
                {
                    succceed = true;
                    ref.add(seq, stoi(d[1]), con, m); break;
                }
                    
                case ID_Mix:
                {
                    succceed = true;
                    ref.add(seq, 0.0, con, m); break;
                }
            }
        }, delim);
        
        return succceed;
    };
    
    if (!f("\t") && !f(","))
    {
        throw std::runtime_error("No sequin is found in the mixture file. Please check and try again.");
    }
}

ChrID Standard::toReverse(const ChrID &cID)
{
    auto x = cID;
    
    // Eg; chr2 to chrev2
    boost::replace_all(x, "chr", "chrev");

    A_ASSERT(!x.empty());
    
    return x;
}

bool Standard::isSynthetic(const ChrID &cID)
{
    assert(!cID.empty());
    
    const std::set<ChrID> sIDs = { "chrT", "chrIS" };

    // Can we match by exact?
    if (sIDs.count(cID))
    {
        return true;
    }

    // Eg; chrev10
    else if (cID.find("rev") != std::string::npos)
    {
        return true;
    }
    else if (cID.find("GS_") != std::string::npos)
    {
        return true;
    }
    else if (cID.find("GI_") != std::string::npos)
    {
        return true;
    }
    
    return false;
}

void Standard::addMMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file. Expected three columns for a single mixture.");
    
    readMixture(Reader(r), r_meta, Mix_1, ID_Length_Mix, 2);
}

void Standard::addVMix(const Reader &r)
{
    readMixture(r, r_var, Mix_1, ID_Length_Mix, 2);
}

void Standard::addRMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 3, "Invalid mixture file. Expected three columns for a single mixture.");

    readMixture(Reader(r), r_rna, Mix_1, ID_Length_Mix, 2);
}

void Standard::addRDMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected four columns for a double mixture.");
    
    readMixture(Reader(r), r_rna, Mix_1, ID_Length_Mix, 2);
    readMixture(Reader(r), r_rna, Mix_2, ID_Length_Mix, 3);
}
