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
#include "parsers/parser_bed.hpp"
#include "parsers/parser_gtf.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

enum MixtureFormat
{
    X_M,
    M_X_M,
    X_X_X_M,
    X_M_X_M,
};

enum TranslateFormat
{
    F_T,
    T_F
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

std::shared_ptr<GTFData> Standard::readGTF(const Reader &r)
{
    return std::shared_ptr<GTFData>(new GTFData(gtfData(r)));
}

template <typename Reference> Translate readTranslate(const Reader &r, Reference &ref, TranslateFormat format, Translate x = Translate())
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
                case F_T: { x.add(d[0], d[1]); break; }
                case T_F: { x.add(d[1], d[0]); break; }
            }
        }, delim);
        
        return x.size();
    };
    
    if (!parse("\t") && !parse(","))
    {
        throw std::runtime_error("No sequin found in the reference file. Please check and try again.");
    }
    
    return x;
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
                case X_M:     { x.add(d[0], m, stof(d[1])); break; }
                case X_M_X_M: { x.add(d[1], m, stof(d[3])); break; }
                case M_X_M:   { x.add(d[0], m, stof(d[2])); break; }
                case X_X_X_M: { x.add(d[0], m, stof(d[3])); break; }
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

Ladder Standard::readLength(const Reader &r)
{
    A_CHECK(countColumns(r) >= 2, "Invalid mixture file. Expected two or more columns.");
    return readLadder(Reader(r), r_rna, Mix_1, X_M);
}

Ladder Standard::addCNV(const Reader &r)
{
    A_CHECK(countColumns(r) == 2, "Invalid mixture file for CNV ladder.");
    return readLadder(Reader(r), r_var, Mix_1, X_M);
}

Ladder Standard::addCon1(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readLadder(Reader(r), r_var, Mix_1, M_X_M);
}

Ladder Standard::addCon2(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readLadder(Reader(r), r_var, Mix_1, X_M_X_M);
}

Translate Standard::addSeq2Unit(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readTranslate(Reader(r), r_var, F_T);
}

Translate Standard::addUnit2Seq(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file for conjoint ladder.");
    return readTranslate(Reader(r), r_var, T_F);
}

Ladder Standard::addAF(const Reader &r)
{
    A_CHECK(countColumns(r) == 2, "Invalid mixture file for allele frequnecy ladder.");
    return readLadder(Reader(r), r_var, Mix_1, X_M);
}

Ladder Standard::addMMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected three or more columns.");
    auto l = readLadder(Reader(r), r_meta, Mix_1, X_M);
    return readLadder(Reader(r), r_meta, Mix_2, M_X_M, l);
}

Ladder Standard::readIsoform(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected three columns.");
    auto l = readLadder(Reader(r), r_rna, Mix_1, M_X_M);
    return readLadder(Reader(r), r_rna, Mix_2, X_X_X_M, l);
}

Ladder Standard::readIDiff(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected three columns.");
    auto l = readIsoform(r);
    Ladder f;
    
    for (const auto &i : l.seqs)
    {
        f.add(i, Mix_1, log2(l.input(i, Mix_2) / l.input(i, Mix_1)));
    }
    
    return f;
}

Ladder Standard::readGDiff(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected three columns.");
    auto l = readGene(r);
    Ladder f;
    
    for (const auto &i : l.seqs)
    {
        f.add(i, Mix_1, log2(l.input(i, Mix_2) / l.input(i, Mix_1)));
    }
    
    return f;
}

template <typename F, typename T = VCFData> T parseVCF2(const Reader &r, F f)
{
    T t;
    
    ParserVCF::parse(r, [&](const Variant &x)
    {
        if (f(x))
        {
            t[x.cID].b2v[x.l.start] = x;
            t[x.cID].m2v[x.type()].insert(x);
        }
    });
    
    return t;
}

VCFLadder Standard::addVCF(const Reader &r, const std::set<SequinVariant::Context> &f)
{
    typedef SequinVariant::Context Context;
    
    VCFLadder v;

    // VCF variants
    v.data = parseVCF2(r, [&](const Variant &x)
    {
        A_ASSERT(x.key());
        
        auto longVar = [&]()
        {
            SequinVariant s;
            
            s.gt = Genotype::Heterzygous;
            
            v.vIDs.insert(x.name);
            v.sVars[x.key()] = s;
            
            return true;
        };
        
        auto shortVar = [&]()
        {
            const auto m1 = std::map<std::string, SequinVariant::Context>
            {
                { "cancer",     Context::Cancer       },
                { "common",     Context::Common       },
                { "high_gc",    Context::HighGC       },
                { "long_di",    Context::LongDinRep   },
                { "long_homo",  Context::LongHompo    },
                { "low_gc",     Context::LowGC        },
                { "long_quad",  Context::LongQuadRep  },
                { "long_tri",   Context::LongTrinRep  },
                { "short_di",   Context::ShortDinRep  },
                { "short_homo", Context::ShortHompo   },
                { "short_quad", Context::ShortQuadRep },
                { "v_high_gc",  Context::VeryHighGC   },
                { "v_low_gc",   Context::VeryLowGC    },
                { "short_tri",  Context::ShortTrinRep }
            };
            
            const auto m2 = std::map<std::string, Genotype>
            {
                { "SOM",     Genotype::Somatic     },
                { "HOM",     Genotype::Homozygous  },
                { "HOM_CNV", Genotype::Homozygous  },
                { "HET",     Genotype::Heterzygous },
            };
            
            auto throwInvalidRef = [&](const std::string &x)
            {
                throw std::runtime_error(r.src() + " doesn't seem to be a valid VCF reference file. Reason: " + x);
            };
            
            if (!x.ifs.count("CX") || !m1.count(x.ifs.at("CX")))
            {
                throwInvalidRef("The CX field is not found or invalid");
            }
            else if (!x.ifs.count("GT") || !m2.count(x.ifs.at("GT")))
            {
                throwInvalidRef("The GT field is not found or invalid");
            }
            
            SequinVariant s;
            
            s.gt   = m2.at(x.ifs.at("GT"));
            s.ctx  = m1.at(x.ifs.at("CX"));
            s.copy = x.iff.at("CP");
            
            if (f.count(s.ctx))
            {
                return false;
            }
            
            v.vIDs.insert(x.name);
            
            Concent af;
            
            switch (s.gt)
            {
                case Genotype::Somatic:     { af = x.allF; break; }
                case Genotype::Homozygous:  { af = 1.0;    break; }
                case Genotype::Heterzygous: { af = 0.5;    break; }
            }
            
            // Update allele frequency ladder
            v.af.add(x.name, Mix_1, af);
            
            v.sVars[x.key()] = s;
            
            return true;
        };
        
        if (x.isSV()) { return longVar();  }
        else          { return shortVar(); }
    });
    
    return v;
}

Ladder Standard::readGeneL(const Reader &r)
{
    auto l = readLength(r);
    
    // Ladder for genes
    Ladder genes;
    
    auto aggregate = [&](Mixture m)
    {
        for (const auto &i : l.seqs)
        {
            const auto gene = isoform2Gene(i);
            
            if (genes.seqs.count(gene))
            {
                genes.add(gene, m, std::max(l.input(i, Mix_1), genes.input(gene, m)));
            }
            else
            {
                genes.add(gene, m, l.input(i, Mix_1));
            }
        }
    };
    
    aggregate(Mix_1);
    return genes;
}

Ladder Standard::readGene(const Reader &r)
{
    auto l = readIsoform(r);
    
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
