/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include <fstream>
#include <algorithm>
#include "VarQuin/v_seq.hpp"
#include "parsers/parser_fq.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

void VSeq::analyze(const FileName &f1,
                   const FileName &f2,
                   const FileName &align,
                   const Options &o)
{    
    std::set<ReadName> names;

    ParserSAM::parse(align, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (!x.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        if (Standard::isSynthetic(x.cID))
        {
            names.insert("@" + x.name);
        }
    });
    
    auto f = [&](const FileName &f, const FileName &output)
    {
        std::ofstream out;
        out.open(output, std::ios_base::app);
        
        ParserFQ::parse(Reader(f), [&](ParserFQ::Data &x, const ParserProgress &)
        {
            if (names.count(x.name))
            {
                std::reverse(x.seq.begin(),  x.seq.end());
                std::reverse(x.qual.begin(), x.qual.end());
            }

            out << x.name << " " << x.info << std::endl;
            out << x.seq  << std::endl;
            out << x.opt  << std::endl;
            out << x.qual << std::endl;
        });
        
        out.close();
    };

    f(f1, o.work + "/VarSeq_1.fq");
    f(f2, o.work + "/VarSeq_2.fq");
}

void VSeq::report(const std::vector<FileName> &files, const Options &o)
{
    analyze(files[0], files[1], files[2], o);
}