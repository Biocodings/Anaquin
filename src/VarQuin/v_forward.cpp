/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include <fstream>
#include <algorithm>
#include "VarQuin/v_forward.hpp"
#include "parsers/parser_fq.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

void VForward::analyze(const FileName &f1,
                       const FileName &f2,
                       const FileName &align,
                       const Options &o)
{    
    std::set<ReadName> names;

    ParserSAM::parse(align, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (!x.i && info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (Standard::isSynthetic(x.cID))
        {
            /*
             * Eg: @GV_IDEL_011863_R-1400/1
             * Eg: @GV_IDEL_011863_R-1400/2
             */
            
            names.insert("@" + x.name + "/1");
            names.insert("@" + x.name + "/2");
        }
    });
    
    auto f = [&](const FileName &f, const FileName &syn, const FileName &gen)
    {
        std::ofstream fs, fg;
 
        fs.open(syn, std::ios_base::app);
        fg.open(gen, std::ios_base::app);
        
        ParserFQ::parse(Reader(f), [&](ParserFQ::Data &x, const ParserProgress &)
        {
            if (names.count(x.name))
            {
                std::reverse(x.seq.begin(),  x.seq.end());
                std::reverse(x.qual.begin(), x.qual.end());
                
                fs << x.name << " " << x.info << std::endl;
                fs << x.seq  << std::endl;
                fs << x.opt  << std::endl;
                fs << x.qual << std::endl;
            }
            else
            {
                fg << x.name << " " << x.info << std::endl;
                fg << x.seq  << std::endl;
                fg << x.opt  << std::endl;
                fg << x.qual << std::endl;
            }
        });
        
        fs.close();
        fg.close();
    };

    o.info("Generating first mates");
    f(f1, o.work + "/VarForward_sequins_2.fq", o.work + "/VarForward_genome_1.fq");

    o.info("Generating second mates");
    f(f2, o.work + "/VarForward_sequins_1.fq", o.work + "/VarForward_genome_2.fq");
}

void VForward::report(const std::vector<FileName> &files, const Options &o)
{
    analyze(files[0], files[1], files[2], o);
}