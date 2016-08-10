/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include <fstream>
#include <algorithm>
#include "data/biology.hpp"
#include "VarQuin/v_flip.hpp"
#include "parsers/parser_fq.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

void VFlip::analyze(const FileName &seq1,
                    const FileName &seq2,
                    const FileName &align,
                    const Options  &o)
{    
    std::set<ReadName> names;

    ParserSAM::parse(align, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
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
                // DNA complement of the sequence
                complement(x.seq);
                
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

    o.info("Generating VarFlip_sequins_1.fq & VarFlip_genome_1.fq");
    f(seq1, o.work + "/VarFlip_sequins_1.fq", o.work + "/VarFlip_genome_1.fq");

    o.info("Generating VarFlip_sequins_2.fq & VarFlip_genome_2.fq");
    f(seq2, o.work + "/VarFlip_sequins_2.fq", o.work + "/VarFlip_genome_2.fq");
}

void VFlip::report(const std::vector<FileName> &files, const Options &o)
{
    analyze(files[0], files[1], files[2], o);
}