/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "tools/samtools.hpp"
#include "VarQuin/v_forward.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/writer_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

VForward::Stats VForward::analyze(const FileName &file,
                                  const FileName &output1,
                                  const FileName &output2,
                                  const Options &o)
{
    assert(!output1.empty());
    assert(!output2.empty());

    // Used for genomic alignments
    std::ofstream out1;
    
    // Used for synthetic alignments
    std::ofstream out2;
    
    ParserSAM::parse(file, [&](ParserSAM::Data &data, const ParserSAM::Info &info)
    {
        if (!data.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        /*
         * Reuse the code for generating headers
         */
        
        if (!info.p.i && !data.i)
        {
            auto f = sam_open(output1.c_str(), "w");
            sam_hdr_write(f, reinterpret_cast<bam_hdr_t *>(info.header));
            sam_close(f);

            f = sam_open(output2.c_str(), "w");
            sam_hdr_write(f, reinterpret_cast<bam_hdr_t *>(info.header));
            sam_close(f);
            
            out1.open(output1, std::ios_base::app);
            out2.open(output2, std::ios_base::app);
        }

        const auto isSyn = Standard::isSynthetic(data.cID);
        
        // Only reverse if it's synthetic (mirror image)
        if (isSyn)
        {
            reverse(data, info);
            replace(data.cID, "rev", "r");
        }

        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";
        const auto str = (boost::format(format) % data.name
                                                % data.flag
                                                % data.cID
                                                % data.l.start
                                                % data.mapq
                                                % data.cigar
                                                % data.rnext
                                                % data.pnext
                                                % data.tlen
                                                % data.seq
                                                % data.qual).str();
        if (isSyn)
        {
            out2 << str << std::endl;
        }
        else
        {
            out1 << str << std::endl;
        }
    });

    if (out1.is_open()) { out1.close(); }
    if (out1.is_open()) { out1.close(); }

    return VForward::Stats();
}

void VForward::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file,
                                o.work + "/VarForward_genome.sam",
                                o.work + "/VarForward_sequins.sam");
    
    /*
     * Generating VarForward_summary.stats
     */    
}