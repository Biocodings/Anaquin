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
                                  const FileName &output,
                                  const Writer &writer,
                                  const Options &o)
{
    assert(!output.empty());

    std::ofstream out;
    //const auto &r = Standard::instance().r_var;

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
            auto f = sam_open(output.c_str(), "w");
            sam_hdr_write(f, reinterpret_cast<bam_hdr_t *>(info.header));
            sam_close(f);
            out.open(output, std::ios_base::app);
        }

        // Reverse the alignment
        reverse(data, info);
        
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
        out << str << std::endl;
    });

    if (out.is_open())
    {
        out.close();
    }

    return VForward::Stats();
}

struct FFileWriter : public VForward::Writer
{
    inline void write(void *data, void *header) const
    {
        const auto *b = reinterpret_cast<bam1_t *>(data);
        const auto *h = reinterpret_cast<bam_hdr_t *>(header);
        writer.write(h, b);
    }

    mutable WriterSAM writer;
};

void VForward::report(const FileName &file, const Options &o)
{
    FFileWriter f;
    f.writer.open(o.work + "/VarForward_forward.sam");

    const auto &stats = analyze(file, o.work + "/VarForward_forward.sam", f);
    
    /*
     * Generating VarForward_summary.stats
     */    
}