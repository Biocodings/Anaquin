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

void VForward::analyze(const FileName &file,
                       const FileName &output1,
                       const FileName &output2,
                       const Options &o)
{
    assert(!output1.empty());
    assert(!output2.empty());

    typedef std::string AlignmentID;
    
    /*
     * Attempt 1: Build up the position for next mates
     */
    
    // Mapping for PNEXT from reverse strand to forward strand
    std::map<Base, Base> r2f;
    
    ParserSAM::parse(file, [&](ParserSAM::Data &data, const ParserSAM::Info &info)
    {
        if (!data.i && info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }

        if (!data.i && Standard::isSynthetic(data.cID))
        {
            const auto t = data.l;
            
            reverse(data, info);
            replace(data.cID, "rev", "r");

            // The two strands can't make it equal
            assert(t.start != data.l.start);
            
            // Now save it for the next attempt
            r2f[t.start] = data.l.start;
        }
    });

    // Used for genomic alignments
    std::ofstream out1;
    
    // Used for synthetic alignments
    std::ofstream out2;

    /*
     * Attempt 2: Repeat but we now have positions for the next mates
     */
    
    ParserSAM::parse(file, [&](ParserSAM::Data &data, const ParserSAM::Info &info)
    {
        if (!data.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
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

        // Only if this is the first block
        if (data.i)
        {
            return;
        }
        
        const auto isSync = Standard::isSynthetic(data.cID);

        if (isSync)
        {
            //std::cout << data.cID << std::endl;
            
            reverse(data, info);
            replace(data.cID, "rev", "r");

            // This is the key, we now know the position of the next mate (TODO: Fix this!)
            //data.pnext = r2f.at(data.pnext + 1);
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
        if (isSync)
        {
            out2 << str << std::endl;
        }
        else
        {
            out1 << str << std::endl;
        }
    });

    if (out1.is_open()) { out1.close(); }
    if (out2.is_open()) { out2.close(); }
}

void VForward::report(const FileName &file, const Options &o)
{
    analyze(file, o.work + "/VarForward_genome.sam", o.work + "/VarForward_sequins.sam");
}