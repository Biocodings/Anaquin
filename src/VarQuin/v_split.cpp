#include "VarQuin.hpp"
#include "VarQuin/v_split.hpp"
#include "writers/bam_writer.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

void VSplit::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.info("Generating sample-derived alignments");
    o.info("Generating sequin-derived alignments");
    o.info("Generating sample-derived alignments within sequin regions");
    
    BAMWriter w1, w2, w3;
    w1.open(o.work + "/VarFlip_sample.bam");
    w2.open(o.work + "/VarFlip_sample_regions.bam");
    w3.open(o.work + "/VarFlip_sequins.bam");

    // Regions without edge effects
    const auto regs = r.regs1();
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        /*
         * Eg: samtools view -h -L hg38.bed A.bam | grep -v chrev | samtools view -bS > VarFlip_sample.bam
         */
        
        // Sample-derived reads
        if (x.cID != "*" && !isRevChr(x.cID) && !isRevChr(x.rnext))
        {
            w1.write(x);
            
            /*
             * Eg: samtools view -b -h -L sequin_regions.hg38.bed sample_normal.bam > sample_normal_regions.bam
             */

            if (regs.count(x.cID) && regs.at(x.cID).contains(x.l))
            {
                w2.write(x);
            }
        }
        
        /*
         * Eg: samtools view -b -L hg38rev.bed normal.bam > sequins.bam
         */
        
        // Sequin-derived reads
        if (isRevChr(x.cID))
        {
            w3.write(x);
        }        
    }, true);

    w1.close();
    w2.close();
    w3.close();
}
