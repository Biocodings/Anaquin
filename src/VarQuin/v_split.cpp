#include "VarQuin.hpp"
#include "VarQuin/v_split.hpp"
#include "parsers/parser_bam.hpp"
#include "writers/bam_writer.hpp"

using namespace Anaquin;

void VSplit::report(const FileName &file, const Options &o)
{
    o.info("Generating sample-derived alignments");
    o.info("Generating sequin-derived alignments");
    o.info("Generating sample-derived regional alignments");
    
    BAMWriter w1, w2, w3;
    w1.open(o.work + "/VarFlip_sample.bam");
    w2.open(o.work + "/VarFlip_sample_regions.bam");
    w3.open(o.work + "/VarFlip_sequins.bam");
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        /*
         * Eg: samtools view -h -L hg38.bed alignments.bam | grep -v chrev | samtools view -bS > sample.bam
         */
        
        // Sample-derived reads
        if (!isRevChr(x.cID) && !isRevChr(x.rnext))
        {
            w1.write(x);
        }
        
        /*
         * Eg: samtools view -b -h -L sequin_regions.hg38.bed sample_normal.bam > sample_normal_regions.bam
         */
        
        // Sample-derived regional reads
        if (!isRevChr(x.cID) && !isRevChr(x.rnext))
        {
            // ....
        }
        
        /*
         * Eg: samtools view -b -L hg38rev.bed alignments.bam > sequins.bam
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
