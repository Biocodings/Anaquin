#include "VarQuin.hpp"
#include "VarQuin/v_split.hpp"
#include "parsers/parser_bam.hpp"
#include "writers/bam_writer.hpp"

using namespace Anaquin;

void VSplit::analyze(const FileName &file, const Options &o)
{
    o.info("Generating sample-derived alignments");
    o.info("Generating sequin-derived alignments");
    o.info("Generating sample-derived regional alignments");
    
    BAMWriter w1, w2, w3;
    w1.open("A.bam");
    w2.open("B.bam");
    w3.open("C.bam");
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        /*
         * Eg: samtools view -h -L hg38.bed normal.bam | grep -v chrev | samtools view -bS > sample_normal.bam
         */
        
        // Sample-derived reads
        if (!isRevChr(x.cID) && !isRevChr(x.rnext))
        {
            w1.write(x);
        }
        
        /*
         * Eg: samtools view -b -L hg38rev.bed normal.bam > sequin_normal_reverse.bam
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
