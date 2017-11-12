#include <fstream>
#include <sstream>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "VarQuin.hpp"
#include "VarQuin/v_split.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

struct BAMWriter
{
    inline void close()
    {
        sam_close(_file);
        sam_close(_null);
    }
    
    inline void open(const FileName &file)
    {
        _file = sam_open(file.c_str(), "wb");
        _null = sam_open("/dev/null", "w");
    }

    template <typename F> void filterH(const ParserBAM::Data &x, F f)
    {
        const auto *h1 = reinterpret_cast<bam_hdr_t *>(x.h());
        const auto x1  = std::string(h1->text);
        
        std::stringstream ss;
        std::istringstream iss(x1);
        
        for (std::string line; std::getline(iss, line); )
        {
            if (f(line))
            {
                ss << line << "\n";
            }
        }
        
        std::ofstream w;
        w.open ("__tmp__.sam");
        w << ss.str();
        w.close();
        
        samFile *x2 = sam_open("__tmp__.sam", "r");
        bam_hdr_t * h2 = sam_hdr_read(x2);
        
        if (sam_hdr_write(_file, h2) == -1)
        {
            throw std::runtime_error("sam_hdr_write failed");
        }
        
        system("rm __tmp__.sam");
    }
    
    inline void fHeader(const ParserBAM::Data &x)
    {
        return filterH(x, [&](const std::string &x)
        {
            return x.find("chrev") == std::string::npos;
        });
    }

    inline void header(const ParserBAM::Data &x)
    {
        if (sam_hdr_write(_file, reinterpret_cast<bam_hdr_t *>(x.h())) == -1)
        {
            throw std::runtime_error("sam_hdr_write failed");
        }
    }
    
    inline const char *record(const ParserBAM::Data &x)
    {
        const auto *b = reinterpret_cast<bam1_t *>(x.b());
        const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());
        
        if (sam_write1(_null, h, b) == -1)
        {
            throw std::runtime_error("sam_write1 failed");
        }
        
        return _null->line.s;
    }
    
    inline void write(const ParserBAM::Data &x)
    {
        const auto *b = reinterpret_cast<bam1_t *>(x.b());
        const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());
        
        if (sam_write1(_file, h, b) == -1)
        {
            throw std::runtime_error("sam_write1 failed");
        }
    }
    
    samFile *_file, *_null;
};

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

    const auto regs = r.regs1();
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        static bool header = false;

        if (!header)
        {
            w3.header(x);
            w1.fHeader(x);
            w2.fHeader(x);
        }

        header = true;
        
        /*
         * ******************** Sample derived reads ********************
         *
         *   Eg: samtools view -h -L hg38.bed A.bam | grep -v chrev | samtools view -bS > VarFlip_sample.bam
         */
        
        if (x.cID != "*" && !isRevChr(x.cID) && !isRevChr(x.rnext))
        {
            if (!strstr(w1.record(x), "chrev"))
            {
                w1.write(x);
                
                /*
                 * ******************** Sample derived regional reads ********************
                 *
                 *   Eg: samtools view -b -h -L sequin_regions.hg38.bed sample_normal.bam > sample_normal_regions.bam
                 */
                
                if (regs.count(x.cID) && regs.at(x.cID).contains(x.l))
                {
                    w2.write(x);
                }
            }
        }
        
        /*
         * ******************** Sequin derived reads ********************
         *
         *   Eg: samtools view -b -L hg38rev.bed normal.bam > sequins.bam
         */
        
        if (isRevChr(x.cID))
        {
            w3.write(x);
        }        
    }, true);

    w1.close();
    w2.close();
    w3.close();
}