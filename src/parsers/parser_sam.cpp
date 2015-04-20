#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <htslib/sam.h>
#include "parsers/parser_sam.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

void ParserSAM::parse(const std::string &file, std::function<void (const Alignment &)> x)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }

    auto h = sam_hdr_read(f);
    auto t = bam_init1();

    Alignment align;

    while (sam_read1(f, h, t) >= 0)
    {
        align.id = std::string(h->target_name[0]);
        align.mapped = !(t->core.flag & BAM_FUNMAP);

        if (align.mapped)
        {
            const auto cigar = bam_get_cigar(t);
            
            // Length of the sequence
            int n = 0;
            
            // It's true only if we have a BAM_CREF_SKIP operation
            align.spliced = false;
            
            for (int k = 0; k < t->core.n_cigar; k++)
            {
                const int op = bam_cigar_op(cigar[k]);
                const int ol = bam_cigar_oplen(cigar[k]);
                
                if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CDEL)
                {
                    n += ol;
                }
                else if (op == BAM_CREF_SKIP)
                {
                    n += ol;
                    align.spliced = true;
                }
            }
            
            // 1-based leftmost coordinate is assumed
            align.l.set(t->core.pos, t->core.pos + n - 1);

            assert(n == align.l.length());
        }

        x(align);
    }

    sam_close(f);
}