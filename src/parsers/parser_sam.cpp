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

            /*
             * What to do with something like "528084 50 79M10250N22M"? We'd split it into three independent
             * blocks.
             */
            
            for (int k = 0; k < t->core.n_cigar; k++)
            {
                const int op = bam_cigar_op(cigar[k]);
                const int ol = bam_cigar_oplen(cigar[k]);

                // 1-based leftmost coordinate is assumed
                align.l.set(n + 1, n + ol);

                // We'll need it for the next operation
                n += ol;
                
                if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CDEL)
                {
                    align.spliced = false;
                }
                else if (op == BAM_CREF_SKIP)
                {
                    align.spliced = true;
                }
                
                x(align);
            }
        }
        else
        {
            x(align);
        }
    }

    sam_close(f);
}