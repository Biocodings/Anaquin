#include <htslib/sam.h>
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

void ParserSAM::parse(const FileName &file, Functor x)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }

    auto t = bam_init1();
    auto h = sam_hdr_read(f);

    Alignment align;
    Info info;
    
    while (sam_read1(f, h, t) >= 0)
    {
        info.p.i++;
        info.length = h->target_len[t->core.tid];

        align.i      = 0;
        align.mapped = false;
        info.data    = t;
        info.header  = h;

        if (t->core.tid < 0)
        {
            x(align, info);
            continue;
        }

        align.cID    = std::string(h->target_name[t->core.tid]);
        align.name   = bam_get_qname(t);
        align.mapped = !(t->core.flag & BAM_FUNMAP);

        if (align.mapped)
        {
            const auto cigar = bam_get_cigar(t);

            auto n = t->core.pos;
            
            // It's true only if we have a BAM_CREF_SKIP operation
            align.spliced = false;

            /*
             * Check if this is a spliced alignment. A spliced alignment is an alignment that has bases skipped.
             */
            
            info.spliced = false;
            
            for (auto i = 0; i < t->core.n_cigar; i++)
            {
                switch (bam_cigar_op(cigar[i]))
                {
                    case BAM_CREF_SKIP: { info.spliced = true; break; }
                    default: { break; }
                }
            }
            
            /*
             * What to do with something like "528084 50 79M10250N22M"? We'd split it into three blocks.
             */

            for (; align.i < t->core.n_cigar; align.i++)
            {
                const auto op = bam_cigar_op(cigar[align.i]);
                const auto ol = bam_cigar_oplen(cigar[align.i]);

                // 1-based leftmost coordinate is assumed
                align.l = Locus(n + 1, n + ol);

                // We'll need it for the next operation
                n += ol;
                
                if (op == BAM_CINS || op == BAM_CDEL)
                {
                    continue;
                }
                else if (op == BAM_CMATCH)
                {
                    align.spliced = false;
                }
                else if (op == BAM_CREF_SKIP)
                {
                    align.spliced = true;
                }

                x(align, info);
            }
        }
        else
        {
            x(align, info);
        }
    }

    sam_close(f);
}