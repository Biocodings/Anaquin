#ifndef HTSLIB_HPP
#define HTSLIB_HPP

#include <map>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <htslib/sam.h>
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    typedef std::string CigarStr;
    
    static std::map<int, char> bam2char =
    {
        { BAM_CMATCH,     'M'  },
        { BAM_CINS,       'I'  },
        { BAM_CDEL,       'D'  },
        { BAM_CREF_SKIP,  'N'  },
        { BAM_CSOFT_CLIP, 'S', },
        { BAM_CHARD_CLIP, 'H', },
        { BAM_CPAD,       'P', },
        { BAM_CEQUAL,     '=', },
        { BAM_CDIFF,      'X', },
    };

    inline std::string bam2rnext(bam_hdr_t *h, bam1_t *b)
    {
        const auto cID = std::string(h->target_name[b->core.tid]);
        
        if (b->core.mtid == -1)
        {
            return "";
        }

        const auto rID = std::string(h->target_name[b->core.mtid]);
        
        if (rID == cID)
        {
            return "=";
        }

        assert(!rID.empty());
        return rID;
    }
    
    inline std::string bam2qual(bam1_t *x)
    {
        std::stringstream buf;
        
        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            buf << (char) (bam_get_qual(x)[i] + 33);
        }
        
        return buf.str();
    }
    
    inline std::string bam2seq(bam1_t *x)
    {
        std::stringstream buf;

        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            buf << seq_nt16_str[bam_seqi(bam_get_seq(x),i)];
        }

        return buf.str();
    }
    
    inline CigarStr bam2cigar(bam1_t *x)
    {
        std::stringstream buf;
        const auto t = bam_get_cigar(x);

        for (auto i = 0; i < x->core.n_cigar; i++)
        {
            buf << std::to_string(bam_cigar_oplen(t[i]));
            buf << bam2char.at(bam_cigar_op(t[i]));
        }
        
        return buf.str();
    }
}

#endif
