#ifndef HTSLIB_HPP
#define HTSLIB_HPP

#include <map>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <htslib/sam.h>
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    typedef std::string CigarStr;
    
    // Avoid excessive heap allocation...
    static std::stringstream __buf__;
    
    #define CLEAR_HTSLIB() { __buf__.str(""); }
    
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
        CLEAR_HTSLIB();
        
        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            __buf__ << (char) (bam_get_qual(x)[i] + 33);
        }
        
        return __buf__.str();
    }
    
    inline std::string bam2seq(bam1_t *x)
    {
        CLEAR_HTSLIB();

        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            __buf__ << seq_nt16_str[bam_seqi(bam_get_seq(x),i)];
        }

        return __buf__.str();
    }
    
    inline CigarStr bam2cigar(bam1_t *x)
    {
        CLEAR_HTSLIB();
        const auto t = bam_get_cigar(x);

        for (auto i = 0; i < x->core.n_cigar; i++)
        {
            __buf__ << std::to_string(bam_cigar_oplen(t[i]));
            __buf__ << bam2char.at(bam_cigar_op(t[i]));
        }
        
        return __buf__.str();
    }
    
    inline CigarStr bam2rcigar(bam1_t *x)
    {
        CLEAR_HTSLIB();
        const auto t = bam_get_cigar(x);

        for (int i = x->core.n_cigar - 1; i >= 0; i--)
        {
            __buf__ << std::to_string(bam_cigar_oplen(t[i]));
            __buf__ << bam2char.at(bam_cigar_op(t[i]));
        }
        
        return __buf__.str();
    }

    inline void reverse(ParserSAM::Data &x, const ParserSAM::Info &i)
    {
        /*
         * The following needs to be reversed:
         *
         *   1. POS
         *   2. CIGAR
         *   3. PNEXT
         *   4. SEQ
         *   5. QUAL
         */

        const auto b = reinterpret_cast<bam1_t *>(i.data);
 
        // Length of the chromosome
        const auto clen = i.length;

        const auto t = x.l;
        
        // New position on the forward strand
        x.l = Locus(clen - x.l.start, clen - x.l.start + x.l.length() - 1);
        
        assert(t.length() == x.l.length());
        
        x.cigar = bam2rcigar(b);
        x.pnext = clen - x.pnext;
        
        std::reverse(x.seq.begin(),  x.seq.end());
        std::reverse(x.qual.begin(), x.qual.end());
    }
}

#endif