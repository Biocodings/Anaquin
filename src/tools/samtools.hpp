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

    inline std::vector<int> bam2delta(bam1_t *x)
    {
        std::vector<int> d;
        const auto t = bam_get_cigar(x);
        
        for (int i = x->core.n_cigar - 1; i >= 0; i--)
        {
            const auto val = bam_cigar_oplen(t[i]);
            
            switch (bam_cigar_op(t[i]))
            {
                case BAM_CPAD:
                case BAM_CDIFF:
                case BAM_CEQUAL:
                case BAM_CMATCH:
                case BAM_CSOFT_CLIP:  { d.push_back(0);    break; }
                case BAM_CINS:        { d.push_back(-val); break; }
                case BAM_CDEL:        { d.push_back(val);  break; }
                case BAM_CREF_SKIP:   { d.push_back(-val); break; }
                case BAM_CHARD_CLIP:  { d.push_back(val);  break; }
            }
        }
        
        return d;
    }
    
    inline Base reversePos(const Locus &l, ParserSAM::Data &x, const ParserSAM::Info &i)
    {
        // Length of the chromosome
        const auto clen = i.length;

        // Length of the sequence
        const auto slen = x.seq.size();
        
        // Insertion, deletion etc
        const auto delta = sum(bam2delta(reinterpret_cast<bam1_t *>(i.data)));
        
        return clen - (l.start + slen + delta) + 2;
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

    // Reverse an alignment
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

        x.cigar = bam2rcigar(b);

        std::reverse(x.seq.begin(),  x.seq.end());
        std::reverse(x.qual.begin(), x.qual.end());
        
        // The left-most position in the forward strand
        const auto rstart = reversePos(x.l, x, i);
        
        // The right-most positiion in the forward strand
        const auto rend = rstart + x.l.length() - 1;

        const auto t = x.l;
        
        // New position on the forward strand
        x.l = Locus(rstart, rend);
        
        assert(t.length() == x.l.length());
    }
}

#endif