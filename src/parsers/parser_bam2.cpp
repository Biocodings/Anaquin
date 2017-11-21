//#include "tools/pool.hpp"
#include "tools/samtools.hpp"
#include "parsers/parser_bam2.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

bool ParserBAM2::Data::nextCigar(Locus &l)
{
    A_ASSERT(_h && _b);
    
    auto t = static_cast<bam1_t *>(_b);    
    const auto cig = bam_get_cigar(t);
    auto _n = t->core.pos;
    
    for (auto i = 0; i < t->core.n_cigar;)
    {
        const auto op = bam_cigar_op(cig[i]);
        const auto ol = bam_cigar_oplen(cig[i]);
        
        switch (op)
        {
            case BAM_CDEL:
            {
                // We'll need it for the next operation
                _n += ol;
                
                continue;
            }
                
            case BAM_CMATCH:
            {
                l.start = _n+1;  // 1-based position
                l.end   = _n+ol; // 1-based position
                
                // We'll need it for the next operation
                _n += ol;
                
                break;
            }
                
            case BAM_CREF_SKIP:
            {
                l.start = _n+1;  // 1-based position
                l.end   = _n+ol; // 1-based position
                
                // We'll need it for the next operation
                _n += ol;
                
                break;
            }
                
            case BAM_CSOFT_CLIP:
            {
                /*
                 *  Eg: 4106431    60    80S35M
                 *
                 *  We should simply ignore the clipping and move to the next cigar
                 */
                
                continue;
            }
                
            case BAM_CINS:
            case BAM_CPAD:
            case BAM_CBACK:
            case BAM_CDIFF:
            case BAM_CEQUAL:
            case BAM_CHARD_CLIP: { continue; }
        }
        
        return true;
    }
    
    return false;
}

std::map<ChrID, Base> ParserBAM2::header(const FileName &file)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }
    
    auto h = sam_hdr_read(f);
    std::map<ChrID, Base> c2b;
    
    for (auto i = 0; i < h->n_targets; i++)
    {
        c2b[std::string(h->target_name[i])] = h->target_len[i];
    }
    
    return c2b;
}

void ParserBAM2::parse(const FileName &file, Functor x)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }

    auto t = bam_init1();
    auto h = sam_hdr_read(f);

    Info info;
    Data align;

    while (sam_read1(f, h, t) >= 0)
    {
        align.mapped = false;
        
        align._b = info.b = t;
        align._h = info.h = h;

        align.mapq = t->core.qual;
        align.flag = t->core.flag;
        
        const auto hasCID = t->core.tid >= 0;
        
        #define isPairedEnd(b)    (((b)->core.flag&0x1)   != 0)
        #define isAllAligned(b)   (((b)->core.flag&0x2)   != 0)
        #define isUnmapped(b)     (((b)->core.flag&0x4)   != 0)
        #define isMateUnmapped(b) (((b)->core.flag&0x8)   != 0)
        #define isMateReversed(b) (((b)->core.flag&0x20)  != 0)
        #define isFirstPair(b)    (((b)->core.flag&0x40)  != 0)
        #define isLastPair(b)     (((b)->core.flag&0x80)  != 0)
        #define isSecondary(b)    (((b)->core.flag&0x100) != 0)
        #define isFailed(b)       (((b)->core.flag&0x200) != 0)
        #define isDuplicate(b)    (((b)->core.flag&0x400) != 0)
        #define isSupplement(b)   (((b)->core.flag&0x800) != 0)
        #define isPrimary(b)      (((b)->core.flag&0x900) == 0)

        align.isPaired      = isPairedEnd(t);
        align.isAllAligned  = isAllAligned(t);
        align.isAligned     = !isUnmapped(t);
        align.isMateAligned = !isMateUnmapped(t);
        align.isMateReverse = isMateReversed(t);
        align.isFirstPair   = isFirstPair(t);
        align.isSecondPair  = isLastPair(t);
        align.isPassed      = !isFailed(t);
        align.isDuplicate   = isDuplicate(t);
        align.isSupplement  = isSupplement(t);
        align.isPrimary     = isPrimary(t);
        align.isSecondary   = isSecondary(t);

        if (hasCID)
        {
            align.cID = std::string(h->target_name[t->core.tid]);
            align.nextCigar(align.l);            
        }
        else
        {
            align.cID = "*";
            align.l.start = align.l.end = 0;
        }

        align.rnext = hasCID ? bam2rnext(h, t) : "*";
        
        if (align.rnext == "=")
        {
            align.rnext = align.cID;
        }

        align.mapped = hasCID && !(t->core.flag & BAM_FUNMAP);

        x(align, info);
        info.p.i++;
    }

    bam_destroy1(t);
    sam_close(f);
}
