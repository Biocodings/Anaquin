#ifndef PARSER_BAM2_HPP
#define PARSER_BAM2_HPP

#include <htslib/sam.h>
#include "tools/samtools.hpp"
#include "data/alignment.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserBAM2
    {
        struct Info
        {
            ParserProgress p;

            void *b;
            void *h;
        };

        class Data
        {
            friend class ParserBAM;
            friend class ParserBAM2;
            
            public:
            
                operator const Locus &() const { return l; }
            
                // Primary alignment
                ChrID cID;
            
                // Location of the alignment
                Locus l;
            
                // If this field is false, no assumption can be made to other fields
                bool mapped;
            
                // Mapping quality
                int mapq;
            
                // Bitwise FLAG
                int flag;
            
                /*
                 * SAM flag fields
                 */
            
                bool isPaired;
                bool isAllAligned;
                bool isAligned;
                bool isMateAligned;
                bool isMateReverse;
                bool isFirstPair;
                bool isSecondPair;
                bool isDuplicate;
                bool isPrimary;
                bool isSupplement;
                bool isPassed;
            
                // Secondary alignment? Typically used for alternative mappings when multiple mappings are presented
                bool isSecondary;
            
                /*
                 * Optional fields
                 */
            
                // Reference sequence name of the primary alignment
                ChrID rnext;

                /*
                 * Optional fields
                 */
            
                // Eg: B7_591:6:155:12:674
                inline ReadName name() { return bam_get_qname(static_cast<bam1_t *>(_b)); }

                // Segment sequence (optional)
                inline std::string seq() { return bam2seq(static_cast<bam1_t *>(_b)); }

                // ASCII of base QUALity (optional)
                inline std::string qual() { return bam2qual(static_cast<bam1_t *>(_b)); }

                bool nextCigar(Locus &);
            
                inline void *b() const { return _b; }
                inline void *h() const { return _h; }

            private:
            
                void *_b;
                void *_h;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        
        static std::map<ChrID, Base> header(const FileName &);
        
        static void parse(const FileName &, Functor);
    };
}

#endif
