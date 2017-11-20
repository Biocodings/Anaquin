#ifndef BAM_WRITER_HPP
#define BAM_WRITER_HPP

#include <htslib/sam.h>
#include "tools/samtools.hpp"

namespace Anaquin
{
    class BAMWriter
    {
        public:

            void close();
            void open(const FileName &);

            template <typename T> void write(const T &x)
            {
                const auto *b = reinterpret_cast<bam1_t *>(x.b());
                const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());
                
                if (!_head && sam_hdr_write(_fp, reinterpret_cast<bam_hdr_t *>(x.h())) == -1)
                {
                    throw std::runtime_error("sam_hdr_write failed");
                }
                
                _head = true;
                
                if (sam_write1(_fp, h, b) == -1)
                {
                    throw std::runtime_error("sam_write1 failed");
                }
            }

        private:
            bool _head;        
            samFile *_fp;
    };
}

#endif
