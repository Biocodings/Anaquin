#ifndef BAM_WRITER_HPP
#define BAM_WRITER_HPP

#include <htslib/sam.h>
#include "tools/samtools.hpp"
#include "writers/writer.hpp"

namespace Anaquin
{
    class BAMWriter : public Writer
    {
        public:

            inline void close() override
            {
                sam_close(_fp);
            }

            inline void open(const FileName &file) override
            {
                _fp = sam_open(file.c_str(), "wb");
            }

            inline void write(const std::string &, bool) override
            {
                throw std::runtime_error("Not implemented");
            }

            inline void write(const ParserBAM::Data &x)
            {
                const auto *b = reinterpret_cast<bam1_t *>(x.b());
                const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());

                if (!_header)
                {
                    std::cout << std::string(h->text);
                }
                
                if (sam_write1(_fp, h, b) == -1)
                {
                    throw std::runtime_error("Failed to SAM record");
                }
                
                std::cout << std::string(_fp->line.s);
                _header = true;
            }

            inline void create(const std::string &) override
            {
                throw std::runtime_error("Not supported in SAMWriter");
            }

        private:

            // Whether the header has written
            bool _header = false;

            // File pointer
            samFile *_fp;
    };
}

#endif
