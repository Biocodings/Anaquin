#ifndef WRITER_SAM_HPP
#define WRITER_SAM_HPP

#include <htslib/sam.h>
#include "writers/writer.hpp"

namespace Anaquin
{
    class WriterSAM : public Writer
    {
        public:

            inline void close() override { sam_close(_fp); }

            inline void open(const FileName &file) override
            {
                _fp = sam_open(file.c_str(), "w");
            }

            inline void write(const std::string &, bool) override
            {
                throw std::runtime_error("Not implemented");
            }

            inline void write(const bam_hdr_t *header, const bam1_t *read)
            {
                if (!_header)
                {
                    sam_hdr_write(_fp, header);
                }
                
                _header = true;
                
                if (sam_write1(_fp, header, read) == -1)
                {
                    throw std::runtime_error("Failed to write");
                }
            }

            inline void create(const std::string &) override
            {
                throw std::runtime_error("Not implemented");
            }

        private:

            // Whether the header has been written
            bool _header;

            // File pointer
            samFile *_fp;
    };
}

#endif