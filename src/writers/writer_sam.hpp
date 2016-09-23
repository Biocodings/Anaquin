#ifndef WRITER_SAM_HPP
#define WRITER_SAM_HPP

#include <htslib/sam.h>
#include "tools/system.hpp"
#include "tools/samtools.hpp"
#include "writers/writer.hpp"

// Defined in the HTSLIB library
extern "C" int sam_hdr_print(htsFile *fp, const bam_hdr_t *h);

// Defined in the HTSLIB library
extern "C" int __NO_SAM_FILE_WRITING__;

namespace Anaquin
{
    class WriterSAM : public Writer
    {
        public:

            inline void close() override
            {
                if (!_term)
                {
                    sam_close(_fp);
                }
            }

            inline void openTerm()
            {
                // Disable file writing in the htslib library
                __NO_SAM_FILE_WRITING__ = 1;
                
                _fp = sam_open(System::tmpFile().c_str(), "w");
                
                _term = true;
            }
        
            inline void open(const FileName &file) override
            {
                _term = false;
                _fp = sam_open(file.c_str(), "w");
                
                if (!_fp)
                {
                    throw std::runtime_error("Failed to open " + file);
                }
            }

            inline void write(const std::string &, bool) override
            {
                throw std::runtime_error("Not implemented");
            }

            inline void write(const ParserSAM::Data &x)
            {
                const auto *b = reinterpret_cast<bam1_t *>(x.b());
                const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());
                
                if (_term)
                {
                    if (!_header)
                    {
                        sam_hdr_print(_fp, h);
                    }
                    
                    // This won't write anything (__NO_SAM_FILE_WRITING__)
                    if (sam_write1(_fp, h, b) == -1)
                    {
                        throw std::runtime_error("Failed to write in write()");
                    }
                    
                    //bam2print(x);
                    std::cout << _fp->line.s; //<< std::endl;
                }
                else
                {
                    if (!_fp)
                    {
                        throw std::runtime_error("Failed to initialize the file pointer");
                    }

                    if (!_header)
                    {
                        sam_hdr_write(_fp, h);
                    }
                    
                    if (sam_write1(_fp, h, b) == -1)
                    {
                        throw std::runtime_error("Failed to write in write()");
                    }
                }

                _header = true;
            }

            inline void create(const std::string &) override
            {
                throw std::runtime_error("Not supported in WriterSAM");
            }

        private:

            bool _term;

            // Whether the header has been written
            bool _header = false;

            // File pointer
            samFile *_fp;
    };
}

#endif