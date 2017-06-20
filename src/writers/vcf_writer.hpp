#ifndef VCF_WRITER_HPP
#define VCF_WRITER_HPP

#include "data/data.hpp"
#include "writers/writer.hpp"

namespace Anaquin
{
    class VCFWriter : public Writer
    {
        public:

            void close() override;

            void open(const FileName &) override;

            void write(const std::string &, bool) override {}

            void write(void *hdr, void *line);
        
            inline void create(const std::string &) override
            {
                throw std::runtime_error("Not supported in SAMWriter");
            }

        private:

            bool _head;
            void *_fp;
    };
}

#endif
