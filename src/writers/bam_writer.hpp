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
            void write(const ParserBAM::Data &);

        private:
            bool _head;        
            samFile *_fp;
    };
}

#endif
