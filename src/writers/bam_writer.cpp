#include "writers/bam_writer.hpp"

using namespace Anaquin;

void BAMWriter::close()
{
    sam_close(_fp);
}
        
void BAMWriter::open(const FileName &file)
{
    _fp = sam_open(file.c_str(), "wb");
}

void BAMWriter::write(const ParserBAM::Data &x)
{
    const auto *b = reinterpret_cast<bam1_t *>(x.b());
    const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());
    
    static bool header = false;
    
    if (!header && sam_hdr_write(_fp, reinterpret_cast<bam_hdr_t *>(x.h())) == -1)
    {
        throw std::runtime_error("sam_hdr_write failed");
    }
    
    header = true;
    
    if (sam_write1(_fp, h, b) == -1)
    {
        throw std::runtime_error("sam_write1 failed");
    }
}
