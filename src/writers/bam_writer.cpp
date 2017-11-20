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
