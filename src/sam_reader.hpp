#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include <string>

/*
 * This class represents a line in a SAM file. Please refer to http://samtools.github.io/hts-specs/SAMv1.pdf
 * for more details.
 *
 */

struct SAMAlignment
{
    // Reference sequence name
    std::string rname;

    // Whether this alignment is mapped
    bool mapped;
};

struct SAMUser
{
    virtual void process(const SAMAlignment &alignment);
};

struct SAMReader
{
    static void read(const std::string &file, SAMUser &user);
};

#endif