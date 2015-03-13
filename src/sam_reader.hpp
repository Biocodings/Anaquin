#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include <string>
#include <functional>

/*
 * This class represents a line in a SAM file. Please refer to http://samtools.github.io/hts-specs/SAMv1.pdf
 * for more details.
 */

struct Alignment
{
    std::string name;

    // Whether this alignment is mapped
    bool mapped;
};

struct SAMReader
{
	static bool read(const std::string &file, std::function<void(const Alignment &alignment)> x);
};

#endif