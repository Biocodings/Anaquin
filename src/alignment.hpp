#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string>

struct Alignment
{
    std::string name;

    // Whether this alignment is mapped
    bool mapped;

	std::string seq;
};

#endif