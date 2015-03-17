#ifndef AS_ALIGNMENT_HPP
#define AS_ALIGNMENT_HPP

#include "Types.hpp"
#include "Locus.hpp"

/*
 * This class represents a sequencing alignment.
 */

struct Alignment
{
    std::string id;

    // If this field is false, no assumption can be made to other fields
    bool mapped;

    Locus loc;
	std::string seq;
};

#endif