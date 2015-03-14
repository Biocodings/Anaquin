#ifndef AS_ALIGNMENT_HPP
#define AS_ALIGNMENT_HPP

#include <string>

/*
 * This class represents a sequencing alignment.
 */

struct Alignment
{
    std::string id;

	/*
	 * It's important to note that not every alignment is mapped. If this field is false, no assumption
	 * can be made to other fields.
	 */

    bool mapped;

	std::string seq;
};

#endif