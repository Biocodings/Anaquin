#ifndef AS_SILLICO_FACTORY_HPP
#define AS_SILLICO_FACTORY_HPP

#include <memory>
#include "Sequence.hpp"

/*
 * This factory class provides support for the in-sillico chromosome developed by
 * Gavian Institute of Medical Research.
 */

struct SillicoFactory
{
	static std::string transGTF();

	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif