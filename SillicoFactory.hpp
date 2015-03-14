#ifndef AS_SILLICO_FACTORY_HPP
#define AS_SILLICO_FACTORY_HPP

#include <memory>
#include "Sequence.hpp"

struct SillicoFactory
{
	static std::string transGTF();

	// Returns a in-sillico sequence
	static std::shared_ptr<Sequence> sequence();
};

#endif