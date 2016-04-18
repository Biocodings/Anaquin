#ifndef DATA_HPP
#define DATA_HPP

#include <string>

namespace Anaquin
{
    /*
     * Represents a mathced element that can be identified
     */
    
    struct Matched
    {
        virtual std::string name() const = 0;
    };
}

#endif