#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

#include <string>
#include "exceptions/invalid_extension.hpp"
#include <boost/algorithm/string/predicate.hpp>

namespace Spike
{
    // Verify the extension of a given file
    inline void verify(const std::string &file, const std::string &ext)
    {
        if (!boost::algorithm::ends_with(file, ext))
        {
           // auto file_extension  = boost::filesystem::extension("");
            //throw InvalidExtension(file, ext, file_extension);
        }
    }
    
    struct ParserOptions
    {
        
    };
}

#endif