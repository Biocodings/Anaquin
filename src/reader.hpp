#ifndef GI_READER_HPP
#define GI_READER_HPP

#include <vector>
#include <string>
#include "parsers/parser.hpp"

namespace Spike
{
    struct ReaderInternal;

    class Reader
    {
        public:
            Reader(const std::string &file, ParserMode mode = File);
            ~Reader();

            // Returns the next line in the file
            bool nextLine(std::string &) const;

            // Returns the next line and parse it into tokens
            bool nextTokens(std::vector<std::string> &, const std::string &c) const;

        private:
            ReaderInternal *_imp;
    };
}

#endif