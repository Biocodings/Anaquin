#ifndef AS_PARSER_CSV_HPP
#define AS_PARSER_CSV_HPP

#include <vector>
#include <functional>

struct ParserCSV
{
	static bool parse(const std::string &file, std::function<void (const std::vector<std::string> &)>);
};

#endif