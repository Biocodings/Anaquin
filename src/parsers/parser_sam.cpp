#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "parser_sam.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

void ParserSAM::parse(const std::string &file, std::function<void (const Alignment &)> x)
{
    std::string line;
    std::ifstream in(file);

    Alignment align;
    
    while (std::getline(in, line))
    {
		if (line.empty() || line[0] == '@')
		{
			continue;
		}

        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("\t"));

        // It's unmapped if the bit is 1
        align.mapped = !(stoi(tokens[1]) & (1 << 2));

        align.id  = tokens[2];
		align.seq = tokens[9];

        const auto cigar = tokens[5];        
        
        // Eg: 100M58378N1M
        align.spliced = (cigar.find('N') != std::string::npos);
        
        if (align.spliced)
        {
            char const* delims = "MN";
            
            std::vector<std::string> result;
            boost::algorithm::split(result, cigar, boost::is_any_of(delims));
            
            const auto r1 = stoi(result[0]);
            const auto r2 = stoi(result[1]);
            const auto r3 = stoi(result[2]);

            const auto start = stoi(tokens[3]);
            
            align.l.set(start, start + r1 + r2 + r3 - 1);
        }
        else
        {
            const auto len = static_cast<BasePair>(tokens[9].size());
            const auto start = stoi(tokens[3]);
            
            align.l.set(start, start + len - 1);

            // Don't forget the minus one, since we're "counting fenceposts"
            assert(len == align.l.length() + 1);
        }
        
        x(align);
    }
}