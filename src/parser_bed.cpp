#include <fstream>
#include <assert.h>
#include "parser_bed.hpp"
#include <boost/algorithm/string.hpp>

bool ParserBED::parse(const std::string &file, std::function<void(const BedFeature &)> x)
{
    std::string line;
    std::ifstream in(file);

    BedFeature f;
    
    /*
     *    1. Name of the chromosome
     *    2. Starting position
     *    3. Ending position
     *    4. Name of the BED line
     *    5. Score
     *    6. Strand
     *    7. Thick starting position
     *    8. Thick ending position
     *    9. An RGB value
     *   10. The number of blocks (exons) in the BED line
     *   11. A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
     *   12. A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart.
     *       The number of items in this list should correspond to blockCount.
     */
    
    std::vector<std::string> sizes, starts, tokens;
    
    while (std::getline(in, line))
    {
        boost::split(tokens, line, boost::is_any_of("\t"));

        if (tokens.size() < 9)
        {
            continue;
        }
        
        // Clear previous blocks
        f.blocks.clear();
        
        // Name of the chromosome
        f.id = tokens[0];
        
        // Name of the BED line (eg: gene)
        f.name = tokens[3];
        
        f.l.set(stod(tokens[1]), stod(tokens[2]));

        boost::split(sizes,  tokens[10], boost::is_any_of(","));
        boost::split(starts, tokens[11], boost::is_any_of(","));
        assert(sizes.size() == starts.size());
        
        // For each block...
        for (auto i = 0; i < stod(tokens[9]); i++)
        {
            const BasePair start = stod(starts[i]);
            const BasePair size  = stod(sizes[i]);
            
            f.blocks.push_back(Locus(f.l.start + start, f.l.start + start + size));
        }

        x(f);
    }
    
    return true;
}