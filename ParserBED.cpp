#include <fstream>
#include <assert.h>
#include "ParserBED.hpp"
#include <boost/algorithm/string.hpp>

using namespace std;

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
    
    std::vector<std::string> tokens;
    std::vector<std::string> options;
    std::vector<std::string> nameValue;
    
    while (std::getline(in, line))
    {
        boost::split(tokens, line, boost::is_any_of("\t"));
        
        f.chromo = tokens[0];
        f.loc.update(stod(tokens[1]), stod(tokens[2]));
        
        // For each block...
        for (auto i = 0; i < stod(tokens[9]); i++)
        {
            if (i)
            {
                
                
            }
        }

        f.chromo = tokens[0];
        f.loc.update(stoi(tokens[3]), stoi(tokens[4]));

        if (tokens[2] == "exon")
        {
            f.type = Exon;
        }
        else if (tokens[2] == "CDS")
        {
            f.type = CDS;
        }
        else if (tokens[2] == "start_codon")
        {
            f.type = StartCodon;
        }
        
        boost::split(options, tokens[8], boost::is_any_of(";"));
        
        /*
         * Eg: "gene_id "R_5_3_R"; transcript_id "R_5_3_R";"
         */
        
        for (auto option : options)
        {
            if (!option.empty())
            {
                boost::trim(option);
                const auto &t = boost::split(nameValue, option, boost::is_any_of(" "));

                if (t.size() == 2)
                {
                    f.options[t[0]] = t[1];
                }
            }
        }
        
        x(f);
    }
    
    return true;
}