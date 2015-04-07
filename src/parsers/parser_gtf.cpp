#include <vector>
#include <fstream>
#include <assert.h>
#include "tokens.hpp"
#include "parser_gtf.hpp"

using namespace Spike;

bool ParserGTF::parse(const std::string &file, std::function<void (const Feature &, ParserProgress &)> x)
{
    std::map<std::string, FeatureType> mapper =
    {
        { "exon", Exon },
        { "CDS",  CDS  },
        { "start_codon", StartCodon },
        { "stop_codon",  StopCodon  },
        { "transcript",  Transcript }
    };

    std::string line;
    std::ifstream in(file);

    Feature f;

	/*
	 * Fields must be tab-separated. Also, all but the final field in each feature line must contain a value;
	 * "empty" columns should be denoted with a '.'. Please refer to the online documentation for more details.
	 *
	 *    1. seqname
	 *    2. sourcre
	 *    3. feature
	 *    4. start
	 *    5. end
	 *    6. score
	 *    7. strand
	 *    8. frame
	 *    9. attribute
	 */
    
    std::vector<std::string> tokens;
    std::vector<std::string> options;
    std::vector<std::string> nameValue;

    ParserProgress p;
    
    while (std::getline(in, line))
    {
        p.i++;
        boost::split(tokens, line, boost::is_any_of("\t"));

        f.id = tokens[0];
        f.l.set(stoi(tokens[3]), stoi(tokens[4]));

        if (!mapper.count(tokens[2]))
        {
            throw std::runtime_error("Unknown feature type: " + tokens[2]);
        }
        
        f.type = mapper[tokens[2]];

        /*
         * Eg: "gene_id "R_5_3"; transcript_id "R_5_3_R";"
         */
        
        boost::split(options, tokens[8], boost::is_any_of(";"));
        
        for (auto option : options)
        {
            if (!option.empty())
            {
                boost::trim(option);
                boost::split(nameValue, option, boost::is_any_of(" "));

                if (nameValue.size() == 2)
                {
                    // Make sure that silly characters are removed
                    nameValue[1].erase(std::remove(nameValue[1].begin(), nameValue[1].end(), '\"'), nameValue[1].end());
                    
                    if (nameValue[0] == "gene_id")
                    {
                        f.geneID = nameValue[1];
                    }
                    else if (nameValue[0] == "transcript_id")
                    {
                        f.iID = nameValue[1];
                    }
                    
                    if (!f.geneID.empty() && !f.iID.empty())
                    {
                        assert(f.geneID != f.iID);
                    }
                    
                    //f.options[nameValue[0]] = nameValue[1];
                }
            }
        }

        x(f, p);
        
        if (p.terminate)
        {
            break;
        }
	}
    
	return true;
}