#include <assert.h>
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Anaquin;

void ParserGTF::parse(const Reader &r, Callback x)
{
    std::map<std::string, RNAFeature> mapper =
    {
        { "exon", Exon },
        { "CDS",  CDS  },
        { "start_codon", StartCodon },
        { "stop_codon",  StopCodon  },
        { "transcript",  Transcript }
    };
    
    std::string line;
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
    
    ParserProgress p;
    
    std::vector<std::string> tokens;
    std::vector<std::string> options;
    std::vector<std::string> nameValue;

    while (r.nextLine(line))
    {
        p.i++;
        boost::split(tokens, line, boost::is_any_of("\t"));

        // Empty line? Unknown feature such as mRNA?
        if (tokens.size() == 1 || !mapper.count(tokens[2]))
        {
            continue;
        }
        
        f.id   = tokens[0];
        f.l    = Locus(stoi(tokens[3]), stoi(tokens[4]));
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
                    else if (nameValue[0] == "nearest_ref" || nameValue[0] == "transcript_id")
                    {
                        f.tID = nameValue[1];
                    }
                }
            }
        }

        x(f, line, p);
    }
}